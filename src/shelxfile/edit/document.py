"""The editing façade over a :class:`Shelxfile`.

``ShelxDocument`` is the **only** object that reaches into the private
``_reslist``.  Views ask it what sits on a line, tell it what the user
did, and render whatever it returns.  It holds no Qt and emits no Qt
signals: observers are plain callables (plan decision D-8).

Cascade direction matters here (**D-9**).  Deleting an *atom* may remove
cards and brackets.  Removing a *card* never deletes atoms, no matter how
the card came to be removed -- losing a restraint is a refinement
decision, not a statement about the atoms it mentioned.  The destructive
variant exists but has to be named explicitly:
:meth:`delete_restraint_with_atoms`.
"""

from __future__ import annotations

import functools
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Iterable, Iterator, NamedTuple, Sequence, Union

from shelxfile.edit import structure_edits
from shelxfile.edit.cascade import CascadeEngine, CascadePlan
from shelxfile.edit.eqiv_cleanup import EqivCleaner, validate_symmetry_arity
from shelxfile.edit.eqiv_factory import EqivFactory
from shelxfile.edit.graph import AtomRestraintGraph
from shelxfile.edit.history import EditHistory, HistoryState, UndoResult
from shelxfile.edit.line_map import RenderedFile, render
from shelxfile.edit.reports import (
    DeletionReport,
    EditReport,
    RemovalReason,
    RenameReport,
)
from shelxfile.edit.token_resolver import (
    LAST_KEYWORD,
    AtomTokenResolver,
    split_range_tokens,
)

if TYPE_CHECKING:
    from shelxfile import Shelxfile
    from shelxfile.atoms.atom import Atom
    from shelxfile.shelx.cards import Command, Restraint

#: An observer is called with the document after every successful edit.
Observer = Callable[['ShelxDocument'], None]

ResListItem = Union['Atom', 'Restraint', 'Command', str]

#: Returns whatever a caller wants stored with every history snapshot.
StateProvider = Callable[[], Any]


def _undoable(label: str):
    """Run the decorated edit as one undo step labelled *label*.

    A call that changes nothing (a refused rename, an empty deletion)
    leaves no step behind, because the step is only kept when the edit
    reached :meth:`ShelxDocument._notify`.
    """
    def decorate(method):
        @functools.wraps(method)
        def wrapper(self: ShelxDocument, *args, **kwargs):
            with self.batch(label):
                return method(self, *args, **kwargs)
        return wrapper
    return decorate


class ParseAttempt(NamedTuple):
    """Outcome of :meth:`ShelxDocument.try_from_string`.

    :param document: the parsed document, or ``None`` on failure.
    :param error: the failure message, or ``None`` on success.
    :param error_line: 0-based line the parser choked on, ``-1`` if
        unknown.
    """

    document: ShelxDocument | None
    error: str | None
    error_line: int

    def __bool__(self) -> bool:
        return self.document is not None


def restraint_keywords() -> list[str]:
    """The restraint instructions :meth:`ShelxDocument.add_restraint` accepts.

    Exposed so a view can offer the vocabulary without importing the
    model's card registry.
    """
    from shelxfile import Shelxfile
    return sorted(Shelxfile.RESTRAINT_CARD_CLASSES)


class ShelxDocument:
    """A :class:`Shelxfile` plus the operations an editor needs.

    :param shx: the model to wrap.  The document does not copy it; both
        refer to the same :class:`Shelxfile`.
    """

    def __init__(self, shx: Shelxfile) -> None:
        self._shx = shx
        self._observers: list[Observer] = []
        self._rendered: RenderedFile | None = None
        self._graph: AtomRestraintGraph | None = None
        self._history = EditHistory()
        self._state_provider: StateProvider | None = None
        self._batch_depth = 0

    # ------------------------------------------------------------- model

    @property
    def graph(self) -> AtomRestraintGraph:
        """Atom/card links, rebuilt on demand when the model changes."""
        if self._graph is None:
            self._graph = AtomRestraintGraph(self._shx)
        self._graph.ensure_current()
        return self._graph

    @property
    def shelxfile(self) -> Shelxfile:
        """The wrapped model, for read-only use by callers."""
        return self._shx

    @classmethod
    def from_file(cls, path: str | Path, **kwargs) -> ShelxDocument:
        from shelxfile import Shelxfile
        shx = Shelxfile(**kwargs)
        shx.read_file(str(path))
        return cls(shx)

    @classmethod
    def from_string(cls, text: str, **kwargs) -> ShelxDocument:
        from shelxfile import Shelxfile
        shx = Shelxfile(**kwargs)
        shx.read_string(text)
        return cls(shx)

    @classmethod
    def try_from_string(cls, text: str) -> ParseAttempt:
        """Parse *text*, reporting failure instead of raising.

        Meant for editors, where the user is mid-edit and a broken file is
        an expected state rather than a programming error.  Parsing runs
        with ``debug=True`` so SHELXL problems surface as messages instead
        of being silently tolerated.

        A file without a ``CELL`` is treated as a failure: nothing
        downstream can work without a unit cell, and accepting it would
        replace a good model with an unusable one.
        """
        from shelxfile import Shelxfile
        shx = Shelxfile(debug=True)
        try:
            shx.read_string(text)
        except Exception as exc:  # noqa: BLE001 - any parse failure is a UI message
            message = str(exc) or type(exc).__name__.replace('_', ' ')
            return ParseAttempt(None, message, shx.error_line_num)
        if not shx.cell:
            return ParseAttempt(
                None,
                'Could not parse SHELX file (missing CELL instruction?).',
                shx.error_line_num,
            )
        return ParseAttempt(cls(shx), None, -1)

    def unused_atom_name(self, element: str) -> str:
        """A free atom name for *element*, e.g. ``'C12'``."""
        return self._shx.unused_atom_name(element)

    # --------------------------------------------------------- rendering

    @property
    def text(self) -> str:
        """The file as text, identical to :meth:`Shelxfile.dumps`."""
        return self._render().text

    @property
    def line_count(self) -> int:
        return len(self._render().line_to_index)

    def _render(self) -> RenderedFile:
        if self._rendered is None:
            self._rendered = render(self._shx)
        return self._rendered

    def invalidate(self) -> None:
        """Drop the cached rendering after an external change to the model."""
        self._rendered = None

    # ------------------------------------------------------ line lookups

    def item_at_line(self, line: int) -> ResListItem | None:
        """The ``_reslist`` entry that produced 0-based text *line*."""
        line_to_index = self._render().line_to_index
        if not 0 <= line < len(line_to_index):
            return None
        return self._shx._reslist[line_to_index[line]]

    def items_in_lines(self, start: int, end: int) -> list[ResListItem]:
        """Distinct entries covered by the inclusive line range.

        Order follows the file.  An entry wrapped over several lines is
        returned once.
        """
        seen: set[int] = set()
        items: list[ResListItem] = []
        for line in range(start, end + 1):
            item = self.item_at_line(line)
            if item is None or id(item) in seen:
                continue
            seen.add(id(item))
            items.append(item)
        return items

    def line_of(self, item: ResListItem) -> int | None:
        """First 0-based text line rendered from *item*, if it is visible."""
        try:
            index = self._shx.index_of(item)
        except ValueError:
            return None
        line_to_index = self._render().line_to_index
        for line, origin in enumerate(line_to_index):
            if origin == index:
                return line
        return None

    def line_of_atom(self, fullname_short: str) -> int | None:
        """Line of the atom named like :attr:`Atom.fullname_short`."""
        wanted = fullname_short.upper()
        for atom in self._shx.atoms:
            if atom.fullname_short.upper() == wanted:
                return self.line_of(atom)
        return None

    # ------------------------------------------------------- observation

    def subscribe(self, callback: Observer) -> None:
        """Register *callback*, invoked after every successful edit."""
        if callback not in self._observers:
            self._observers.append(callback)

    def unsubscribe(self, callback: Observer) -> None:
        if callback in self._observers:
            self._observers.remove(callback)

    def _notify(self) -> None:
        if self._batch_depth:
            # Inside a batch: the observers are told once, when the
            # outermost batch ends.
            return
        self.invalidate()
        for callback in list(self._observers):
            callback(self)

    # ----------------------------------------------------- undo and redo

    @property
    def history(self) -> EditHistory:
        """The undo/redo stacks, for views that list the steps."""
        return self._history

    def set_state_provider(self, provider: StateProvider | None) -> None:
        """Store *provider()*'s result with every history snapshot.

        A caller that keeps state of its own next to the document registers
        a provider; :meth:`undo` and :meth:`redo` hand the value stored with
        the restored state back in :attr:`UndoResult.payload`.
        """
        self._state_provider = provider

    def _capture(self) -> HistoryState:
        payload = self._state_provider() if self._state_provider is not None else None
        return self._history.capture(self._shx.dumps(), payload)

    @contextmanager
    def batch(self, label: str = 'Edit') -> Iterator[ShelxDocument]:
        """Group every edit made inside the block into one undo step.

        Observers are told once, when the outermost block ends.  Nested
        blocks fold into the outer one.  If the block raises, the model is
        restored to its state before the block and the exception
        propagates, so a composite edit is never left half-applied.

        Objects taken from the model before a failed block are stale
        afterwards, because the restore re-reads the file text.
        """
        if self._batch_depth:
            self._batch_depth += 1
            try:
                yield self
            finally:
                self._batch_depth -= 1
            return

        before = self._capture()
        self._batch_depth = 1
        try:
            yield self
        except BaseException:
            self._batch_depth = 0
            self._restore(before)
            self._notify()
            raise
        self._batch_depth = 0
        if self._shx.dumps() != before.text:
            self._history.push(label, before)
            self._notify()

    def _restore(self, state: HistoryState) -> None:
        """Replace the model with the one serialised in *state*."""
        from shelxfile import Shelxfile

        old = self._shx
        shx = Shelxfile(debug=old.debug, verbose=old.verbose)
        shx.read_string(state.text)
        shx.resfile = old.resfile
        shx.encoding = old.encoding
        self._shx = shx
        self._graph = None
        self.invalidate()

    @property
    def can_undo(self) -> bool:
        return self._history.can_undo

    @property
    def can_redo(self) -> bool:
        return self._history.can_redo

    @property
    def is_modified(self) -> bool:
        """Whether the model differs from the state last written or loaded."""
        return self._history.is_modified

    def mark_saved(self) -> None:
        """Declare the current state to be the one on disk."""
        self._history.mark_saved()

    def undo(self) -> UndoResult | None:
        """Revert the most recent step; ``None`` when there is none.

        The model is re-read from the stored text, so :attr:`shelxfile` is
        a **new** object afterwards and any atom or card taken from the old
        one must be looked up again (by name).
        """
        if self._batch_depth:
            raise RuntimeError('Cannot undo inside a batch')
        if not self._history.can_undo:
            return None
        label, state = self._history.pop_undo(self._capture())
        self._restore(state)
        self._notify()
        return UndoResult(label, state.payload)

    def redo(self) -> UndoResult | None:
        """Re-apply the most recently undone step; ``None`` when there is none."""
        if self._batch_depth:
            raise RuntimeError('Cannot redo inside a batch')
        if not self._history.can_redo:
            return None
        label, state = self._history.pop_redo(self._capture())
        self._restore(state)
        self._notify()
        return UndoResult(label, state.payload)

    # ---------------------------------------------------------- deletion

    @_undoable('Delete atoms')
    def delete_atoms(self, atoms: Iterable[Atom]) -> DeletionReport:
        """Delete *atoms* and everything that cannot survive without them.

        The consequences are worked out in full before anything is
        touched, so a cascade is never half-applied.  Restraints that lose
        an atom are trimmed, or removed when too little is left to mean
        anything; ``AFIX`` groups that fall below their required member
        count go whole, taking their atoms with them.

        Symmetry-generated atoms are skipped: they are products of
        :meth:`Shelxfile.grow`/:meth:`~Shelxfile.pack`, not lines in the
        file.
        """
        engine = CascadeEngine(self._shx, self.graph)
        plan = engine.plan(atoms)
        if plan.is_empty:
            return DeletionReport()
        report = engine.apply(plan)
        self._collect_orphaned_eqivs(report)
        self._notify()
        return report

    def plan_deletion(self, atoms: Iterable[Atom]) -> CascadePlan:
        """Work out what deleting *atoms* would do, without doing it.

        Useful for confirming a destructive edit before it happens.
        """
        return CascadeEngine(self._shx, self.graph).plan(atoms)

    # ---------------------------------------------------------- renaming

    @_undoable('Rename atom')
    def rename_atom(self, atom: Atom, new_name: str) -> RenameReport:
        """Rename *atom* and follow it through the instructions.

        A token is rewritten only when it refers to this atom **and
        nothing else**.  On a residue-class-scoped card such as
        ``SADI_TOL C1 C2``, the token ``C1`` stands for that atom in every
        ``TOL`` residue, so rewriting it would silently redirect all of
        them.  Those references are left alone and listed in
        :attr:`RenameReport.skipped`; SHELXL then ignores the instruction
        for this residue, which is its normal behaviour for a name it
        cannot find.

        Atoms in the middle of a range need no change at all -- a range is
        positional -- but an endpoint does.
        """
        report = RenameReport(atom=atom, old_name=atom.fullname)
        if atom.symmgen:
            report.error = 'symmetry-generated atoms are not part of the file'
            return report
        problem = atom.name_problem(new_name)
        if problem:
            report.error = problem
            return report
        if new_name.strip().upper() == atom.name.upper():
            report.new_name = atom.fullname
            return report

        old_fullname = atom.fullname
        affected = self._rename_targets(atom, old_fullname)

        atom.name = new_name.strip()
        report.new_name = atom.fullname

        resolver = AtomTokenResolver(self._shx)
        for card, exclusive, shared in affected:
            before = str(card)
            if exclusive:
                self._retarget(card, exclusive, atom.name)
                report.add_update(card, before)
            for token in shared:
                if self._still_covers(resolver, card, token, atom.fullname):
                    # An element reference such as '$H' names a class of
                    # atoms, so it follows the rename by itself.
                    continue
                report.add_skip(
                    card, token,
                    'token also names atoms in other residues, so '
                    'rewriting it would redirect them too',
                )
        self._notify()
        return report

    @staticmethod
    def _still_covers(resolver, card, token: str, fullname: str) -> bool:
        """Whether *token* reaches *fullname* after the rename."""
        scope = getattr(card, 'residue_number', None)
        if isinstance(scope, int):
            scope = [scope]
        return fullname in resolver.resolve([token], scope).fullnames

    def _rename_targets(self, atom: Atom, old_fullname: str):
        """Split each referencing card's tokens into exclusive and shared.

        *exclusive* tokens name this atom and nothing else, so they may be
        rewritten; *shared* ones also reach other atoms and must be left
        as they are.
        """
        resolver = AtomTokenResolver(self._shx)
        found = []
        for card in self.graph.cards_for_atom(atom):
            tokens = split_range_tokens(card.referenced_atoms)
            scope = getattr(card, 'residue_number', None)
            if isinstance(scope, int):
                scope = [scope]
            exclusive: list[str] = []
            shared: list[str] = []
            for token in tokens:
                if token in ('>', '<', '='):
                    continue
                names = resolver.resolve([token], scope).fullnames
                if old_fullname not in names:
                    continue
                if _names_a_single_atom(token) and len(names) == 1:
                    exclusive.append(token)
                else:
                    shared.append(token)
            if exclusive or shared:
                found.append((card, exclusive, shared))
        return found

    @staticmethod
    def _retarget(card, tokens: list[str], new_name: str) -> None:
        """Point *tokens* on *card* at *new_name*, keeping any suffix."""
        replacements = {}
        for token in tokens:
            suffix = token.split('_', 1)[1] if '_' in token else ''
            replacements[token] = f'{new_name}_{suffix}' if suffix else new_name
        if hasattr(card, 'atoms') and isinstance(card.atoms, list):
            card.atoms[:] = [replacements.get(t, t) for t in card.atoms]
            return
        for attribute in ('atom1', 'atom2', 'donor', 'acceptor'):
            current = getattr(card, attribute, None)
            if current in replacements:
                setattr(card, attribute, replacements[current])

    def _collect_orphaned_eqivs(self, report: DeletionReport) -> None:
        """Drop ``EQIV`` cards left with nothing referencing them.

        Ids are never reassigned: a surviving ``_$n`` elsewhere would
        otherwise start pointing at a different operation.
        """
        self.graph.rebuild()
        EqivCleaner(self._shx, self.graph).remove_orphans(report)

    def delete_atom(self, atom: Atom) -> DeletionReport:
        return self.delete_atoms([atom])

    @_undoable('Remove restraint')
    def remove_restraint(self, restraint: Restraint) -> DeletionReport:
        """Remove the card and nothing else.

        The atoms it referenced stay (**D-9b**).  This is also the path
        taken when a restraint is removed as a side effect, so no implicit
        route can ever reach an atom.
        """
        report = DeletionReport()
        report.add_card(restraint, RemovalReason.REQUESTED)
        restraint.delete()
        self._collect_orphaned_eqivs(report)
        self._notify()
        return report

    @_undoable('Delete restraint with atoms')
    def delete_restraint_with_atoms(self, restraint: Restraint) -> DeletionReport:
        """Remove the card **and** the atoms it names.

        Destructive, and deliberately separate from
        :meth:`remove_restraint` rather than hidden behind a flag.
        """
        named = self._referenced_atoms(restraint)
        report = DeletionReport()
        report.add_card(restraint, RemovalReason.REQUESTED)
        restraint.delete()
        for atom in named:
            report.add_atom(atom, RemovalReason.REQUESTED)
            atom.delete()
        self._notify()
        return report

    def _referenced_atoms(self, restraint: Restraint) -> list[Atom]:
        """Atoms a restraint names, as far as they can be resolved today.

        Range and wildcard tokens are not expanded yet; that arrives with
        the token resolver.  Unresolvable tokens are skipped rather than
        guessed at.
        """
        found: list[Atom] = []
        for token in restraint.atoms:
            if token in ('>', '<', '=') or token.startswith('$'):
                continue
            name = token if '_' in token else f'{token}_0'
            atom = self._shx.atoms.get_atom_by_name(name.upper())
            if atom is not None and atom not in found:
                found.append(atom)
        return found

    # ---------------------------------------------------------- addition

    #: Cards that name atoms and can be added directly. Kept apart from
    #: ``RESTRAINT_CARD_CLASSES``, which is documented as restraints only.
    ATOM_CARD_CLASSES: dict[str, str] = {
        'BIND': 'bind',
        'FREE': 'free',
        'CONN': 'conn',
        'HTAB': 'htab',
    }

    def add_bind(self, atom1, atom2) -> EditReport:
        """Add a bond between two atoms.

        Either operand may be a symmetry image; the ``EQIV`` it needs is
        created if this is the first reference to that operation.  Only
        one of the two may be an image: *"Only one of the two atoms may
        be an equivalent atom"*.
        """
        return self._add_atom_pair_card('BIND', atom1, atom2)

    def add_free(self, atom1, atom2) -> EditReport:
        """Remove a bond between two atoms, with the same rules as
        :meth:`add_bind`."""
        return self._add_atom_pair_card('FREE', atom1, atom2)

    def add_htab(self, donor, acceptor) -> EditReport:
        """Record a hydrogen bond.

        Only the acceptor may be a symmetry image: *"Only the acceptor
        atom may specify a symmetry operation (_$n) because CIF requires
        this"*.
        """
        return self._add_atom_pair_card('HTAB', donor, acceptor)

    @_undoable('Add instruction')
    def _add_atom_pair_card(self, keyword: str, first, second) -> EditReport:
        names = []
        for operand in (first, second):
            name = self.name_for(operand)
            if name is None:
                raise ValueError(
                    f'Cannot express {operand!r} as an atom reference; its '
                    f'symmetry operation could not be resolved.'
                )
            names.append(name)
        card = self._insert_card(keyword, names)
        problem = validate_symmetry_arity(card)
        if problem:
            self._undo_insert(card, keyword)
            raise ValueError(problem)
        self._notify()
        return EditReport(added=[card])

    def name_for(self, target) -> str | None:
        """How to refer to *target* in an instruction.

        Accepts an :class:`Atom`, a :class:`SymmetryMate`, a
        :class:`SymBond` neighbour, or a name that is already a string,
        so a viewer can pass whatever its click handler produced.
        """
        if isinstance(target, str):
            return target
        return EqivFactory(self._shx).symmetry_atom_name(target)

    def _insert_card(self, keyword: str, names: list[str]):
        from shelxfile.shelx import cards as card_module

        card_class = getattr(card_module, keyword)
        card = card_class(self._shx, [keyword, *names])
        position = self._card_insert_position(names)
        self._shx.insert_into_reslist(position, card)
        collection = getattr(self._shx, self.ATOM_CARD_CLASSES[keyword], None)
        if isinstance(collection, list):
            collection.append(card)
        else:
            setattr(self._shx, self.ATOM_CARD_CLASSES[keyword], card)
        self._shx.touch()
        return card

    def _undo_insert(self, card, keyword: str) -> None:
        try:
            self._shx.remove_from_reslist(card)
        except ValueError:
            pass
        collection = getattr(self._shx, self.ATOM_CARD_CLASSES[keyword], None)
        if isinstance(collection, list) and card in collection:
            collection.remove(card)

    def _card_insert_position(self, names: list[str]) -> int:
        """Below any ``EQIV`` the card references, above the atoms.

        An ``EQIV`` has to be defined before it is used, so a card that
        names one cannot sit above it.
        """
        from shelxfile.shelx.cards import EQIV as EqivCard

        referenced = {n.rsplit('_', 1)[1] for n in names if '_$' in n}
        lowest_allowed = 0
        for index, item in enumerate(self._shx._reslist):
            if isinstance(item, EqivCard) and item.id in referenced:
                lowest_allowed = max(lowest_allowed, index + 1)
        return max(lowest_allowed, self._shx._header_insert_position())

    @_undoable('Remove instruction')
    def remove_card(self, card) -> DeletionReport:
        """Remove an instruction, leaving its atoms alone (**D-9**).

        Accepts any card the document knows how to detach, restraints
        included, so a caller never has to work out which family it holds.
        An ``EQIV`` left with nothing referencing it is collected too.
        """
        from shelxfile.shelx.cards import Restraint as RestraintCard
        if isinstance(card, RestraintCard):
            return self.remove_restraint(card)
        report = DeletionReport()
        report.add_card(card, RemovalReason.REQUESTED)
        try:
            self._shx.remove_from_reslist(card)
        except ValueError:
            return DeletionReport()
        for name in self.ATOM_CARD_CLASSES.values():
            collection = getattr(self._shx, name, None)
            if isinstance(collection, list) and card in collection:
                collection.remove(card)
            elif collection is card:
                setattr(self._shx, name, None)
        self._collect_orphaned_eqivs(report)
        self._notify()
        return report

    @_undoable('Add restraint')
    def add_restraint(self, text: str, *, header: bool = False) -> EditReport:
        """Parse and insert a restraint instruction line.

        :param header: Put it just before the first atom, outside every
            ``PART``/``AFIX``/``RESI`` scope, instead of after the last
            existing restraint.  Inside a ``RESI`` block an unqualified name
            refers to that residue, so generated restraints that name atoms
            by their full name belong in the header.
        """
        restraint = self._shx.add_restraint(text, header=header)
        self._notify()
        return EditReport(added=[restraint])

    @_undoable('Add atom')
    def add_atom(self, **kwargs) -> EditReport:
        atom = self._shx.add_atom(**kwargs)
        self._notify()
        return EditReport(added=[atom])

    # ------------------------------------------------- structural edits

    def _to_fractional(self, xyz: Sequence[float], cartesian: bool) -> tuple[float, float, float]:
        from shelxfile.misc.misc import cart_to_frac

        if not cartesian:
            return float(xyz[0]), float(xyz[1]), float(xyz[2])
        return cart_to_frac([float(v) for v in xyz], list(self._shx.cell))

    @_undoable('Move atoms')
    def move_atoms(self, moves: Iterable[tuple[Atom, Sequence[float]]], *,
                   cartesian: bool = True) -> EditReport:
        """Put each atom of *moves* (``(atom, xyz)`` pairs) at a new position.

        :param cartesian: Whether *xyz* is Cartesian (Å, the default) or
            fractional.
        :returns: The moved atoms in :attr:`EditReport.changed`, and a message
            for every coordinate whose fixing code was released because it
            moved (see :func:`~shelxfile.edit.structure_edits.move_atom`).
        """
        report = EditReport()
        for atom, xyz in moves:
            report.messages.extend(
                structure_edits.move_atom(self._shx, atom, self._to_fractional(xyz, cartesian)))
            report.changed.append(atom)
        if report.changed:
            self._notify()
        return report

    @_undoable('Add free variable')
    def add_free_variable(self, value: float = 0.5) -> int:
        """Append a free variable starting at *value* and return its number."""
        number = structure_edits.add_free_variable(self._shx, value)
        self._notify()
        return number

    @_undoable('Assign PART')
    def assign_part(self, atoms: Iterable[Atom], part: int,
                    sof: float | Sequence[float]) -> EditReport:
        """Put *atoms* (all in ``PART 0``) into ``PART part`` with raw *sof*.

        See :func:`~shelxfile.edit.structure_edits.assign_part`.  The new
        ``PART`` cards are listed in :attr:`EditReport.added`.
        """
        atoms = list(atoms)
        cards = structure_edits.assign_part(self._shx, atoms, part, sof)
        if cards:
            self._notify()
        return EditReport(added=list(cards), changed=atoms)

    @_undoable('Duplicate atoms')
    def duplicate_atoms(
        self,
        atoms: Sequence[Atom],
        names: Sequence[str],
        part: int,
        sof: float | Sequence[float],
        coordinates: Sequence[Sequence[float] | None] | None = None,
        *,
        cartesian: bool = True,
        uvals: Sequence[Sequence[float] | None] | None = None,
    ) -> EditReport:
        """Copy *atoms* into a new ``PART part`` block.

        See :func:`~shelxfile.edit.structure_edits.duplicate_atoms`.  All
        per-atom arguments are parallel to *atoms*; the copies are returned
        in :attr:`EditReport.added` in the same order.
        """
        frac = None
        if coordinates is not None:
            frac = [None if xyz is None else self._to_fractional(xyz, cartesian)
                    for xyz in coordinates]
        copies = structure_edits.duplicate_atoms(
            self._shx, atoms, names, part, sof, frac_coordinates=frac, uvals=uvals)
        if copies:
            self._notify()
        return EditReport(added=list(copies))

    def split_name(self, atom: Atom, suffix: str, taken: Iterable[str] = (), *,
                   keep_own_name: bool = False) -> str:
        """``C1`` → ``C1A``/``C1B``, with fallbacks when that does not fit.

        See :func:`~shelxfile.edit.structure_edits.split_name`.
        """
        return structure_edits.split_name(self._shx, atom, suffix, taken,
                                          keep_own_name=keep_own_name)

    def free_atom_name(self, element: str, resinum: int = 0,
                       taken: Iterable[str] = ()) -> str:
        """The first unused ``<element><number>`` name in residue *resinum*."""
        return structure_edits.free_atom_name(self._shx, element, resinum, taken)

    def atom_problem_for_part(self, atoms: Iterable[Atom]) -> str | None:
        """Why *atoms* cannot be bracketed or duplicated as a unit, or ``None``."""
        return structure_edits.afix_closure_problem(self._shx, atoms)

    @_undoable('Set U values')
    def set_uvals(self, changes: Iterable[tuple[Atom, Sequence[float]]]) -> EditReport:
        """Replace the raw U values of atoms (``(atom, uvals)`` pairs).

        One value makes the atom isotropic, six set ``U11 U22 U33 U23 U13
        U12``.  Values are written as given, so riding codes such as ``-1.2``
        and free-variable codes are accepted.
        """
        report = EditReport()
        pending = []
        for atom, uvals in changes:
            values = [float(u) for u in uvals]
            if len(values) == 1:
                values += [0.0] * 5
            if len(values) != 6:
                raise ValueError(f'{atom.name}: expected 1 or 6 U values, got {len(values)}')
            pending.append((atom, values))
        for atom, values in pending:
            atom.set_uvals(values)
            atom.uvals_orig = list(values)
            report.changed.append(atom)
        if report.changed:
            self._shx.touch()
            self._notify()
        return report

    @_undoable('Add EQIV')
    def name_for_operation(self, atom: Atom, symmop: str) -> str:
        """How to write *atom* transformed by *symmop* in an instruction.

        ``'C2'`` when *symmop* is the identity, otherwise ``'C2_$n'`` with
        an ``EQIV`` for *exactly* this operation reused or created (e.g.
        ``'-X+1, Y, -Z+1/2'``).

        Unlike :meth:`EqivFactory.find`, lattice translations count here: a
        distance restraint to ``C2`` one cell further along is a different
        restraint.

        :raises ValueError: when *symmop* cannot be parsed or every ``$n`` is
            taken.
        """
        from shelxfile.edit.eqiv_factory import IDENTITY_SYMMOP, canonical_symmop

        wanted = canonical_symmop(symmop, reduce_translation=False)
        if wanted is None:
            raise ValueError(f'Cannot parse symmetry operation {symmop!r}')
        if wanted == IDENTITY_SYMMOP:
            return atom.fullname_short
        factory = EqivFactory(self._shx)
        existing = factory.find(symmop, reduce_translation=False)
        if existing is not None:
            return f'{atom.fullname_short}_{existing.id}'
        number = factory.next_free_number()
        if number is None:
            raise ValueError('No free EQIV number left')
        card = factory._insert(number, symmop)
        self._notify()
        return f'{atom.fullname_short}_{card.id}'

    # ------------------------------------------------------------ output

    def write(self, path: str | Path | None = None) -> None:
        """Write the file and remember this state as the saved one."""
        self._shx.write_shelx_file(path)
        self._history.mark_saved()

    def __repr__(self) -> str:
        return (f'<ShelxDocument {len(self._shx.atoms)} atoms, '
                f'{self.line_count} lines>')


def _names_a_single_atom(token: str) -> bool:
    """Whether *token* is a literal atom name rather than a class of them.

    Class references must never be rewritten by a rename, however few
    atoms they happen to match today.  ``BOND $H`` in a structure with
    one hydrogen still means "every hydrogen", and turning it into
    ``BOND H9`` would quietly change the instruction.

    Excluded: ``$element`` references, ``_*`` wildcards, and ``LAST``.
    """
    if not token or token.startswith('$'):
        return False
    if token.upper() == LAST_KEYWORD:
        return False
    return not token.endswith('_*')
