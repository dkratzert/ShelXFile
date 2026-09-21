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

from pathlib import Path
from typing import TYPE_CHECKING, Callable, Iterable, Union

from shelxfile.edit.line_map import RenderedFile, render
from shelxfile.edit.reports import DeletionReport, EditReport, RemovalReason

if TYPE_CHECKING:
    from shelxfile import Shelxfile
    from shelxfile.atoms.atom import Atom
    from shelxfile.shelx.cards import Command, Restraint

#: An observer is called with the document after every successful edit.
Observer = Callable[['ShelxDocument'], None]

ResListItem = Union['Atom', 'Restraint', 'Command', str]


class ShelxDocument:
    """A :class:`Shelxfile` plus the operations an editor needs.

    :param shx: the model to wrap.  The document does not copy it; both
        refer to the same :class:`Shelxfile`.
    """

    def __init__(self, shx: Shelxfile) -> None:
        self._shx = shx
        self._observers: list[Observer] = []
        self._rendered: RenderedFile | None = None

    # ------------------------------------------------------------- model

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
        self.invalidate()
        for callback in list(self._observers):
            callback(self)

    # ---------------------------------------------------------- deletion

    def delete_atoms(self, atoms: Iterable[Atom]) -> DeletionReport:
        """Delete *atoms*, reporting everything that went with them.

        Symmetry-generated atoms are skipped: they are products of
        :meth:`Shelxfile.grow`/:meth:`~Shelxfile.pack`, not lines in the
        file.
        """
        report = DeletionReport()
        for atom in list(atoms):
            if atom.symmgen:
                continue
            report.add_atom(atom, RemovalReason.REQUESTED)
            atom.delete()
        if not report.is_empty:
            self._notify()
        return report

    def delete_atom(self, atom: Atom) -> DeletionReport:
        return self.delete_atoms([atom])

    def remove_restraint(self, restraint: Restraint) -> DeletionReport:
        """Remove the card and nothing else.

        The atoms it referenced stay (**D-9b**).  This is also the path
        taken when a restraint is removed as a side effect, so no implicit
        route can ever reach an atom.
        """
        report = DeletionReport()
        report.add_card(restraint, RemovalReason.REQUESTED)
        restraint.delete()
        self._notify()
        return report

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

    def add_restraint(self, text: str) -> EditReport:
        """Parse and insert a restraint instruction line."""
        restraint = self._shx.add_restraint(text)
        self._notify()
        return EditReport(added=[restraint])

    def add_atom(self, **kwargs) -> EditReport:
        atom = self._shx.add_atom(**kwargs)
        self._notify()
        return EditReport(added=[atom])

    # ------------------------------------------------------------ output

    def write(self, path: str | Path | None = None) -> None:
        self._shx.write_shelx_file(path)

    def __repr__(self) -> str:
        return (f'<ShelxDocument {len(self._shx.atoms)} atoms, '
                f'{self.line_count} lines>')
