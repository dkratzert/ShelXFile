"""Working out everything a deletion drags with it, then doing it.

Deleting an atom is rarely a local change.  It can empty a restraint,
break a rigid group, orphan a symmetry definition -- and removing those
can in turn strand more atoms.  The work is therefore split in two:

**plan**
    walk the consequences to a fixpoint, touching nothing.
**apply**
    carry the whole plan out at once.

Computing first means a cascade can never be half-applied: if planning
goes wrong the model is still exactly as it was.

Direction matters (**D-9**).  Deleting an *atom* may remove cards and,
for ``AFIX`` groups, further atoms.  Removing a *card* never deletes
atoms -- losing a restraint is a refinement decision, not a statement
about the atoms it mentioned.

One manual rule does most of the work in keeping edits small: *"If ...
some or all of the named atoms cannot be found for a particular residue,
the instruction is simply ignored for that residue."*  So a token that
still resolves somewhere is left alone, and a card scoped over several
residues survives losing one of them untouched.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Iterable

from shelxfile.edit.brackets import BracketResolver
from shelxfile.edit.card_meta import AtomGrouping, AtomListSemantics
from shelxfile.edit.reports import DeletionReport, RemovalReason
from shelxfile.edit.same_links import SameResolver
from shelxfile.edit.token_resolver import AtomTokenResolver, split_range_tokens

if TYPE_CHECKING:
    from shelxfile import Shelxfile
    from shelxfile.atoms.atom import Atom
    from shelxfile.edit.card_meta import AtomReferencingCard
    from shelxfile.edit.graph import AtomRestraintGraph

RANGE_MARKERS = ('>', '<')


@dataclass
class CascadePlan:
    """Everything a deletion will do, before any of it happens."""

    atoms: dict[int, tuple] = field(default_factory=dict)
    cards: dict[int, tuple] = field(default_factory=dict)
    #: card id -> (card, surviving token list)
    edits: dict[int, tuple] = field(default_factory=dict)

    def wants_atom(self, atom: Atom) -> bool:
        return id(atom) in self.atoms

    def wants_card(self, card) -> bool:
        return id(card) in self.cards

    def add_atom(self, atom: Atom, reason: RemovalReason) -> bool:
        """Record an atom for deletion. ``False`` if already planned."""
        if id(atom) in self.atoms:
            return False
        self.atoms[id(atom)] = (atom, reason)
        return True

    def add_card(self, card, reason: RemovalReason) -> bool:
        if id(card) in self.cards:
            return False
        self.cards[id(card)] = (card, reason)
        self.edits.pop(id(card), None)
        return True

    def edit_card(self, card, tokens: list[str]) -> None:
        if id(card) in self.cards:
            return
        self.edits[id(card)] = (card, tokens)

    @property
    def is_empty(self) -> bool:
        return not (self.atoms or self.cards or self.edits)


class CascadeEngine:
    """Plans and applies the fallout of deleting atoms."""

    def __init__(self, shx: Shelxfile, graph: AtomRestraintGraph) -> None:
        self.shx = shx
        self.graph = graph
        self._tokens = AtomTokenResolver(shx)
        self._same = SameResolver(shx)

    # --------------------------------------------------------------- plan

    def plan(self, atoms: Iterable[Atom]) -> CascadePlan:
        """Walk the consequences of deleting *atoms* to a fixpoint.

        Every requested atom is registered before any consequence is
        examined, so a batch is judged as a batch: deleting both atoms of
        a two-atom restraint is recognised as emptying it, rather than
        looking like two separate near-misses.
        """
        plan = CascadePlan()
        requested = [a for a in atoms if not a.symmgen]
        for atom in requested:
            plan.add_atom(atom, RemovalReason.REQUESTED)

        queue = list(requested)
        processed: set[int] = set()
        while queue:
            atom = queue.pop(0)
            if id(atom) in processed:
                continue
            processed.add(id(atom))
            self._plan_card_fallout(atom, plan)
            for extra, reason in self._plan_afix_fallout(atom, plan):
                plan.add_atom(extra, reason)
                if id(extra) not in processed:
                    queue.append(extra)
        return plan

    def _plan_card_fallout(self, atom: Atom, plan: CascadePlan) -> None:
        """Cards that reference *atom*: trim their tokens, or drop them."""
        for card in self.graph.cards_for_atom(atom):
            if plan.wants_card(card):
                continue
            survivors = self._surviving_tokens(card, atom, plan)
            if survivors is None:
                continue  # nothing about this card actually changes
            reason = self._removal_reason(card, survivors)
            if reason is not None:
                plan.add_card(card, reason)
            else:
                plan.edit_card(card, survivors)

    def _surviving_tokens(
        self,
        card: AtomReferencingCard,
        atom: Atom,
        plan: CascadePlan,
    ) -> list[str] | None:
        """Tokens left on *card* once everything planned is gone.

        ``None`` means the card is unaffected, which is the common case
        for a residue-scoped instruction losing one of its residues.
        """
        tokens = split_range_tokens(card.referenced_atoms)
        doomed = {id(a) for a in [atom]} | set(plan.atoms)
        dead_positions = {
            index for index, token in enumerate(tokens)
            if token not in RANGE_MARKERS
            and self._token_is_dead(card, token, doomed)
        }
        if not dead_positions:
            return None
        if card.atom_grouping is AtomGrouping.PAIRS:
            dead_positions |= self._partner_positions(tokens, dead_positions)
        survivors = [t for i, t in enumerate(tokens) if i not in dead_positions]
        return _strip_dangling_markers(survivors)

    def _token_is_dead(
        self,
        card: AtomReferencingCard,
        token: str,
        doomed: set[int],
    ) -> bool:
        """Whether *token* will resolve to nothing once *doomed* are gone.

        A token that still names a surviving atom stays: SHELXL ignores an
        instruction for the residues where the atoms are missing rather
        than rejecting it outright.
        """
        scope = getattr(card, 'residue_number', None)
        if isinstance(scope, int):
            scope = [scope]
        resolved = self._tokens.resolve([token], scope)
        names = resolved.fullnames
        if not names:
            return False  # already unresolvable; not ours to judge
        for name in names:
            found = self.shx.atoms.get_atom_by_name(name)
            if found is not None and id(found) not in doomed:
                return False
        return True

    @staticmethod
    def _partner_positions(tokens: list[str], dead: set[int]) -> set[int]:
        """The other half of each pair, for ``DFIX``/``DANG``/``SADI``.

        "The distances between the first and second named atoms, the third
        and fourth ...", so a surviving partner has nothing to pair with.
        """
        atom_positions = [i for i, t in enumerate(tokens)
                          if t not in RANGE_MARKERS]
        partners: set[int] = set()
        for slot, position in enumerate(atom_positions):
            if position not in dead:
                continue
            twin = slot + 1 if slot % 2 == 0 else slot - 1
            if 0 <= twin < len(atom_positions):
                partners.add(atom_positions[twin])
        return partners

    def _removal_reason(
        self,
        card: AtomReferencingCard,
        survivors: list[str],
    ) -> RemovalReason | None:
        """Why *card* cannot survive on *survivors*, if it cannot."""
        remaining = [t for t in survivors if t not in RANGE_MARKERS]
        if not remaining:
            if card.EMPTY_MEANS is not AtomListSemantics.EXPLICIT:
                # Emptying it would silently turn a targeted instruction
                # into one that applies to the whole structure.
                return RemovalReason.WOULD_BECOME_GLOBAL
            return RemovalReason.BELOW_MIN_ATOMS
        if len(remaining) < card.MIN_ATOMS:
            return RemovalReason.BELOW_MIN_ATOMS
        return None

    def _plan_afix_fallout(
        self,
        atom: Atom,
        plan: CascadePlan,
    ) -> list[tuple[Atom, RemovalReason]]:
        """``AFIX`` groups invalidated by losing *atom*.

        A group below its required member count, or one whose outside
        pivot has gone, is an instruction SHELXL cannot carry out -- so it
        goes whole, taking its members with it.
        """
        resolver = BracketResolver(self.shx)
        stranded: list[tuple[Atom, RemovalReason]] = []
        for scope in resolver.afix_scopes():
            if scope.is_reset or plan.wants_card(scope.card):
                continue
            members = [m for m in scope.members if not plan.wants_atom(m)]
            if any(m is atom for m in scope.members):
                if scope.expected_count is None:
                    continue
                # *atom* is already planned, so *members* is the count
                # that would remain.
                if len(members) >= scope.expected_count:
                    continue
                reason = RemovalReason.AFIX_UNDER_POPULATED
            elif resolver.pivot_of(scope) is atom:
                reason = RemovalReason.AFIX_PIVOT_DELETED
            else:
                continue
            plan.add_card(scope.card, reason)
            if scope.closer is not None:
                plan.add_card(scope.closer, reason)
            for member in scope.members:
                if member is not atom:
                    stranded.append((member, reason))
        return stranded

    # -------------------------------------------------------------- apply

    def apply(self, plan: CascadePlan) -> DeletionReport:
        """Carry out *plan* in one go."""
        report = DeletionReport()
        for card, tokens in plan.edits.values():
            before = str(card)
            self._rewrite(card, tokens)
            report.add_edit(card, before)
        for card, reason in plan.cards.values():
            if self._remove_card(card):
                report.add_card(card, reason)
        for atom, reason in plan.atoms.values():
            report.add_atom(atom, reason)
            atom.delete()
        return report

    @staticmethod
    def _rewrite(card: AtomReferencingCard, tokens: list[str]) -> None:
        """Replace a card's atom list in place, keeping its identity.

        Cards whose grammar couples a count to the list -- only ``MPLA``
        so far -- get that count brought back into range, since SHELXL
        cannot satisfy an ``na`` larger than the list it indexes.
        """
        if hasattr(card, 'atoms') and isinstance(card.atoms, list):
            card.atoms[:] = tokens
            clamp = getattr(card, 'clamp_na', None)
            if callable(clamp):
                clamp()
            return
        # FREE/HTAB keep their operands in named fields.
        names = [t for t in tokens if t not in RANGE_MARKERS]
        for attribute, value in zip(
            ('atom1', 'atom2') if hasattr(card, 'atom1') else ('donor', 'acceptor'),
            names + [None, None],
        ):
            setattr(card, attribute, value)

    def _remove_card(self, card) -> bool:
        from shelxfile.shelx.cards import Restraint

        try:
            if isinstance(card, Restraint):
                self.shx.restraints.remove(card)
            else:
                self.shx.remove_from_reslist(card)
        except ValueError:
            return False
        return True


def _strip_dangling_markers(tokens: list[str]) -> list[str]:
    """Drop range markers left without an atom on both sides."""
    cleaned: list[str] = []
    for index, token in enumerate(tokens):
        if token not in RANGE_MARKERS:
            cleaned.append(token)
            continue
        has_left = bool(cleaned) and cleaned[-1] not in RANGE_MARKERS
        has_right = any(t not in RANGE_MARKERS for t in tokens[index + 1:])
        if has_left and has_right:
            cleaned.append(token)
    return cleaned
