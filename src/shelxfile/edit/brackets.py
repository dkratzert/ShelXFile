"""``RESI``, ``AFIX`` and ``PART`` as brackets around runs of atoms.

All three scope the atoms that follow them, but they end differently, and
the original edit plan's single shared "is this bracket empty?" helper
would have papered over that:

``AFIX``
    *"applies constraints ... for all atoms until the next AFIX
    instruction"*, except that a rotating group (``n`` = 7) runs *"until
    the next AFIX with n not equal to 7"*.
``PART``
    *"remains in force until a further PART instruction is read"*.  Note
    ``PART -n`` is a real group; only ``PART 0`` resets.
``RESI``
    *"Until the next RESI instruction"*, or the end of the atom list.

``AFIX`` carries a second constraint the plan missed entirely: most codes
require a **specific number** of atoms -- *"Each AFIX instruction must be
followed by the required number of hydrogen or other atoms"*.  Deleting
one carbon from an ``AFIX 66`` hexagon leaves five atoms for a six-atom
fit, which SHELXL cannot satisfy.  That is not a weakened restraint but
an invalid instruction, which is why such a group is removed whole.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from shelxfile.edit.card_meta import AfixDependency

if TYPE_CHECKING:
    from shelxfile import Shelxfile
    from shelxfile.atoms.atom import Atom

#: Atoms each ``AFIX`` *geometry* code expects to follow it.
#:
#: ``m`` selects the geometry, ``n`` the constraint, so this is keyed on
#: ``mn // 10``.  Codes above 16 fit a ``FRAG``/``FEND`` block whose size
#: is only known from that block, and ``m`` = 0 means "no action"; both
#: map to ``None``.
AFIX_EXPECTED_ATOMS: dict[int, int] = {
    1: 1,    # tertiary C-H
    2: 2,    # secondary CH2
    3: 3,    # CH3
    4: 1,    # aromatic C-H or amide N-H
    5: 5,    # regular pentagon
    6: 6,    # regular hexagon
    7: 6,    # currently identical to m=6
    8: 1,    # idealised OH
    9: 2,    # terminal X=CH2 / X=NH2+
    10: 10,  # pentamethylcyclopentadienyl
    11: 8,   # naphthalene, numbered as a figure of eight
    12: 6,   # disordered methyl, two staggered positions
    13: 3,   # CH3 with torsion from the difference map
    14: 1,   # idealised OH
    15: 1,   # BH in a polyhedral fragment
    16: 1,   # acetylenic C-H
}


def expected_atom_count(mn: int | None) -> int | None:
    """How many atoms an ``AFIX mn`` group must be followed by.

    ``None`` when the code sets no requirement: ``m`` = 0 (no action) or
    ``m`` > 16, where the count comes from the matching ``FRAG`` block.
    """
    if mn is None:
        return None
    return AFIX_EXPECTED_ATOMS.get(abs(mn) // 10)


@dataclass
class BracketScope:
    """One bracket card and the atoms it governs."""

    card: object
    members: list[Atom] = field(default_factory=list)
    #: The card that ends this scope, if any.
    closer: object | None = None

    @property
    def is_empty(self) -> bool:
        return not self.members

    @property
    def is_reset(self) -> bool:
        """``AFIX 0`` / ``PART 0`` / ``RESI 0`` -- never removable.

        These are not groups but the instruction that ends one, so
        deleting them would silently extend the preceding group over
        everything that follows.
        """
        card = self.card
        if hasattr(card, 'mn'):
            return not card.mn
        if hasattr(card, 'n'):
            return card.n == 0
        if hasattr(card, 'residue_number'):
            return card.residue_number == 0
        return False

    @property
    def expected_count(self) -> int | None:
        """Required member count, for ``AFIX`` groups that have one.

        ``None`` for an ``n`` = 5 card.  Those are *dependent* atoms
        continuing the enclosing rigid group -- "automatically generated
        for the atoms following an n=6 or n=9 atom, so does not need to be
        included specifically unless m has to be changed" -- so the
        geometry's member count belongs to the group that opened it, not
        to this continuation.
        """
        if self.dependency is AfixDependency.RIGID_DEPENDENT:
            return None
        return expected_atom_count(getattr(self.card, 'mn', None))

    @property
    def is_under_populated(self) -> bool:
        """Fewer atoms than the ``AFIX`` geometry needs."""
        required = self.expected_count
        if required is None or self.is_reset:
            return False
        return len(self.members) < required

    @property
    def dependency(self) -> AfixDependency:
        card = self.card
        if hasattr(card, 'dependency'):
            return card.dependency
        return AfixDependency.NONE


class BracketResolver:
    """Finds bracket scopes and their members in a parsed file."""

    def __init__(self, shx: Shelxfile) -> None:
        self.shx = shx

    # ------------------------------------------------------------- scopes

    def afix_scopes(self) -> list[BracketScope]:
        return self._scopes('AFIX')

    def part_scopes(self) -> list[BracketScope]:
        return self._scopes('PART')

    def resi_scopes(self) -> list[BracketScope]:
        return self._scopes('RESI')

    def scope_of(self, card: object) -> BracketScope | None:
        """The scope opened by *card*, if it is a bracket card."""
        kind = type(card).__name__
        if kind not in ('AFIX', 'PART', 'RESI'):
            return None
        for scope in self._scopes(kind):
            if scope.card is card:
                return scope
        return None

    def _scopes(self, kind: str) -> list[BracketScope]:
        from shelxfile.atoms.atom import Atom as AtomClass

        scopes: list[BracketScope] = []
        current: BracketScope | None = None
        for item in self.shx._reslist:
            if type(item).__name__ == kind:
                if current is not None and self._continues(kind, current.card, item):
                    # A rotating AFIX group runs on through further n=7
                    # cards rather than being closed by them.
                    continue
                if current is not None:
                    current.closer = item
                current = BracketScope(card=item)
                scopes.append(current)
                continue
            if current is None:
                continue
            if isinstance(item, AtomClass) and not item.qpeak:
                current.members.append(item)
        return scopes

    @staticmethod
    def _continues(kind: str, open_card: object, next_card: object) -> bool:
        """Whether *next_card* extends rather than closes the open scope.

        Only ``AFIX`` does this: *"The following atoms (until the next
        AFIX with n not equal to 7) are allowed to ride ... and rotate"*.
        """
        if kind != 'AFIX':
            return False
        open_mn = getattr(open_card, 'mn', None)
        next_mn = getattr(next_card, 'mn', None)
        if open_mn is None or next_mn is None:
            return False
        return abs(open_mn) % 10 == 7 and abs(next_mn) % 10 == 7

    # ------------------------------------------------------------- pivots

    def pivot_of(self, scope: BracketScope) -> Atom | None:
        """The atom a riding or rotating group hangs from.

        Riding groups *"ride on the previous non-riding atom"*, which sits
        **before** the ``AFIX`` card rather than inside its bracket, so
        deleting it orphans an otherwise complete group.  Rigid groups
        keep their pivot inside the bracket and are covered by the member
        count instead.
        """
        from shelxfile.atoms.atom import Atom as AtomClass

        if not scope.dependency.pivot_is_outside_bracket:
            return None
        try:
            start = self.shx.index_of(scope.card)
        except ValueError:
            return None
        for item in reversed(self.shx._reslist[:start]):
            if not isinstance(item, AtomClass) or item.qpeak:
                continue
            # "The atom on which riding is performed may not itself be a
            # riding atom", so skip over other riders.
            rider = AfixDependency.from_code(getattr(item.afix, 'mn', None))
            if rider.pivot_is_outside_bracket:
                continue
            return item
        return None

    def scopes_pivoted_on(self, atom: Atom) -> list[BracketScope]:
        """Riding/rotating groups that would be orphaned by losing *atom*."""
        return [scope for scope in self.afix_scopes()
                if self.pivot_of(scope) is atom]
