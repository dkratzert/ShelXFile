"""Finding or creating the ``EQIV`` a symmetry reference needs.

To write ``BIND C1 C2_$3`` ShelXFile first has to know that ``$3`` is the
operation the user actually clicked on.  This module turns a symmetry
operation into an ``EQIV`` card, reusing an existing definition whenever
one matches so a file does not accumulate duplicates of the same
operation under different numbers.

Matching is done on the parsed rotation and translation rather than the
text, because ``1-x, y, 1-z`` and ``-x+1, +y, -z+1`` are the same
operation written two ways.

Two manual rules constrain what may be done:

* *"Such a symmetry operation must be defined before it is used"* -- a
  new ``EQIV`` goes above the instructions that will reference it.
* *"The same $n may not appear on two separate EQIV instructions"* --
  numbers are handed out from the unused ones and never reassigned.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from shelxfile.misc.dsrmath import SymmetryElement

if TYPE_CHECKING:
    from shelxfile import Shelxfile
    from shelxfile.atoms.atom import Atom
    from shelxfile.atoms.symmetry_mate import SymmetryMate
    from shelxfile.shelx.cards import EQIV


def canonical_symmop(symmop: str) -> tuple | None:
    """A comparable form of *symmop*, or ``None`` if it will not parse.

    Two spellings of one operation give the same value, so this is what
    lookups compare on rather than the original text.  Translations are
    reduced into ``[0, 1)`` because a lattice translation does not make a
    different operation.
    """
    parts = [p.strip() for p in symmop.split(',') if p.strip()]
    if len(parts) != 3:
        return None
    try:
        element = SymmetryElement(parts)
    except Exception:  # noqa: BLE001 - malformed input is the caller's problem
        return None
    matrix = tuple(tuple(round(float(v), 6) for v in row) for row in element.matrix)
    trans = tuple(round(float(v) % 1.0, 6) for v in element.trans)
    return matrix, trans


class EqivFactory:
    """Reuses or mints the ``EQIV`` cards symmetry references need."""

    def __init__(self, shx: Shelxfile) -> None:
        self.shx = shx

    # ----------------------------------------------------------- lookup

    def find(self, symmop: str) -> EQIV | None:
        """An existing live ``EQIV`` for *symmop*, if there is one.

        Post-``END`` definitions are ignored: that block is a record of a
        finished refinement and must not be referenced by new input.
        """
        wanted = canonical_symmop(symmop)
        if wanted is None:
            return None
        for card in self.shx.eqiv:
            if not card.lifetime or card.number is None:
                continue
            if canonical_symmop(card.symmop) == wanted:
                return card
        return None

    def next_free_number(self) -> int | None:
        """The lowest ``$n`` not yet in use, or ``None`` if exhausted."""
        from shelxfile.shelx.cards import EQIV as EqivCard

        taken = {c.number for c in self.shx.eqiv if c.number is not None}
        for candidate in range(EqivCard.MIN_NUMBER, EqivCard.MAX_NUMBER + 1):
            if candidate not in taken:
                return candidate
        return None

    # ------------------------------------------------------------ mint

    def eqiv_for(self, symmop: str, *, create: bool = True) -> EQIV | None:
        """The ``EQIV`` describing *symmop*, creating one if needed.

        :returns: the matching card, or ``None`` when *symmop* cannot be
            parsed, ``create`` is false and none exists, or every number
            is already taken.
        """
        existing = self.find(symmop)
        if existing is not None or not create:
            return existing
        if canonical_symmop(symmop) is None:
            return None
        number = self.next_free_number()
        if number is None:
            return None
        return self._insert(number, symmop)

    def _insert(self, number: int, symmop: str) -> EQIV:
        from shelxfile.shelx.cards import EQIV as EqivCard

        card = EqivCard(self.shx, ['EQIV', f'${number}', *symmop.split()])
        position = self.insert_position()
        self.shx._reslist.insert(position, card)
        self.shx.eqiv.append(card)
        self.shx.touch()
        return card

    def insert_position(self) -> int:
        """Where a new ``EQIV`` belongs in ``_reslist``.

        After the last existing ``EQIV`` so the block stays together, and
        otherwise before the first atom -- an ``EQIV`` has to be defined
        ahead of anything that uses it, and references live among the
        instructions that follow.
        """
        from shelxfile.atoms.atom import Atom as AtomClass
        from shelxfile.shelx.cards import EQIV as EqivCard

        last_eqiv = None
        first_atom = None
        for index, item in enumerate(self.shx._reslist):
            if isinstance(item, EqivCard) and item.lifetime:
                last_eqiv = index
            elif first_atom is None and isinstance(item, AtomClass):
                first_atom = index
        if last_eqiv is not None:
            return last_eqiv + 1
        if first_atom is not None:
            return first_atom
        return max(len(self.shx._reslist) - 1, 0)

    # ------------------------------------------------------- references

    def symmetry_atom_name(self, target: Atom | SymmetryMate,
                           *, create: bool = True) -> str | None:
        """How to write *target* in an instruction.

        ``'C2'`` for an atom of the asymmetric unit, ``'C2_$3'`` for a
        symmetry image -- creating the ``EQIV`` if the operation has not
        been named yet.

        :returns: ``None`` when the image's operation cannot be resolved.
        """
        mate = getattr(target, 'symm_mate', None)
        if mate is None and hasattr(target, 'parent'):
            mate = target
        if mate is None:
            return target.fullname_short
        if mate.is_identity:
            return mate.parent.fullname_short
        label = mate.symm_label
        if not label:
            return None
        card = self.eqiv_for(label, create=create)
        if card is None:
            return None
        return f'{mate.parent.fullname_short}_{card.id}'
