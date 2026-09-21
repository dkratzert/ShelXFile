"""Where a symmetry-generated atom came from.

:meth:`Shelxfile.grow` and :meth:`~Shelxfile.pack` invent a fresh name
for every image they produce, which is enough to draw a picture but not
to edit one.  A viewer that lets the user click a symmetry mate and ask
for a bond needs to say *which* operation produced it, so ShelXFile can
write ``BIND C1 C2_$3`` against the right ``EQIV``.

:class:`SymmetryMate` records that provenance.  The label it produces
matches :class:`~shelxfile.atoms.pairs.SymBond`'s, so an atom click and a
bond click describe the same operation the same way.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from shelxfile.atoms.atom import Atom


@dataclass(frozen=True)
class SymmetryMate:
    """Identifies one image of an asymmetric-unit atom.

    :param parent: the atom in the asymmetric unit this was generated
        from.
    :param symm_number: index of the operation in
        ``Shelxfile.symmcards``; 0 is the identity.
    :param hkl_shift: extra integer lattice translation applied on top of
        the symmetry operation.
    """

    parent: Atom
    symm_number: int = 0
    hkl_shift: tuple[int, int, int] = (0, 0, 0)

    @property
    def is_identity(self) -> bool:
        """``True`` when this is the parent atom itself, unmoved."""
        return self.symm_number == 0 and self.hkl_shift == (0, 0, 0)

    @property
    def symm_label(self) -> str:
        """The operation in SHELXL/CIF form, e.g. ``'-x+1, y+1/2, -z'``.

        Empty for the identity, matching
        :attr:`~shelxfile.atoms.pairs.SymBond.symm_label`.
        """
        if self.is_identity:
            return ''
        from shelxfile.atoms.atoms import _build_symm_label

        shx = self.parent.shx
        try:
            element = shx.symmcards[self.symm_number]
        except (IndexError, TypeError):
            return ''
        return _build_symm_label(element, *self.hkl_shift)

    def __str__(self) -> str:
        label = self.symm_label
        if not label:
            return self.parent.fullname_short
        return f'{self.parent.fullname_short} [{label}]'
