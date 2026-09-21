"""Working out what a ``SAME`` instruction actually links.

``SAME`` is the one restraint whose partners are not all written on the
card, which makes it the easiest to corrupt by editing.  The manual
describes two distinct modes:

**Plain ``SAME`` and ``SAME_<n>``**
    *"The list of atoms ... is compared with the same number of atoms
    which follow the SAME instruction."*  The counterpart fragment is the
    run of atoms after the card in the file, so deleting an atom that is
    never named on the card can still break it.  *"The position of a SAME
    instruction in the input file is critical."*

**``SAME_<class>``**
    *"SAME_XYZ no longer uses the following atoms but is applied to all
    residues with the name XYZ, so the atoms must have the same names in
    the same order in each of these residues."*

Both reduce to the same shape: a set of **positionally aligned
fragments** that SHELXL restrains to be alike.  Position *i* of every
fragment describes the same site, so if one is lost the rest of that
column is meaningless too.

Hydrogens are excluded throughout: *"Since hydrogen atoms are ignored by
SAME"*.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from shelxfile.edit.card_meta import AtomGrouping
from shelxfile.edit.token_resolver import AtomTokenResolver

if TYPE_CHECKING:
    from shelxfile import Shelxfile
    from shelxfile.atoms.atom import Atom
    from shelxfile.shelx.cards import SAME


@dataclass
class SameFragments:
    """The aligned fragments one ``SAME`` instruction ties together.

    :param fragments: one list of atoms per fragment.  Entry *i* of each
        fragment refers to the same site.
    :param mode: which of the two ``SAME`` readings produced these.
    :param truncated: ``True`` when the fragments could not be filled to
        equal length, e.g. because the file ends before enough following
        atoms were found.
    """

    fragments: list[list[Atom]]
    mode: AtomGrouping
    truncated: bool = False

    @property
    def width(self) -> int:
        """Number of aligned positions."""
        return min((len(f) for f in self.fragments), default=0)

    def column(self, index: int) -> list[Atom]:
        """Every atom describing site *index*."""
        return [f[index] for f in self.fragments if index < len(f)]

    def position_of(self, atom: Atom) -> int | None:
        """Which aligned site *atom* occupies, if any."""
        for fragment in self.fragments:
            for index, candidate in enumerate(fragment):
                if candidate is atom:
                    return index
        return None

    def counterparts(self, atom: Atom) -> list[Atom]:
        """The other atoms restrained to match *atom*.

        These are what an edit has to account for: dropping *atom* without
        them leaves the remaining fragments misaligned, and SHELXL would
        silently pair up the wrong sites.
        """
        index = self.position_of(atom)
        if index is None:
            return []
        return [a for a in self.column(index) if a is not atom]

    @property
    def all_atoms(self) -> list[Atom]:
        seen: set[int] = set()
        result: list[Atom] = []
        for fragment in self.fragments:
            for atom in fragment:
                if id(atom) not in seen:
                    seen.add(id(atom))
                    result.append(atom)
        return result


class SameResolver:
    """Expands ``SAME`` instructions into aligned fragments."""

    def __init__(self, shx: Shelxfile) -> None:
        self.shx = shx
        self._tokens = AtomTokenResolver(shx)

    def fragments(self, card: SAME) -> SameFragments:
        if card.atom_grouping is AtomGrouping.SAME_RESIDUE_CLASS:
            return self._residue_class_fragments(card)
        return self._following_atom_fragments(card)

    # ------------------------------------------------ plain SAME / SAME_n

    def _following_atom_fragments(self, card: SAME) -> SameFragments:
        template = self._named_atoms(card)
        following = self._atoms_after(card, len(template))
        return SameFragments(
            fragments=[template, following],
            mode=AtomGrouping.SAME_FOLLOWING_ATOMS,
            truncated=len(following) < len(template),
        )

    def _named_atoms(self, card: SAME) -> list[Atom]:
        """Atoms written on the card, ranges expanded, hydrogens dropped."""
        scope = card.residue_number if card.residue_number else [0]
        resolved = self._tokens.resolve(list(card.atoms), scope)
        atoms = []
        for name in resolved.fullnames:
            atom = self.shx.atoms.get_atom_by_name(name)
            if atom is not None and not atom.is_hydrogen and not atom.qpeak:
                atoms.append(atom)
        return atoms

    def _atoms_after(self, card: SAME, count: int) -> list[Atom]:
        """The next *count* non-hydrogen atoms following the card.

        Walks the file rather than the atom list, because ``SAME`` is
        defined by its position among the instructions.
        """
        from shelxfile.atoms.atom import Atom as AtomClass

        try:
            start = self.shx.index_of(card)
        except ValueError:
            return []
        found: list[AtomClass] = []
        for item in self.shx._reslist[start + 1:]:
            if len(found) >= count:
                break
            if not isinstance(item, AtomClass):
                continue
            if item.is_hydrogen or item.qpeak:
                continue
            found.append(item)
        return found

    # ------------------------------------------------------- SAME_<class>

    def _residue_class_fragments(self, card: SAME) -> SameFragments:
        """One fragment per residue of the card's class.

        The same atom names are looked up in each residue, so a name
        missing from one residue shortens that fragment and the alignment
        is reported as truncated.
        """
        numbers = self.shx.residues.residue_classes.get(card.residue_class, [])
        fragments: list[list[Atom]] = []
        for number in numbers:
            resolved = self._tokens.resolve(list(card.atoms), [number])
            atoms = []
            for name in resolved.fullnames:
                atom = self.shx.atoms.get_atom_by_name(name)
                if atom is not None and not atom.is_hydrogen and not atom.qpeak:
                    atoms.append(atom)
            if atoms:
                fragments.append(atoms)
        widths = {len(f) for f in fragments}
        return SameFragments(
            fragments=fragments,
            mode=AtomGrouping.SAME_RESIDUE_CLASS,
            truncated=len(widths) > 1,
        )
