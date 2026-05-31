"""Tests for build_conntable (shelxfile.misc.misc) and Atoms.conntable property."""
from unittest import TestCase

import numpy as np

from shelxfile.atoms.pairs import Bond, SymBond
from shelxfile.misc.misc import build_conntable
from shelxfile.shelx.shelx import Shelxfile

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_CO_DIST = 1.4   # typical C-O single bond (Å)


def _two_coords(dist: float = _CO_DIST) -> np.ndarray:
    return np.array([[0.0, 0.0, 0.0], [dist, 0.0, 0.0]])


# ---------------------------------------------------------------------------
# Unit tests for build_conntable
# ---------------------------------------------------------------------------

class TestBuildConntableBasic(TestCase):

    def test_two_bonded_atoms(self):
        """C and O at 1.4 Å should be bonded."""
        bonds = build_conntable(_two_coords(), ['C', 'O'], [0, 0])
        self.assertIn((0, 1), bonds)

    def test_two_far_atoms_not_bonded(self):
        """C and O at 5 Å should NOT be bonded (exceeds 4 Å pre-filter)."""
        bonds = build_conntable(np.array([[0, 0, 0], [5, 0, 0]], dtype=float),
                                ['C', 'O'], [0, 0])
        self.assertNotIn((0, 1), bonds)

    def test_same_atom_not_self_bonded(self):
        """A single atom should produce no bonds."""
        bonds = build_conntable(np.array([[0, 0, 0]], dtype=float), ['C'], [0])
        self.assertEqual(bonds, ())

    def test_empty_input(self):
        """Empty atom list returns empty tuple."""
        bonds = build_conntable(np.zeros((0, 3)), [], [])
        self.assertEqual(bonds, ())

    def test_result_i_less_than_j(self):
        """All returned pairs should have i < j."""
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        for i, j in shx.atoms.conntable:
            self.assertLess(i, j)

    def test_hh_not_bonded(self):
        """H–H pair at bonding distance should NOT produce a bond."""
        coords = np.array([[0, 0, 0], [0.74, 0, 0]], dtype=float)  # H2 bond length
        bonds = build_conntable(coords, ['H', 'H'], [0, 0])
        self.assertNotIn((0, 1), bonds)

    def test_ch_bonded(self):
        """C–H pair at 1.09 Å should be bonded."""
        coords = np.array([[0, 0, 0], [1.09, 0, 0]], dtype=float)
        bonds = build_conntable(coords, ['C', 'H'], [0, 0])
        self.assertIn((0, 1), bonds)


class TestBuildConntableDisorder(TestCase):

    def test_different_nonzero_parts_not_bonded(self):
        """Atoms in different non-zero disorder parts should NOT be bonded."""
        bonds = build_conntable(_two_coords(), ['C', 'O'], [1, 2])
        self.assertNotIn((0, 1), bonds)

    def test_same_nonzero_parts_bonded(self):
        """Atoms in the same non-zero disorder part should be bonded."""
        bonds = build_conntable(_two_coords(), ['C', 'O'], [1, 1])
        self.assertIn((0, 1), bonds)

    def test_nonzero_to_zero_part_bonded(self):
        """Atom in a non-zero part is allowed to bond to an atom in part 0."""
        bonds = build_conntable(_two_coords(), ['C', 'O'], [1, 0])
        self.assertIn((0, 1), bonds)


class TestBuildConntableNegativePart(TestCase):

    def test_neg_part_both_base_bonded(self):
        """Two base atoms with the same negative part should bond."""
        bonds = build_conntable(_two_coords(), ['C', 'O'], [-1, -1],
                                symmgen=[False, False])
        self.assertIn((0, 1), bonds)

    def test_neg_part_both_symmgen_bonded(self):
        """Two symmgen atoms with the same negative part should bond."""
        bonds = build_conntable(_two_coords(), ['C', 'O'], [-1, -1],
                                symmgen=[True, True])
        self.assertIn((0, 1), bonds)

    def test_neg_part_cross_boundary_excluded(self):
        """Base + symmgen with negative part should NOT be bonded."""
        bonds = build_conntable(_two_coords(), ['C', 'O'], [-1, -1],
                                symmgen=[False, True])
        self.assertNotIn((0, 1), bonds)

    def test_positive_parts_cross_boundary_allowed(self):
        """Normal atoms that happen to cross base/symmgen boundary are still bonded."""
        bonds = build_conntable(_two_coords(), ['C', 'O'], [0, 0],
                                symmgen=[False, True])
        self.assertIn((0, 1), bonds)

    def test_no_symmgen_no_filter(self):
        """When symmgen is None, the negative-part filter is not applied."""
        bonds = build_conntable(_two_coords(), ['C', 'O'], [-1, -1],
                                symmgen=None)
        self.assertIn((0, 1), bonds)


# ---------------------------------------------------------------------------
# Integration tests via Atoms.conntable property
# ---------------------------------------------------------------------------

class TestAtomConntableProperty(TestCase):

    def setUp(self):
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_returns_tuple(self):
        """conntable must return a tuple."""
        self.assertIsInstance(self.shx.atoms.conntable, tuple)

    def test_bonds_are_pairs(self):
        """Every entry in conntable must be a (int, int) pair."""
        for bond in self.shx.atoms.conntable:
            self.assertEqual(len(bond), 2)
            self.assertIsInstance(bond[0], int)
            self.assertIsInstance(bond[1], int)

    def test_nonzero_bond_count(self):
        """A real structure should have at least one bond."""
        self.assertGreater(len(self.shx.atoms.conntable), 0)

    def test_bond_distances_sensible(self):
        """All detected bonds between non-q-peak atoms should be between 0.5 and 4.0 Å."""
        atoms = self.shx.atoms.all_atoms
        for i, j in self.shx.atoms.conntable:
            ai, aj = atoms[i], atoms[j]
            if ai.qpeak or aj.qpeak:
                continue  # Q-peaks can sit on top of atoms; skip distance check
            d = np.linalg.norm([ai.xc - aj.xc, ai.yc - aj.yc, ai.zc - aj.zc])
            self.assertGreater(d, 0.5, f"Bond {ai.name}-{aj.name} too short: {d:.3f} Å")
            self.assertLess(d, 4.0, f"Bond {ai.name}-{aj.name} too long: {d:.3f} Å")

    def test_o1_c1_bonded(self):
        """O1 and C1 (indices 0, 1) in p21c.res should be bonded."""
        atoms = self.shx.atoms.all_atoms
        idx_o1 = next(i for i, a in enumerate(atoms) if a.name == 'O1')
        idx_c1 = next(i for i, a in enumerate(atoms) if a.name == 'C1')
        pair = tuple(sorted((idx_o1, idx_c1)))
        self.assertIn(pair, self.shx.atoms.conntable)

    def test_no_hh_bonds(self):
        """The connectivity table must not contain any H–H bonds."""
        atoms = self.shx.atoms.all_atoms
        for i, j in self.shx.atoms.conntable:
            both_h = atoms[i].is_hydrogen and atoms[j].is_hydrogen
            self.assertFalse(both_h,
                             f"Unexpected H–H bond: {atoms[i].name}–{atoms[j].name}")

    def test_symmetry_consistency_with_build_conntable(self):
        """conntable matches a manual call to build_conntable with same inputs."""
        atoms = self.shx.atoms.all_atoms
        coords = np.array([[at.xc, at.yc, at.zc] for at in atoms])
        types = [at.element for at in atoms]
        parts = [at.part.n for at in atoms]
        radii = np.array([at.radius for at in atoms])
        symmgen = [at.symmgen for at in atoms]
        expected = build_conntable(coords, types, parts, radii=radii, symmgen=symmgen)
        self.assertEqual(self.shx.atoms.conntable, expected)


# ---------------------------------------------------------------------------
# Tests for Atoms.bonds (human-friendly Bond list)
# ---------------------------------------------------------------------------

class TestAtomsBonds(TestCase):

    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_bonds_returns_list(self):
        self.assertIsInstance(self.shx.atoms.bonds, list)

    def test_bonds_items_are_bond_objects(self):
        for b in self.shx.atoms.bonds:
            self.assertIsInstance(b, Bond)

    def test_bond_count_matches_conntable(self):
        self.assertEqual(len(self.shx.atoms.bonds),
                         len(self.shx.atoms.conntable))

    def test_bond_distance_is_positive(self):
        for b in self.shx.atoms.bonds:
            self.assertGreater(b.distance, 0.0)

    def test_bond_distance_is_realistic(self):
        """All bonds should be between 0.01 Å (build_conntable threshold) and 4.0 Å."""
        for b in self.shx.atoms.bonds:
            self.assertGreater(b.distance, 0.01,
                               f"Bond {b!r} is unrealistically short")
            self.assertLess(b.distance, 4.0,
                            f"Bond {b!r} is unrealistically long")

    def test_bond_atoms_are_atom_objects(self):
        from shelxfile.atoms.atom import Atom
        for b in self.shx.atoms.bonds:
            self.assertIsInstance(b.atom1, Atom)
            self.assertIsInstance(b.atom2, Atom)

    def test_bond_str_contains_atom_names(self):
        b = self.shx.atoms.bonds[0]
        s = str(b)
        self.assertIn(b.atom1.fullname_short, s)
        self.assertIn(b.atom2.fullname_short, s)

    def test_bond_str_contains_angstrom(self):
        b = self.shx.atoms.bonds[0]
        self.assertIn('Å', str(b))

    def test_bond_repr(self):
        b = self.shx.atoms.bonds[0]
        r = repr(b)
        self.assertTrue(r.startswith('Bond('))
        self.assertIn('Å', r)

    def test_bond_unpack(self):
        """Bond supports tuple-style unpacking into (atom1, atom2, distance)."""
        b = self.shx.atoms.bonds[0]
        a1, a2, dist = b
        self.assertIs(a1, b.atom1)
        self.assertIs(a2, b.atom2)
        self.assertEqual(dist, b.distance)

    def test_bonds_sorted_by_name(self):
        """Bond list is sorted by atom1 name then atom2 name."""
        bonds = self.shx.atoms.bonds
        names = [(b.atom1.fullname_short, b.atom2.fullname_short) for b in bonds]
        self.assertEqual(names, sorted(names))

    def test_no_hh_bonds(self):
        for b in self.shx.atoms.bonds:
            self.assertFalse(b.atom1.is_hydrogen and b.atom2.is_hydrogen,
                             f"Unexpected H–H bond: {b!r}")


# ---------------------------------------------------------------------------
# Tests for Atoms.full_bond_list() — asymmetric unit + symmetry neighbors
# ---------------------------------------------------------------------------

class TestFullBondList(TestCase):
    """Uses I-43d.res which has genuine symmetry bonds."""

    def setUp(self) -> None:
        # I-43d has many SYMM cards → guaranteed symmetry bonds
        self.shx_symm = Shelxfile()
        self.shx_symm.read_file('tests/resources/I-43d.res')
        # p21c.res: all bonds are within the asymmetric unit (full molecule in AU)
        self.shx_plain = Shelxfile()
        self.shx_plain.read_file('tests/resources/p21c.res')

    # --- type and structure ---

    def test_returns_list(self):
        self.assertIsInstance(self.shx_symm.atoms.full_bond_list(), list)

    def test_items_are_symbond(self):
        for b in self.shx_symm.atoms.full_bond_list():
            self.assertIsInstance(b, SymBond)

    def test_symbond_is_subclass_of_bond(self):
        for b in self.shx_symm.atoms.full_bond_list():
            self.assertIsInstance(b, Bond)

    # --- counts ---

    def test_plain_bonds_match_conntable(self):
        """Plain (no-symm) bonds must be a subset of the conntable."""
        plain = [b for b in self.shx_plain.atoms.full_bond_list()
                 if not b.is_symmetry_bond]
        # full_bond_list excludes Q-peaks by default; conntable includes them.
        # The plain bonds should be a non-empty subset.
        self.assertGreater(len(plain), 0)
        self.assertLessEqual(len(plain), len(self.shx_plain.atoms.conntable))

    def test_has_symmetry_bonds(self):
        symm = [b for b in self.shx_symm.atoms.full_bond_list()
                if b.is_symmetry_bond]
        self.assertGreater(len(symm), 0)

    # --- symmetry label format ---

    def test_symm_label_is_string(self):
        for b in self.shx_symm.atoms.full_bond_list():
            self.assertIsInstance(b.symm_label, str)

    def test_plain_bond_has_empty_label(self):
        for b in self.shx_plain.atoms.full_bond_list():
            self.assertEqual('', b.symm_label)

    def test_symm_bond_label_contains_axis(self):
        """Symmetry labels must contain at least one axis letter."""
        symm = [b for b in self.shx_symm.atoms.full_bond_list()
                if b.is_symmetry_bond]
        for b in symm:
            self.assertTrue(
                any(ax in b.symm_label for ax in ('x', 'y', 'z')),
                f"Label {b.symm_label!r} contains no axis"
            )

    def test_symm_str_shows_label_in_brackets(self):
        symm = [b for b in self.shx_symm.atoms.full_bond_list()
                if b.is_symmetry_bond]
        for b in symm[:5]:
            s = str(b)
            self.assertIn('[', s)
            self.assertIn(b.symm_label, s)

    def test_plain_str_has_no_brackets(self):
        for b in self.shx_plain.atoms.full_bond_list():
            self.assertNotIn('[', str(b))

    # --- distance ---

    def test_distances_positive(self):
        for b in self.shx_symm.atoms.full_bond_list():
            self.assertGreater(b.distance, 0.0)

    def test_angstrom_in_str(self):
        b = self.shx_symm.atoms.full_bond_list()[0]
        self.assertIn('Å', str(b))

    # --- unpack ---

    def test_unpack_four_fields(self):
        b = self.shx_symm.atoms.full_bond_list()[0]
        a1, a2, dist, label = b
        self.assertIs(a1, b.atom1)
        self.assertIs(a2, b.atom2)
        self.assertEqual(dist, b.distance)
        self.assertEqual(label, b.symm_label)

    # --- sorting ---

    def test_sorted_by_name(self):
        bl = self.shx_symm.atoms.full_bond_list()
        names = [(b.atom1.fullname_short, b.atom2.fullname_short) for b in bl]
        self.assertEqual(names, sorted(names))

    # --- repr ---

    def test_repr_plain(self):
        b = [b for b in self.shx_plain.atoms.full_bond_list()][0]
        self.assertTrue(repr(b).startswith('SymBond('))

    def test_repr_symm_contains_bracket(self):
        b = [b for b in self.shx_symm.atoms.full_bond_list()
             if b.is_symmetry_bond][0]
        self.assertIn('[', repr(b))

    # --- with_qpeaks=False (default) ---

    def test_no_qpeaks_by_default(self):
        for b in self.shx_symm.atoms.full_bond_list():
            self.assertFalse(b.atom1.qpeak or b.atom2.qpeak,
                             f"Q-peak appeared in default bond list: {b!r}")


