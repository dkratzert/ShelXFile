"""
Additional tests to improve coverage for atoms, shelx, cards, and misc modules.
"""
import tempfile
import unittest
from pathlib import Path
from unittest import TestCase

from shelxfile import Shelxfile
from shelxfile.atoms.pairs import AtomPair
from shelxfile.misc.misc import (
    split_fvar_and_parameter, flatten, remove_file, multiline_test, wrap_line
)
from shelxfile.shelx.cards import (
    DAMP, HKLF, SUMP, BLOC, BOND, BIND, REM, FMAP, MOVE, MERG, SIZE,
    HTAB, GRID, SHEL, STIR, TWIN, SWAT, WIGL, WPDB, PRIG, CONF,
    XNPD, FVAR, ACTA, Restraints, SAME, RIGU, SIMU, DELU, ISOR, FLAT, BUMP, CHIV, EADP, EXYZ,
    RTAB, HFIX, MORE, FREE,
)


# =============================================================================
# Tests for shelxfile/atoms/atom.py
# =============================================================================

class TestAtomProperties(TestCase):
    """Tests for uncovered Atom properties."""

    def setUp(self):
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_atomid_regular_atom(self):
        a = self.shx.atoms.get_atom_by_name('O1_4')
        self.assertGreater(a.atomid, 0)

    def test_fullname_with_residue(self):
        a = self.shx.atoms.get_atom_by_name('F1_2')
        self.assertEqual('F1_2', a.fullname)

    def test_fullname_short_residue_zero(self):
        a = self.shx.atoms.get_atom_by_name('Ga1')
        self.assertEqual('GA1', a.fullname_short)

    def test_fullname_short_with_residue(self):
        a = self.shx.atoms.get_atom_by_name('F1_2')
        self.assertEqual('F1_2', a.fullname_short)

    def test_chain_id(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        # chain_id should be None for most atoms in this file
        self.assertIsNone(a.chain_id)

    def test_fvar(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        # sof is -31.0 => fvar is -3
        self.assertEqual(-3, a.fvar)

    def test_fvar_regular(self):
        a = self.shx.atoms.get_atom_by_name('Ga1')
        # sof is 11.0 => fvar = 1
        self.assertEqual(1, a.fvar)

    def test_occupancy_positive_fvar(self):
        a = self.shx.atoms.get_atom_by_name('C1_1')
        # This atom has positive fvar occupancy
        occ = a.occupancy
        self.assertIsInstance(occ, float)
        self.assertGreater(occ, 0)

    def test_occupancy_negative_fvar(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        # negative sof => should use negative fvar calc
        occ = a.occupancy
        self.assertIsInstance(occ, float)
        self.assertGreater(occ, 0)

    def test_occupancy_setter(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        a.occupancy = 0.75
        self.assertEqual(0.75, a._occupancy)

    def test_is_hydrogen_true(self):
        a = self.shx.atoms.get_atom_by_name('H34')
        self.assertTrue(a.is_hydrogen)

    def test_is_hydrogen_false(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        self.assertFalse(a.is_hydrogen)

    def test_ishydrogen_alias(self):
        a = self.shx.atoms.get_atom_by_name('H34')
        self.assertTrue(a.ishydrogen)

    def test_is_isotropic_true(self):
        a = self.shx.atoms.get_atom_by_name('H34')
        self.assertTrue(a.is_isotropic)

    def test_is_isotropic_false(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        self.assertFalse(a.is_isotropic)

    def test_uvals_as_list(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        ulist = a.uvals_as_list
        self.assertEqual(3, len(ulist))
        self.assertEqual(3, len(ulist[0]))
        # Check symmetry: U12 == U21 etc.
        self.assertEqual(ulist[0][1], ulist[1][0])  # U12 == U21
        self.assertEqual(ulist[0][2], ulist[2][0])  # U13 == U31
        self.assertEqual(ulist[1][2], ulist[2][1])  # U23 == U32

    def test_atom_iter(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        parts = list(a)
        self.assertIn('Atom', parts)

    def test_atom_repr(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        self.assertTrue(repr(a).startswith('Atom ID:'))

    def test_frac_coords_setter(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        a.frac_coords = [0.1, 0.2, 0.3]
        self.assertEqual((0.1, 0.2, 0.3), a.frac_coords)

    def test_element_property(self):
        a = self.shx.atoms.get_atom_by_name('F1_2')
        self.assertEqual('F', a.element)

    def test_an_property(self):
        a = self.shx.atoms.get_atom_by_name('O1')
        self.assertEqual(8, a.an)

    def test_radius(self):
        a = self.shx.atoms.get_atom_by_name('O1')
        self.assertEqual(0.73, a.radius)

    def test_ueq(self):
        a = self.shx.atoms.get_atom_by_name('Ga1')
        ueq = a.ueq
        self.assertIsInstance(ueq, float)
        self.assertGreater(ueq, 0)

    def test_ucif(self):
        a = self.shx.atoms.get_atom_by_name('Ga1')
        ucif = a.ucif
        self.assertIsNotNone(ucif)

    def test_qpeak_atom(self):
        qpeaks = self.shx.atoms.q_peaks
        self.assertGreater(len(qpeaks), 0)
        q = qpeaks[0]
        self.assertTrue(q.qpeak)
        self.assertIn('Q', str(q))


class TestAtomNameProperty(TestCase):
    def setUp(self):
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_name_getter(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        self.assertEqual('C1', a.name)

    def test_name_setter_valid(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        a.name = 'C99'
        self.assertEqual('C99', a.name)

    def test_name_setter_invalid(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        a.name = 'C99_3'  # Invalid: contains underscore
        self.assertEqual('C1', a.name)  # Should stay C1


# =============================================================================
# Tests for shelxfile/atoms/atoms.py
# =============================================================================

class TestAtomsCollection(TestCase):
    """Tests for Atoms class uncovered methods."""

    def setUp(self):
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_repr_non_empty(self):
        r = repr(self.shx.atoms)
        self.assertIn('O1', r)

    def test_repr_empty(self):
        shx = Shelxfile()
        r = repr(shx.atoms)
        self.assertEqual('No Atoms in file.', r)

    def test_getitem(self):
        a = self.shx.atoms.get_atom_by_name('C1_4')
        atomid = a.atomid
        same_atom = self.shx.atoms[atomid]
        self.assertEqual(str(a), str(same_atom))

    def test_len(self):
        self.assertEqual(148, len(self.shx.atoms))

    def test_get_atom_by_name_none_for_gt(self):
        self.assertIsNone(self.shx.atoms.get_atom_by_name('>'))

    def test_get_atom_by_name_none_for_lt(self):
        self.assertIsNone(self.shx.atoms.get_atom_by_name('<'))

    def test_get_atom_by_name_not_found(self):
        self.assertIsNone(self.shx.atoms.get_atom_by_name('ZZZZZ'))

    def test_get_atom_by_id_not_found(self):
        self.assertIsNone(self.shx.atoms.get_atom_by_id(99999))

    def test_nameslist(self):
        names = self.shx.atoms.nameslist
        self.assertIsInstance(names, tuple)
        self.assertIn('O1_4', names)
        self.assertIn('GA1_0', names)

    def test_distance_not_found(self):
        d = self.shx.atoms.distance('NONEXISTENT1', 'NONEXISTENT2')
        self.assertEqual(0.0, d)

    def test_has_atom_false(self):
        self.assertFalse(self.shx.atoms.has_atom('XXXXXXX'))


# =============================================================================
# Tests for shelxfile/atoms/pairs.py
# =============================================================================

class TestAtomPair(TestCase):
    def test_repr_both_atoms(self):
        p = AtomPair('C1', 'C2')
        self.assertEqual('C1 C2', repr(p))

    def test_repr_missing_atom(self):
        p = AtomPair(None, 'C2')
        self.assertEqual('', repr(p))

    def test_len_both(self):
        p = AtomPair('C1', 'C2')
        self.assertEqual(2, len(p))

    def test_len_none(self):
        p = AtomPair(None, 'C2')
        self.assertEqual(0, len(p))

    def test_len_both_none(self):
        p = AtomPair(None, None)
        self.assertEqual(0, len(p))


# =============================================================================
# Tests for shelxfile/shelx/shelx.py
# =============================================================================

class TestShelxfileOperations(TestCase):
    """Tests for Shelxfile methods with uncovered lines."""

    def setUp(self):
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_verbose_and_debug_raises(self):
        with self.assertRaises(ValueError):
            Shelxfile(verbose=True, debug=True)

    def test_add_atom(self):
        count_before = len(self.shx.atoms)
        self.shx.add_atom(name='C99', coordinates=[0.1, 0.2, 0.3], element='C',
                          uvals=[0.04, 0.0, 0.0, 0.0, 0.0, 0.0])
        # The atom list should be repopulated on next access:
        # Since add_atom uses _append_card at position 0, let's check it was added
        found = self.shx.atoms.get_atom_by_name('C99')
        self.assertIsNotNone(found)
        self.assertEqual('C99', found.name)
        self.assertGreater(len(self.shx.atoms), count_before)

    def test_add_atom_with_part(self):
        self.shx.add_atom(name='N1', coordinates=[0.5, 0.5, 0.5], element='N', part=1, sof=21.0,
                          uvals=[0.04, 0.0, 0.0, 0.0, 0.0, 0.0])
        found = self.shx.atoms.get_atom_by_name('N1')
        self.assertIsNotNone(found)

    def test_frac_to_cart(self):
        result = self.shx.frac_to_cart([0.1, 0.2, 0.3])
        self.assertEqual(3, len(result))

    def test_repr(self):
        r = repr(self.shx)
        self.assertIn('TITL', r)
        self.assertIn('CELL', r)

    def test_sum_formula(self):
        sf = self.shx.sum_formula
        self.assertIn('C', sf)
        self.assertIsNotNone(self.shx.formula_weight)

    def test_sum_formula_exact(self):
        sf = self.shx.sum_formula_exact
        self.assertIn('C', sf)
        self.assertIn('H', sf)

    def test_sum_formula_exact_as_dict(self):
        d = self.shx.sum_formula_exact_as_dict()
        self.assertIsInstance(d, dict)
        self.assertIn('C', d)

    def test_insert_anis_default(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        shx.insert_anis()
        # Check that ANIS was inserted
        found = any('ANIS' in str(x) for x in shx._reslist)
        self.assertTrue(found)

    def test_insert_anis_with_atoms(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        shx.insert_anis(atoms='C1 C2')
        found = any('ANIS C1 C2' in str(x) for x in shx._reslist)
        self.assertTrue(found)

    def test_insert_anis_with_residue(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        shx.insert_anis(atoms='$C', residue='CCF3')
        found = any('ANIS_CCF3 $C' in str(x) for x in shx._reslist)
        self.assertTrue(found)

    def test_add_line(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        old_len = len(shx._reslist)
        shx.add_line(5, 'REM Test Line')
        self.assertEqual(old_len + 1, len(shx._reslist))
        self.assertEqual('REM Test Line', shx._reslist[6])

    def test_replace_line(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        old_plan = shx.plan
        shx.replace_line(old_plan, 'PLAN 99')
        idx = shx.index_of('PLAN 99')
        self.assertGreater(idx, 0)

    def test_index_of(self):
        idx = self.shx.index_of(self.shx.plan)
        self.assertGreater(idx, 0)

    def test_weight_converged_true(self):
        self.assertTrue(self.shx._weight_converged([0.0, 0.0]))

    def test_weight_converged_false(self):
        self.assertFalse(self.shx._weight_converged([0.1, 0.0]))

    def test_update_weight(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        # No suggested weight, should just return
        shx.wght_suggested = None
        shx.update_weight()  # Should not raise

    def test_is_atom_coordinates_unrealistic(self):
        self.assertTrue(Shelxfile._coordinates_are_unrealistic(['O1', '3', '5.0', '0.5', '0.5']))
        self.assertFalse(Shelxfile._coordinates_are_unrealistic(['O1', '3', '0.5', '0.5', '0.5']))

    def test_write_shelx_file_no_reslist(self):
        shx = Shelxfile()
        # No file loaded, _reslist is empty
        result = shx.write_shelx_file('nonexistent.ins')
        self.assertIsNone(result)

    def test_write_shelx_file(self):
        with tempfile.NamedTemporaryFile(suffix='.ins', delete=False) as f:
            tmppath = Path(f.name)
        try:
            self.shx.write_shelx_file(str(tmppath))
            self.assertTrue(tmppath.exists())
            content = tmppath.read_text()
            self.assertIn('TITL', content)
            self.assertIn('CELL', content)
        finally:
            tmppath.unlink(missing_ok=True)

    def test_is_atom_blank_line(self):
        self.assertFalse(Shelxfile.is_atom(''))

    def test_is_atom_with_exclamation(self):
        self.assertFalse(Shelxfile.is_atom('O1  3  0.1  0.2  0.3  !'))


class TestShelxfileResiduals(TestCase):

    def setUp(self):
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_r1(self):
        self.assertEqual(0.04, self.shx.R1)

    def test_wr2(self):
        self.assertEqual(0.1005, self.shx.wr2)

    def test_goof(self):
        self.assertEqual(1.016, self.shx.goof)

    def test_space_group(self):
        self.assertEqual('P2(1)/c', self.shx.space_group)

    def test_wavelength(self):
        self.assertEqual(0.71073, self.shx.wavelength)

    def test_temp(self):
        self.assertAlmostEqual(-173.18, self.shx.temp, places=2)

    def test_temp_in_kelvin(self):
        self.assertAlmostEqual(99.97, self.shx.temp_in_kelvin, places=2)

    def test_data(self):
        self.assertEqual(10786, self.shx.data)

    def test_parameters(self):
        self.assertEqual(945, self.shx.parameters)

    def test_dat_to_param(self):
        self.assertAlmostEqual(11.4138, self.shx.dat_to_param, places=3)

    def test_num_restraints(self):
        self.assertEqual(1842, self.shx.num_restraints)

    def test_highest_peak(self):
        self.assertAlmostEqual(0.407, self.shx.highest_peak, places=3)

    def test_deepest_hole(self):
        self.assertAlmostEqual(-0.691, self.shx.deepest_hole, places=3)


# =============================================================================
# Tests for shelxfile/shelx/cards.py
# =============================================================================

class TestRestraintsClass(TestCase):
    def test_empty_restraints(self):
        r = Restraints()
        self.assertEqual('No Restraints in file.', repr(r))

    def test_restraints_getitem(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        r = shx.restraints[0]
        self.assertIsNotNone(r)

    def test_restraints_iter(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        count = sum(1 for _ in shx.restraints)
        self.assertGreater(count, 0)


class TestDAMP(TestCase):
    def test_damp_with_limse(self):
        shx = Shelxfile()
        d = DAMP(shx, 'DAMP 0.5 20'.split())
        self.assertEqual(0.5, d.damp)
        self.assertEqual(20, d.limse)
        self.assertEqual('DAMP  0.5 20', repr(d))

    def test_damp_without_limse(self):
        shx = Shelxfile()
        d = DAMP(shx, 'DAMP 0.7'.split())
        self.assertEqual(0.7, d.damp)
        self.assertEqual(0, d.limse)
        self.assertEqual('DAMP  0.7', repr(d))


class TestHKLF(TestCase):
    def test_hklf_basic(self):
        shx = Shelxfile()
        h = HKLF(shx, 'HKLF 4'.split())
        self.assertEqual(4, h.n)
        self.assertEqual(1, h.s)
        self.assertEqual([1, 0, 0, 0, 1, 0, 0, 0, 1], h.matrix)

    def test_hklf_repr(self):
        shx = Shelxfile()
        h = HKLF(shx, 'HKLF 4'.split())
        self.assertIn('HKLF 4', repr(h))


class TestSUMP(TestCase):
    def test_sump_basic(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        s = SUMP(shx, 'SUMP 1 0.001 1 2 1 3'.split())
        self.assertEqual(1.0, s.c)
        self.assertEqual(0.001, s.sigma)
        self.assertEqual([1.0, 2], s[0])
        self.assertEqual([1.0, 3], s[1])


class TestBLOC(TestCase):
    def test_bloc(self):
        shx = Shelxfile()
        b = BLOC(shx, 'BLOC 1 5 C1 C2'.split())
        self.assertEqual(1, b.n1)
        self.assertEqual(5, b.n2)
        self.assertEqual(['C1', 'C2'], b.atoms)


class TestBOND(TestCase):
    def test_bond(self):
        shx = Shelxfile()
        b = BOND(shx, 'BOND C1 C2 C3'.split())
        self.assertEqual(['C1', 'C2', 'C3'], b.atoms)


class TestBIND(TestCase):
    def test_bind(self):
        shx = Shelxfile()
        b = BIND(shx, 'BIND C1 C2'.split())
        self.assertEqual(['C1', 'C2'], b.atoms)


class TestREM(TestCase):
    def test_rem(self):
        shx = Shelxfile()
        r = REM(shx, 'REM This is a remark'.split())
        self.assertEqual('REM This is a remark', str(r))


class TestFMAP(TestCase):
    def test_fmap_all_params(self):
        shx = Shelxfile()
        f = FMAP(shx, 'FMAP 2 1 53'.split())
        self.assertEqual(2, f.code)
        self.assertEqual(1, f.axis)
        self.assertEqual(53, f.nl)


class TestMOVE(TestCase):
    def test_move(self):
        shx = Shelxfile()
        m = MOVE(shx, 'MOVE 0.1 0.2 0.3 1'.split())
        self.assertEqual([0.1, 0.2, 0.3], m.dxdydz)
        self.assertEqual(1, m.sign)


class TestMERG(TestCase):
    def test_merg(self):
        shx = Shelxfile()
        m = MERG(shx, 'MERG 3'.split())
        self.assertEqual(3, m.n)


class TestSIZE(TestCase):
    def test_size_all_params(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        self.assertIsNotNone(shx.size)
        self.assertEqual(0.12, shx.size.dx)
        self.assertEqual(0.23, shx.size.dy)
        self.assertEqual(0.33, shx.size.dz)
        self.assertEqual(0.33, shx.size.max)
        self.assertEqual(0.23, shx.size.mid)
        self.assertEqual(0.12, shx.size.min)

    def test_size_as_text(self):
        shx = Shelxfile()
        s = SIZE(shx, 'SIZE 0.1 0.2 0.3'.split())
        self.assertEqual('SIZE 0.1 0.2 0.3', s._as_text())


class TestHTAB(TestCase):
    def test_htab_dh(self):
        shx = Shelxfile()
        h = HTAB(shx, 'HTAB 2.5'.split())
        self.assertEqual(2.5, h.dh)

    def test_htab_atoms(self):
        shx = Shelxfile()
        h = HTAB(shx, 'HTAB O1 N1'.split())
        self.assertEqual('O1', h.donor)
        self.assertEqual('N1', h.acceptor)


class TestGRID(TestCase):
    def test_grid(self):
        shx = Shelxfile()
        g = GRID(shx, 'GRID 1.0 2.0 3.0 4.0 5.0 6.0'.split())
        self.assertEqual(1.0, g.sl)
        self.assertEqual(2.0, g.sa)
        self.assertEqual(3.0, g.sd)
        self.assertEqual(4.0, g.dl)
        self.assertEqual(5.0, g.da)
        self.assertEqual(6.0, g.dd)


class TestSHEL(TestCase):
    def test_shel(self):
        shx = Shelxfile()
        s = SHEL(shx, 'SHEL 99.0 0.5'.split())
        self.assertEqual(99.0, s.lowres)
        self.assertEqual(0.5, s.highres)


class TestSTIR(TestCase):
    def test_stir(self):
        shx = Shelxfile()
        s = STIR(shx, 'STIR 0.8 0.02'.split())
        self.assertEqual(0.8, s.sres)
        self.assertEqual(0.02, s.step)


class TestTWIN(TestCase):
    def test_twin(self):
        shx = Shelxfile()
        t = TWIN(shx, 'TWIN -1 0 0 0 -1 0 0 0 -1 2'.split())
        self.assertEqual([-1, 0, 0, 0, -1, 0, 0, 0, -1], t.matrix)
        self.assertEqual(2, t.n_value)
        self.assertEqual(1, t.allowed_N)


class TestSWAT(TestCase):
    def test_swat(self):
        shx = Shelxfile()
        s = SWAT(shx, 'SWAT 0.5 3.0'.split())
        self.assertEqual(0.5, s.g)
        self.assertEqual(3.0, s.U)


class TestWIGL(TestCase):
    def test_wigl(self):
        shx = Shelxfile()
        w = WIGL(shx, 'WIGL 0.3 0.4'.split())
        self.assertEqual(0.3, w.d)
        self.assertEqual(0.4, w.dU)


class TestWPDB(TestCase):
    def test_wpdb(self):
        shx = Shelxfile()
        w = WPDB(shx, 'WPDB 2'.split())
        self.assertEqual(2, w.n)


class TestPRIG(TestCase):
    def test_prig(self):
        shx = Shelxfile()
        p = PRIG(shx, 'PRIG 5'.split())
        self.assertEqual(5, p.p)


class TestCONF(TestCase):
    def test_conf(self):
        shx = Shelxfile()
        c = CONF(shx, 'CONF'.split())
        self.assertEqual('CONF', str(c))


class TestXNPD(TestCase):
    def test_xnpd_default(self):
        shx = Shelxfile()
        x = XNPD(shx, 'XNPD'.split())
        self.assertEqual(-0.001, x.Umin)

    def test_xnpd_custom(self):
        shx = Shelxfile()
        x = XNPD(shx, 'XNPD -0.005'.split())
        self.assertEqual(-0.005, x.Umin)


class TestFVAR(TestCase):
    def test_fvar_str(self):
        f = FVAR(1, 0.5)
        self.assertEqual('0.5', str(f))

    def test_fvar_repr(self):
        f = FVAR(2, 0.753)
        self.assertEqual('0.753', repr(f))


class TestFVARs(TestCase):
    def setUp(self):
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_len(self):
        self.assertGreater(len(self.shx.fvars), 0)

    def test_getitem(self):
        # SHELXL counts fvars from 1
        val = self.shx.fvars[1]
        self.assertIsNotNone(val)

    def test_str(self):
        s = str(self.shx.fvars)
        self.assertIn('FVAR', s)

    def test_repr(self):
        r = repr(self.shx.fvars)
        self.assertIsNotNone(r)

    def test_fvars_used(self):
        result = self.shx.fvars.fvars_used()
        self.assertIsInstance(result, dict)


class TestACTAClass(TestCase):
    def test_acta_with_twotheta(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        a = ACTA(shx, 'ACTA 45'.split())
        self.assertEqual('ACTA 45', repr(a))
        self.assertEqual('ACTA 45', str(a))

    def test_acta_without_twotheta(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        a = ACTA(shx, 'ACTA'.split())
        self.assertEqual('ACTA', repr(a))


class TestWGHTClass(TestCase):
    def setUp(self):
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_wght_a(self):
        self.assertEqual(0.049, self.shx.wght.a)

    def test_wght_b(self):
        self.assertEqual(0.0, self.shx.wght.b)


class TestDEFS(TestCase):
    def test_defs(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        self.assertIsNotNone(shx.defs)
        self.assertEqual(0.0234, shx.defs.sd)
        self.assertIsNotNone(shx.defs.all)


class TestRestraintTypes(TestCase):
    """Tests for various restraint card types."""

    def setUp(self):
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_same(self):
        s = SAME(self.shx, 'SAME 0.03 0.05 C1 C2 C3'.split())
        self.assertEqual(0.03, s.s1)
        self.assertEqual(0.05, s.s2)
        self.assertEqual(['C1', 'C2', 'C3'], s.atoms)

    def test_rigu(self):
        r = RIGU(self.shx, 'RIGU 0.005 0.006 C1 C2'.split())
        self.assertEqual(0.005, r.s1)
        self.assertEqual(0.006, r.s2)

    def test_simu(self):
        s = SIMU(self.shx, 'SIMU 0.05 0.09 3.0 C1 C2'.split())
        self.assertEqual(0.05, s.s)
        self.assertEqual(0.09, s.st)
        self.assertEqual(3.0, s.dmax)

    def test_delu(self):
        d = DELU(self.shx, 'DELU 0.02 0.03 C1 C2'.split())
        self.assertEqual(0.02, d.s1)
        self.assertEqual(0.03, d.s2)

    def test_isor(self):
        i = ISOR(self.shx, 'ISOR 0.15 0.25 C1 C2'.split())
        self.assertEqual(0.15, i.s)
        self.assertEqual(0.25, i.st)

    def test_flat(self):
        f = FLAT(self.shx, 'FLAT 0.2 C1 C2 C3 C4'.split())
        self.assertEqual(0.2, f.s)
        self.assertEqual(['C1', 'C2', 'C3', 'C4'], f.atoms)

    def test_bump(self):
        b = BUMP(self.shx, 'BUMP 0.03'.split())
        self.assertEqual(0.03, b.s)

    def test_chiv(self):
        c = CHIV(self.shx, 'CHIV 0 0.2 C1'.split())
        self.assertIsNotNone(c)

    def test_eadp(self):
        e = EADP(self.shx, 'EADP C1 C2'.split())
        self.assertEqual(['C1', 'C2'], e.atoms)

    def test_exyz(self):
        e = EXYZ(self.shx, 'EXYZ C1 C2'.split())
        self.assertEqual(['C1', 'C2'], e.atoms)


class TestRTAB(TestCase):
    def test_rtab(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        r = RTAB(shx, 'RTAB omeg C1 C2 C3 C4'.split())
        self.assertEqual('omeg', r.code)
        self.assertEqual(['C1', 'C2', 'C3', 'C4'], r.atoms)


class TestLSCycles(TestCase):
    def test_ls_cycles(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        self.assertIsNotNone(shx.cycles)
        self.assertIsNotNone(shx.cycles.number)


class TestCELLCard(TestCase):
    def test_cell(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        c = shx.cell
        self.assertEqual(10.5086, c.a)
        self.assertEqual(20.9035, c.b)
        self.assertEqual(20.5072, c.c)
        self.assertEqual(90.0, c.alpha)
        self.assertEqual(94.13, c.beta)
        self.assertEqual(90.0, c.gamma)
        self.assertAlmostEqual(4493.0474, c.volume, places=4)

    def test_cell_iter(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        cell_list = list(shx.cell)
        self.assertEqual(6, len(cell_list))
        self.assertEqual(10.5086, cell_list[0])


class TestZERR(TestCase):
    def test_zerr(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        self.assertIsNotNone(shx.zerr)
        self.assertEqual(4, shx.Z)


class TestLATT(TestCase):
    def test_latt_centric(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        self.assertTrue(shx.latt.centric)


class TestSFACTable(TestCase):
    def test_sfac_elements(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        self.assertEqual('SFAC C  H  O  F  Al  Ga', str(shx.sfac_table))
        elements = shx.sfac_table.elements_list
        self.assertEqual(['C', 'H', 'O', 'F', 'AL', 'GA'], elements)


class TestUNIT(TestCase):
    def test_unit(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        self.assertEqual('UNIT 1  2  3  4  5  6', str(shx.unit))
        self.assertEqual([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], shx.unit.values)


class TestHFIX(TestCase):
    def test_hfix(self):
        shx = Shelxfile()
        h = HFIX(shx, 'HFIX 23 C1 C2'.split())
        self.assertEqual([23], h.params)
        self.assertEqual(['C1', 'C2'], h.atoms)
        self.assertIn('HFIX', repr(h))


class TestMORE(TestCase):
    def test_more(self):
        shx = Shelxfile()
        m = MORE(shx, 'MORE 2'.split())
        self.assertEqual(2, m.m)


class TestFREE(TestCase):
    def test_free(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        f = FREE(shx, 'FREE C1 C2'.split())
        self.assertEqual('C1', f.atom1)
        self.assertEqual('C2', f.atom2)


class TestCommandBase(TestCase):
    """Test Command base class methods."""

    def test_set(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        shx.plan.set('PLAN 50')
        self.assertEqual('PLAN 50', repr(shx.plan))

    def test_split(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        parts = shx.plan.split()
        self.assertEqual(['PLAN', '20'], parts)

    def test_iter(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        parts = list(shx.plan)
        self.assertIn('PLAN', parts)


class TestRestraintBase(TestCase):
    """Test Restraint base class methods."""

    def test_iter(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        r = shx.restraints[0]
        parts = list(r)
        self.assertGreater(len(parts), 0)

    def test_split(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        r = shx.restraints[0]
        parts = r.split()
        self.assertIsInstance(parts, list)


# =============================================================================
# Tests for shelxfile/misc/misc.py
# =============================================================================

class TestSplitFvarAndParameter(TestCase):
    def test_positive(self):
        self.assertEqual((3, 0.5), split_fvar_and_parameter(30.5))

    def test_negative(self):
        self.assertEqual((-3, -0.5), split_fvar_and_parameter(-30.5))

    def test_eleven(self):
        self.assertEqual((1, 1.0), split_fvar_and_parameter(11.0))

    def test_negative_eleven(self):
        self.assertEqual((-1, -1.0), split_fvar_and_parameter(-11.0))


class TestFlatten(TestCase):
    def test_flatten_nested(self):
        self.assertEqual([1, 2, 3, 4], flatten([[1, 2], [3, [4]]]))

    def test_flatten_empty(self):
        self.assertEqual([], flatten([]))

    def test_flatten_no_nesting(self):
        self.assertEqual([1, 2, 3], flatten([1, 2, 3]))


class TestRemoveFile(TestCase):
    def test_remove_existing_file(self):
        with tempfile.NamedTemporaryFile(delete=False) as f:
            path = f.name
        self.assertTrue(remove_file(path))

    def test_remove_nonexistent_file(self):
        result = remove_file(str(Path(tempfile.gettempdir()) / 'nonexistent_file_xyz123.txt'))
        self.assertFalse(result)


class TestWrapLine(TestCase):
    def test_short_line(self):
        self.assertEqual('SHORT LINE', wrap_line('SHORT LINE'))

    def test_long_line(self):
        line = 'A' * 100
        wrapped = wrap_line(line)
        self.assertIn(' =\n', wrapped)


class TestMultilineTest(TestCase):
    def test_multiline_true(self):
        self.assertTrue(multiline_test('C1  1  0.1  0.2  0.3  11.0  0.05 ='))

    def test_multiline_false(self):
        self.assertFalse(multiline_test('C1  1  0.1  0.2  0.3  11.0  0.05'))

    def test_multiline_rem_not_dsr(self):
        self.assertFalse(multiline_test('REM this has = sign'))


# =============================================================================
# Tests for edge cases and integration
# =============================================================================

class TestReadString(TestCase):
    def test_read_string(self):
        shx = Shelxfile()
        content = Path('tests/resources/p21c.res').read_text(encoding='latin1')
        shx.read_string(content)
        self.assertEqual(148, len(shx.atoms))
        self.assertIsNotNone(shx.cell)


class TestSymmCards(TestCase):
    def test_symmcards(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        self.assertGreater(len(shx.symmcards), 0)
        s = str(shx.symmcards)
        self.assertIn('1  0  0', s)


class TestDifferentFiles(TestCase):
    """Test parsing different res files for broader coverage."""

    def test_I43d_file(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/I-43d.res')
        self.assertGreater(len(shx.atoms), 0)

    def test_sad_final_file(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/sad-final.res')
        self.assertGreater(len(shx.atoms), 0)

    def test_2240189_file(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/2240189.res')
        self.assertGreater(len(shx.atoms), 0)


if __name__ == '__main__':
    unittest.main()
