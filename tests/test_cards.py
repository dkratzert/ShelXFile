from pathlib import Path
from unittest import TestCase

from shelxfile.refine.refine import ShelxlRefine
from shelxfile.shelx.cards import RESI, ANIS, ABIN, SADI, DFIX, AFIX
from shelxfile.shelx.shelx import Shelxfile


class TestCELL(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_volume(self):
        self.assertEqual(4493.0474, round(self.shx.cell.volume, 4))


class TestRESI(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_get_resi_1_tol(self):
        r = RESI(None, 'RESI 1 TOL'.split())
        self.assertEqual(('TOL', 1, None, None), (r.residue_class, r.residue_number, r.chain_id, r.alias))

    def test_tol_1(self):
        r = RESI(None, 'RESI TOL 1'.split())
        self.assertEqual(('TOL', 1, None, None), (r.residue_class, r.residue_number, r.chain_id, r.alias))

    def test_resi_chainid_tol(self):
        r = RESI(None, 'RESI A:100 TOL'.split())
        self.assertEqual(('TOL', 100, 'A', None), (r.residue_class, r.residue_number, r.chain_id, r.alias))

    def test_resi__negative_chain_tol(self):
        r = RESI(None, 'RESI -10 TOL'.split())
        self.assertEqual(('TOL', -10, None, None), (r.residue_class, r.residue_number, r.chain_id, r.alias))

    def test_resi_b_chain_tol(self):
        r = RESI(None, 'RESI b:-10 TOL'.split())
        self.assertEqual(('TOL', -10, 'b', None), (r.residue_class, r.residue_number, r.chain_id, r.alias))


class TestFMAP(TestCase):

    def test_fmap_number(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        self.assertEqual(2.0, shx.fmap.code)


class TestACTA(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')
        self.ref = ShelxlRefine(self.shx, Path('./tests/resources/p21c.res'))

    def test_acta(self):
        self.assertEqual('ACTA 45', self.shx._reslist[12].__repr__())
        self.assertEqual('ACTA 45', self.shx.acta.__repr__())

    def test_remove_acta(self):
        self.ref.remove_acta_card(self.shx.acta)
        self.assertEqual('SIZE 0.12 0.23 0.33', self.shx._reslist[12].__repr__())

    def test_restore_acta(self):
        self.ref.remove_acta_card(self.shx.acta)
        self.ref.restore_acta_card()
        self.assertEqual(8, self.shx.index_of(self.shx.acta))
        self.assertEqual(
            ['UNIT 1  2  3  4  5  6', 'ACTA 45', 'LIST 4 ! automatically inserted. Change 6 to 4 for CHECKCIF!!'],
            [str(x) for x in self.shx._reslist[7:10]])


class TestLSCycles(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')
        self.ref = ShelxlRefine(self.shx, Path('./resources/p21c.res'))

    def test_set_refine_cycles(self):
        self.shx.cycles.set_refine_cycles(44)
        self.assertEqual('L.S. 44', self.shx._reslist[self.shx.cycles.index].__repr__())


class TestSFACTable(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_sfac_table(self):
        self.assertEqual('SFAC C  H  O  F  Al  Ga', str(self.shx.sfac_table))

    def test_unit_table(self):
        self.assertEqual('UNIT 1  2  3  4  5  6', str(self.shx.unit))

    def test_sfac_add_element(self):
        self.shx.sfac_table.add_element('Au')
        self.assertEqual('SFAC C  H  O  F  Al  Ga  Au', str(self.shx.sfac_table))

    def test_unit_add_element(self):
        self.shx.sfac_table.add_element('Au')
        self.assertEqual('UNIT 1  2  3  4  5  6  1', self.shx.unit.__repr__())

    def test_unit_add_element_string_repr(self):
        self.shx.sfac_table.add_element('Au')
        self.assertEqual('UNIT 1  2  3  4  5  6  1', str(self.shx.unit))


class TestWGHT(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_wght(self):
        self.assertEqual('WGHT   0.049 0.0', self.shx.wght.__repr__())

    def test_wght_differene(self):
        self.assertEqual([0.0, 0.0], self.shx.wght.difference())

    def test_wght_set(self):
        self.shx.wght.set('WGHT 0.1 99')
        self.assertEqual('WGHT   0.1 99.0', str(self.shx.wght))
        self.assertEqual(0.1, self.shx.wght.a)
        self.assertEqual(99.0, self.shx.wght.b)


class TestANIS(TestCase):

    def test_anis_n(self):
        a = ANIS(Shelxfile(), 'ANIS 12'.split())
        self.assertEqual([], a.atoms)
        self.assertEqual(12, a.n)

    def test_anis_atoms(self):
        a = ANIS(Shelxfile(), 'ANIS C12 C23 Al3a'.split())
        self.assertEqual(['C12', 'C23', 'Al3a'], a.atoms)
        self.assertEqual(False, hasattr(a, 'n'))


class TestABIN(TestCase):
    def test_abin(self):
        a = ABIN(Shelxfile(), 'ABIN 1.234 4'.split())
        self.assertEqual(1.234, a.n1)
        self.assertEqual(4, a.n2)


class TestSADI(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/I-43d.res')

    def test_sadi_normal(self):
        a = SADI(self.shx, 'SADI C1 C2 C2 C3 C3 C4'.split())
        self.assertEqual(['C1', 'C2', 'C2', 'C3', 'C3', 'C4'], a.atoms)
        self.assertEqual('', a.residue_class)
        self.assertEqual([0], a.residue_number)

    def test_sadi_normal_with_resinum(self):
        a = SADI(self.shx, 'SADI_1 C1 C2 C2 C3 C3 C4'.split())
        self.assertEqual(['C1', 'C2', 'C2', 'C3', 'C3', 'C4'], a.atoms)
        self.assertEqual('', a.residue_class)
        self.assertEqual([1], a.residue_number)

    def test_sadi_normal_with_resinum_on_atoms(self):
        a = SADI(self.shx, 'SADI C1_1 C2_1 C2_2 C3_2 C3_2 C4_2'.split())
        self.assertEqual(['C1_1', 'C2_1', 'C2_2', 'C3_2', 'C3_2', 'C4_2'], a.atoms)
        self.assertEqual('', a.residue_class)
        self.assertEqual([0], a.residue_number)

    def test_sadi_normal_with_name(self):
        a = SADI(self.shx, 'SADI_CCF3 C1 C2 C2 C3 C3 C4'.split())
        self.assertEqual(['C1', 'C2', 'C2', 'C3', 'C3', 'C4'], a.atoms)

    def test_sadi_with_name_2(self):
        a = SADI(self.shx, 'SADI_CCF3 C1 C2 C2 C3 C3 C4'.split())
        self.assertEqual('CCF3', a.residue_class)
        self.assertEqual([0], a.residue_number)

    def test_sadi_with_name_3(self):
        a = SADI(self.shx, 'SADI_CCF3 C1 C2 C2 C3 C3 C4'.split())
        self.assertEqual(0.02, a.s)

    def test_sadi_with_name_4(self):
        a = SADI(self.shx, 'SADI_CCF3 0.001 C1 C2 C2 C3 C3 C4'.split())
        self.assertEqual(0.001, a.s)


class TestwithDEFS(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_sadi_with_s(self):
        a = SADI(self.shx, 'SADI_CCF3 0.001 C1 C2 C2 C3 C3 C4'.split())
        self.assertEqual(0.001, a.s)

    def test_sadi_without_s(self):
        a = SADI(self.shx, 'SADI C1 C2 C2 C3 C3 C4'.split())
        self.assertEqual(0.0234, a.s)

    def test_DFIX_without_s(self):
        a = DFIX(self.shx, 'DFIX 1.45 C1 C2 C2 C3 C3 C4'.split())
        self.assertEqual(1.45, a.d)
        self.assertEqual(0.0234, a.s)
        self.assertEqual(0.0234, self.shx.defs.sd)

    def test_DFIX_with_s(self):
        a = DFIX(self.shx, 'DFIX 1.45 0.022 C1 C2 C2 C3 C3 C4'.split())
        self.assertEqual(1.45, a.d)
        self.assertEqual(0.022, a.s)


class TestPLAN(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_plan__repr__(self):
        self.assertEqual('PLAN 20', self.shx.plan.__repr__())

    def test_set_plan_value(self):
        self.shx.plan.set('PLAN 23')
        self.assertEqual('PLAN 23', self.shx.plan.__repr__())

    def test_set_plan_all_values(self):
        self.shx.plan.set('PLAN 23 1.34 1.1')
        self.assertEqual('PLAN 23 1.34 1.1', self.shx.plan.__repr__())


class TestAFIX(TestCase):

    def test_afix_n(self):
        a = AFIX(Shelxfile(), 'AFIX 33'.split())
        self.assertEqual(33, a.mn)
        self.assertEqual(None, a.d)
        self.assertEqual(10.08, a.U)
        self.assertEqual(11.0, a.sof)
        self.assertEqual(True, a.__bool__())

    def test_afix_zero(self):
        a = AFIX(Shelxfile(), 'AFIX 0'.split())
        self.assertEqual(0, a.mn)
        self.assertEqual(False, a.__bool__())

    def test_afix_empty(self):
        a = AFIX(Shelxfile(), 'AFIX'.split())
        self.assertEqual(None, a.mn)
        self.assertEqual(False, a.__bool__())
