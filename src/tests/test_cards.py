from unittest import TestCase

from src.shelxfile.cards import RESI
from src.shelxfile.refine.refine import ShelxlRefine
from src.shelxfile.shelx import Shelxfile


class TestCELL(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile('./resources/p21c.res')

    def test_volume(self):
        self.assertEqual(4493.0474, round(self.shx.cell.volume, 4))


class TestRESI(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile('./resources/p21c.res')

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
        shx = Shelxfile('./resources/p21c.res')
        self.assertEqual(2.0, shx.fmap.code)


"""
    >>> shx._reslist[12]
    ACTA 45
    >>> shx.acta
    ACTA 45
    >>> ref.remove_acta_card(shx.acta)
    >>> shx._reslist[12]
    SIZE 0.12 0.23 0.33
    >>> ref.restore_acta_card()
    >>> shx.index_of(shx.acta)
    8
    >>> shx._reslist[7:10]
    [UNIT 1  2  3  4  5  6, ACTA 45, 'LIST 4 ! automatically inserted. Change 6 to 4 for CHECKCIF!!']
"""


class TestACTA(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile('./resources/p21c.res')
        self.ref = ShelxlRefine(self.shx, './resources/p21c.res')

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
            ['UNIT 1 2 3 4 5 6', 'ACTA 45', 'LIST 4 ! automatically inserted. Change 6 to 4 for CHECKCIF!!'],
            [str(x) for x in self.shx._reslist[7:10]])


class TestLSCycles(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile('./resources/p21c.res')
        self.ref = ShelxlRefine(self.shx, './resources/p21c.res')

    def test_set_refine_cycles(self):
        self.shx.cycles.set_refine_cycles(44)
        self.assertEqual('L.S. 44', self.shx._reslist[self.shx.cycles.index].__repr__())


class TestSFACTable(TestCase):
    def test_parse_element_line(self):
        self.fail()
        """
        >>> from src.shelxfile import ShelXFile
        >>> shx = ShelXFile('./tests/p21c.res')
        >>> shx.sfac_table
        SFAC C  H  O  F  Al  Ga
        >>> shx.unit
        UNIT 1  2  3  4  5  6
        >>> shx.sfac_table.add_element('Au')
        >>> shx.sfac_table
        SFAC C  H  O  F  Al  Ga  Au
        >>> shx.unit
        UNIT 1  2  3  4  5  6  1
        """
