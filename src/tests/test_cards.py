from unittest import TestCase

from src.shelxfile.cards import RESI
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
        self.assertEqual(('TOL', 1, None, None),(r.residue_class, r.residue_number, r.chain_id, r.alias))

    def test_resi_chainid_tol(self):
        r = RESI(None, 'RESI A:100 TOL'.split())
        self.assertEqual(('TOL', 100, 'A', None),(r.residue_class, r.residue_number, r.chain_id, r.alias))

    def test_resi__negative_chain_tol(self):
        r = RESI(None, 'RESI -10 TOL'.split())
        self.assertEqual(('TOL', -10, None, None), (r.residue_class, r.residue_number, r.chain_id, r.alias))

    def test_resi_b_chain_tol(self):
        r = RESI(None, 'RESI b:-10 TOL'.split())
        self.assertEqual(('TOL', -10, 'b', None),(r.residue_class, r.residue_number, r.chain_id, r.alias))

