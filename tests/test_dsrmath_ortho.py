from unittest import TestCase

from shelxfile import Shelxfile
from shelxfile.misc.dsrmath import OrthogonalMatrix, Array


class TestOrthogonalMatrix(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/jkd77.res')

    def test_frac_to_cart(self):
        o = OrthogonalMatrix(*self.shx.cell)
        self.assertListEqual(list(self.shx.atoms.all_atoms[2].cart_coords),
                             list(o * Array(self.shx.atoms.all_atoms[2].frac_coords)))

    def test_cart_to_frac(self):
        o = OrthogonalMatrix(*self.shx.cell)
        self.assertListEqual(list(self.shx.atoms.all_atoms[2].frac_coords),
                             [round(x, 8) for x in o.m.inversed * Array(self.shx.atoms.all_atoms[2].cart_coords)])
