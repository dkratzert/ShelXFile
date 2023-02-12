from unittest import TestCase

from shelxfile import Shelxfile
from shelxfile.misc.dsrmath import Array


class TestOrthogonalMatrix(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile(debug=True)
        self.shx.read_file('tests/resources/jkd77.res')

    def test_frac_to_cart(self):
        self.assertListEqual(list(self.shx.atoms.all_atoms[2].cart_coords),
                             list(self.shx.cell.o * Array(self.shx.atoms.all_atoms[2].frac_coords)))

    def test_cart_to_frac(self):
        self.assertListEqual(list(self.shx.atoms.all_atoms[2].frac_coords),
                             [round(x, 8) for x in
                              self.shx.cell.o.m.inversed * Array(self.shx.atoms.all_atoms[2].cart_coords)])

    def test_frac_to_cart_all(self):
        for atom in self.shx.atoms.all_atoms:
            self.assertListEqual(list(atom.cart_coords), list(self.shx.cell.o * Array(atom.frac_coords)))

    def test_cart_to_frac_all(self):
        for atom in self.shx.atoms.all_atoms:
            self.assertListEqual(list(atom.frac_coords),
                                 [round(x, 9) for x in list(self.shx.cell.o.m.inversed * Array(atom.cart_coords))])
