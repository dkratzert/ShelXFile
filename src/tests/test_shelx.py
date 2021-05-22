from pathlib import Path
from unittest import TestCase

from src.shelxfile.shelx import Shelxfile


class TestShelxFileIsAtom(TestCase):

    def test_is_atom_long_line_true(self):
        """
        """
        is_an_atom = Shelxfile.is_atom(atomline='O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        self.assertEqual(True, is_an_atom)

    def test_is_no_atoms_because_type_is_missing(self):
        is_an_atom = Shelxfile.is_atom(atomline='O1    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        self.assertEqual(False, is_an_atom)

    def test_is_no_atom_because_coordinate_is_missing(self):
        is_an_atom = Shelxfile.is_atom(atomline='O1  4  0.120080    0.494426  11.00000   0.01445 ...')
        self.assertEqual(False, is_an_atom)

    def test_is_no_atom_because_its_a_command(self):
        is_atom = Shelxfile.is_atom("AFIX 123")
        self.assertEqual(False, is_atom)

    def test_is_no_atom_because_its_a_command2(self):
        is_atom = Shelxfile.is_atom("AFIX")
        self.assertEqual(False, is_atom)

    def test_is_an_atom_short(self):
        is_atom = Shelxfile.is_atom('O1    3    0.120080   0.336659   0.494426')
        self.assertEqual(True, is_atom)


class TestShelxfileElementToSfac(TestCase):
    """
    SFAC C H O F Al Ga
    """

    def setUp(self) -> None:
        self.shx = Shelxfile(Path('resources/p21c.res').resolve())

    def test_elem2sfac_oxygen(self):
        self.assertEqual(3, self.shx.elem2sfac('O'))

    def test_elem2sfac_carbon(self):
        self.assertEqual(1, self.shx.elem2sfac('c'))

    def test_elem2sfac_Argon(self):
        self.assertEqual(None, self.shx.elem2sfac('Ar'))


class TestShelxfile(TestCase):

    def setUp(self) -> None:
        self.shx = Shelxfile(Path('resources/p21c.res').resolve())

    def test_sfac2elem_C(self):
        self.assertEqual('C', self.shx.sfac2elem(1))

    def test_sfac2elem_H(self):
        self.assertEqual('H', self.shx.sfac2elem(2))

    def test_sfac2elem_O(self):
        self.assertEqual('O', self.shx.sfac2elem(3))

    def test_sfac2elem_Al(self):
        self.assertEqual('Al', self.shx.sfac2elem(5))

    def test_sfac2elem_not_existent(self):
        self.assertEqual('', self.shx.sfac2elem(8))

    def test_sfac2elem_zero(self):
        self.assertEqual('', self.shx.sfac2elem(0))

    def test_sum_formula(self):
        self.assertEqual('C0.25 H0.5 O0.75 F1 AL1.25 GA1.5', self.shx.sum_formula)
