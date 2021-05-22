from unittest import TestCase

from src.shelxfile.shelx import Shelxfile


class TestAtoms(TestCase):

    def setUp(self) -> None:
        self.shx = Shelxfile('resources/p21c.res')

    def test_number(self):
        self.assertEqual(148, self.shx.atoms.number)

    def test_has_atom_al1(self):
        self.assertEqual(True, self.shx.atoms.has_atom('Al1'))

    def test_has_atom_al1_0(self):
        self.assertEqual(True, self.shx.atoms.has_atom('Al1_0'))

    def test_has_atom_al2_0(self):
        self.assertEqual(False, self.shx.atoms.has_atom('Al2_0'))

    def test_get_atom_by_name(self):
        self.assertEqual("Atom ID: 73", self.shx.atoms.get_atom_by_name('Al1').__repr__())
