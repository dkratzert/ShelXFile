import unittest

from shelxfile import Shelxfile


class TestRestraintsWarnings(unittest.TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile(debug=True)

    def test_class_and_number_on_restr_and_atoms_but_wrong(self):
        self.shx.read_file('tests/resources/restraint_tests/class_and_number_on_restr_and_atoms_but_wrong.res')
        result = self.shx._assign_atoms_to_restraints()
        expected = ['*** Unknown atoms in restraint: SADI_CCF3 N1 N2_3 C1_0 C2, line 13 ***',
                    '*** Atom list has no --> C2_1, N2_3 ***']
        self.assertListEqual(expected, result)

    @unittest.skip('Unfinished')
    def test_wildcard_on_atom_names(self):
        self.shx.read_file('tests/resources/restraint_tests/wildcard_on_atom.res')
        result = self.shx._assign_atoms_to_restraints()
        expected = ['*** Unknown atoms in restraint: RIGU $N_* $C_*, line 13 ***',
                    '*** Atom list has no --> $C_*, $N_* ***']
        self.assertListEqual(expected, result)

    def test_class_and_numbers(self):
        self.shx.read_file('tests/resources/restraint_tests/class_and_numbers.res')
        result = self.shx._assign_atoms_to_restraints()
        expected = []
        self.assertListEqual(expected, result)

    def test_class_but_missing_atoms(self):
        self.shx.read_file('tests/resources/restraint_tests/class_but_missing_atoms.res')
        result = self.shx._assign_atoms_to_restraints()
        expected = ['*** Unknown atoms in restraint: SADI_CCF3 N1 N2 C1 C2, line 13 ***',
                    '*** Atom list has no --> C1_1, C2_1 ***']
        self.assertListEqual(expected, result)

    def test_no_class_or_number(self):
        self.shx.read_file('tests/resources/restraint_tests/no_class_or_number.res')
        result = self.shx._assign_atoms_to_restraints()
        expected = ['*** Unknown atoms in restraint: SADI N1 N2 C1 C2, line 13 ***',
                    '*** Atom list has no --> N1, N2 ***']
        self.assertListEqual(expected, result)

    def test_number_but_missing_atoms(self):
        self.shx.read_file('tests/resources/restraint_tests/number_but_missing_atoms.res')
        result = self.shx._assign_atoms_to_restraints()
        expected = ['*** Unknown atoms in restraint: SADI_1 N1 N2 C1 C2, line 13 ***',
                    '*** Atom list has no --> C1_1, C2_1 ***']
        self.assertListEqual(expected, result)

    def test_number_on_restr_and_atoms(self):
        self.shx.read_file('tests/resources/restraint_tests/number_on_restr_and_atoms.res')
        result = self.shx._assign_atoms_to_restraints()
        expected = []
        self.assertListEqual(expected, result)

    def test_number_on_restr_and_atoms_but_wrong(self):
        self.shx.read_file('tests/resources/restraint_tests/number_on_restr_and_atoms_but_wrong.res')
        result = self.shx._assign_atoms_to_restraints()
        expected = ['*** Unknown atoms in restraint: SADI_2 N1 N2_3 C1_0 C2_0, line 13 ***',
                    '*** Atom list has no --> N1_2, N2_3 ***']
        self.assertListEqual(expected, result)

    def test_only_number_0_on_atoms(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/restraint_tests/only_number_0_on_atoms.res')
        result = shx._assign_atoms_to_restraints()
        expected = ['*** Unknown atoms in restraint: SADI N1 N2 C1_0 C2_0, line 13 ***',
                    '*** Atom list has no --> N1, N2 ***']
        self.assertListEqual(expected, result)
