from unittest import TestCase

from shelxfile import Shelxfile


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
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

    def test_elem2sfac_oxygen(self):
        self.assertEqual(3, self.shx.elem2sfac('O'))

    def test_elem2sfac_carbon(self):
        self.assertEqual(1, self.shx.elem2sfac('c'))

    def test_elem2sfac_Argon(self):
        self.assertEqual(None, self.shx.elem2sfac('Ar'))


class TestShelxfile(TestCase):

    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/p21c.res')

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

    def test_read_file_to_list(self):
        shx = Shelxfile()
        shx.read_file('tests/resources/p21c.res')
        self.assertEqual(['TITL p21c in P2(1)/c',
                          '    created by SHELXL-2018/3 at 16:18:25 on 03-May-2018',
                          'CELL 0.71073 10.5086 20.9035 20.5072 90 94.13 90',
                          'ZERR 4 0.0003 0.0005 0.0005 0 0.001 0',
                          'LATT 1'],
                         [str(x) for x in shx._reslist[:5]])

    def test_read_file(self):
        shx = Shelxfile()
        with self.assertRaises(FileNotFoundError):
            shx.read_file('tests/resources/foobar.res')


class TestShelxfileGoodModel(TestCase):

    def setUp(self) -> None:
        self.shx = Shelxfile()
        self.shx.read_file('tests/resources/model_finished/p21c.res')

    def test_atom_representation(self):
        a = self.shx._reslist[39]
        self.assertEqual(
            'C1     1   -0.187504    0.282782    0.527531    11.00000    0.01807    0.02353      0.01797    0.00008   -0.00179    0.00036',
            str(a))

    def test_atom_object(self):
        a = self.shx._reslist[39]
        self.assertEqual(
            'C1     1   -0.187504    0.282782    0.527531    11.00000    0.01807    0.02353      0.01797    0.00008   -0.00179    0.00036',
            a)
