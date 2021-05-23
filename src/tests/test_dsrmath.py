from unittest import TestCase

from src.shelxfile.dsrmath import vol_unitcell, distance, dice_coefficient, levenshtein, dice_coefficient2, \
    subtract_vect, SymmetryElement
from src.shelxfile.shelx import Shelxfile


class Testdsrmath(TestCase):
    def test_vol_unitcell(self):
        volume = vol_unitcell(2, 2, 2, 90, 90, 90)
        self.assertEqual(8.0, volume)

    def test_distance_1(self):
        d = distance(1, 1, 1, 2, 2, 2, 4)
        self.assertEqual(1.7321, d)

    def test_distance_2(self):
        d = distance(1, 0, 0, 2, 0, 0, 4)
        self.assertEqual(1.0, d)

    def test_levenshtein(self):
        l = levenshtein('hallo', 'holla')
        self.assertEqual(2, l)

    def test_dice(self):
        d = dice_coefficient('hallo', 'holla')
        self.assertEqual(0.25, d)

    def test_dice2(self):
        self.assertEqual(0.75, dice_coefficient2('hallo', 'holla'))

    def test_dice3(self):
        self.assertEqual(0.6, dice_coefficient2('Banze', 'Benzene'))

    def test_dice4(self):
        self.assertEqual(0.333333, dice_coefficient2('halo', 'Haaallo'))

    def test_dice5(self):
        self.assertEqual(0.2, dice_coefficient2('hallo', 'Haaallo'))

    def test_dice6(self):
        self.assertEqual(0.0, dice_coefficient2('hallo', 'Hallo'))

    def test_dice_7(self):
        self.assertEqual(1.0, dice_coefficient2('aaa', 'BBBBB'))

    def test_dice_8(self):
        self.assertEqual(1.0, dice_coefficient2('', ''))

    def test_subtract_vect(self):
        self.assertEqual((-2, 0, 1), subtract_vect([1, 2, 3], [3, 2, 2]))


class TestSymmetryElement(TestCase):
    def setUp(self) -> None:
        self.shx = Shelxfile('resources/p21c.res')

    def test_to_shelxl(self):
        self.assertEqual('[+X, +Y, +Z, -X, 0.5+Y, 0.5-Z, +X, -0.5-Y, -0.5+Z, -X, -Y, -Z]',
                         self.shx.symmcards._symmcards.__repr__())

    def test_repr(self):
        self.assertEqual(SymmetryElement(['-X', '-Y', '-Z']), self.shx.symmcards[3])

    def test_string(self):
        self.assertEqual("|-1  0  0|   | 0.0|\n"
                         "| 0  1  0| + | 0.5|\n"
                         "| 0  0 -1|   | 0.5|\n", self.shx.symmcards[1].__str__())

    def test_equals_false(self):
        self.assertEqual(False, self.shx.symmcards[0] == self.shx.symmcards[1])

    def test_equals_True(self):
        self.assertEqual(True, self.shx.symmcards[1] == self.shx.symmcards[1])


    def test_s12_equals(self):
        s1 = SymmetryElement(['0.5', '0.5', '0.5'])
        s2 = SymmetryElement(['0.5', '0.5', '0.5'])
        self.assertEqual(True, s1 == s2)

    def test_s12_equals2(self):
        s1 = SymmetryElement(['1.5', '1.5', '1.5'])
        s2 = SymmetryElement(['0.5', '0.5', '0.5'])
        self.assertEqual(True, s1 == s2)

    def test_s34_not_equals(self):
        s3 = SymmetryElement(['1', '0.5', '0.5'])
        s4 = SymmetryElement(['0.5', '0.5', '0.5'])
        self.assertEqual(False, s3 == s4)