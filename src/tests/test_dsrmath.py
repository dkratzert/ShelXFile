from unittest import TestCase

from src.shelxfile.dsrmath import vol_unitcell, distance, dice_coefficient, levenshtein, dice_coefficient2, \
    subtract_vect


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
