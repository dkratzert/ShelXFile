from unittest import TestCase

from src.shelxfile.dsrmath import vol_unitcell, distance, dice_coefficient, levenshtein, dice_coefficient2, \
    SymmetryElement, Matrix, Array, mean, median, std_dev, nalimov_test, id_generator, atomic_distance, \
    almost_equal
from src.shelxfile.misc import flatten, frac_to_cart, cart_to_frac, determinante, subtract_vect
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


class TestMatrix(TestCase):
    def test_det(self):
        m1 = Matrix([[2, 0, 0], [0, 2, 0], [0, 0, 2]])
        self.assertEqual(8, m1.det)

    def test_zero(self):
        self.assertEqual("| 0.0000  0.0000  0.0000|\n"
                         "| 0.0000  0.0000  0.0000|\n"
                         "| 0.0000  0.0000  0.0000|\n"
                         "| 0.0000  0.0000  0.0000|\n"
                         "| 0.0000  0.0000  0.0000|\n", Matrix.zero(5, 3).__repr__())

    def test_equal(self):
        m1 = Matrix([(1., 2., 3.), (1., 2., 3.), (1., 2., 3.)])
        m2 = Matrix([(1, 2, 3), (1, 2, 3), (1, 2, 3)])
        self.assertEqual(True, m1 == m2)

    def test_equal2(self):
        m1 = Matrix([(1, 2, 3), (1, 2, 3), (1, 2, 3)])
        m2 = Matrix([(1, 2, 3), (3, 2, 3), (1, 2, 3)])
        self.assertEqual(False, m1 == m2)

    def test_subtract_matrix_from_matrix(self):
        self.assertEqual([[2, 0, -2], [2, 0, -2], [2, 0, 0]],
                         Matrix([[3, 2, 1], [3, 2, 1], [3, 2, 3]]) - Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]]))

    def test_transposed(self):
        m = Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
        self.assertEqual([(1, 1, 1), (2, 2, 2), (3, 3, 3)], m.transposed.values)

    def transpose_alt(self):
        m = Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
        self.assertEqual([[1, 1, 1], [2, 2, 2], [3, 3, 3]], m.transpose_alt().values)

    def test_cholesky(self):
        m = Matrix([[25, 15, -5], [15, 18, 0], [-5, 0, 11]])
        self.assertEqual(
            "| 5.0000  0.0000  0.0000|\n"
            "| 3.0000  3.0000  0.0000|\n"
            "|-1.0000  1.0000  3.0000|\n", m.cholesky().__repr__())

    def test_inversed(self):
        self.assertEqual("|-0.8125  0.1250  0.1875|\n"
                         "| 0.1250 -0.2500  0.1250|\n"
                         "| 0.5208  0.1250 -0.1458|\n",
                         Matrix([[1, 2, 3], [4, 1, 6], [7, 8, 9]]).inversed.__repr__())

    def test_foo(self):
        m = Matrix([[1, 2, 300], [4.1, 4.2, 4.3], [5, 6, 7]])
        x = m + m
        self.assertEqual("| 2.0000  4.0000 600.0000|\n"
                         "| 8.2000  8.4000  8.6000|\n"
                         "|10.0000 12.0000 14.0000|\n", x.__repr__())
        m *= 3
        self.assertEqual("| 3.0000  6.0000 900.0000|\n"
                         "|12.3000 12.6000 12.9000|\n"
                         "|15.0000 18.0000 21.0000|\n", m.__repr__())

    def test_getitem(self):
        m = Matrix([[2., 2., 3.], [1., 2.2, 3.], [1., 2., 3.]])
        self.assertEqual(2.2, m[1, 1])
        self.assertEqual(2.2, m[1][1])

    def test_matrix_add_matrix(self):
        m1 = Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
        m2 = Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
        t1 = m1 + m2
        self.assertEqual(
            "| 2.0000  2.0000  2.0000|\n"
            "| 2.0000  2.0000  2.0000|\n"
            "| 2.0000  2.0000  2.0000|\n", t1.__repr__())
        t2 = m1 + 0.5
        self.assertEqual(
            "| 1.5000  1.5000  1.5000|\n"
            "| 1.5000  1.5000  1.5000|\n"
            "| 1.5000  1.5000  1.5000|\n", t2.__repr__())

    def test_matrix_multiply1(self):
        m = Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]]) * 2
        self.assertEqual("| 2.0000  2.0000  2.0000|\n"
                         "| 2.0000  2.0000  2.0000|\n"
                         "| 2.0000  2.0000  2.0000|\n", m.__repr__())

    def test_matrix_multiply2(self):
        m = Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
        m2 = Matrix([[2, 2, 2], [0.5, 0.5, 0.5], [2, 2, 1]])
        x = m * m2
        self.assertEqual("| 6.0000  1.5000  5.0000|\n"
                         "| 6.0000  1.5000  5.0000|\n"
                         "| 6.0000  1.5000  5.0000|\n", x.__repr__())
        self.assertEqual("| 1.0000  1.0000  1.0000|\n"
                         "| 1.0000  1.0000  1.0000|\n"
                         "| 1.0000  1.0000  1.0000|\n", m.__repr__())

    def test_matrix_multiply_array(self):
        m = Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
        self.assertEqual(Array([6, 6, 6]), m * Array([2, 2, 2]))

    def test_matrix_multiply4(self):
        m = Matrix([(0, 1, 0), (-1, -1, 0), (0, 0, 1)])
        self.assertEqual(Array([-0.666667, -0.333334, 0.45191]), Array([0.333333, 0.666667, 0.45191]) * m)


class TestArray(TestCase):

    def test_add_array_to_array_inplace(self):
        a = Array([1, 2, 3, 4.1])
        a += a
        self.assertEqual(Array([2, 4, 6, 8.2]), a)

    def test_add_constant_to_array(self):
        a = Array([1, 2, 3, 4.1])
        self.assertEqual(Array([5, 6, 7, 8.1]), a + 4)

    def test_add_array_to_array(self):
        a = Array([1, 2, 3, 4.1])
        self.assertEqual(Array([2, 4, 6, 8.2]), a + a)

    def test_getitem(self):
        a = Array([2, 4, 6, 8.2])
        self.assertEqual(4, a[1])

    def test_get_item(self):
        self.assertEqual(2, Array([1, 2, 3])[1])

    def test_multiply_constant_inplace(self):
        a = Array([1, 2, 3, 4.1])
        a *= 3
        self.assertEqual(Array([3, 6, 9, 12.299999999999999]), a)

    def test_dot_multiply(self):
        a = Array([1, 2, 3, 4.1])
        self.assertEqual(30.81, a.dot(a))

    def test_multiply_array_with_array(self):
        a = Array([1, 2, 3, 4.1])
        self.assertEqual(30.81, a * a)

    def test_norm(self):
        a = Array([1, 2, 3, 4])
        self.assertEqual(30, a.norm())

    def test_normalized(self):
        a = Array([2, 2, 1])
        self.assertEqual(3.0, a.normalized())

    def test_zero(self):
        self.assertEqual(Array([0.0, 0.0, 0.0, 0.0, 0.0]), Array.zero(5))

    def test_floor(self):
        self.assertEqual(Array([3, 2, 0]), Array([3.634, 2, 0.345]).floor)

    def test_cross(self):
        a = Array([1, 2, 3])
        b = Array([-7, 8, 9])
        self.assertEqual(Array([-6, -30, 22]), a.cross(b))

    def test_angle(self):
        a = Array([1, 0, 1])
        b = Array([1, 0, 0])
        self.assertEqual(45.0, a.angle(b))

    def test_angle_2(self):
        va = Array([0.03562, 0.14298, 0.24008]) - Array([0.04402, 0.16614, 0.22275])
        vb = Array([0.07078, 0.17382, 0.22106]) - Array([0.04402, 0.16614, 0.22275])
        self.assertEqual(120.9401, round(va.angle(vb), 4))

    def test_add(self):
        a = Array([1, 2, 3])
        b = Array([1, 1, 1])
        self.assertEqual(Array([2, 3, 4]), a + b)

    def test_setitem(self):
        a = Array([0, 0, 0])
        a[1] = 5
        self.assertEqual(Array([0, 5, 0]), a)

    def test_multiply(self):
        a1 = Array([1, 2, 3])
        a2 = Array([1, 2, 3])
        self.assertEqual(14, a1 * a2)

    def test_subtract(self):
        a = Array([1, 2, 3])
        b = Array([1, 1, 1])
        self.assertEqual(Array([0, 1, 2]), a - b)
        self.assertEqual(Array([0, -1, -2]), b - a)

    def test_equality(self):
        a1 = Array([1, 2, 3, 4])
        a2 = Array([1, 2, 3.0, 4.0])
        self.assertEqual(True, a1 == a2)

    def test_equality_false(self):
        a1 = Array([1, 2, 3, 4])
        a2 = Array([2, 2, 3.0, 4.0])
        self.assertEqual(False, a1 == a2)


class TestMisc(TestCase):
    def test_mean(self):
        self.assertEqual(2.5, mean([1, 2, 3, 4, 1, 2, 3, 4]))

    def test_median1(self):
        self.assertEqual(2, median([2]))

    def test_median2(self):
        self.assertEqual(2.5, median([1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]))

    def test_median3(self):
        self.assertEqual(3, median([1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4.1, 1000000]))

    def test_median4(self):
        with self.assertRaises(ValueError):
            median([])

    def test_std_dev1(self):
        l1 = [1.334, 1.322, 1.345, 1.451, 1.000, 1.434, 1.321, 1.322]
        self.assertEqual(0.13797871, round(std_dev(l1), 8))

    def test_std_dev2(self):
        l2 = [1.234, 1.222, 1.345, 1.451, 2.500, 1.234, 1.321, 1.222]
        self.assertEqual(0.43536797, round(std_dev(l2), 8))

    def test_std_dev3(self):
        l1 = [1.334, 1.322, 1.345, 1.451, 1.000, 1.434, 1.321, 1.322]
        self.assertEqual(1.328, median(l1))
        self.assertEqual(1.316125, mean(l1))

    def test_nalimov_test(self):
        data = [1.120, 1.234, 1.224, 1.469, 1.145, 1.222, 1.123, 1.223, 1.2654, 1.221, 1.215]
        self.assertEqual([3], nalimov_test(data))

    def test_flatten_list(self):
        self.assertEqual(['wer', 234, 'brdt5', 'dfg', 21, 34, 5, 'fhg', 4],
                         flatten([['wer', 234, 'brdt5'], ['dfg'], [[21, 34, 5], ['fhg', 4]]]))

    def test_id_generator(self):
        self.assertEqual('a', id_generator(1, 'a'))

    def test_atomic_distance(self):
        cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
        coord1 = [-0.186843, 0.282708, 0.526803]
        coord2 = [-0.155278, 0.264593, 0.600644]
        self.assertEqual(1.5729229943265979, atomic_distance(coord1, coord2, cell))

    def test_determinante(self):
        m1 = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
        self.assertEqual(8, determinante(m1))

    def test_almost_equal(self):
        self.assertEqual(True, almost_equal(1.0001, 1.0005))
        self.assertEqual(False, almost_equal(1.1, 1.0005))
        self.assertEqual(False, almost_equal(2, 1))

    def test_fractional_to_cartesian(self):
        cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
        coord1 = [-0.186843, 0.282708, 0.526803]
        self.assertEqual([-2.741505423999065, 5.909586678000002, 10.775200700893734], frac_to_cart(coord1, cell))

    def test_cart_to_frac(self):
        cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
        coords = [-2.74150542399906, 5.909586678, 10.7752007008937]
        self.assertEqual((-0.1868429999999998, 0.28270799999999996, 0.5268029999999984), cart_to_frac(coords, cell))
