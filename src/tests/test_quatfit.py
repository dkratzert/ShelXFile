from unittest import TestCase

from src.shelxfile.fit.quatfit import matrix_minus_vect, matrix_plus_vect, transpose, centroid


class Test(TestCase):
    def test_matrix_minus_vect(self):
        source = [[-0.01453, 1.6659, 0.10966], [-0.00146, 0.26814, 0.06351], [-0.27813, -0.21605, 1.52795]]
        cent = [-0.09804, 0.57266333, 0.56704]

        self.assertEqual([[0.08351, 1.09323667, -0.45738], [0.09658, -0.30452333000000004, -0.50353],
                          [-0.18008999999999997, -0.78871333, 0.9609099999999999]], matrix_minus_vect(source, cent))
        self.assertEqual([[0, 0, 0], [1, 1, 1], [2, 2, 2]],
                         matrix_minus_vect([[1, 1, 1], [2, 2, 2], [3, 3, 3]], [1, 1, 1]))

    def test_matrix_plus_vect(self):
        source = [[-0.01453, 1.6659, 0.10966], [-0.00146, 0.26814, 0.06351], [-0.27813, -0.21605, 1.52795]]
        cent = [-0.09804, 0.57266333, 0.56704]
        self.assertEqual([[-0.11257, 2.23856333, 0.6767], [-0.0995, 0.84080333, 0.6305499999999999],
                          [-0.37617, 0.35661333000000006, 2.09499]], matrix_plus_vect(source, cent))
        self.assertEqual([[2, 2, 2], [3, 3, 3], [4, 4, 4]],
                         matrix_plus_vect([[1, 1, 1], [2, 2, 2], [3, 3, 3]], [1, 1, 1]))

    def test_transpose(self):
        m = [[1, 2, 3], [1, 2, 3], [1, 2, 3]]
        self.assertEqual([(1, 1, 1), (2, 2, 2), (3, 3, 3)], transpose(m))

    def test_centroid(self):
        self.assertEqual((-0.09804, 0.5726633333333333, 0.56704), centroid(
            [[-1.45300e-02, 1.66590e+00, 1.09660e-01], [-1.46000e-03, 2.68140e-01, 6.35100e-02],
             [-2.78130e-01, -2.16050e-01, 1.52795e+00]]))
