from unittest import TestCase

from shelxfile.misc.elements import get_radius, get_radius_from_element, get_atomic_number, get_element, get_atomlabel


class TestElements(TestCase):
    def test_get_radius_6(self):
        self.assertEqual(0.77, get_radius(6))

    def test_get_radius_7(self):
        self.assertEqual(0.74, get_radius(7))

    def test_get_radius_8(self):
        self.assertEqual(0.71, get_radius(8))

    def test_get_radius_1(self):
        self.assertEqual(0.45, get_radius(1))

    def test_get_radius_2(self):
        self.assertEqual(0.45, get_radius(2))

    def test_get_radius_3(self):
        self.assertEqual(1.23, get_radius(3))

    def test_get_radius_from_element_F(self):
        self.assertEqual(0.72, get_radius_from_element('F'))

    def test_get_radius_from_element_C(self):
        self.assertEqual(0.77, get_radius_from_element('C'))

    def test_get_radius_from_element_O(self):
        self.assertEqual(0.73, get_radius_from_element('O'))

    def test_get_radius_from_element_H(self):
        self.assertEqual(0.45, get_radius_from_element('H'))

    def test_get_radius_from_element_N(self):
        self.assertEqual(0.75, get_radius_from_element('N'))

    def test_get_atomic_number(self):
        self.assertEqual(9, get_atomic_number('F'))

    def test_get_element(self):
        self.assertEqual('N', get_element(7))

    def test_get_atomlabel(self):
        self.assertEqual('C', get_atomlabel('C12'))
