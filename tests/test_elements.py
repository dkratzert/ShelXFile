from unittest import TestCase

from shelxfile.elements import get_radius, get_radius_from_element, get_atomic_number, get_element, get_atomlabel


class Test(TestCase):
    def test_get_radius(self):
        self.assertEqual(0.77, get_radius(6))

    def test_get_radius_from_element(self):
        self.assertEqual(0.72, get_radius_from_element('F'))

    def test_get_atomic_number(self):
        self.assertEqual(9, get_atomic_number('F'))

    def test_get_element(self):
        self.assertEqual('N', get_element(7))

    def test_get_atomlabel(self):
        self.assertEqual('C', get_atomlabel('C12'))
