from unittest import TestCase

from shelxfile.refine.refine import get_xl_version_string


class TestRefine(TestCase):
    def test_get_xl_version_string_with_no_path(self):
        self.assertEqual('', get_xl_version_string(''))

    def test_get_xl_version_string_with_real_path(self):
        self.assertEqual('2018/3', get_xl_version_string('/usr/local/bin/shelxl'))


    #def test_reload(self):
    #    self.assertEqual('', )