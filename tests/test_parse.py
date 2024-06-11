import unittest

from shelxfile import Shelxfile


class ParseTestCase(unittest.TestCase):
    def test_parse_hklf_after_uncloded_resi(self):
        shx = Shelxfile(debug=True)
        shx.read_file(r'tests/resources/unclosed_resi.res')
        self.assertEqual('HKLF 4 1  1 0 0 0 1 0 0 0 1  1 0', shx.hklf.__repr__())  # add assertion here

    def test_parse_only_isotropic_after_solution(self):
        shx = Shelxfile(debug=True)
        shx.read_file(r'tests/resources/p21c_a_only_isotropic.res')
        self.assertEqual(146, len(shx.atoms))
        self.assertEqual({}, shx.fvars.fvars_used())


if __name__ == '__main__':
    unittest.main()
