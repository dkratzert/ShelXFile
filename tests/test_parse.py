import unittest

from shelxfile import Shelxfile


class ParseTestCase(unittest.TestCase):
    def test_parse_hklf_after_uncloded_resi(self):
        shx = Shelxfile(debug=True)
        shx.read_file(r'tests/resources/unclosed_resi.res')
        self.assertEqual('HKLF 4 1  1 0 0 0 1 0 0 0 1  1 0', shx.hklf.__repr__())  # add assertion here

    def test_closing_at_hklf_keeps_atoms_in_their_resi_part_and_afix(self):
        shx = Shelxfile(debug=True)
        shx.read_string('TITL t\n'
                        'CELL 0.71 10 10 10 90 90 90\n'
                        'ZERR 1 0 0 0 0 0 0\n'
                        'LATT 1\n'
                        'SFAC C H\n'
                        'UNIT 1 1\n'
                        'FVAR 1\n'
                        'RESI TOL 1\n'
                        'PART 1\n'
                        'C1 1 0.1 0.1 0.1 11.0 0.05\n'
                        'AFIX 43\n'
                        'H1 2 0.2 0.1 0.1 11.0 -1.2\n'
                        'HKLF 4\n'
                        'END\n')
        carbon, hydrogen = list(shx.atoms)
        self.assertEqual('C1_1', carbon.fullname)
        self.assertEqual(1, carbon.part.n)
        self.assertEqual(43, hydrogen.afix.mn)
        self.assertEqual(['RESI TOL 1', 'PART 1', 'AFIX 43'],
                         [line for line in shx.dumps().splitlines() if line.startswith(('RESI', 'PART', 'AFIX'))])

    def test_parse_only_isotropic_after_solution(self):
        shx = Shelxfile(debug=True)
        shx.read_file(r'tests/resources/p21c_a_only_isotropic.res')
        self.assertEqual(146, len(shx.atoms))
        self.assertEqual({}, shx.fvars.fvars_used())


if __name__ == '__main__':
    unittest.main()
