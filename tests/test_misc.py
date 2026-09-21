from unittest import TestCase

from shelxfile.misc.misc import range_resolver, wrap_line, multiline_test, chunks, find_line, eigenvals


class Test(TestCase):
    def test_range_resolver_1(self):
        r = "C2 > C5".split()
        atlist = 'C1 C2 C3 C4 C5'.split()
        self.assertEqual(['C2', 'C3', 'C4', 'C5'], range_resolver(r, atlist))

    def test_range_resolver_2(self):
        r = "C2_2 > C5_2".split()
        atlist = 'C1_1 C1_2 C2_2 C3_2 C4_2 C5_2'.split()
        self.assertEqual(['C2_2', 'C3_2', 'C4_2', 'C5_2'], range_resolver(r, atlist))

    def test_range_resolver_3(self):
        r = "C2_1 > C5_1".split()
        atlist = 'C1_1 C1_2 C2_2 C3_2 C4_2 C5_2'.split()
        with self.assertRaises(ValueError):
            range_resolver(r, atlist)

    def test_wrap_line(self):
        self.assertEqual('This is a really long line with over 79 characters. Shelxl wants it to be =\n   wrapped.',
                         wrap_line(
                             "This is a really long line with over 79 characters. Shelxl wants it to be wrapped."))

    def test_wrap_line_is_stable(self):
        """A wrapped line must survive being read back and wrapped again.

        The continuation marker used to be appended without stripping the
        chunk first, emitting '  ='. Re-reading collapses that run of
        spaces, so the file changed on every round-trip even though
        nothing had been edited.
        """
        long_line = ("REM wR2 = 0.627712, GooF = S = 5.48377, "
                     "Restrained GooF = 103.69011 for all data and more text")
        once = wrap_line(long_line)
        rejoined = ' '.join(once.replace('=\n', '').split())
        self.assertEqual(once, wrap_line(rejoined))

    def test_wrap_line_leaves_blank_lines_alone(self):
        """A blank line is a legal comment and must not gain a '='.

        Wrapping one would append a continuation marker, and SHELXL would
        then glue the following instruction onto the empty line.
        """
        blank = ' ' * 100
        self.assertEqual(blank, wrap_line(blank))

    def test_wrap_line_ignores_trailing_padding(self):
        """A line padded to a fixed column width is not a long line.

        Counting the padding wrapped the instruction and left the
        remainder as a whitespace continuation, so the '=' it gained made
        the next line a continuation on the following read.
        """
        padded = 'TITL p-1 in P-1'.ljust(80)
        self.assertEqual('TITL p-1 in P-1', wrap_line(padded))

    def test_multiline_test1(self):
        line = 'C1    1    0.278062    0.552051    0.832431    11.00000    0.02895    0.02285 ='
        self.assertEqual(True, multiline_test(line))

    def test_multiline_test2(self):
        line = 'C1    1    0.278062    0.552051    0.832431    11.00000    0.05 '
        self.assertEqual(False, multiline_test(line))

    def test_chunks(self):
        alist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e', 'f']
        self.assertEqual([[1, 2, 3, 4, 5], [6, 7, 8, 9, 0], ['a', 'b', 'c', 'd', 'e'], ['f']], chunks(alist, 5))
        self.assertEqual([[1], [2], [3], [4], [5], [6], [7], [8], [9], [0], ['a'], ['b'], ['c'], ['d'], ['e'], ['f']],
                         chunks(alist, 1))
        self.assertEqual([[1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e', 'f']], chunks(alist, 50))


class Testfind_line(TestCase):
    def test_find_line_found_something(self):
        inp = ['Hallo blub', 'foo bar blub', '123', '1 blub 2 3 4']
        self.assertEqual(0, find_line(inp, '.*blub.*'))

    def test_dont_find_something(self):
        inp = [['foo'], ['bar']]

        with self.assertRaises(TypeError):
            find_line(inp, '.*blub.*')


class TestEigenvalues(TestCase):
    def test_eigenvalues_1(self):
        # [0.08422976 0.15727835 0.20849189]
        matrix = [
            [0.1, 0.02, 0.03],
            [0.02, 0.2, 0.01],
            [0.03, 0.01, 0.15]
        ]
        eigenvalues = eigenvals(matrix)
        self.assertAlmostEqual(0.08422976, eigenvalues[0])
        self.assertAlmostEqual(0.15727835, eigenvalues[1])
        self.assertAlmostEqual(0.20849189, eigenvalues[2])
