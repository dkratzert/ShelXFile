from unittest import TestCase

from shelxfile.misc.misc import range_resolver, wrap_line, multiline_test, chunks, find_line


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
        self.assertEqual('This is a really long line with over 79 characters. Shelxl wants it to be  =\n   wrapped.',
                         wrap_line(
                             "This is a really long line with over 79 characters. Shelxl wants it to be wrapped."))

    def test_multiline_test1(self):
        line = 'C1    1    0.278062    0.552051    0.832431    11.00000    0.02895    0.02285 ='
        self.assertEqual(True, multiline_test(line))

    def test_multiline_test2(self):
        line = 'C1    1    0.278062    0.552051    0.832431    11.00000    0.05 '
        self.assertEqual(False, multiline_test(line))

    def test_chunks(self):
        l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e', 'f']
        self.assertEqual([[1, 2, 3, 4, 5], [6, 7, 8, 9, 0], ['a', 'b', 'c', 'd', 'e'], ['f']], chunks(l, 5))
        self.assertEqual([[1], [2], [3], [4], [5], [6], [7], [8], [9], [0], ['a'], ['b'], ['c'], ['d'], ['e'], ['f']],
                         chunks(l, 1))
        self.assertEqual([[1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e', 'f']], chunks(l, 50))


class Testfind_line(TestCase):
    def test_find_line_found_something(self):
        inp = ['Hallo blub', 'foo bar blub', '123', '1 blub 2 3 4']
        self.assertEqual(0, find_line(inp, '.*blub.*'))

    def test_dont_find_something(self):
        inp = [['foo'], ['bar']]

        with self.assertRaises(TypeError) as e:
            find_line(inp, '.*blub.*')
