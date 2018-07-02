# -*- encoding: utf-8 -*-
# möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#

import random
import string
from math import sqrt, radians, cos, sin


class Array(object):
    """
    MIT License

    Copyright (c) 2018 Jens Luebben

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    >>> a = Array([1, 2, 3, 4.1])
    >>> a+a
    Array([2, 4, 6, 8.2])
    >>> a+4
    Array([5, 6, 7, 8.1])
    >>> a[1]
    2
    >>> a*=3
    >>> a
    Array([3, 6, 9, 12.299999999999999])
    >>> a = Array([1, 2, 3, 4.1])
    >>> a.dot(a)
    30.81
    >>> a*a
    30.81
    """

    def __init__(self, values):
        self.values = values

    def __iter__(self):
        for v in self.values:
            yield v

    def __len__(self):
        return len(self.values)

    def __add__(self, other):
        if isinstance(other, Array):
            if not len(self) == len(other):
                raise ValueError('Arrays are not of equal length.')
            return Array([i + j for i, j in zip(self, other)])
        elif type(other) == float or type(other) == int:
            return Array([i + other for i in self.values])
        else:
            raise TypeError('Cannot add type Array to type {}.'.format(str(type(other))))

    def __imul__(self, other):
        if isinstance(other, int):
            self.values = [v * other for v in self.values]
            return self
        else:
            raise TypeError('Unsupported operation.')

    def __mul__(self, other):
        """
        a*b = axbx + ayby + azbz
        """
        return self.dot(other)

    def __repr__(self):
        return 'Array({})'.format(str(self.values))

    def __getitem__(self, val):
        #try:
        #    start, stop, step = val.start, val.stop, val.step
        #except AttributeError:
        #    pass
        #else:
        #    return self.values[val]
        #finally:
        return self.values[val]

    def dot(self, other):
        return sum([i * j for i, j in zip(self, other)])

    def cross(self, other):
        """
        Cross product of the Array
        M = |a| * |b| * sin((a, b) or
        M = (ay * bz - az * by)*i + (az * bx - ax bz)j + (ax * by - ay bx)k
        """
        pass


class Matrix(object):
    """
    MIT License

    Copyright (c) 2018 Jens Luebben

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    DK: Missing: inverse, eigenvalues

    >>> m = Matrix([[1, 2, 3], [4.1, 4.2, 4.3], [5, 6, 7]])

    #>>> m+m
    #Array([2, 4, 6, 8.2])
    #>>> m+4
    #Array([5, 6, 7, 8.1])
    >>> m*=3
    >>> m.values
    [[3, 6, 9], [12.299999999999999, 12.600000000000001, 12.899999999999999], [15, 18, 21]]
    >>> m = Matrix([[1, 2, 3], [4.1, 4.2, 4.3], [5, 6, 7]])
    >>> m.det(m)
    """

    def __init__(self, values):
        self.shape = (len(values[0]), len(values))
        self.values = values

    def __getitem__(self, val):
        if isinstance(val, tuple):
            if len(val) == 1:
                return self.values[val[0]]
            return self.values[val[1]][val[0]]

    def __str__(self):
        print(self.values)
        return '\n'.join([str(row) for row in self.values])

    def __imul__(self, other):
        if isinstance(other, int):
            self.values = [[v * other for v in row] for row in self.values]
            return self
        else:
            raise TypeError('Unsupported operation.')

    def __mul__(self, other):
        return self.dot(other)

    def __len__(self):
        return self.shape[1]

    @property
    def t(self):
        return self.transpose()

    def transpose(self):
        """
        maybe: return list(zip(*m))
        """
        rows = []
        for i in range(self.shape[0]):
            rows.append([r[i] for r in self.values])
        return Matrix(rows)

    def dot(self, other):
        newA = []
        for i, col in enumerate(self.transpose().values):
            s = sum([v * o for v, o in zip(col, other)])
            newA.append(s)
        return Array(newA)

    def det(self, a):
        """
        #>>> a = Matrix([[1, 2], [3, 4]])
        #>>> Matrix.det(a)
        """
        pass


class SymmetryElement(object):
    """
    Class representing a symmetry operation.
    """
    symm_ID = 1

    def __init__(self, symms, centric=False):
        """
        Constructor.
        """
        self.centric = centric
        self.symms = symms
        self.ID = SymmetryElement.symm_ID
        SymmetryElement.symm_ID += 1
        lines = []
        trans = []
        for symm in self.symms:
            line, t = self._parse_line(symm)
            lines.append(line)
            trans.append(t)
        self.matrix = Matrix(lines).transpose()
        self.trans = Array(trans)
        if centric:
            self.matrix *= -1
            self.trans *= -1

    def __str__(self):
        string = r'''
        |{aa:2} {ab:2} {ac:2}|   |{v:>4.2}|
        |{ba:2} {bb:2} {bc:2}| + |{vv:>4.2}|
        |{ca:2} {cb:2} {cc:2}|   |{vvv:>4.2}|'''.format(aa=self.matrix[0, 0],
                                                        ab=self.matrix[0, 1],
                                                        ac=self.matrix[0, 2],
                                                        ba=self.matrix[1, 0],
                                                        bb=self.matrix[1, 1],
                                                        bc=self.matrix[1, 2],
                                                        ca=self.matrix[2, 0],
                                                        cb=self.matrix[2, 1],
                                                        cc=self.matrix[2, 2],
                                                        v=float(self.trans[0]),
                                                        vv=float(self.trans[1]),
                                                        vvv=float(self.trans[2]))
        return string

    def __eq__(self, other):
        """
        Check two SymmetryElement instances for equivalence.
        Note that differences in lattice translation are ignored.
        :param other: SymmetryElement instance
        :return: True/False
        """
        m = (self.matrix == other.matrix).all()
        t1 = Array([v % 1 for v in self.trans])
        t2 = Array([v % 1 for v in other.trans])
        t = (t1 == t2).all()
        return m and t

    def __sub__(self, other):
        """
        Computes and returns the translational difference between two SymmetryElements. Returns 999.0 if the elements
        cannot be superimposed via an integer shift of the translational parts.
        :param other: SymmetryElement instance
        :return: float
        """
        if not self == other:
            return 999.
        return self.trans - other.trans

    def applyLattSymm(self, lattSymm):
        """
        Copies SymmetryElement instance and returns the copy after applying the translational part of 'lattSymm'.
        :param lattSymm: SymmetryElement.
        :return: SymmetryElement.
        """
        # newSymm = deepcopy(self)
        newSymm = SymmetryElement(self.toShelxl().split(','))
        newSymm.trans = Array([(self.trans[0] + lattSymm.trans[0]) / 1,
                               (self.trans[1] + lattSymm.trans[1]) / 1,
                               (self.trans[2] + lattSymm.trans[2]) / 1])
        newSymm.centric = self.centric
        return newSymm

    def toShelxl(self):
        """
        Generate and return string representation of Symmetry Operation in Shelxl syntax.
        :return: string.
        """
        axes = ['X', 'Y', 'Z']
        lines = []
        for i in range(3):
            text = str(self.trans[i]) if self.trans[i] else ''
            for j in range(3):
                s = '' if not self.matrix[i, j] else axes[j]
                if self.matrix[i, j] < 0:
                    s = '-' + s
                elif s:
                    s = '+' + s
                text += s
            lines.append(text)
        return ', '.join(lines)

    def _parse_line(self, symm):
        symm = symm.upper().replace(' ', '')
        chars = ['X', 'Y', 'Z']
        line = []
        for char in chars:
            element, symm = self._partition(symm, char)
            line.append(element)
        if symm:
            trans = self._float(symm)
        else:
            trans = 0
        return line, trans

    def _float(self, string):
        try:
            return float(string)
        except ValueError:
            if '/' in string:
                string = string.replace('/', './') + '.'
                return eval('{}'.format(string))

    def _partition(self, symm, char):
        parts = symm.partition(char)
        if parts[1]:
            if parts[0]:
                sign = parts[0][-1]
            else:
                sign = '+'
            if sign is '-':
                return -1, ''.join((parts[0][:-1], parts[2]))
            else:
                return 1, ''.join((parts[0], parts[2])).replace('+', '')
        else:
            return 0, symm


##### End of work by Jens Lübben #############


def my_isnumeric(value: str):
    """
    Determines if a string can be converted to a number.
    """
    try:
        float(value)
    except ValueError:
        return False
    return True


def mean(values):
    """
    returns mean value of a list of numbers

    >>> mean([1, 2, 3, 4, 1, 2, 3, 4])
    2.5
    >>> round(mean([1, 2, 3, 4, 1, 2, 3, 4.1, 1000000]), 4)
    111113.3444
    """
    return sum(values) / float(len(values))


def median(nums):
    """
    calculates the median of a list of numbers
    >>> median([2])
    2
    >>> median([1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4])
    2.5
    >>> median([1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4.1, 1000000])
    3
    >>> median([])
    Traceback (most recent call last):
    ...
    ValueError: Need a non-empty iterable
    """
    ls = sorted(nums)
    n = len(ls)
    if n == 0:
        raise ValueError("Need a non-empty iterable")
    # for uneven list length:
    elif n % 2 == 1:
        # // is floordiv:
        return ls[n // 2]
    else:
        return sum(ls[int(int(n) / 2 - 1):int(int(n) / 2 + 1)]) / 2.0


def std_dev(data):
    """
    returns standard deviation of values rounded to pl decimal places
    S = sqrt( (sum(x-xm)^2) / n-1 )
    xm = sum(x)/n
    :param data: list with integer or float values
    :type data: list
    >>> l1 = [1.334, 1.322, 1.345, 1.451, 1.000, 1.434, 1.321, 1.322]
    >>> l2 = [1.234, 1.222, 1.345, 1.451, 2.500, 1.234, 1.321, 1.222]
    >>> round(std_dev(l1), 8)
    0.13797871
    >>> round(std_dev(l2), 8)
    0.43536797
    >>> median(l1)
    1.328
    >>> mean(l1)
    1.316125
    """
    if len(data) == 0:
        return 0
    K = data[0]
    n = 0
    Sum = 0
    Sum_sqr = 0
    for x in data:
        n += 1
        Sum += x - K
        Sum_sqr += (x - K) * (x - K)
    variance = (Sum_sqr - (Sum * Sum) / n) / (n - 1)
    # use n instead of (n-1) if want to compute the exact variance of the given data
    # use (n-1) if data are samples of a larger population
    return sqrt(variance)


def nalimov_test(data):
    """
    returns a index list of outliers base on the Nalimov test for data.
    Modified implementation of:
    "R. Kaiser, G. Gottschalk, Elementare Tests zur Beurteilung von Messdaten
    Bibliographisches Institut, Mannheim 1972."

    >>> data = [1.120, 1.234, 1.224, 1.469, 1.145, 1.222, 1.123, 1.223, 1.2654, 1.221, 1.215]
    >>> nalimov_test(data)
    [3]
    """
    # q-values for degrees of freedom:
    f = {1: 1.409, 2: 1.645, 3: 1.757, 4: 1.814, 5: 1.848, 6: 1.870, 7: 1.885, 8: 1.895,
         9: 1.903, 10: 1.910, 11: 1.916, 12: 1.920, 13: 1.923, 14: 1.926, 15: 1.928,
         16: 1.931, 17: 1.933, 18: 1.935, 19: 1.936, 20: 1.937, 30: 1.945}
    fact = sqrt(float(len(data)) / (len(data) - 1))
    fval = len(data) - 2
    if fval < 2:
        return []
    outliers = []
    if fval in f:
        # less strict than the original:
        q_crit = f[fval]
    else:
        q_crit = 1.95
    for num, i in enumerate(data):
        q = abs(((i - median(data)) / std_dev(data)) * fact)
        if q > q_crit:
            outliers.append(num)
    return outliers


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    """
    returns a random ID like 'L5J74W'
    :param size: length of the string
    :type size: integer
    :param chars: characters used for the ID
    :type chars: string

    >>> id_generator(1, 'a')
    'a'
    """
    return ''.join(random.choice(chars) for _ in range(size))


def atomic_distance(p1: list, p2: list, cell=None, shortest_dist=False):
    """
    p1 and p2 are x, y , z coordinates as list ['x', 'y', 'z']
    cell are the cell parameters as list: ['a', 'b', 'c', 'alpha', 'beta', 'gamma']

    Returns the distance between the two points (Atoms). If shortest_dist is True, the
    shortest distance ignoring translation is computed.

    >>> cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
    >>> coord1 = [-0.186843,   0.282708,   0.526803]
    >>> coord2 = [-0.155278,   0.264593,   0.600644]
    >>> atomic_distance(coord1, coord2, cell)
    1.5729229943265979
    """
    if cell:
        a, b, c = cell[:3]
        al = radians(cell[3])
        be = radians(cell[4])
        ga = radians(cell[5])
    if shortest_dist:
        x1, y1, z1 = [x + 99.5 for x in p1]
        x2, y2, z2 = [x + 99.5 for x in p2]
        dx = (x1 - x2) % 1 - 0.5
        dy = (y1 - y2) % 1 - 0.5
        dz = (z1 - z2) % 1 - 0.5
    else:
        x1, y1, z1 = p1
        x2, y2, z2 = p2
        dx = (x1 - x2)
        dy = (y1 - y2)
        dz = (z1 - z2)
    if cell:
        return sqrt((a * dx) ** 2 + (b * dy) ** 2 + (c * dz) ** 2 + 2 * b * c * cos(al) * dy * dz + \
                    2 * dx * dz * a * c * cos(be) + 2 * dx * dy * a * b * cos(ga))
    else:
        return sqrt(dx ** 2 + dy ** 2 + dz ** 2)


def zero(m, n):
    """
    Create zero matrix of dimension m,n
    :param m: integer
    :param n: integer

    >>> zero(5, 3)
    [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    """
    new_matrix = [[0 for row in range(n)] for col in range(m)]  # @UnusedVariable
    return new_matrix


def determinante(a):
    """
    return determinant of 3x3 matrix

    >>> m1 = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    >>> determinante(m1)
    8
    """
    return (a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
            - a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2])
            + a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]))


def subtract_vect(a, b):
    """
    subtract vector b from vector a
    Deprecated, use mpmath instead!!!
    :param a: [float, float, float]
    :param b: [float, float, float]

    >>> subtract_vect([1, 2, 3], [3, 2, 2])
    (-2, 0, 1)
    """
    return (a[0] - b[0],
            a[1] - b[1],
            a[2] - b[2])


def transpose(a):
    """
    transposes a matrix

    >>> m = [[1, 2, 3], [1, 2, 3], [1, 2, 3]]
    >>> transpose(m)
    [(1, 1, 1), (2, 2, 2), (3, 3, 3)]
    """
    return list(zip(*a))


def norm_vec(a):
    """
    returns a normalized vector

    >>> norm_vec([1, 2, 1])
    (0.4082482904638631, 0.8164965809277261, 0.4082482904638631)
    """
    l = sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2)
    return a[0] / l, a[1] / l, a[2] / l


def dice_coefficient(a, b, case_insens=True):
    """
    :type a: str
    :type b: str
    :type case_insens: bool
    dice coefficient 2nt/na + nb.
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Dice%27s_coefficient#Python
    >>> dice_coefficient('hallo', 'holla')
    0.25
    >>> dice_coefficient('Banze', 'Benzene')
    0.444444
    >>> dice_coefficient('halo', 'Haaallo')
    0.75
    >>> dice_coefficient('hallo', 'Haaallo')
    0.888889
    >>> dice_coefficient('hallo', 'Hallo')
    1.0
    >>> dice_coefficient('aaa', 'BBBBB')
    0.0
    """
    if case_insens:
        a = a.lower()
        b = b.lower()
    if not len(a) or not len(b):
        return 0.0
    if len(a) == 1:
        a = a + u'.'
    if len(b) == 1:
        b = b + u'.'
    a_bigram_list = []
    for i in range(len(a) - 1):
        a_bigram_list.append(a[i:i + 2])
    b_bigram_list = []
    for i in range(len(b) - 1):
        b_bigram_list.append(b[i:i + 2])
    a_bigrams = set(a_bigram_list)
    b_bigrams = set(b_bigram_list)
    overlap = len(a_bigrams & b_bigrams)
    dice_coeff = overlap * 2.0 / (len(a_bigrams) + len(b_bigrams))
    return round(dice_coeff, 6)


def dice_coefficient2(a, b, case_insens=True):
    """
    :type a: str
    :type b: str
    :type case_insens: bool
    duplicate bigrams in a word should be counted distinctly
    (per discussion), otherwise 'AA' and 'AAAA' would have a
    dice coefficient of 1...
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Dice%27s_coefficient#Python

    This implementation is reverse. 1 means not hit, 0 means best match
    >>> dice_coefficient2('hallo', 'holla')
    0.75
    >>> dice_coefficient2('Banze', 'Benzene')
    0.6
    >>> dice_coefficient2('halo', 'Haaallo')
    0.333333
    >>> dice_coefficient2('hallo', 'Haaallo')
    0.2
    >>> dice_coefficient2('hallo', 'Hallo')
    0.0
    >>> dice_coefficient2('aaa', 'BBBBB')
    1.0
    >>> dice_coefficient2('', '')
    1.0
    """
    if case_insens:
        a = a.lower()
        b = b.lower()
    if not len(a) or not len(b):
        return 1.0
    # quick case for true duplicates
    if a == b:
        return 0.0
    # if a != b, and a or b are single chars, then they can't possibly match
    if len(a) == 1 or len(b) == 1:
        return 1.0
    # use python list comprehension, preferred over list.append()
    a_bigram_list = [a[i:i + 2] for i in range(len(a) - 1)]
    b_bigram_list = [b[i:i + 2] for i in range(len(b) - 1)]
    a_bigram_list.sort()
    b_bigram_list.sort()
    # assignments to save function calls
    lena = len(a_bigram_list)
    lenb = len(b_bigram_list)
    # initialize match counters
    matches = i = j = 0
    while i < lena and j < lenb:
        if a_bigram_list[i] == b_bigram_list[j]:
            matches += 2
            i += 1
            j += 1
        elif a_bigram_list[i] < b_bigram_list[j]:
            i += 1
        else:
            j += 1
    score = float(matches) / float(lena + lenb)
    score = 1 - score
    return round(score, 6)


def fft(x):
    """
    fft implementation from rosettacode.
    The purpose of this task is to calculate the FFT (Fast Fourier Transform) of an input sequence.
    The most general case allows for complex numbers at the input and results in a sequence of
    equal length, again of complex numbers. If you need to restrict yourself to real numbers,
    the output should be the magnitude (i.e. sqrt(re²+im²)) of the complex result.
    :param x:
    :type x:

    >>> print( ' '.join("%5.3f" % abs(f) for f in fft([1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])) )
    4.000 2.613 0.000 1.082 0.000 1.082 0.000 2.613
    """
    from cmath import exp, pi
    N = len(x)
    if N <= 1: return x
    even = fft(x[0::2])
    odd = fft(x[1::2])
    T = [exp(-2j * pi * k / N) * odd[k] for k in range(int(N / 2))]
    return [even[k] + T[k] for k in range(int(N / 2))] + \
           [even[k] - T[k] for k in range(int(N / 2))]


def levenshtein(s1, s2):
    """
    >>> levenshtein('hallo', 'holla')
    2
    >>> dice_coefficient('hallo', 'holla')
    0.25
    """
    s1 = s1.lower()
    s2 = s2.lower()
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            # j+1 instead of j since previous_row and current_row are one character longer:
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1  # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def distance(x1, y1, z1, x2, y2, z2, round_out=False):
    """
    distance between two points in space for orthogonal axes.
    >>> distance(1, 1, 1, 2, 2, 2, 4)
    1.7321
    >>> distance(1, 0, 0, 2, 0, 0, 4)
    1.0
    """
    import math as m
    d = m.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
    if round_out:
        return round(d, round_out)
    else:
        return d


def vol_unitcell(a, b, c, al, be, ga):
    """
    calculates the volume of a unit cell
    >>> v = vol_unitcell(2, 2, 2, 90, 90, 90)
    >>> print(v)
    8.0
    """
    ca, cb, cg = cos(radians(al)), cos(radians(be)), cos(radians(ga))
    v = a * b * c * sqrt(1 + 2 * ca * cb * cg - ca ** 2 - cb ** 2 - cg ** 2)
    return v


class A(object):
    """
    orthogonalization matrix
    e.g. converts fractional coordinates to cartesian coodinates
    >>> import mpmath as mpm
    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    >>> coord = (-0.186843,   0.282708,   0.526803)
    >>> A = A(cell).orthogonal_matrix
    >>> print(mpm.nstr(A*mpm.matrix(coord)))
    [-2.74151]
    [ 5.90959]
    [ 10.7752]
    >>> cartcoord = mpm.matrix([['-2.74150542399906'], ['5.909586678'], ['10.7752007008937']])
    >>> print(mpm.nstr(A**-1*cartcoord))
    [-0.186843]
    [ 0.282708]
    [ 0.526803]
    """

    def __init__(self, cell):
        self.a, self.b, self.c, alpha, beta, gamma = cell
        self.V = vol_unitcell(self.a, self.b, self.c, alpha, beta, gamma)
        self.alpha = radians(alpha)
        self.beta = radians(beta)
        self.gamma = radians(gamma)

    @property
    def orthogonal_matrix(self):
        """
        Converts von fractional to cartesian.
        Invert the matrix to do the opposite.
        """
        import mpmath as mpm
        Am = mpm.matrix([[self.a, self.b * cos(self.gamma), self.c * cos(self.beta)],
                         [0, self.b * sin(self.gamma),
                          (self.c * (cos(self.alpha) - cos(self.beta) * cos(self.gamma)) / sin(self.gamma))],
                         [0, 0, self.V / (self.a * self.b * sin(self.gamma))]])
        return Am


def calc_ellipsoid_axes(coords: list, uvals: list, cell: list, probability: float = 0.5, longest: bool = True):
    """
    This method calculates the principal axes of an ellipsoid as list of two
    fractional coordinate triples.
    Many thanks to R. W. Grosse-Kunstleve and P. D. Adams
    for their great publication on the handling of atomic anisotropic displacement
    parameters:
    R. W. Grosse-Kunstleve, P. D. Adams, J Appl Crystallogr 2002, 35, 477–480.

    F = ... * exp ( -2π²[ h²(a*)²U11 + k²(b*)²U22 + ... + 2hka*b*U12 ] )

    SHELXL atom:
    Name type  x      y      z    occ     U11 U22 U33 U23 U13 U12
    F3    4    0.210835   0.104067   0.437922  21.00000   0.07243   0.03058 =
       0.03216  -0.01057  -0.01708   0.03014
    >>> import mpmath as mpm
    >>> cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
    >>> coords = [0.210835,   0.104067,   0.437922]
    >>> uvals = [0.07243, 0.03058, 0.03216, -0.01057, -0.01708, 0.03014]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=True)
    [(mpf('0.24765096088767435'), mpf('0.11383281312627748'), mpf('0.43064756017994188')), (mpf('0.17401903911232561'), mpf('0.094301186873722534'), mpf('0.44519643982005791'))]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=False)
    [[(mpf('0.24765096088767435'), mpf('0.11383281312627748'), mpf('0.43064756017994188')), (mpf('0.21840599615600714'), mpf('0.096261419454711866'), mpf('0.4374612735649967')), (mpf('0.21924358264443167'), mpf('0.10514684301214304'), mpf('0.44886867568495714'))], [(mpf('0.17401903911232561'), mpf('0.094301186873722534'), mpf('0.44519643982005791')), (mpf('0.20326400384399282'), mpf('0.11187258054528816'), mpf('0.43838272643500309')), (mpf('0.20242641735556829'), mpf('0.10298715698785697'), mpf('0.4269753243150427'))]]
    >>> cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
    >>> coords = [0.210835,   0.104067,   0.437922]
    >>> uvals = [0.07243, -0.03058, 0.03216, -0.01057, -0.01708, 0.03014]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=True)
    *** Ellipsoid is non positive definite! ***
    (False, False)

    >>> uvals = [0.07243, 0.03058, 0.03216, -0.01057, -0.01708]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=False)
    Traceback (most recent call last):
    ...
    Exception: 6 Uij values have to be supplied!

    >>> cell = [10.5086, 20.9035, 90, 94.13, 90]
    >>> coords = [0.210835,   0.104067,   0.437922]
    >>> uvals = [0.07243, 0.03058, 0.03216, -0.01057, -0.01708, 0.03014]
    >>> calc_ellipsoid_axes(coords, uvals, cell, longest=True)
    Traceback (most recent call last):
    ...
    Exception: cell needs six parameters!

    :param coords: coordinates of the respective atom in fractional coordinates
    :type coords: list
    :param uvals: Uij valiues of the respective ellipsoid on fractional
                  basis like in cif and SHELXL format
    :type uvals: list
    :param cell: unit cell of the structure: a, b, c, alpha, beta, gamma
    :type cell:  list
    :param probability: thermal probability of the ellipsoid
    :type probability: float or int
    :param longest: not always the length is important. make to False to
                    get all three coordiantes of the ellipsoid axes.
    :type longest: boolean

    """
    import mpmath as mpm
    probability += 1
    # Uij is symmetric:
    if len(uvals) != 6:
        raise Exception('6 Uij values have to be supplied!')
    if len(cell) != 6:
        raise Exception('cell needs six parameters!')
    U11, U22, U33, U23, U13, U12 = uvals
    U21 = U12
    U32 = U23
    U31 = U13
    Uij = mpm.matrix([[U11, U12, U13], [U21, U22, U23], [U31, U32, U33]])
    a, b, c, alpha, beta, gamma = cell
    V = vol_unitcell(*cell)
    # calculate reciprocal lattice vectors:
    astar = (b * c * sin(radians(alpha))) / V
    bstar = (c * a * sin(radians(beta))) / V
    cstar = (a * b * sin(radians(gamma))) / V
    # orthogonalization matrix that transforms the fractional coordinates
    # with respect to a crystallographic basis system to coordinates
    # with respect to a Cartesian basis:
    amatrix = A(cell).orthogonal_matrix
    # matrix with the reciprocal lattice vectors:
    N = mpm.matrix([[astar, 0, 0],
                    [0, bstar, 0],
                    [0, 0, cstar]])
    # Finally transform Uij values from fractional to cartesian axis system:
    Ucart = amatrix * N * Uij * N.T * amatrix.T
    # E => eigenvalues, Q => eigenvectors:
    E, Q = mpm.eig(Ucart)
    # calculate vectors of ellipsoid axes
    try:
        sqrt(E[0])
        sqrt(E[1])
        sqrt(E[2])
    except ValueError:
        print('*** Ellipsoid is non positive definite! ***')
        return (False, False)
    v1 = mpm.matrix([Q[0, 0], Q[1, 0], Q[2, 0]])
    v2 = mpm.matrix([Q[0, 1], Q[1, 1], Q[2, 1]])
    v3 = mpm.matrix([Q[0, 2], Q[1, 2], Q[2, 2]])
    v1i = v1 * (-1)
    v2i = v2 * (-1)
    v3i = v3 * (-1)
    # multiply probability (usually 50%)
    e1 = sqrt(E[0]) * probability
    e2 = sqrt(E[1]) * probability
    e3 = sqrt(E[2]) * probability
    # scale axis vectors to eigenvalues
    v1, v2, v3, v1i, v2i, v3i = v1 * e1, v2 * e2, v3 * e3, v1i * e1, v2i * e2, v3i * e3
    # find out which vector is the longest:
    length = mpm.norm(v1)
    v = 0
    if mpm.norm(v2) > length:
        length = mpm.norm(v2)
        v = 1
    elif mpm.norm(v3) > length:
        length = mpm.norm(v3)
        v = 2
    # move vectors back to atomic position
    atom = amatrix * mpm.matrix(coords)
    v1, v1i = v1 + atom, v1i + atom
    v2, v2i = v2 + atom, v2i + atom
    v3, v3i = v3 + atom, v3i + atom
    # go back into fractional coordinates:
    a1 = cart_to_frac(v1, cell)
    a2 = cart_to_frac(v2, cell)
    a3 = cart_to_frac(v3, cell)
    a1i = cart_to_frac(v1i, cell)
    a2i = cart_to_frac(v2i, cell)
    a3i = cart_to_frac(v3i, cell)
    allvec = [[a1, a2, a3], [a1i, a2i, a3i]]
    if longest:
        # only the longest vector
        return [allvec[0][v], allvec[1][v]]
    else:
        # all vectors:
        return allvec


def almost_equal(a, b, places=3):
    """
    Returns True or False if the number a and b are are equal inside the
    decimal places "places".
    :param a: a real number
    :type a: int/float
    :param b: a real number
    :type b: int/float
    :param places: number of decimal places
    :type places: int

    >>> almost_equal(1.0001, 1.0005)
    True
    >>> almost_equal(1.1, 1.0005)
    False
    >>> almost_equal(2, 1)
    False
    """
    return round(abs(a - b), places) == 0


def frac_to_cart(frac_coord: list, cell: list) -> tuple:
    """
    Converts fractional coordinates to cartesian coodinates
    :param frac_coord: [float, float, float]
    :param cell:       [float, float, float, float, float, float]

    >>> cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
    >>> coord1 = [-0.186843,   0.282708,   0.526803]
    >>> print(frac_to_cart(coord1, cell))
    (-2.741505423999065, 5.909586678000002, 10.775200700893734)
    """
    a, b, c, alpha, beta, gamma = cell
    x, y, z = frac_coord
    alpha = radians(alpha)
    beta = radians(beta)
    gamma = radians(gamma)
    cosastar = (cos(beta) * cos(gamma) - cos(alpha)) / (sin(beta) * sin(gamma))
    sinastar = sqrt(1 - cosastar ** 2)
    xc = a * x + (b * cos(gamma)) * y + (c * cos(beta)) * z
    yc = 0 + (b * sin(gamma)) * y + (-c * sin(beta) * cosastar) * z
    zc = 0 + 0 + (c * sin(beta) * sinastar) * z
    return xc, yc, zc


def cart_to_frac(cart_coord: list, cell: list) -> tuple:
    """
    converts cartesian coordinates to fractional coordinates
    :param cart_coord: [float, float, float]
    :param cell:       [float, float, float, float, float, float]
    >>> cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
    >>> coords = [-2.74150542399906, 5.909586678, 10.7752007008937]
    >>> cart_to_frac(coords, cell)
    (-0.1868429999999998, 0.28270799999999996, 0.5268029999999984)
    """
    a, b, c, alpha, beta, gamma = cell
    xc, yc, zc = cart_coord
    alpha = radians(alpha)
    beta = radians(beta)
    gamma = radians(gamma)
    cosastar = (cos(beta) * cos(gamma) - cos(alpha)) / (sin(beta) * sin(gamma))
    sinastar = sqrt(1 - cosastar ** 2)
    z = zc / (c * sin(beta) * sinastar)
    y = (yc - (-c * sin(beta) * cosastar) * z) / (b * sin(gamma))
    x = (xc - (b * cos(gamma)) * y - (c * cos(beta)) * z) / a
    return x, y, z
