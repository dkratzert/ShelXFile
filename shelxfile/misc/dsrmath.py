# -*- encoding: utf-8 -*-
# möp
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <dkratzert@gmx.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
import random
import string
from math import sqrt, radians, cos, sin, acos, degrees, floor
from operator import sub, add
from typing import List, Union, Optional, Iterable, Tuple

from shelxfile.misc.misc import flatten, determinante


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
    """
    __slots__ = ['values']

    def __init__(self, values: Union[list, tuple]):
        self.values = values

    def __iter__(self) -> Iterable[Union[int, float]]:
        for v in self.values:
            yield v

    def __len__(self) -> int:
        return len(self.values)

    def __hash__(self):
        return hash(self.values)

    def __add__(self, other: (list, 'Array')) -> 'Array':
        """
        This method is optimized for speed.
        """
        if isinstance(other, Array):
            return Array(tuple(map(add, self.values, other)))
        elif isinstance(other, (float, int)):
            # return Array( list(map(lambda x: x - other, self.values)) )
            return Array([i + other for i in self.values])
        else:
            raise TypeError(f'Cannot add type Array to type {str(type(other))}.')

    def __iadd__(self, other):
        return self.__add__(other)

    def __eq__(self, other):
        """
        Test for equality.
        """
        return all([a == b for (a, b) in zip(self.values, other.values)])

    def __sub__(self, other):
        """
        Subtracts eiter an Array or a value from the self Array.
        This method is optimized for speed.
        """
        if isinstance(other, Array):
            return Array(list(map(sub, self.values, other)))  # slightly faster
        elif isinstance(other, float) or isinstance(other, int):
            return Array([i - other for i in self.values])
        else:
            raise TypeError('Cannot add type Array to type {}.'.format(str(type(other))))

    def __imul__(self, other):
        """
        Currently supports multiplication by a number.
        __imul__ means a *= b
        """
        if isinstance(other, (int, float)):
            self.values = [v * other for v in self.values]
            return self
        else:
            raise TypeError('Unsupported operation.')

    def __mul__(self, other: ('Array', 'Matrix')) -> (float, 'Array'):
        """
        Calculates: a * b = axbx + ayby + azbz
        """
        if isinstance(other, Matrix):
            # Array() * Matrix()
            a = other[0]
            b = other[1]
            c = other[2]
            x = self.values[0] * a[0] + self.values[1] * b[0] + self.values[2] * c[0]
            y = self.values[0] * a[1] + self.values[1] * b[1] + self.values[2] * c[1]
            z = self.values[0] * a[2] + self.values[1] * b[2] + self.values[2] * c[2]
            return Array([x, y, z])
        else:
            return self.dot(other)

    def __repr__(self):
        return 'Array({})'.format(str(self.values))

    def __getitem__(self, val):
        """
        Get one item from the array.
        """
        return self.values[val]

    def __setitem__(self, pos, val):
        """
        Get one item from the array.
        """
        self.values[pos] = val

    def norm(self):
        """
        The squared lenght of an array
        """
        return sum([n ** 2 for n in self.values])

    def normalized(self):
        """
        Euclidean norm (straight-line distance) of a vector array.
        """
        return sqrt(self.norm())

    @staticmethod
    def zero(m: int) -> 'Array':
        """
        Create zero Array of dimension m
        """
        return Array([0.0 for _ in range(m)])

    @staticmethod
    def randarray(m: int) -> 'Array':
        """
        Create random Array of dimension m
        """
        return Array([random.randint(1, 99) for _ in range(m)])

    @property
    def floor(self):
        return Array(list(map(floor, self.values)))

    def dot(self, other: 'Array') -> float:
        """
        Dot product of an array in kartesian space.
        """
        if len(self) != len(other):
            raise ValueError('Vector sizes must match')
        return sum([i * j for i, j in zip(self, other)])

    def cross(self, other: 'Array') -> 'Array':
        """
        Cross product of the Array (currently only for 3D vectors).
        """
        if len(self) != len(other) != 3:
            raise ValueError('For 3D vectors only')
        a1, a2, a3 = self
        b1, b2, b3 = other
        return Array([(a2 * b3 - a3 * b2), (a3 * b1 - a1 * b3), (a1 * b2 - a2 * b1)])

    def angle(self, other: 'Array') -> float:
        """
        Calculates the angle between two vectors.
        """
        return round(degrees(acos(self.dot(other) / (self.normalized() * other.normalized()))), 9)


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
    """
    __slots__ = ['values', 'shape', 'rows', 'columns']

    def __init__(self, values):
        self.shape = (len(values[0]), len(values))
        self.rows = len(values[0])
        self.columns = len(values)
        self.values = values

    def __getitem__(self, val):
        """
        """
        if isinstance(val, (tuple, list)):
            if len(val) == 1:
                return self.values[val[0]]
            return self.values[val[1]][val[0]]
        else:
            return self.values[val]

    def __repr__(self):
        rows = ''
        for row in self.values:
            rows += '|' + ' '.join(['{:>7.4f}'.format(float(x)) for x in row]) + '|' + '\n'
        return rows

    def __add__(self, other: (list, 'Matrix')) -> 'Matrix':
        """
        Matrix addition
        """
        if isinstance(other, Array):
            if len(self) != len(other):
                raise ValueError('Matrix and Array are not of equal length.')
            return Matrix(
                [[sum(e1, e2) for e1, e2 in zip(row1, row2)] for row1, row2 in zip(self.values, other.values)])
        elif isinstance(other, (float, int)):
            return Matrix([[x + other for x in i] for i in self.values])
        elif isinstance(other, Matrix):
            return Matrix([[e1 + e2 for e1, e2 in zip(row1, row2)] for row1, row2 in zip(self.values, other.values)])
        else:
            raise TypeError('Cannot add type {} Array to Matrix.'.format(str(type(other))))

    def __mul__(self, other: ('Matrix', 'Array', int, float)) -> ('Matrix', 'Array'):
        """
        a * b operation
        """
        if isinstance(other, (int, float)):
            return Matrix([[v * other for v in row] for row in self.values])
        elif isinstance(other, (Matrix, OrthogonalMatrix)):
            return Matrix([[sum(ea * eb for ea, eb in zip(a, b)) for b in other.values] for a in self.values])
        elif isinstance(other, Array):
            return Array([sum([b * x for (b, x) in zip(other.values, row)]) for row in self.values])
        else:
            raise TypeError('Cannot add type {} to Matrix.'.format(str(type(other))))

    def __len__(self):
        return self.shape[1]

    def __iter__(self) -> Iterable[List[float]]:
        for n in self.values:
            yield n

    def __eq__(self, other):
        """
        Test for equality.
        """
        return all([b == x for (b, x) in zip(other.values, self.values)])

    def __sub__(self, other):
        """
        Substract two matrices.
        """
        output = []
        for idx in range(len(self)):
            tmp = []
            for val_a, val_b in zip(self[idx], other[idx]):
                tmp.append(val_a - val_b)
            output.append(tmp[:])
        return output[:]

    def __truediv__(self, other: Union['Matrix', float, int]):
        """
        #>>> Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]]) / Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
        A / B = A * B^-1
        """
        if isinstance(other, Matrix):
            return self * other.inversed
        elif isinstance(other, (int, float)):
            return Matrix([[v / other for v in row] for row in self.values])
        else:
            raise NotImplementedError

    def __setitem__(self, key, value):
        # TODO: Implement setitem
        pass

    @property
    def T(self):
        return self.transposed

    @property
    def transposed(self) -> 'Matrix':
        """
        transposes a matrix
        """
        return Matrix(list(zip(*self.values)))

    @property
    def transposed_alt(self) -> 'Matrix':
        """
        Transposes the current matrix.
        """
        rows = []
        for i in range(self.shape[0]):
            rows.append([r[i] for r in self.values])
        return Matrix(rows)

    @property
    def trace(self):
        return self.values[0][0] + self.values[1][1] + self.values[2][2]

    def dot(self, other):
        """
        Dot product of two matrices.
        """
        result = [[sum(a * b for a, b in zip(x_row, y_col)) for y_col in zip(*other)] for x_row in self]
        return Matrix(result)

    @staticmethod
    def zero(m: int, n: int) -> 'Matrix':
        """
        Create zero matrix of dimension m,n
        """
        return Matrix([[0.0 for _ in range(n)] for _ in range(m)])

    @staticmethod
    def randmat(m: int, n: int) -> 'Matrix':
        """
        Create random matrix of dimension m, n (rows, columns)

        #>>> Matrix.randmat(5, 3)
        |  2.000   1.000   2.000|
        |  8.000   1.000   2.000|
        |  3.000   5.000   1.000|
        |  4.000   1.000   9.000|
        |  2.000   8.000   4.000|
        <BLANKLINE>
        """
        return Matrix([[random.randint(1, 99) for _ in range(n)] for _ in range(m)])

    def cholesky(self) -> 'Matrix':
        """
        """
        L = Matrix.zero(*self.shape)
        for i, (ai, li) in enumerate(zip(self.values, L.values)):
            for j, lj in enumerate(L.values[:i + 1]):
                s = sum(li[k] * lj[k] for k in range(j))
                li[j] = sqrt(ai[i] - s) if (i == j) else (1.0 / lj[j] * (ai[j] - s))
        return L

    @property
    def inversed(self) -> 'Matrix':
        """
        Inversion of 3 × 3 matrices
        """
        if self.shape != (3, 3):
            raise ValueError('Inversion is only valid for 3x3 Matrix.')
        d = self.det
        m1, m2, m3, m4, m5, m6, m7, m8, m9 = flatten(self.values)
        inv = Matrix([[(m5 * m9 - m6 * m8) / d, (m3 * m8 - m2 * m9) / d, (m2 * m6 - m3 * m5) / d],
                      [(m6 * m7 - m4 * m9) / d, (m1 * m9 - m3 * m7) / d, (m3 * m4 - m1 * m6) / d],
                      [(m4 * m8 - m5 * m7) / d, (m2 * m7 - m1 * m8) / d, (m1 * m5 - m2 * m4) / d]])
        return inv

    @property
    def det(self):
        """
        Return determinant of 3x3 matrix.
        """
        return determinante(self.values)

    @property
    def norm(self):
        return self.frobenius_norm()

    def frobenius_norm(self):
        # To store the sum of squares of the
        # elements of the given matrix
        sumSq = 0
        for i in range(self.rows):
            for j in range(self.columns):
                sumSq += pow(self[i][j], 2)
        # Return the square root of
        # the sum of squares
        return sqrt(sumSq)

    def power_iteration(self, num_simulations=10):
        """
        Eigenvalue algorythm from https://en.wikipedia.org/wiki/Power_iteration

        This does not work. DK
        """
        # Ideally choose a random vector
        # To decrease the chance that our vector
        # Is orthogonal to the eigenvector
        b_k = Matrix.randmat(self.rows, self.columns)
        for _ in range(num_simulations):
            # calculate the matrix-by-vector product Ab
            b_k1 = self.dot(b_k)
            # calculate the norm
            b_k1_norm = b_k1.norm
            # re normalize the vector
            b_k = b_k1 / b_k1_norm
        return b_k


class SymmetryElement(object):
    """
    Class representing a symmetry operation.
    """
    symm_id = 1
    __slots__ = ['centric', 'symms', 'ID', 'matrix', 'trans']

    def __init__(self, symms, centric=False):
        """
        Constructor.
        """
        self.centric = centric
        self.symms = symms
        self.ID = SymmetryElement.symm_id
        SymmetryElement.symm_id += 1
        lines = []
        trans = []
        for symm in self.symms:
            line, t = self._parse_line(symm)
            lines.append(line)
            trans.append(t)
        self.matrix = Matrix(lines).transposed
        self.trans = Array(trans)
        if centric:
            self.matrix *= -1
            self.trans *= -1

    def __str__(self):
        string = "|{aa:2} {ab:2} {ac:2}|   |{v:>4.2}|\n" \
                 "|{ba:2} {bb:2} {bc:2}| + |{vv:>4.2}|\n" \
                 "|{ca:2} {cb:2} {cc:2}|   |{vvv:>4.2}|\n".format(aa=self.matrix[0, 0],
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

    def __repr__(self):
        return self.to_shelxl()

    def __eq__(self, other):
        """
        Check two SymmetryElement instances for equivalence.
        Note that differences in lattice translation are ignored.
        :param other: SymmetryElement instance
        :return: True/False
        """
        m = (self.matrix == other.matrix)
        t1 = Array([v % 1 for v in self.trans])
        t2 = Array([v % 1 for v in other.trans])
        t = (t1 == t2)
        return m and t

    def __sub__(self, other):
        """
        Computes and returns the translational difference between two SymmetryElements. Returns 999.0 if the elements
        cannot be superimposed via an integer shift of the translational parts.
        :param other: SymmetryElement instance
        :return: float
        """
        if self != other:
            return 999.0
        return self.trans - other.trans

    def apply_latt_symm(self, latt_symm):
        """
        Copies SymmetryElement instance and returns the copy after applying the translational part of 'lattSymm'.
        :param latt_symm: SymmetryElement.
        :return: SymmetryElement.
        """
        new_symm = SymmetryElement(self.to_shelxl().split(','))
        new_symm.trans = Array([(self.trans[0] + latt_symm.trans[0]) / 1,
                                (self.trans[1] + latt_symm.trans[1]) / 1,
                                (self.trans[2] + latt_symm.trans[2]) / 1])
        new_symm.centric = self.centric
        return new_symm

    def to_shelxl(self):
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

    def to_cif(self) -> str:
        return self._replace_float_values(self.to_shelxl()).lower()

    def _replace_float_values(self, val: str) -> str:
        val = val.replace('1.25', '5/4')
        val = val.replace('0.75', '3/4')
        val = val.replace('0.5', '1/2')
        val = val.replace('0.33', '1/3')
        val = val.replace('0.25', '1/4')
        val = val.replace('0.125', '1/6')
        return val

    def _parse_line(self, symm: str) -> Tuple[List[int], float]:
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

    def _float(self, string: str) -> float:
        try:
            return float(string)
        except ValueError:
            if '/' in string:
                string = string.replace('/', './') + '.'
                return eval('{}'.format(string))

    def _partition(self, symm: str, char: str) -> Tuple[int, str]:
        parts = symm.partition(char)
        if parts[1]:
            if parts[0]:
                sign = parts[0][-1]
            else:
                sign = '+'
            if sign == '-':
                return -1, ''.join((parts[0][:-1], parts[2]))
            else:
                return 1, ''.join((parts[0], parts[2])).replace('+', '')
        else:
            return 0, symm


# End of work by Jens Lübben #############


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
    """
    return sum(values) / float(len(values))


def median(nums):
    """
    calculates the median of a list of numbers
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


def std_dev(data: List) -> float:
    """
    returns standard deviation of values rounded to pl decimal places
    S = sqrt( (sum(x-xm)^2) / n-1 )
    xm = sum(x)/n
    :param data: list with integer or float values
    """
    if len(data) == 0:
        return 0
    K = data[0]
    n = 0
    summ = 0
    sum_sqr = 0
    for x in data:
        n += 1
        summ += x - K
        sum_sqr += (x - K) * (x - K)
    variance = (sum_sqr - (summ * summ) / n) / (n - 1)
    # use n instead of (n-1) if want to compute the exact variance of the given data
    # use (n-1) if data are samples of a larger population
    return sqrt(variance)


def nalimov_test(data):
    """
    returns a index list of outliers base on the Nalimov test for data.
    Modified implementation of:
    "R. Kaiser, G. Gottschalk, Elementare Tests zur Beurteilung von Messdaten
    Bibliographisches Institut, Mannheim 1972."
    """
    # q-values for degrees of freedom:
    f = {1 : 1.409, 2: 1.645, 3: 1.757, 4: 1.814, 5: 1.848, 6: 1.870, 7: 1.885, 8: 1.895,
         9 : 1.903, 10: 1.910, 11: 1.916, 12: 1.920, 13: 1.923, 14: 1.926, 15: 1.928,
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
    """
    return ''.join(random.choice(chars) for _ in range(size))


def atomic_distance(p1: List, p2: List, cell=None, shortest_dist=False):
    """
    p1 and p2 are x, y , z coordinates as list ['x', 'y', 'z']
    cell are the cell parameters as list: ['a', 'b', 'c', 'alpha', 'beta', 'gamma']

    Returns the distance between the two points (Atoms). If shortest_dist is True, the
    shortest distance ignoring translation is computed.
    """
    a, b, c, al, be, ga = 1, 1, 1, 1, 1, 1
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
        return sqrt((a * dx) ** 2 + (b * dy) ** 2 + (c * dz) ** 2 + 2 * b * c * cos(al) * dy * dz +
                    2 * dx * dz * a * c * cos(be) + 2 * dx * dy * a * b * cos(ga))
    else:
        return sqrt(dx ** 2 + dy ** 2 + dz ** 2)


def dice_coefficient(a, b, case_insens=True):
    """
    :type a: str
    :type b: str
    :type case_insens: bool
    dice coefficient 2nt/na + nb.
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Dice%27s_coefficient#Python
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


def levenshtein(s1, s2):
    """
    The levensteins distance of two strings.
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
    """
    ca, cb, cg = cos(radians(al)), cos(radians(be)), cos(radians(ga))
    v = a * b * c * sqrt(1 + 2 * ca * cb * cg - ca ** 2 - cb ** 2 - cg ** 2)
    return v


class OrthogonalMatrix():
    """
    Orthogonalization matrix used to convert fractional coordinates to cartesian.
    """

    def __init__(self, a, b, c, alpha, beta, gamma):
        self.a, self.b, self.c = a, b, c
        self.V = vol_unitcell(a, b, c, alpha, beta, gamma)
        self.alpha = radians(alpha)
        self.beta = radians(beta)
        self.gamma = radians(gamma)
        self.m = Matrix(((self.a, self.b * cos(self.gamma), self.c * cos(self.beta)),
                         (0, self.b * sin(self.gamma),
                          (self.c * (cos(self.alpha) - cos(self.beta) * cos(self.gamma)) / sin(self.gamma))),
                         (0, 0, self.V / (self.a * self.b * sin(self.gamma)))))
        self.metric_matrix = self.transposed.dot(self.m)
        self._inversed: Optional[Matrix] = None

    def __mul__(self, other: Array) -> Array:
        """
        To convert from fractional to cartesian.
        """
        return self.m * other

    @property
    def inversed(self):
        """
        To convert from cartesian to fractional.
        """
        if not self._inversed:
            self._inversed = self.m.inversed
            return self._inversed
        else:
            return self._inversed

    @property
    def transposed(self):
        return self.m.transposed

    # noinspection PyPep8Naming
    @property
    def T(self):
        return self.m.transposed

    @property
    def values(self):
        return self.m.values


def almost_equal(a: Union[int, float], b: Union[int, float], places=3) -> float:
    """
    Returns True or False if the numbers a and b are equal inside the
    decimal places "places".
    """
    return round(abs(a - b), places) == 0
