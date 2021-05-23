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
from math import sqrt, radians, cos, sin, acos, degrees, floor, tan
from operator import sub, add


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
    >>> a += a
    >>> a
    Array([2, 4, 6, 8.2])
    >>> a = Array([1, 2, 3, 4.1])
    >>> a + 4
    Array([5, 6, 7, 8.1])
    >>> a + a
    Array([2, 4, 6, 8.2])
    >>> a[1]
    2
    >>> a *= 3
    >>> a
    Array([3, 6, 9, 12.299999999999999])
    >>> a = Array([1, 2, 3, 4.1])
    >>> a.dot(a)
    30.81
    >>> a * a
    30.81
    """
    __slots__ = ['values']

    def __init__(self, values: list):
        self.values = values

    def __iter__(self) -> iter:
        for v in self.values:
            yield v

    def __len__(self) -> int:
        return len(self.values)

    def __hash__(self):
        return hash(self.values)

    def __add__(self, other: (list, 'Array')) -> 'Array':
        """
        This method is optimized for speed.
        >>> a = Array([1, 2, 3])
        >>> b = Array([1, 1, 1])
        >>> a+b
        Array([2, 3, 4])
        """
        if isinstance(other, Array):
            return Array(list(map(add, self.values, other)))
        elif type(other) == float or type(other) == int:
            # return Array( list(map(lambda x: x - other, self.values)) )
            return Array([i + other for i in self.values])
        else:
            raise TypeError('Cannot add type Array to type {}.'.format(str(type(other))))

    def __iadd__(self, other):
        return self.__add__(other)

    def __eq__(self, other):
        """
        >>> a1 = Array([1, 2, 3, 4])
        >>> a2 = Array([1, 2, 3.0, 4.0])
        >>> a1 == a2
        True
        >>> a1 = Array([1, 2, 3, 4])
        >>> a2 = Array([2, 2, 3.0, 4.0])
        >>> a1 == a2
        False
        """
        return all([a == b for (a, b) in zip(self.values, other.values)])

    def __sub__(self, other):
        """
        Subtracts eiter an Array or a value from the self Array.
        This method is optimized for speed.
        >>> a = Array([1, 2, 3])
        >>> b = Array([1, 1, 1])
        >>> a-b
        Array([0, 1, 2])
        >>> b-a
        Array([0, -1, -2])
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
        a * b = axbx + ayby + azbz

        >>> a1 = Array([1, 2, 3])
        >>> a2 = Array([1, 2, 3])
        >>> a1 * a2
        14
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
        >>> Array([1, 2, 3])[1]
        2
        """
        return self.values[val]

    def __setitem__(self, pos, val):
        """
        Get one item from the array.
        >>> a = Array([0, 0, 0])
        >>> a[1] = 5
        >>> a
        Array([0, 5, 0])
        """
        self.values[pos] = val

    def norm(self):
        """
        The squared lenght of an array

        >>> a = Array([1, 2, 3, 4])
        >>> a.norm()
        30
        """
        return sum([n ** 2 for n in self.values])

    def normalized(self):
        """
        Euclidean norm (straight-line distance) of a vector array.
        >>> a = Array([2, 2, 1])
        >>> a.normalized()
        3.0
        """
        return sqrt(self.norm())

    @staticmethod
    def zero(m: int) -> 'Array':
        """
        Create zero Array of dimension m

        >>> Array.zero(5)
        Array([0.0, 0.0, 0.0, 0.0, 0.0])
        """
        return Array([0.0 for row in range(m)])

    @staticmethod
    def randarray(m: int) -> 'Array':
        """
        Create zero Array of dimension m

        #>>> Array.randarray(5)
        #Array([1, 4, 6, 4, 8])
        """
        return Array([random.randint(1, 9) for row in range(m)])

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

        >>> a = Array([1, 2, 3])
        >>> b = Array([-7, 8, 9])
        >>> a.cross(b)
        Array([-6, -30, 22])
        """
        if len(self) != len(other) != 3:
            raise ValueError('For 3D vectors only')
        a1, a2, a3 = self
        b1, b2, b3 = other
        return Array([(a2 * b3 - a3 * b2), (a3 * b1 - a1 * b3), (a1 * b2 - a2 * b1)])

    def angle(self, other: 'Array') -> float:
        """
        Calculates the angle between two vectors.
        >>> a = Array([1, 0, 1])
        >>> b = Array([1, 0, 0])
        >>> a.angle(b)
        45.0
        >>> va = Array([0.03562, 0.14298, 0.24008]) - Array([0.04402, 0.16614, 0.22275])
        >>> vb = Array([0.07078, 0.17382, 0.22106]) - Array([0.04402, 0.16614, 0.22275])
        >>> round(va.angle(vb), 4)
        120.9401
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

    >>> m = Matrix([[1, 2, 300], [4.1, 4.2, 4.3], [5, 6, 7]])

    #>>> m+m
    #Array([2, 4, 6, 8.2])
    #>>> m+4
    #Array([5, 6, 7, 8.1])
    >>> m*=3
    >>> m
    | 3.0000  6.0000 900.0000|
    |12.3000 12.6000 12.9000|
    |15.0000 18.0000 21.0000|
    <BLANKLINE>
    """
    __slots__ = ['values', 'shape']

    def __init__(self, values):
        self.shape = (len(values[0]), len(values))
        self.values = values

    def __getitem__(self, val):
        """
        >>> m = Matrix([[2., 2., 3.], [1., 2.2, 3.], [1., 2., 3.]])
        >>> m[1,1]
        2.2
        >>> m[1][1]
        2.2
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

        >>> m1 = Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
        >>> m2 = Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
        >>> m1 + m2
        | 2.0000  2.0000  2.0000|
        | 2.0000  2.0000  2.0000|
        | 2.0000  2.0000  2.0000|
        <BLANKLINE>
        >>> m1 + 0.5
        | 1.5000  1.5000  1.5000|
        | 1.5000  1.5000  1.5000|
        | 1.5000  1.5000  1.5000|
        <BLANKLINE>
        """
        if isinstance(other, Array):
            if not len(self) == len(other):
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
        >>> m = Matrix([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
        >>> m * 2
        | 2.0000  2.0000  2.0000|
        | 2.0000  2.0000  2.0000|
        | 2.0000  2.0000  2.0000|
        <BLANKLINE>
        >>> m2 = Matrix([[2, 2, 2], [0.5, 0.5, 0.5], [2, 2, 1]])
        >>> m * m2
        | 6.0000  1.5000  5.0000|
        | 6.0000  1.5000  5.0000|
        | 6.0000  1.5000  5.0000|
        <BLANKLINE>
        >>> m
        | 1.0000  1.0000  1.0000|
        | 1.0000  1.0000  1.0000|
        | 1.0000  1.0000  1.0000|
        <BLANKLINE>
        >>> m * Array([2, 2, 2])
        Array([6, 6, 6])
        >>> m = Matrix([(0, 1, 0), (-1, -1, 0), (0, 0, 1)])
        >>> Array([0.333333, 0.666667, 0.45191]) * m
        Array([-0.666667, -0.333334, 0.45191])
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

    def __iter__(self) -> list:
        return [n for n in self.values]

    def __eq__(self, other):
        """
        >>> m1 = Matrix([(1., 2., 3.), (1., 2., 3.), (1., 2., 3.)])
        >>> m2 = Matrix([(1, 2, 3), (1, 2, 3), (1, 2, 3)])
        >>> m1 == m2
        True
        >>> m1 = Matrix([(1, 2, 3), (1, 2, 3), (1, 2, 3)])
        >>> m2 = Matrix([(1, 2, 3), (3, 2, 3), (1, 2, 3)])
        >>> m1 == m2
        False
        """
        return all([b == x for (b, x) in zip(other.values, self.values)])

    def __sub__(self, other):
        """
        Substract two matrices.
        >>> Matrix([[3, 2, 1], [3, 2, 1], [3, 2, 3]]) - Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
        [[2, 0, -2], [2, 0, -2], [2, 0, 0]]
        """
        output = []
        for idx in range(len(self)):
            tmp = []
            for val_a, val_b in zip(self[idx], other[idx]):
                tmp.append(val_a - val_b)
            output.append(tmp[:])
        return output[:]

    def __truediv__(self, other):
        """
        #>>> Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]]) / Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
        A / B = A * A^-1
        """
        return NotImplementedError

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

        >>> m = Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
        >>> m.transposed.values
        [(1, 1, 1), (2, 2, 2), (3, 3, 3)]
        """
        return Matrix(list(zip(*self.values)))

    def transpose_alt(self):
        """
        Transposes the current matrix.
        >>> m = Matrix([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
        >>> m.transpose_alt().values
        [[1, 1, 1], [2, 2, 2], [3, 3, 3]]
        """
        rows = []
        for i in range(self.shape[0]):
            rows.append([r[i] for r in self.values])
        return Matrix(rows)

    def dot(self, other):
        """
        Dot product of two matrices.
        """
        new_a = []
        for i, col in enumerate(self.transposed.values):
            s = sum([v * o for v, o in zip(col, other)])
            new_a.append(s)
        return Matrix(new_a)

    @staticmethod
    def zero(m: int, n: int) -> 'Matrix':
        """
        Create zero matrix of dimension m,n

        >>> Matrix.zero(5, 3)
        | 0.0000  0.0000  0.0000|
        | 0.0000  0.0000  0.0000|
        | 0.0000  0.0000  0.0000|
        | 0.0000  0.0000  0.0000|
        | 0.0000  0.0000  0.0000|
        <BLANKLINE>
        """
        return Matrix([[0.0 for row in range(n)] for col in range(m)])

    @staticmethod
    def randmat(m: int, n: int) -> 'Matrix':
        """
        Create random matrix of dimension m, n

        #>>> Matrix.randmat(5, 3)
        |  2.000   1.000   2.000|
        |  8.000   1.000   2.000|
        |  3.000   5.000   1.000|
        |  4.000   1.000   9.000|
        |  2.000   8.000   4.000|
        <BLANKLINE>
        """
        return Matrix([[random.randint(1, 9) for y in range(n)] for x in range(m)])

    def cholesky(self) -> 'Matrix':
        """
        >>> m = Matrix([[25, 15, -5], [15, 18,  0], [-5,  0, 11]])
        >>> m.cholesky()
        | 5.0000  0.0000  0.0000|
        | 3.0000  3.0000  0.0000|
        |-1.0000  1.0000  3.0000|
        <BLANKLINE>
        """
        L = Matrix.zero(*self.shape)
        for i, (Ai, Li) in enumerate(zip(self.values, L.values)):
            for j, Lj in enumerate(L.values[:i + 1]):
                s = sum(Li[k] * Lj[k] for k in range(j))
                Li[j] = sqrt(Ai[i] - s) if (i == j) else (1.0 / Lj[j] * (Ai[j] - s))
        return L

    def norm(self):
        return self.det

    @property
    def inversed(self) -> 'Matrix':
        """
        Inversion of 3 × 3 matrices
         -0.812500    0.125000    0.187500
          0.125000   -0.250000    0.125000
          0.520833    0.125000   -0.145833

        >>> Matrix([ [1, 2, 3], [4, 1, 6], [7, 8, 9] ]).inversed
        |-0.8125  0.1250  0.1875|
        | 0.1250 -0.2500  0.1250|
        | 0.5208  0.1250 -0.1458|
        <BLANKLINE>
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
        a = self.values
        return (a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
                - a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2])
                + a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]))

    def power_iteration(self, num_simulations=10):
        """
        Eigenvalue algorythm from https://en.wikipedia.org/wiki/Power_iteration
        """
        # Ideally choose a random vector
        # To decrease the chance that our vector
        # Is orthogonal to the eigenvector
        b_k = Array.randarray(self.shape[1])

        for _ in range(num_simulations):
            # calculate the matrix-by-vector product Ab
            b_k1 = self.dot(b_k)

            # calculate the norm
            b_k1_norm = b_k1.det

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
        if not self == other:
            return 999.
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
            if sign == '-':
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


def flatten(lis):
    """
    Given a list, possibly nested to any level, return it flattened.
    From: http://code.activestate.com/recipes/578948-flattening-an-arbitrarily-nested-list-in-python/

    >>> flatten([['wer', 234, 'brdt5'], ['dfg'], [[21, 34,5], ['fhg', 4]]])
    ['wer', 234, 'brdt5', 'dfg', 21, 34, 5, 'fhg', 4]
    """
    new_lis = []
    for item in lis:
        if isinstance(item, list):
            new_lis.extend(flatten(item))
        else:
            new_lis.append(item)
    return new_lis


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
        return sqrt((a * dx) ** 2 + (b * dy) ** 2 + (c * dz) ** 2 + 2 * b * c * cos(al) * dy * dz + \
                    2 * dx * dz * a * c * cos(be) + 2 * dx * dy * a * b * cos(ga))
    else:
        return sqrt(dx ** 2 + dy ** 2 + dz ** 2)


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
    """
    return (a[0] - b[0],
            a[1] - b[1],
            a[2] - b[2])


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

    def __init__(self, a, b, c, alpha, beta, gamma):
        """
        Converts von fractional to cartesian and vice versa.

        >>> cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
        >>> coord = [-0.186843,   0.282708,   0.526803]
        >>> ort = OrthogonalMatrix(*cell)
        >>> [ round(x, 6) for x in (ort * Array(coord)).values ]
        [-2.741505, 5.909587, 10.775201]
        >>> c_coord = [-2.741505423999065, 5.909586678000002, 10.775200700893732]
        >>> [ round(x, 6) for x in (ort.inversed * Array(c_coord)).values]
        [-0.186843, 0.282708, 0.526803]
        """
        self.a, self.b, self.c = a, b, c
        self.V = vol_unitcell(a, b, c, alpha, beta, gamma)
        self.alpha = radians(alpha)
        self.beta = radians(beta)
        self.gamma = radians(gamma)
        phi = sqrt(1 - cos(self.alpha) ** 2 - cos(self.beta) ** 2 - cos(self.gamma) ** 2 +
                   2 * cos(self.alpha) * cos(self.beta) * cos(self.gamma))
        self.m = Matrix([[self.a, self.b * cos(self.gamma), self.c * cos(self.beta)],
                         [0, self.b * sin(self.gamma),
                          (self.c * (cos(self.alpha) - cos(self.beta) * cos(self.gamma)) / sin(self.gamma))],
                         [0, 0, self.V / (self.a * self.b * sin(self.gamma))]])

        # The inverted matrix:
        self.mi = \
            Matrix([[1.0 / self.a, -1.0 / (self.a * tan(self.gamma)),
                     (cos(self.alpha) * cos(self.gamma) - cos(self.beta)) / (self.a * phi * sin(self.gamma))],
                    [0.0, 1 / (self.b * sin(self.gamma)), (cos(self.beta) * cos(self.gamma) - cos(self.alpha)) /
                     self.b * phi * sin(self.gamma)],
                    [0.0, 0.0, sin(self.gamma) / (self.c * phi)]])

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
        return self.mi

    @property
    def transposed(self):
        return self.m.transposed

    @property
    def T(self):
        return self.m.transposed


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


def frac_to_cart(frac_coord: (list, tuple), cell: list) -> list:
    """
    Converts fractional coordinates to cartesian coodinates
    :param frac_coord: [float, float, float]
    :param cell:       [float, float, float, float, float, float]

    >>> cell = [10.5086, 20.9035, 20.5072, 90, 94.13, 90]
    >>> coord1 = [-0.186843,   0.282708,   0.526803]
    >>> print(frac_to_cart(coord1, cell))
    [-2.741505423999065, 5.909586678000002, 10.775200700893734]
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
    return [xc, yc, zc]


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
