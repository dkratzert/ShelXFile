# -*- encoding: utf-8 -*-
# möpß
#
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#
"""
This is a full implementation of the SHELXL file syntax. Additionally it is able to edit SHELX properties with Python.
The implementation is Python3-only and supports SHELXL after 2017 (You should not use old versions anyway).
"""
from __future__ import print_function

import os
import re
import sys
import textwrap
import time
from math import radians, cos, sin, sqrt

from shelxfile.dsrmath import atomic_distance, Matrix, frac_to_cart, subtract_vect, determinante, my_isnumeric

PROFILE = False
DEBUG = False

"""
TODO:

Needed for DSR:
- fit fragment without shelxl
------------------------------------------------------------
- make all __repr__() unwrapped strings and wrap lines during write_res_file()
- add remove_hydrogen_atoms(atom) method.
- shx.remove_all_H([list of atoms], or all)
- check if atoms in restraints are also in structure
- deleting atoms should also remove them from restraints  
- shx.update_weight
- shx.weight_difference
- bond list
- shx.atoms.angle(at1, at2, at3)
- shx.atoms.tors(at1, at2, at3, at4)
- shx.atom.change_type('xx')
- restraints involved with an atom should also be part of the atoms properties
- implement an add_afix(afixm, afixn, atoms, frag_fend=False, position=None, afix_options=None)
  default position is directly behind FVAR or FRAG/FEND if enabled
- will read in lst file after refinement to fill shx.lst_file properties.
- shx.move_in_part([list of atoms])
- shx.move_in_resi([list of atoms])
- shx.sort.residue(3) (low priority)
- shx.lst_file.residuals -> density, r-values, goof, movements, bad reflections
   bad restraints, shelx *** errors
- shx.unused_atom_name('C') -> get a carbon atom with unused number
- shx.sort -> sort file (low priority)
"""

SHX_CARDS = ('TITL', 'CELL', 'ZERR', 'LATT', 'SYMM', 'SFAC', 'UNIT', 'LIST', 'L.S.', 'CGLS',
             'BOND', 'FMAP', 'PLAN', 'TEMP', 'ACTA', 'CONF', 'SIMU', 'RIGU', 'WGHT', 'FVAR',
             'DELU', 'SAME', 'DISP', 'LAUE', 'REM ', 'MORE', 'TIME', 'END ', 'HKLF', 'OMIT',
             'SHEL', 'BASF', 'TWIN', 'EXTI', 'SWAT', 'HOPE', 'MERG', 'SPEC', 'RESI', 'MOVE',
             'ANIS', 'AFIX', 'HFIX', 'FRAG', 'FEND', 'EXYZ', 'EADP', 'EQIV', 'CONN', 'BIND',
             'FREE', 'DFIX', 'BUMP', 'SADI', 'CHIV', 'FLAT', 'DEFS', 'ISOR', 'NCSY', 'SUMP',
             'BLOC', 'DAMP', 'STIR', 'MPLA', 'RTAB', 'HTAB', 'SIZE', 'WPDB', 'GRID', 'MOLE',
             'XNPD', 'REST', 'CHAN', 'FLAP', 'RNUM', 'SOCC', 'PRIG', 'WIGL', 'RANG', 'TANG',
             'ADDA', 'STAG', 'NEUT', 'ABIN', 'ANSC', 'ANSR', 'NOTR', 'TWST', 'PART', 'DANG',
             'BEDE', 'LONE', 'REM', 'END')

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


def remove_line(reslist, linenum, rem=False, remove=False, frontspace=False):
    """
    removes a single line from the res file with tree different methods.
    The default is a space character in front of the line (frontspace).
    This removes the line in the next refinement cycle. "rem" writes rem
    in front of the line and "remove" clears the line.
    :param reslist: .res file list
    :param linenum: integer, line number
    :param rem:     True/False, activate comment with 'REM' in front
    :param remove:  True/False, remove the line
    :param frontspace: True/False, activate removing with a front space

    >>> reslist = ['foo', 'bar', 'baz']
    >>> remove_line(reslist, 1, rem=True)
    ['foo', 'rem bar', 'baz']
    >>> remove_line(reslist, 2, remove=True)
    ['foo', 'rem bar', '']
    """
    line = reslist[linenum]
    if rem:  # comment out with 'rem ' in front
        reslist[linenum] = 'rem ' + line
        if multiline_test(line):
            reslist[linenum + 1] = 'rem ' + reslist[linenum + 1]
    elif remove:  # really delete the line "linenum"
        if multiline_test(line):
            reslist[linenum] = ''
            reslist[linenum + 1] = ''
        else:
            reslist[linenum] = ''
    if frontspace:  # only put a space in front
        reslist[linenum] = ' ' + line
        if multiline_test(line):
            reslist[linenum + 1] = ' ' + reslist[linenum + 1]
    return reslist


def vol_tetrahedron(a, b, c, d, cell=None):
    """
    Returns the volume of a terahedron spanned by four points.

    No cell is needed for orthogonal coordinates.

    e.g. A = (3, 2, 1), B = (1, 2, 4), C = (4, 0, 3), D = (1, 1, 7)
            |u1 u2 u3|
    v = 1/6*|v1 v2 v3|
            |w1 w2 w3|
    AB = (1-3, 2-2, 4-1) = (-2, 0, 3)
    AC = ...
    AD = ...
    V = 1/6[u,v,w]
              |-2,  0, 3|
    [u,v,w] = | 1, -2, 2| = 24-3-12 = 5
              |-2, -1, 6|
    V = 1/6*5
    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    >>> a = (0.838817,   0.484526,   0.190081) # a ist um 0.01 ausgelenkt
    >>> b = (0.875251,   0.478410,   0.256955)
    >>> c = (0.789290,   0.456520,   0.301616)
    >>> d = (0.674054,   0.430194,   0.280727)
    >>> print('volume of Benzene ring atoms:')
    volume of Benzene ring atoms:
    >>> print(round(vol_tetrahedron(a, b, c, d, cell), 7))
    0.0633528
    """
    A = [float(i) for i in a]
    B = [float(i) for i in b]
    C = [float(i) for i in c]
    D = [float(i) for i in d]
    if cell:
        A = frac_to_cart(a, cell)
        B = frac_to_cart(b, cell)
        C = frac_to_cart(c, cell)
        D = frac_to_cart(d, cell)
    AB = subtract_vect(A, B)
    AC = subtract_vect(A, C)
    AD = subtract_vect(A, D)
    D = determinante([AB, AC, AD])
    return abs((D / 6))


def time_this_method(f):
    """
    Rather promitive way of timing a method. More advanced would be the profilehooks module.
    """
    if PROFILE:
        from functools import wraps

        @wraps(f)
        def wrapper(*args, **kwargs):
            t1 = time.clock()
            result = f(*args, **kwargs)
            t2 = time.clock()
            if PROFILE:
                print('Time for "{}": {:5.3} ms\n'.format(f.__name__ + '()', (t2 - t1) * 1000))
            return result
    else:
        wrapper = f
    return wrapper


def chunks(l: list, n: int) -> list:
    """
    returns successive n-sized chunks from l.
    >>> l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e', 'f']
    >>> chunks(l, 5)
    [[1, 2, 3, 4, 5], [6, 7, 8, 9, 0], ['a', 'b', 'c', 'd', 'e'], ['f']]
    >>> chunks(l, 1)
    [[1], [2], [3], [4], [5], [6], [7], [8], [9], [0], ['a'], ['b'], ['c'], ['d'], ['e'], ['f']]
    >>> chunks(l, 50)
    [[1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 'a', 'b', 'c', 'd', 'e', 'f']]
    """
    return [l[i:i + n] for i in range(0, len(l), n)]


def multiline_test(line: str) -> bool:
    """
    test if the current line is a multiline with "=" at the end
    :param line: 'O1 3 -0.01453 1.66590 0.10966 11.00 0.05 ='
    :type line: string
    >>> line = 'C1    1    0.278062    0.552051    0.832431    11.00000    0.02895    0.02285 ='
    >>> multiline_test(line)
    True
    >>> line = 'C1    1    0.278062    0.552051    0.832431    11.00000    0.05 '
    >>> multiline_test(line)
    False
    """
    if line.rfind('=') > -1:
        # A '=' character in a rem line is not a line break!
        if line.startswith("REM") and not ShelXlFile.dsr_regex.match(line):
            return False
        return True
    else:
        return False


def split_fvar_and_parameter(parameter: float) -> tuple:
    """
    Returns the free variable and value of a given parameter e.g. 30.5 for the occupancy.
    :return (fvar: int, value: float)

    >>> split_fvar_and_parameter(30.5)
    (3, 0.5)
    >>> split_fvar_and_parameter(31.0)
    (3, 1.0)
    >>> split_fvar_and_parameter(-30.5)
    (-3, 0.5)
    >>> split_fvar_and_parameter(11.0)
    (1, 1.0)
    >>> split_fvar_and_parameter(-11.0)
    (-1, 1.0)
    >>> split_fvar_and_parameter(-10.33333333)
    (-1, 0.33333333)
    """
    fvar = abs(int(str(parameter).split('.')[0])) // 10  # The free variable number e.g. 2
    value = abs(float(parameter)) % 10  # The value with which the free variable was multiplied e.g. 0.5
    if parameter < 0:
        fvar *= -1
    return fvar, round(value, 8)


class ParseOrderError(Exception):
    def __init__(self, arg=None):
        if DEBUG:
            if arg:
                print(arg)
            else:
                print("*** WRONG ODER of INSTRUCTIONS ***")


class ParseNumError(Exception):
    def __init__(self, arg=None):
        if DEBUG:
            if arg:
                print(arg)
            print("*** WRONG NUMBER OF NUMERICAL PARAMETERS ***")


class ParseParamError(Exception):
    def __init__(self, arg=None):
        if DEBUG:
            if arg:
                print(arg)
            print("*** WRONG NUMBER OF PARAMETERS ***")


class ParseUnknownParam(Exception):
    def __init__(self):
        if DEBUG:
            print("*** UNKNOWN PARAMETER ***")


def range_resolver(atoms_range: list, atom_names: list) -> list:
    """
    Resolves the atom names of ranges like "C1 > C5"
    and works for each restraint line separately.
    :param atoms_range: atoms with a range definition
    :param atom_names: names of atoms in the fragment
    >>> r = "C2 > C5".split()
    >>> atlist = 'C1 C2 C3 C4 C5'.split()
    >>> range_resolver(r, atlist)
    ['C2', 'C3', 'C4', 'C5']
    >>> r = "C2_2 > C5_2".split()
    >>> atlist = 'C1_1 C1_2 C2_2 C3_2 C4_2 C5_2'.split()
    >>> range_resolver(r, atlist)
    ['C2_2', 'C3_2', 'C4_2', 'C5_2']
    >>> r = "C2_1 > C5_1".split()
    >>> atlist = 'C1_1 C1_2 C2_2 C3_2 C4_2 C5_2'.split()
    >>> range_resolver(r, atlist) # doctest +ELLIPSIS
    Traceback (most recent call last):
     ...
    ValueError: 'C2_1' is not in list
    """
    # dict with lists of positions of the > or < sign:
    rightleft = {'>': [], '<': []}
    for rl in rightleft:
        for num, i in enumerate(atoms_range):
            i = i.upper()
            if rl == i:
                # fill the dictionary:
                rightleft[rl].append(num)
    for rl in rightleft:
        # for each sign:
        for i in rightleft[rl]:
            # for each position of < or >:
            if rl == '>':
                # forward range
                left = atom_names.index(atoms_range[i - 1]) + 1
                right = atom_names.index(atoms_range[i + 1])
                atoms_range[i:i + 1] = atom_names[left:right]
            else:
                # backward range
                left = atom_names.index(atoms_range[i - 1])
                right = atom_names.index(atoms_range[i + 1]) + 1
                names = atom_names[right:left]
                names.reverse()  # counting backwards
                atoms_range[i:i + 1] = names
    return atoms_range


class Restraints():
    """
    Base class for the list of restraints.
    """

    def __init__(self):
        self.restraints = []

    def append(self, restr):
        self.restraints.append(restr)

    def __iter__(self):
        for x in self.restraints:
            yield x

    def __getitem__(self, item):
        return self.restraints[item]

    def __repr__(self):
        if self.restraints:
            for x in self.restraints:
                return x
        else:
            return 'No Restraints in file.'


class Restraint():

    def __init__(self, spline: list, line_nums: list):
        """
        Base class for parsing restraints.
        TODO: resolve ranges like SADI_CCF3 O1 > F9
        Therefore, make method to get atoms of residue CCF3 and residue x<-number
        and between C21 ans C25 for atoms outside residues.
        """
        self.line_numbers = line_nums
        self.residue_class = ''
        self.residue_number = 0
        self.textline = ' '.join(spline)
        self.name = None
        self.atoms = []
        self.atoms_involved = []

    def parse_line(self, spline, pairs=False):
        self.spline = spline
        if '_' in spline[0]:
            self.name, suffix = spline[0].upper().split('_')
            if any([x.isalpha() for x in suffix]):
                self.residue_class = suffix
            else:
                # TODO: implement _+, _- and _*
                if '*' in suffix:
                    self.residue_number = suffix
                else:
                    self.residue_number = int(suffix)
        else:
            self.name = spline[0].upper()
        # Beware! DEFS changes only the non-defined default values:
        # DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
        if DEFS.active:
            if self.name == 'DFIX':
                self.s = DEFS.sd
            if self.name == 'SAME':
                self.s1 = DEFS.sd
                self.s2 = DEFS.sd * 2
            if self.name == 'SADI':
                self.s = DEFS.sd
            if self.name == 'CHIV':
                self.s = DEFS.sf
            if self.name == 'FLAT':
                self.s = DEFS.sf
            if self.name == 'DELU':
                self.s1 = DEFS.su
                self.s2 = DEFS.su
            if self.name == 'SIMU':
                self.s = DEFS.ss
                self.st = DEFS.ss * 2
        params = []
        atoms = []
        for x in spline[1:]:
            if my_isnumeric(x):
                params.append(float(x))
            else:
                atoms.append(x)
        if pairs:
            return params, chunks(atoms, 2)
        else:
            return params, atoms

    def paircheck(self):
        if not self.atoms:
            return
        if len(self.atoms[-1]) != 2:
            print('*** Wrong number of numerical parameters ***')
            print('Instruction: {}'.format(self.textline))
            # raise ParseNumError

    def __iter__(self):
        for x in self.textline.split():
            yield x

    def __repr__(self):
        return self.textline  # + ':' + '\n' + str(self.__dict__)

    def __str__(self):
        return self.textline

    def split(self):
        return self.textline.split()


class DEFS(Restraint):
    """
    DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
    Changes the *default* effective standard deviations for the following
    DFIX, SAME, SADI, CHIV, FLAT, DELU and SIMU restraints.
    """
    # keeps track if DEFS was previously activated:
    active = False
    sd = 0.02
    sf = 0.1
    su = 0.01
    ss = 0.04
    maxsof = 1

    def __init__(self, spline: list, line_nums: list):
        super(DEFS, self).__init__(spline, line_nums)
        DEFS.active = True
        p, _ = self.parse_line(spline)
        if _:
            raise ParseParamError
        if len(p) > 0:
            DEFS.sd = p[0]
        if len(p) > 1:
            DEFS.sf = p[1]
        if len(p) > 2:
            DEFS.su = p[2]
        if len(p) > 3:
            DEFS.ss = p[3]
        if len(p) > 4:
            DEFS.maxsof = p[4]

    @property
    def all(self):
        return DEFS.sd, DEFS.sf, DEFS.su, DEFS.ss, DEFS.maxsof


class NCSY(Restraint):
    """
    NCSY DN sd[0.1] su[0.05] atoms
    """

    def __init__(self, spline: list, line_nums: list):
        super(NCSY, self).__init__(spline, line_nums)
        self.sd = 0.1
        self.su = 0.05
        self.DN = None
        p, self.atoms = self.parse_line(spline, pairs=False)
        if len(p) > 0:
            self.DN = p[0]
        if len(p) > 1:
            self.sd = p[1]
        if len(p) > 2:
            self.su = p[2]
        if not self.DN:
            raise ParseNumError


class ISOR(Restraint):
    """
    ISOR s[0.1] st[0.2] atomnames
    """

    def __init__(self, spline: list, line_nums: list):
        super(ISOR, self).__init__(spline, line_nums)
        self.s = 0.1
        self.st = 0.2
        p, self.atoms = self.parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s = p[0]
        if len(p) > 1:
            self.st = p[1]


class FLAT(Restraint):
    """
    FLAT s[0.1] four or more atoms
    """

    def __init__(self, spline: list, line_nums: list):
        super(FLAT, self).__init__(spline, line_nums)
        self.s = 0.1
        p, self.atoms = self.parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s = p[0]
        # TODO: Have to resolve ranges first:
        # if len(self.atoms) < 4:
        #    raise ParseParamError


class BUMP(Restraint):
    """
    BUMP s [0.02]
    """

    def __init__(self, spline, line_nums):
        super(BUMP, self).__init__(spline, line_nums)
        self.s = 0.02
        p, _ = self.parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s = p[0]
        if _:
            raise ParseParamError


class DFIX(Restraint):
    """
    DFIX d s[0.02] atom pairs
    """

    def __init__(self, spline, line_nums):
        super(DFIX, self).__init__(spline, line_nums)
        self.s = 0.02
        p, self.atoms = self.parse_line(spline, pairs=True)
        if len(p) > 0:
            self.d = p[0]
        if len(p) > 1:
            self.s = p[1]
        self.paircheck()
        if not self.d:
            raise ParseNumError
        if 0.0001 < self.d <= self.s:  # Raise exception if d is smaller than s
            print('*** WRONG ODER of INSTRUCTIONS. d is smaller than s ***')
            print("{}".format(self.textline))


class DANG(Restraint):
    """
    DANG d s[0.04] atom pairs
    """

    def __init__(self, spline, line_nums):
        super(DANG, self).__init__(spline, line_nums)
        self.s = 0.04
        p, self.atoms = self.parse_line(spline, pairs=True)
        if len(p) > 0:
            self.d = p[0]
        if len(p) > 1:
            self.s = p[1]
        self.paircheck()
        if not self.d:
            raise ParseNumError
        if 0.0001 < self.d <= self.s:  # Raise exception if d is smaller than s
            raise ParseOrderError


class SADI(Restraint):
    """
    SADI s[0.02] pairs of atoms
    """

    def __init__(self, spline, line_nums):
        super(SADI, self).__init__(spline, line_nums)
        self.s = 0.02
        p, self.atoms = self.parse_line(spline, pairs=True)
        if len(p) > 0:
            self.s = p[0]
        self.paircheck()


class SAME(Restraint):
    """
    SAME s1[0.02] s2[0.04] atomnames
    """

    def __init__(self, spline, line_nums):
        super(SAME, self).__init__(spline, line_nums)
        self.s1 = 0.02
        self.s2 = 0.04
        p, self.atoms = self.parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s1 = p[0]
        if len(p) > 1:
            self.s2 = p[1]


class RIGU(Restraint):
    """
    RIGU s1[0.004] s2[0.004] atomnames
    """

    def __init__(self, spline: list, line_nums: list):
        super(RIGU, self).__init__(spline, line_nums)
        self.s1 = 0.004
        self.s2 = 0.004
        p, self.atoms = self.parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s1 = p[0]
        if len(p) > 1:
            self.s2 = p[1]


class SIMU(Restraint):
    """
    SIMU s[0.04] st[0.08] dmax[2.0] atomnames
    """

    def __init__(self, spline: list, line_nums: list):
        super(SIMU, self).__init__(spline, line_nums)
        self.s = 0.04
        self.st = 0.08
        self.dmax = 2.0
        p, self.atoms = self.parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s = p[0]
        if len(p) > 1:
            self.st = p[1]
        if len(p) > 2:
            self.dmax = p[2]


class DELU(Restraint):
    """
    DELU s1[0.01] s2[0.01] atomnames
    """

    def __init__(self, spline: list, line_nums: list):
        super(DELU, self).__init__(spline, line_nums)
        self.s1 = 0.01
        self.s2 = 0.01
        p, self.atoms = self.parse_line(spline, pairs=False)
        if len(p) > 0:
            self.s1 = p[0]
        if len(p) > 1:
            self.s2 = p[1]


class CHIV(Restraint):
    """
    CHIV V[0] s[0.1] atomnames
    """

    def __init__(self, spline: list, line_nums: list):
        super(CHIV, self).__init__(spline, line_nums)
        self.s = 0.1
        self.V = 0.0
        p, self.atoms = self.parse_line(spline, pairs=False)
        if len(p) > 0:
            self.V = p[0]
        if len(p) > 1:
            self.s = p[1]


class EADP(Restraint):
    """
    EADP atomnames
    """

    def __init__(self, spline: list, line_nums: list) -> None:
        super(EADP, self).__init__(spline, line_nums)
        _, self.atoms = self.parse_line(spline, pairs=False)


class EXYZ(Restraint):
    """
    EADP atomnames
    """

    def __init__(self, spline: list, line_nums: list) -> None:
        super(EXYZ, self).__init__(spline, line_nums)
        _, self.atoms = self.parse_line(spline, pairs=False)


class Command():
    """
    A class to parse all general commands except restraints.
    """

    def __init__(self, spline: list, line_nums: list):
        self.line_numbers = line_nums
        self.residue_class = ''
        self.residue_number = 0
        self.textline = ' '.join(spline)
        self.name = None
        self.atoms = []

    def parse_line(self, spline, intnums=False):
        """
        :param spline: Splitted shelxl line
        :param intnums: if numerical parameters should be integer
        :return: numerical parameters and words
        """
        if '_' in spline[0]:
            self.name, suffix = spline[0].upper().split('_')
            if any([x.isalpha() for x in suffix]):
                self.residue_class = suffix
            else:
                # TODO: implement _+, _- and _*
                self.residue_number = int(suffix)
        else:
            self.name = spline[0].upper()
        numparams = []
        words = []
        for x in spline[1:]:
            if str.isdigit(x[0]) or x[0] in '+-':
                if intnums:
                    numparams.append(int(x))
                else:
                    numparams.append(float(x))
            else:
                words.append(x)
        return numparams, words

    def split(self):
        return self.textline.split()

    def __str__(self):
        return self.textline


class DAMP(Command):
    """
    DAMP damp[0.7] limse[15]
    """
    def __init__(self, spline, line_nums):
        super(DAMP, self).__init__(spline, line_nums)
        values, _ = self.parse_line(spline, intnums=False)
        self.damp, self.limse = 0, 0
        if len(values) > 0:
            self.damp = values[0]
        if len(values) > 1:
            self.damp, self.limse = values

    def __repr__(self) -> str:
        if self.limse == 0:
            return "DAMP  {:,g}".format(self.damp)
        else:
            return "DAMP  {:,g} {:,g}".format(self.damp, self.limse)


class HFIX(Command):
    """
    HFIX mn U[#] d[#] atomnames
    """
    def __init__(self, spline: list, line_nums: list):
        super(HFIX, self).__init__(spline, line_nums)
        self.params, self.atoms = self.parse_line(spline, intnums=True)

    def __repr__(self):
        return "HFIX {} {}".format(" ".join([str(x) for x in self.params]) if self.params else '',
                                   " ".join(self.atoms) if self.atoms else '')


class ACTA(Command):
    """
    ACTA 2θfull[#]
    """
    def __init__(self, shx: 'ShelXlFile', spline: list, line_nums: list):
        super(ACTA, self).__init__(spline, line_nums)
        self.twotheta, _ = self.parse_line(spline)
        self.shx = shx

    def remove_acta_card(self):
        self.shx.delete_on_write.update([self.shx._reslist.index(self)])
        return self.textline.strip('\r\n')

    def __repr__(self):
        return "ACTA {}".format(self.twotheta if self.twotheta else '')


class HKLF(Command):
    """
    HKLF N[0] S[1] r11...r33[1 0 0 0 1 0 0 0 1] sm[1] m[0]
    """
    def __init__(self, spline: list, line_nums: list):
        super(HKLF, self).__init__(spline, line_nums)
        p, _ = self.parse_line(spline)
        self.n = 0
        self.s = 1
        self.matrix = [1, 0, 0, 0, 1, 0, 0, 0, 1]
        self.sm = 1
        self.m = 0
        if len(p) > 1:
            self.n = p[0]
        if len(p) > 2:
            self.s = p[1]
        if len(p) > 10:
            self.matrix = p[3:11]
        if len(p) > 11:
            self.sm = p[12]
        if len(p) > 12:
            self.m = p[13]

    def __repr__(self):
        return "HKLF {} {}  {}  {} {}".format(self.n, self.s, ' '.join([str(i) for i in self.matrix]), self.sm, self.m)


class SUMP(Command):
    """
    SUMP for linear equation eypressions with free variables.
    SUMP c sigma c1 m1 c2 m2 ...
    """

    def __init__(self, spline, line_nums):
        super(SUMP, self).__init__(spline, line_nums)
        p, _ = self.parse_line(spline)
        self.c = p.pop(0)
        self.fvars = {}
        self.sigma = p.pop(0)
        # this is to have integer free variables
        fvars = [int(x) for x in p[1::2]]
        times = [x for x in p[0::2]]
        self.fvars = [[x, y] for x, y in zip(times, fvars)]

    def __getitem__(self, item):
        return self.fvars[item]

    def __repr__(self):
        return self.textline


class SYMM(Command):
    """
    Container for a symm card.
    """

    def __init__(self, spline, line_nums):
        super(SYMM, self).__init__(spline, line_nums)
        self.symmcards = self.parse_line(spline)

    def parse_line(self, spline, intnums=False):
        symmcards = []
        line = ''.join(spline[1:])  # removes whitespace
        symmcards.append(line.split(','))
        return symmcards

    def __repr__(self):
        return "\n".join(["SYMM  " + "  ".join(x) for x in self.symmcards])
    
    def __str__(self):
        return "\n".join(["SYMM  " + "  ".join(x) for x in self.symmcards])


class FVAR():
    def __init__(self, number=1, value=0.0):
        """
        FVAR osf[1] free variables
        """
        self.fvar_value = value  # value
        self.number = number  # occurence inside of FVAR instructions
        self.usage = 1

    def __str__(self):
        return str(float(self.fvar_value))


class FVARs():
    def __init__(self, shx: 'ShelXlFile'):
        super(FVARs, self).__init__()
        self.fvars = []  # free variables
        self.shx = shx
        self._fvarline = 0

    def __iter__(self):
        """
        Must be defined for __repr__() to work.
        """
        for x in self.fvars:
            yield x

    def __getitem__(self, item: int) -> str:
        # SHELXL counts fvars from 1 to x:
        item = item - 1
        if item < 0:
            raise IndexError("*** Illegal free variable number ***")
        return self.fvars[item].fvar_value

    def __setitem__(self, key, fvar_value):
        self.fvars[key] = fvar_value

    def __len__(self) -> int:
        return len(self.fvars)

    def __str__(self) -> str:
        lines = chunks(self.as_stringlist, 7)
        fvars = [' '.join(i) for i in lines]
        fvars = ['FVAR ' + i for i in fvars]
        return "\n".join(fvars)

    @property
    def fvarline(self) -> int:
        self._fvarline = self.shx._reslist.index(self)
        return self._fvarline

    def set_free_variables(self, fvar: int, dummy_fvar: float = 0.5):
        """
        Inserts additional free variables according to the fvar number.
        #TODO: make this in FVARS()
        """
        if fvar > 99:
            print('*** SHELXL allows only 99 free variables! ***')
            raise ParseParamError
        varlen = len(self.fvars)
        difference = (abs(fvar) - varlen)
        if difference > 0:
            for n in range(int(difference)):
                fv = FVAR(varlen + n, dummy_fvar)
                self.fvars.append(fv)

    def append(self, fvar) -> None:
        self.fvars.append(fvar)

    def set_fvar_usage(self, fvarnum: int, times: int = 1) -> None:
        """
        Incerements the usage count of a free variable by times.
        """
        fvarnum = abs(fvarnum)
        if len(self.fvars) >= abs(fvarnum):
            self.fvars[fvarnum - 1].usage += times
        elif fvarnum > 1:
            print('*** Free variable {} is not defined but used! ***'.format(fvarnum))
            # raise Exception

    def get_fvar_usage(self, fvarnum):
        """
        Returns the usage (count) of a certain free variable.
        """
        try:
            usage = self.fvars[fvarnum - 1].usage
        except IndexError:
            return 0
        return usage

    def fvars_used(self):
        """
        Retruns a dictionary with the usage of all free variables.
        """
        used = {}
        for num, fv in enumerate(self.fvars):
            used[num + 1] = fv.usage
        return used

    @property
    def as_stringlist(self):
        return [str(x.fvar_value) for x in self.fvars]

    @property
    def line_number(self) -> list:
        return self.fvarline


class LSCycles():
    def __init__(self, shx: 'ShelXlFile', spline: list, line_number: int = 0):
        """
        L.S. nls[0] nrf[0] nextra[0]
        If nrf is positive, it is the number of these cycles that should be performed before applying ANIS.
        Negative nrf indicates which reflections should be ignored during the refinement but used instead for
        the calculation of free R-factors in the final structure factor summation.
        nextra is the number of additional parameters that were derived from the data when 'squeezing' the
        structure etc.
        """
        self.shx = shx
        self.cgls = False
        self.cycles = 0
        self.nrf = ''
        self.nextra = ''
        self.line_number = line_number  # line number in res file
        try:
            self.cycles = int(spline[1])
        except (IndexError, NameError):
            raise ParseNumError
        try:
            self.nrf = spline[2]
        except IndexError:
            pass
        try:
            self.nextra = spline[3]
        except IndexError:
            pass
        if spline[0].upper() == 'CGLS':
            self.cgls = True

    @property
    def number(self):
        return self.cycles

    def set_refine_cycles(self, number: int):
        """
        Sets the number of refinement cycles for the current res file.
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.cycles.set_refine_cycles(44)
        >>> shx._reslist[shx.cycles.line_number]
        L.S. 44
        """
        self.cycles = number
        #self.shx.reslist[self.shx.reslist.index(self)] = self.text.strip('\r\n')

    @property
    def text(self):
        """
        'CGLS 10 2 '
        """
        return self.__repr__()

    def __repr__(self):
        return '{} {} {} {}'.format('CGLS' if self.cgls else 'L.S.', self.cycles,
                                    self.nrf if self.nrf else '', self.nextra if self.nextra else '').strip()



class SFACTable():
    def __init__(self, shx: 'ShelXlFile'):
        """
        Holds the information of SFAC instructions. Either with default values and only elements
        SFAC elements
        or as explicit scattering factor in the form of an exponential series, followed by real and
        imaginary dispersion terms, linear absorption coefficient, covalent radius and atomic weight.

        SFAC elements  or  SFAC E a1 b1 a2 b2 a3 b3 a4 b4 c f' f" mu r wt
        """
        self.sfac_table = []
        self.shx = shx
        self.elements_list = []

    def __iter__(self):
        for x in self.sfac_table:
            yield x['element'].capitalize()

    def __repr__(self):
        regular = []
        exponential = []
        for sf in self.sfac_table:
            if not self.is_exp(sf):
                regular.append(sf['element'].capitalize())
            else:
                values = []
                for x in ['element', 'a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'a4', 'b4', 'c',
                          'fprime', 'fdprime', 'mu', 'r', 'wt']:
                    values.append(sf[x])
                exponential.append("SFAC " + "  ".join(values))
        sfac = "SFAC {:<4s}".format("  ".join(regular))
        if exponential:
            return "\n".join(exponential) + "\n" + sfac
        else:
            return sfac

    def __getitem__(self, index: int):
        """
        Returns the n-th element in the sfac table, beginning with 1.
        """
        if index == 0:
            raise IndexError
        if index < 0:
            index = len(self.sfac_table) + index + 1
        return self.sfac_table[index - 1]['element'].capitalize()


    def parse_element_line(self, spline: list):
        """
        Adds a new SFAC card to the list of cards.
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.sfac_table
        SFAC C  H  O  F  Al  Ga
        >>> shx.unit
        UNIT 1  2  3  4  5  6
        >>> shx.sfac_table.add_element('Au')
        >>> shx.sfac_table
        SFAC C  H  O  F  Al  Ga  Au
        >>> shx.unit
        UNIT 1  2  3  4  5  6  1
        """
        if not ''.join(spline[1:]).isalpha():  # joining, because space is not alphabetical
            # Excplicit with all values
            sfdic = {}
            for n, x in enumerate(['element', 'a1', 'b1', 'a2', 'b2', 'a3', 'b3', 'a4', 'b4', 'c',
                                   'fprime', 'fdprime', 'mu', 'r', 'wt']):
                if n == 0:
                    self.elements_list.append(spline[n + 1].upper())
                try:
                    sfdic[x] = spline[n + 1]
                except IndexError:
                    raise ParseNumError()
            self.sfac_table.append(sfdic)
        else:
            # Just the elements
            for x in spline[1:]:
                self.elements_list.append(x.upper())
                self.sfac_table.append({'element': x.upper()})

    def has_element(self, element):
        return element.upper() in self.elements_list

    @staticmethod
    def is_exp(item):
        return 'a1' in item

    def add_element(self, element: str):
        """
        Adds an element to the SFAC list. 
        """
        if self.has_element(element):
            return
        self.elements_list.append(element.upper())
        self.sfac_table.append({'element': element.upper(), 'line_number': None})
        self.shx.unit.add_number(1)

    def remove_element(self, element: str):
        del self.sfac_table[self.shx.elem2sfac(element.upper()) - 1]
        del self.elements_list[self.elements_list.index(element.upper())]


class UNIT(Command):
    """
    UNIT n1 n2 ...
    """
    def __init__(self, spline: list, line_nums: list):
        super(UNIT, self).__init__(spline, line_nums)
        self.values, _ = self.parse_line(spline)

    def add_number(self, number: float):
        self.values.append(number)

    def __iter__(self):
        yield [x for x in self.values]

    def __repr__(self):
        return "UNIT " + "  ".join(["{:,g}".format(x) for x in self.values])

    def __setitem__(self, key, value):
        self.values[key] = value

    def __getitem__(self, item):
        return self.values[item]

    def __add__(self, other):
        self.values.append(other)


class BASF(Command):
    """
    BASF scale factors
    BASF can occour in multiple lines.
    """

    def __init__(self, spline: list, line_numbers: list):
        super(BASF, self).__init__(spline, line_numbers)
        self.scale_factors, _ = self.parse_line(spline)
        del self.atoms

    def __iter__(self):
        yield self.scale_factors


class TWIN(Command):
    """
    TWIN 3x3 matrix [-1 0 0 0 -1 0 0 0 -1] N[2]
    +N     -N  m = |N|
    m-1 to 2m-1
    m-1   (2*abs(m)/2)-1
    """

    def __init__(self, spline: list, line_nums: list):
        super(TWIN, self).__init__(spline, line_nums)
        self.matrix = [-1, 0, 0, 0, -1, 0, 0, 0, -1]
        self.allowed_N = 2
        self.n_value = 2
        if len(spline) > 1:
            p, _ = self.parse_line(spline, intnums=False)
            if len(p) == 9:
                self.matrix = p
            elif len(p) == 10:
                self.matrix = p[:9]
                self.n_value = int(p[9])
            else:
                raise ParseNumError("*** Check TWIN instruction. ***")
        m = abs(self.n_value) / 2
        if self.n_value > 0:
            self.allowed_N = abs(self.n_value) - 1
        else:
            self.allowed_N = (2 * m) - 1


class WGHT(Command):
    """
    The weighting scheme is defined as follows:
    w = q / [ σ²(Fo²) + (a*P)² + b*P + d + e*sin(θ)/$lambda; ]

    WGHT a[0.1] b[0] c[0] d[0] e[0] f[.33333]
    Usually only WGHT a b
    """

    def __init__(self, spline: list, line_nums: list):
        super(WGHT, self).__init__(spline, line_nums)
        self.a = 0.1
        self.b = 0.0
        self.c = 0.0
        self.d = 0.0
        self.e = 0.0
        self.f = 0.33333
        self.line_numbers = line_nums[0]
        p, _ = self.parse_line(spline)
        if len(p) > 0:
            self.a = p[0]
        if len(p) > 1:
            self.b = p[1]
        if len(p) > 2:
            self.c = p[2]
        if len(p) > 3:
            self.d = p[3]
        if len(p) > 4:
            self.e = p[4]
        if len(p) > 5:
            self.f = p[5]

    def __repr__(self):
        wght = 'WGHT {} {}'.format(self.a, self.b)
        # It is very unlikely that someone changes other parameter than a and b:
        if (self.c + self.d + self.e + self.f) != 0.33333:
            wght += ' {} {} {} {}'.format(self.c, self.d, self.e, self.f)
        return wght


class REM(Command):
    """
    Parses REM lines
    """

    def __init__(self, spline: list, line_nums: list) -> None:
        super(REM, self).__init__(spline, line_nums)

    def __repr__(self):
        return self.textline

    def __str__(self):
        return self.textline


class BOND(Command):
    """
    BOND atomnames
    """

    def __init__(self, spline: list, line_nums: list) -> None:
        super(BOND, self).__init__(spline, line_nums)
        _, self.atoms = self.parse_line(spline)

    def __repr__(self) -> str:
        return self.textline.strip('\n\r')


class Atoms():
    """
    All atoms from a SHELXL file with their properties.
    """
    def __init__(self, shx: 'ShelXlFile'):
        self.shx = shx
        self.atoms = []
        self.atomsdict = {}
        self.nameslist = []

    def append(self, atom: 'Atom') -> None:
        """
        Adds a new atom to the list of atoms. Using append is essential.
        """
        self.atoms.append(atom)
        name = atom.name + '_{}'.format(atom.resinum)
        self.atomsdict[name] = atom
        self.nameslist.append(name.upper())

    def __repr__(self):
        if self.atoms:
            for x in self.atoms:
                return x
        else:
            return 'No Atoms in file.'

    def __iter__(self):
        for x in self.atoms:
            yield x

    def __getitem__(self, item: int) -> None:
        return self.atoms[item]

    def __len__(self) -> int:
        return len(self.atoms)

    def __delitem__(self, key):
        """
        Delete an atom by its atomid:
        del atoms[4]
        """
        for n, at in enumerate(self.atoms):
            if key == at.atomid:
                if DEBUG:
                    print("deleting atom", at.name)
                del self.atoms[n]
                del self.atomsdict[at.name+'_{}'.format(at.resinum)]
                del self.nameslist[self.nameslist.index(at.fullname.upper())]
                self.shx.delete_on_write.update(at.line_numbers)

    @property
    def number(self) -> int:
        """
        The number of atoms in the current SHELX file.
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.atoms.number
        148
        """
        return len(self.atoms)

    def get_atom_by_id(self, aid: int) -> 'Atom':
        """
        Returns the atom objext with atomId id.
        """
        for a in self.atoms:
            if aid == a.atomid:
                return a

    def has_atom(self, atom_name: str) -> bool:
        """
        Returns true if shelx file has atom.
        
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.atoms.has_atom('Al1')
        True
        >>> shx.atoms.has_atom('Al1_0')
        True
        >>> shx.atoms.has_atom('Al2_0')
        False
        """
        if '_' not in atom_name:
            atom_name += '_0'
        if atom_name.upper() in self.nameslist:
            return True
        else:
            return False

    def get_atom_by_name(self, atom_name: str) -> 'Atom' or None:
        """
        Returns an Atom object using an atom name with residue number like C1, C1_0, F2_4, etc.
        C1 means atom C1 in residue 0.
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.atoms.get_atom_by_name('Al1')
        ID:70
        """
        if '_' not in atom_name:
            atom_name += '_0'
        try:
            at = self.atomsdict[atom_name.upper()]
        except KeyError:
            print("Atom {} not found in atom list.".format(atom_name))
            return None
        return at

    def get_all_atomcoordinates(self) -> dict:
        """
        Returns a dictionary {'C1': ['1.123', '0.7456', '3.245'], 'C2_2': ...}
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.atoms.get_all_atomcoordinates() # doctest: +ELLIPSIS
        {'O1_4': [0.074835, 0.238436, 0.402457], 'C1_4': [0.028576, 0.234542, 0.337234], ...}
        """
        atdict = {}
        for at in self.atoms:
            #if at.qpeak:
            #    atdict[at.name] = at.frac_coords
            #else:
            atdict[at.name.upper() + '_' + str(at.resinum)] = at.frac_coords
        return atdict

    def get_frag_fend_atoms(self) -> list:
        """
        Returns a list of atoms with cartesian coordinates. Atom names and sfac are ignored. They come from AFIX 17x.
        [[0.5316439256202359, 7.037351406500001, 10.112963255220803],
        [-1.7511017452002604, 5.461541059000001, 10.01187984858907]]
        """
        atoms = []
        for at in self.atoms:
            if at.frag_atom:
                atoms.append([at.xc, at.yc, at.zc])
        return atoms

    @property
    def residues(self) -> list:
        """
        Returns a list of the residue numbers in the shelx file.
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.atoms.residues
        [0, 1, 2, 3, 4]
        """
        return list({x.resinum for x in self.atoms})

    @property
    def q_peaks(self) -> list:
        r"""
        Returns a list of q-peaks in the file.
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.atoms.q_peaks # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
        [ID:330, ID:331, ID:332, ID:333, ID:334, ID:335, ID:336, ID:337, ID:338, ID:339, ID:340, ID:341, ID:342, ID:343, ID:344, ID:345, ID:346, ID:347, ID:348, ID:349]
        """
        return [x for x in self.atoms if x.qpeak]

    def distance(self, atom1: str, atom2: str) -> float:
        """
        Calculates the (shortest) distance of two atoms given as text names e.g. C1_3.
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.atoms.distance('C22', 'C23')
        1.4016798832482247
        """
        a1 = self.get_atom_by_name(atom1)
        a2 = self.get_atom_by_name(atom2)
        try:
            return atomic_distance([a1.xc, a1.yc, a1.zc], [a2.xc, a2.yc, a2.zc])
        except AttributeError:
            return 0.0

    def atoms_in_class(self, name: str) -> list:
        """
        Returns a list of atoms in residue class 'name'
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.atoms.atoms_in_class('CCF3')
        ['O1', 'C1', 'C2', 'F1', 'F2', 'F3', 'C3', 'F4', 'F5', 'F6', 'C4', 'F7', 'F8', 'F9']
        """
        atoms = []
        for x in self.atoms:
            if x.resiclass == name:
                if x.name not in atoms:
                    atoms.append(x.name)
        return atoms


class Atom(Atoms):
    """
    An Opbect holding all Properties of a shelxl atom plus some extra information like
    kartesian coordinates and element type.
    """
    #                name    sfac     x         y        z       occ      u11      u12 ...
    anisatomstr = '{:<4.4s}{:>3}{:>12.6f}{:>12.6f}{:>12.6f}{:>12.5f}{:>11.5f}{:>11.5f}' \
                  ' {:>12.5f}{:>11.5f}{:>11.5f}{:>11.5f}'
    #               name    sfac     x         y         z         occ      u11
    isoatomstr = '{:<5.5s} {:<3}{:>10.6f}  {:>10.6f}  {:>9.6f}  {:>9.5f}  {:>9.5f}'
    fragatomstr = '{:<5.5s} {:>10.6f}  {:>10.6f}  {:>9.6f}'

    def __init__(self, shelx: 'ShelXlFile', spline: list, line_nums: list, line_number: int, part: int = 0,
                 afix: int = 0, residict: dict = None, sof: float = 0) -> None:
        super(Atom, self).__init__(shelx)
        self._line_number = line_number
        self._lines = line_nums
        self.sfac_num = None
        self.name = None      # Name without residue number like "C1"
        self.fullname = None  # Name including residue nimber like "C1_2"
        # Site occupation factor including free variable like 31.0
        self.sof = None
        self.atomid = line_number
        self.shx = shelx
        self.element = None
        # fractional coordinates:
        self.x = None
        self.y = None
        self.z = None
        # cartesian coordinates:
        self.xc = None
        self.yc = None
        self.zc = None
        self.frag_atom = False
        self.previous_non_h = self.shx.non_h
        if self.element not in ['H', 'D']:
            self.shx.non_h = line_number
        else:
            self.shx.non_h = None
        if sof:
            # sof defined from outside e.g. by PART 1 31
            self.sof = float(sof)
        elif len(spline) > 5:
            self.sof = float(spline[5])
        else:
            self.sof = 11.0
        # Only the occupancy of the atom *without* the free variable like 0.5
        fvar, self.occupancy = split_fvar_and_parameter(self.sof)
        self.shx.fvars.set_fvar_usage(fvar)
        self.uvals = [0.04]  # [u11 u12 u13 u21 u22 u23]
        self.resiclass = residict['class']
        if not residict['number']:
            self.resinum = 0  # all other atoms are residue 0
        else:
            self.resinum = residict['number']
        self.chain_id = residict['ID']
        self.part = part
        self.afix = afix
        self.qpeak = False
        self.peak_height = 0.0
        self.cell = shelx.cell
        self.parse_line(spline)
        if self.shx.frag:
            self.afix = int(self.shx.frag[0])  # The FRAG AFIX fit code like 176 (must be greater 16)
        if self.shx.anis:
            self.parse_anis()
        for n, u in enumerate(self.uvals):
            if abs(u) > 4.0:
                fvar, uval = split_fvar_and_parameter(u)
                self.uvals[n] = uval
                self.shx.fvars.set_fvar_usage(fvar)

    def parse_anis(self):
        """
        Parses the ANIS card. It can be either ANIS, ANIS name(s) or ANIS number.
        # TODO: Test if ANIS $CL works and if ANIS_* $C works
        """
        try:
            # ANIS with a number as parameter
            if int(self.shx.anis[1]) > 0:
                self.shx.anis[1] -= 1
                if len(self.uvals) < 6:
                    self.uvals = [0.04, 0.0, 0.0, 0.0, 0.0, 0.0]
        except (TypeError, KeyError, ValueError, IndexError):
            # ANIS with a list of atoms
            if len(self.shx.anis) > 1:
                # if '_' in self.shx.anis[0]:
                #    resinum = self.shx.anis[0].upper().split('_')[1]
                for x in self.shx.anis[1:]:
                    if '_' in x:
                        name, resinum = x.upper().split('_')
                    else:
                        name = x.upper()
                        resinum = 0
                    if self.name == name and (int(self.resinum) == int(resinum) or resinum == '*'):
                        self.uvals = [0.04, 0.0, 0.0, 0.0, 0.0, 0.0]
                        self.shx.anis.pop()
                        if self.shx.anis == ['ANIS']:
                            # ANIS finished, deactivating again:
                            self.shx.anis = None
                    if x.startswith('$'):
                        if name[1:].upper() == self.element \
                                and (int(self.resinum) == int(resinum) or resinum == '*'):
                            self.uvals = [0.04, 0.0, 0.0, 0.0, 0.0, 0.0]
                            self.shx.anis.pop()
                            # TODO: This is a mess. Test and fix all sorts of ANIS possibilities.
                            if self.shx.anis == ['ANIS']:
                                # ANIS finished, deactivating again:
                                self.shx.anis = None
            # ANIS for all atoms
            else:
                if len(self.uvals) < 6:
                    self.uvals = [0.04, 0.0, 0.0, 0.0, 0.0, 0.0]

    def parse_line(self, line):
        self.name = line[0][:4]
        self.fullname = self.name+'_{}'.format(self.resinum)
        uvals = [float(x) for x in line[6:12]]
        try:
            x, y, z = [float(x) for x in line[2:5]]
        except ValueError as e:
            if DEBUG:
                print(e, 'Line:', self.line_numbers[-1])
            raise ParseUnknownParam
        if abs(x) > 4:
            fvar, x = split_fvar_and_parameter(x)
            self.shx.fvars.set_fvar_usage(fvar)
        if abs(y) > 4:
            fvar, x = split_fvar_and_parameter(y)
            self.shx.fvars.set_fvar_usage(fvar)
        if abs(z) > 4:
            fvar, x = split_fvar_and_parameter(z)
            self.shx.fvars.set_fvar_usage(fvar)
        self.x = x
        self.y = y
        self.z = z
        self.uvals = uvals
        if len(self.uvals) == 2:
            self.peak_height = uvals.pop()
        if self.shx.end:  # After 'END' can only be Q-peaks!
            self.qpeak = True
        self.sfac_num = int(line[1])
        self.element = self.shx.sfac2elem(self.sfac_num).upper()
        self.xc, self.yc, self.zc = frac_to_cart([self.x, self.y, self.z], self.cell)

    def __repr__(self) -> str:
        return 'ID:' + str(self.atomid)

    def __str__(self) -> str:
        """
        Returns a text line of the Atom with SHELXL syntax.
        :return: SHELX-formated atom string
        """
        if self.afix and self.shx.frag:
            # An atom from a FRAG/FEND instruction
            return Atom.fragatomstr.format(self.name, self.x, self.y, self.z)
        else:
            if len(self.uvals) > 2:
                # anisotropic atom
                try:
                    return Atom.anisatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof, *self.uvals)
                except(IndexError):
                    return 'REM Error in U values.'
            else:
                # isotropic atom
                try:
                    return Atom.isoatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof, *self.uvals)
                except(IndexError):
                    return Atom.isoatomstr.format(self.name, self.sfac_num, self.x, self.y, self.z, self.sof, 0.04)

    @property
    def line_numbers(self) -> list:
        return self._lines

    @line_numbers.setter
    def line_numbers(self, value: list):
        self._lines = value

    @property
    def frac_coords(self):
        return [self.x, self.y, self.z]

    @property
    def cart_coords(self):
        return [self.xc, self.yc, self.zc]

    def delete(self):
        del self.shx.atoms[self.atomid]

    def to_isotropic(self) -> None:
        """
        Makes the current atom isotropic.
        """
        self.uvals = [0.04]
        # TODO: Check if this works
        #if len(self.line_numbers) > 1:
        #    self.shx.delete_on_write.update([self.line_numbers[-1]])
        #    self.shx.reslist[self.line_numbers[0]] = self.shx.reslist[self.line_numbers[0]].strip('=')

    '''
    def __eq__(self, other) -> bool:
        """
        Returns True if two atoms are of same name, part and residue.
        """
        if isinstance(other, str):
            return self.__repr__() == other
        if self.atomid == other.atomid:
            return True
        else:
            return False
    '''

    def find_atoms_around(self, dist=1.2, only_part=0) -> list:
        """
        Finds atoms around the current atom.

        >>> file = './p21c.res'
        >>> shx = ShelXlFile(os.path.normpath(file))
        >>> at = shx.atoms.get_atom_by_name('Al1')
        >>> found = at.find_atoms_around(2)
        >>> [n.atomid for n in found]
        [72, 74]
        >>> del shx.atoms[17]

        """
        found = []
        for at in self.shx.atoms:
            if atomic_distance([self.x, self.y, self.z], [at.x, at.y, at.z], self.cell) < dist:
                # Not the atom itselv:
                if not self == at:
                    # only in special part and no q-peaks:
                    if at.part == only_part and not at.qpeak:
                        found.append(at)
        return found


class TextLine:
    def __init__(self, initdata):
        """
        #>>> t = TextLine('foo')
        #>>> t.id
        """
        self.data = initdata
        self.next = None
        self.id = time.time()

    def get_data(self):
        return self.data

    def get_next(self):
        return self.next

    def set_data(self, newdata):
        self.data = newdata

    def set_next(self, newnext):
        self.next = newnext


class ResList():
    """
    Contains the lines of the res file as unordered list.
    """
    def __init__(self):
        """
        >>> res = ResList()
        >>> res
        <BLANKLINE>
        >>> res.append('zweiter')
        >>> res.add('erster')
        >>> res.add('dritter')
        >>> res.append('vierter')
        >>> res.size
        4
        >>> res.tail.get_data()
        'vierter'
        >>> res
        dritter
        erster
        zweiter
        vierter
        """
        self.head = None
        self.tail = None
        self.size = 0

    def __repr__(self):
        current = self.head
        datalist = []
        if not self.head:
            return ''
        while current != self.tail:
            datalist.append(current.get_data())
            current = current.get_next()
        datalist.append(self.tail.get_data())
        return "\n".join(datalist)

    def is_empty(self):
        return self.head is None

    def add(self, item):
        temp = TextLine(item)
        temp.set_next(self.head)
        self.head = temp
        self.size += 1

    def search(self,item):
        current = self.head
        found = False
        while current is not None and not found:
            if current.get_data() == item:
                found = True
            else:
                current = current.get_next()
        return found

    def remove(self, item):
        current = self.head
        previous = None
        found = False
        while not found:
            if current.get_data() == item:
                found = True
            else:
                previous = current
                current = current.get_next()
        if previous is None:
            self.head = current.get_next()
        else:
            previous.setNext(current.get_next())
        self.size -= 1

    def append(self, item):
        temp = TextLine(item)
        if not self.head:
            self.head = temp
        else:
            last = self.tail
            last.set_next(temp)
        self.tail = temp
        self.size += 1


class ShelXlFile():
    """
    Class for data from a SHELXL res file. Includes Atoms, cards and unit cell.
    """
    delete_on_write = None
    atoms = None
    sump = []
    dsr_regex = re.compile(r'^rem\s+DSR\s+(PUT|REPLACE).*', re.IGNORECASE)
    fvars = None
    sfac_table = None
    _reslist = None
    restraints = None
    dsrlines = None
    symmcards = None

    def __init__(self: 'ShelXlFile', resfile: str):
        """
        Reads the shelx file and extracts information.
        #TODO: Exchange reslist items with SHELX objects like SFAC or ATOM. Add Reslist() class.
        #TODO: line number of objects are then inside each object self.sxh.index(self)
        :param resfile: file path
        """
        self.shelx_max_line_length = 79  # maximum character lenth per line in SHELXL
        self.nohkl = False
        self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.V = None, None, None, None, None, None, None
        self.ansc = None
        self.abin = None
        self.acta = None
        self.fmap = None
        self.xnpd = None
        self.wpdb = None
        self.wigl = None
        self.temp = 20
        self.swat = None
        self.stir = None
        self.spec = None
        self.twst = None
        self.plan = None
        self.prig = None
        self.merg = None
        self.more = None
        self.move = None
        self.defs = None
        self.zerr = None
        self.wght = None
        self.frag = None
        self.twin = None
        self.basf = None
        self.latt = None
        self.anis = None
        self.damp = None
        self.unit = None
        self.sump = []
        self.end = False
        self.maxsof = 1.0
        self.commands = []
        self.size = {}
        self.htab = []
        self.shel = []
        self.mpla = []
        self.rtab = []
        self.omit = []
        self.hklf = []
        self.grid = []
        self.free = []
        self.titl = ""
        self.exti = 0
        self.eqiv = []
        self.disp = []
        self.conn = []
        self.conv = []
        self.bind = []
        self.ansr = 0.001
        self.bloc = []
        self.cell = []
        self.dsrlines = []
        self.dsrline_nums = []
        self.symmcards = []
        self.hfixes = []
        self.Z = 1
        self.rem = []
        self.indexes = {}
        self.atoms = Atoms(self)
        self.fvars = FVARs(self)
        self.restraints = Restraints()
        self.sfac_table = SFACTable(self)
        self.delete_on_write = set()
        self.wavelen = None
        self.global_sadi = None
        self.cycles = None
        self.list = 0
        self.theta_full = 0
        self.non_h = None
        self.error_line_num = -1  # Only used to tell the line number during an exception.
        self.restrdict = {}
        self.resfile = os.path.abspath(resfile)
        if DEBUG:
            print('Resfile is:', self.resfile)
        try:
            self._reslist = self.read_file_to_list(self.resfile)
        except UnicodeDecodeError:
            if DEBUG:
                print('*** Unable to read file', self.resfile, '***')
            return
        try:
            self.parse_shx_file()
            pass
        except Exception as e:
            if DEBUG:
                print(e)
                print("*** Syntax error found in file {}, line {} ***".format(self.resfile, self.error_line_num + 1))
            if DEBUG:
                raise
            else:
                return
        else:
            self.run_after_parse()


    @time_this_method
    def parse_shx_file(self):
        """
        Extracts the atoms and other information from the res file.

        line is upper() after multiline_test()
        spline is as in .res file.
        """
        lastcard = ''
        part = False
        partnum = 0
        resi = False
        residict = {'class': '', 'number': 0, 'ID': ''}
        sof = 0
        afix = False
        afixnum = 0
        fvarnum = 1
        resinull = re.compile(r'^RESI\s+0')
        partnull = re.compile(r'^PART\s+0')
        afixnull = re.compile(r'^AFIX\s+0')
        for line_num, line in enumerate(self._reslist):
            self.error_line_num = line_num  # For exception during parsing.
            list_of_lines = [line_num]  # list of lines where a card appears, e.g. for atoms with two lines
            if line[:1] == ' ' or line == '':
                continue
            if not self.titl and line[:4] == 'TITL':
                # TITL[]  ->  = and ! can be part of the TITL!
                self.titl = line[5:76]
                lastcard = 'TITL'
                continue
            wrapindex = 0
            # This while loop makes wrapped lines look like they are not wrapped. The following lines are then
            # beginning with a space character and thus are ignored. The 'lines' list holds the line nnumbers where
            # 'line' is located ([line_num]) plus the wrapped lines.
            while multiline_test(self._reslist[line_num + wrapindex]):
                # Glue together the two lines wrapped with "=":
                wrapindex += 1
                line = line.rpartition('=')[0] + self._reslist[line_num + wrapindex]
                self.delete_on_write.update([line_num + wrapindex])
                list_of_lines.append(line_num + wrapindex)  # list containing the lines of a multiline command
            # The current line splitted:
            spline = line.split('!')[0].split()  # Ignore comments with "!", see how this performes
            # The current line as string:
            line = line.upper().split('!')[0]  # Ignore comments with "!", see how this performes
            # RESI class[ ] number[0] alias
            if resinull.match(line) and resi:  # RESI 0
                resi = False
                residict['number'] = 0
                residict['class'] = ''
                residict['ID'] = ''
                continue
            if resinull.match(line) and not resi:
                # A second RESI 0
                continue
            if line.startswith(('END', 'HKLF', 'RESI')) and resi:
                self._reslist.insert(line_num, "RESI 0")
                if DEBUG:
                    print('RESI in line {} was not closed'.format(line_num + 1))
                resi = False
                continue
            if line.startswith('RESI') and not resinull.match(line):
                resi = True
                residict = self.get_resi_definition_dict(spline[1:])
                continue
            # Now collect the part:
            # PART n sof
            if partnull.match(line) and part:  # PART 0
                part = False
                partnum = 0
                sof = 0
                continue
            if partnull.match(line) and not part:
                # A second PART 0
                continue
            if line.startswith(('END', 'HKLF', 'PART')) and part:
                self._reslist.insert(line_num, "PART 0")
                if DEBUG:
                    print('PART in line {} was not closed'.format(line_num + 1))
                part = False
                continue
            if line.startswith('PART') and not partnull.match(line):
                part = True
                partnum, sof = self.get_partnumber(line)
                continue
            # Now collect the AFIXes:
            # AFIX mn d[#] sof[11] U[10.08]
            if afixnull.match(line) and afix:  # AFIX 0
                afix = False
                afixnum = 0
                continue
            if afixnull.match(line) and not afix:
                # A second AFIX 0
                continue
            if line.startswith(('END', 'HKLF', 'AFIX')) and afix:
                self._reslist.insert(line_num, "AFIX 0")
                if DEBUG:
                    print('AFIX in line {} was not closed'.format(line_num + 1))
                afix = False
                continue
            elif line.startswith('AFIX') and not afixnull.match(line):
                # TODO: if afixnum > 17x: compare atoms count in frag and afix
                afix = True
                # TODO: only afixnum and sof are used at the moment. Use also u and d for AFIX():
                afixnum, d, sof, u = self.get_afix_numbers(spline, line_num)
                continue
            elif self.is_atom(line):
                # A SHELXL atom:
                # F9    4    0.395366   0.177026   0.601546  21.00000   0.03231  ( 0.03248 =
                #            0.03649  -0.00522  -0.01212   0.00157 )
                a = Atom(self, spline, list_of_lines, line_num, part=partnum, afix=afixnum, residict=residict, sof=sof)
                self.append_card(self.atoms, a, line_num)
                continue
            elif self.end:
                # Prevents e.g. parsing of second WGHT after END:
                continue
            elif line[:4] == 'SADI':
                # SADI s[0.02] pairs of atoms
                # or SADI
                if len(spline) == 1:
                    self.global_sadi = line_num
                self.append_card(self.restraints, SADI(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'DFIX':
                # DFIX d s[0.02] atom pairs
                self.append_card(self.restraints, DFIX(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'SIMU':
                # SIMU s[0.04] st[0.08] dmax[2.0] atomnames
                self.append_card(self.restraints, SIMU(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'DELU':
                # DELU s1[0.01] s2[0.01] atomnames
                self.append_card(self.restraints, DELU(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'RIGU':
                # RIGU s1[0.004] s2[0.004] atomnames
                self.append_card(self.restraints, RIGU(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'BASF':
                # BASF scale factors
                self.append_card(self.restraints, BASF(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'HFIX':
                # HFIX mn U[#] d[#] atomnames
                self.append_card(self.hfixes, HFIX(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'DANG':
                # DANG d s[0.04] atom pairs
                self.append_card(self.restraints, DANG(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'EADP':
                self.append_card(self.restraints, EADP(spline, list_of_lines), line_num)
                continue
            elif line[:3] == 'REM':
                if ShelXlFile.dsr_regex.match(line):
                    self.dsrlines.append(" ".join(spline))
                    self.dsrline_nums.extend(list_of_lines)
                self.append_card(self.rem, REM(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'AFIX':
                pass
            elif line[:4] == 'CELL':
                # CELL λ a b c α β γ
                if not lastcard == 'TITL':
                    if DEBUG:
                        print('TITL is missing.')
                    # raise ParseOrderError
                if len(spline) >= 8:
                    self.cell = [float(x) for x in spline[2:8]]
                    self.a, self.b, self.c, self.alpha, self.beta, self.gamma = self.cell
                    self.V = self.vol_unitcell(self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
                    # self.A = self.orthogonal_matrix()
                self.wavelen = float(spline[1])
                lastcard = 'CELL'
                continue
            elif line[:4] == "ZERR":
                # ZERR Z esd(a) esd(b) esd(c) esd(α) esd(β) esd(γ)
                if not lastcard == 'CELL':
                    if DEBUG:
                        print('*** Invalid SHELX file!')
                    raise ParseOrderError
                if not self.cell:
                    raise ParseOrderError('*** Cell parameters missing! ***')
                if len(spline) >= 8:
                    self.Z = spline[1]
                    self.zerr = [float(x) for x in spline[2:8]]
                lastcard = 'ZERR'
                continue
            elif line[:4] == "SYMM":
                # SYMM symmetry operation
                #  Being more greedy, because many files do this wrong:
                # if not lastcard == 'ZERR':
                #    raise ParseOrderError
                # if not self.zerr:
                #    raise ParseOrderError
                self.append_card(self.symmcards, SYMM(spline, list_of_lines), line_num)
                lastcard = 'SYMM'
                continue
            elif line[:4] == 'SFAC':
                # SFAC elements or
                # SFAC E a1 b1 a2 b2 a3 b3 a4 b4 c f' f" mu r wt
                # Being less strict to be able to parse files without cell errors:
                # if not (lastcard == 'LATT' or lastcard == 'ZERR'):
                #    raise ParseOrderError
                # if not self.symmcards:
                #    raise ParseOrderError
                if len(spline) <= 1:
                    continue
                self.sfac_table.parse_element_line(spline)
                if self.sfac_table not in self._reslist:
                    self._reslist[line_num] = self.sfac_table
                else:
                    self.delete_on_write.update([line_num])
                lastcard = 'SFAC'
                continue
            elif line[:4] == 'UNIT':
                # UNIT n1 n2 ...
                # Number of atoms of each type in the unit-cell, in SFAC order.
                if not lastcard == 'SFAC':
                    raise ParseOrderError
                if self.sfac_table:
                    try:
                        self.unit = self.assign_card(UNIT(spline, list_of_lines), line_num)
                    except ValueError:
                        if DEBUG:
                            print('*** Non-numeric value in SFAC instruction! ***')
                        raise
                else:
                    raise ParseOrderError
                if len(self.unit.values) != len(self.sfac_table.elements_list):
                    raise ParseNumError
                lastcard = 'UNIT'
                continue
            elif line[:4] == "LATT":
                # LATT N[1]
                # 1=P, 2=I, 3=rhombohedral obverse on hexagonal axes, 4=F, 5=A, 6=B, 7=C.
                # negative is non-centrosymmetric
                self.latt = spline[1]
                if not lastcard == 'ZERR':
                    if DEBUG:
                        print('*** ZERR instruction is missing! ***')
                    # raise ParseOrderError
                continue
            elif line[:4] in ['L.S.', 'CGLS']:
                # CGLS nls[0] nrf[0] nextra[0]
                # L.S. nls[0] nrf[0] nextra[0]
                self.cycles = self.assign_card(LSCycles(self, spline, line_num), line_num)
                continue
            elif line[:4] == "LIST":
                # LIST m[#] mult[1] (mult is for list 4 only)
                self.list = int(spline[1])
                continue
            elif line[:4] == "FVAR":
                # FVAR osf[1] free variables
                # TODO: assign value
                for fvvalue in spline[1:]:
                    fv = FVAR(fvarnum, fvvalue)
                    fvarnum += 1
                    self.fvars.append(fv)
                    if self.fvars not in self._reslist:
                        self._reslist[line_num] = self.fvars
                    else:
                        self.delete_on_write.update([line_num])
            elif line[:4] == 'ANIS':
                # ANIS n or ANIS names
                # Must be before Atom(), to know which atom is anis.
                self.anis = spline
                continue
            elif line[:4] == 'WGHT':
                # WGHT a[0.1] b[0] c[0] d[0] e[0] f[.33333]
                self.wght = self.assign_card(WGHT(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'ACTA':
                # ACTA 2θfull[#] -> optional parameter NOHKL
                self.acta = self.assign_card(ACTA(self, spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'DAMP':
                # DAMP damp[0.7] limse[15]
                self.damp = self.assign_card(DAMP(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'ABIN':
                # ABIN n1 n2   ->   Reads h, k, l, A and B from the file name.fab
                self.abin = [float(x) for x in spline[1:]]
                continue
            elif line[:4] == 'ANSC':
                # ANSC six coefficients
                if len(spline) == 7:
                    self.ansc = [float(x) for x in spline[:1]]
                continue
            elif line[:4] == 'ANSR':
                # ANSR anres[0.001]
                if len(spline) == 2:
                    self.ansr = float(spline[1])
                continue
            elif line[:4] == 'BIND':
                # BIND atom1 atom2
                if len(spline) == 3:
                    self.bind.append(spline[1:])
                continue
            elif line[:4] == 'BLOC':
                # BLOC n1 n2 atomnames
                # TODO: Make class that resolves atomnames and cycles
                self.bloc.append(spline[1:])
                continue
            elif line[:4] == 'BOND':
                # BOND atomnames
                self.append_card(self.commands, BOND(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'BUMP':
                # BUMP s [0.02]
                self.append_card(self.restraints, BUMP(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'CHIV':
                # CHIV V[0] s[0.1] atomnames
                self.append_card(self.restraints, CHIV(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'CONF':
                # CONF atomnames max_d[1.9] max_a[170]
                self.conv.append(spline[1:])
                continue
            elif line[:4] == 'CONN':
                # CONN bmax[12] r[#] atomnames or CONN bmax[12]
                # bonded are d < (r1 + r2 + 0.5) Å
                self.conn.append(spline[1:])
                continue
            elif line[:4] == 'DEFS':
                # DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
                self.defs = DEFS(spline, list_of_lines)
                self.assign_card(self.defs, line_num)
                continue
            elif line[:4] == 'DISP':
                # DISP E f' f"[#] mu[#]
                if not lastcard == 'SFAC':
                    raise ParseOrderError
                self.disp.append(spline[1:])
                continue
            elif line[:4] == 'EQIV':
                # EQIV $n symmetry operation
                if len(spline) > 1:
                    if spline[1].startswith('$'):
                        self.eqiv.append(spline[1:])
                continue
            elif line[:4] == 'EXTI':
                # EXTI x[0]
                self.exti = float(spline[1])
                continue
            elif line[:4] == 'EXYZ':
                # EXYZ atomnames
                self.append_card(self.restraints, EXYZ(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'FRAG':
                # FRAG code[17] a[1] b[1] c[1] α[90] β[90] γ[90]
                if len(spline) == 8:
                    self.frag = spline[1:]
                continue
            elif line[:4] == 'FEND':
                # FEND (must follow FRAG)
                if not self.frag:
                    raise ParseOrderError
                self.frag = None  # Turns frag mode off.
                continue
            elif line[:4] == 'FLAT':
                # FLAT s[0.1] four or more atoms
                self.append_card(self.restraints, FLAT(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'FREE':
                # FREE atom1 atom2
                self.free = spline[1:]
                continue
            elif line[:4] == 'GRID':
                # GRID sl[#] sa[#] sd[#] dl[#] da[#] dd[#]
                self.grid = spline[1:]
                continue
            elif line[:4] == 'HKLF':
                # HKLF N[0] S[1] r11...r33[1 0 0 0 1 0 0 0 1] sm[1] m[0]
                self.hklf = HKLF(spline, list_of_lines)
                self.assign_card(self.hklf, line_num)
                continue
            elif line.startswith('END'):
                # END (after HKLF or ends an include file)
                # TODO: run sanity checks after END like checking if EXYZ and
                # anisotropy fit togeter
                self.end = True
                continue
            elif line[:4] == 'HTAB':
                # HTAB dh[2.0]  or  HTAB donor-atom acceptor-atom
                self.htab = spline[1:]
                continue
            elif line[:4] == 'ISOR':
                # ISOR s[0.1] st[0.2] atomnames
                self.append_card(self.restraints, ISOR(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'LAUE':
                # LAUE E
                # I completely do not understand the LAUE instruction description in the manual!
                continue
            elif line[:4] == 'MERG':
                # MERG n[2]
                self.merg = spline[1:]
                continue
            elif line[:4] == 'MORE':
                # MORE m[1]
                self.more = spline[1]
                continue
            elif line[:4] == 'FMAP':
                # FMAP code[2] axis[#] nl[53]
                self.fmap = spline[1:]
                continue
            elif line[:4] == 'MOVE':
                # MOVE dx[0] dy[0] dz[0] sign[1]
                self.move = spline[1:]
                continue
            elif line[:4] == 'MPLA':
                # MPLA na atomnames
                self.mpla.append(spline[1:])
                continue
            elif line[:4] == 'NCSY':
                # NCSY DN sd[0.1] su[0.05] atoms
                self.append_card(self.restraints, NCSY(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'NEUT':
                # NEUT
                if not lastcard == 'SYMM':
                    raise ParseOrderError
                continue
            elif line[:4] == 'OMIT':
                # OMIT atomnames  or  OMIT s[-2] 2θ(lim)[180]  or  OMIT h k l
                self.omit.append(spline[1:])
                continue
            elif line[:4] == 'PLAN':
                # PLAN npeaks[20] d1[#] d2[#]
                self.plan = spline[1:]
                continue
            elif line[:4] == 'PRIG':
                # PRIG p[#]
                self.prig = spline[1:]
                continue
            elif line[:4] == 'RTAB':
                # RTAB codename atomnames  -->  codename: e.g. 'omeg' gets tabualted in the lst
                self.rtab.append(spline[1:])
                continue
            elif line[:4] == 'SAME':
                # SAME s1[0.02] s2[0.04] atomnames
                self.append_card(self.restraints, SAME(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'SHEL':
                # SHEL lowres[infinite] highres[0]
                self.shel = spline[1:]
                continue
            elif line[:4] == 'SIZE':
                # SIZE dx dy dz
                if len(spline) == 4:
                    self.size['x'] = spline[1]
                    self.size['y'] = spline[2]
                    self.size['z'] = spline[3]
                continue
            elif line[:4] == 'SPEC':
                # SPEC del[0.2]
                # This implementation is not enough, but more is maybe only needed for a
                # refinement program:
                if len(spline) > 1:
                    self.spec = spline[1]
                continue
            elif line[:4] == 'STIR':
                # STIR sres step[0.01]   -> stepwise improvement in the resolution sres
                self.stir = spline[1:]
                continue
            elif line[:4] == 'SUMP':
                # SUMP c sigma c1 m1 c2 m2 ...
                self.append_card(self.sump, SUMP(spline, list_of_lines), line_num)
                continue
            elif line[:4] == 'SWAT':
                # SWAT g[0] U[2]
                self.swat = spline[1:]
                continue
            elif line[:4] == 'TEMP':
                # TEMP T[20]  -> in Celsius
                self.temp = spline[1]
                continue
            elif line[:4] == 'TWIN':
                # TWIN 3x3 matrix [-1 0 0 0 -1 0 0 0 -1] N[2]
                self.twin = TWIN(spline, list_of_lines)
                self.assign_card(self.twin, line_num)
                continue
            elif line[:4] == 'TWST':
                # TWST N[0] (N[1] after SHELXL-2018/3)
                if len(spline) > 1:
                    self.twst = spline[1]
                continue
            elif line[:4] == 'WIGL':
                # WIGL del[0.2] dU[0.2]
                if len(spline) == 1:
                    self.wigl = True
                if len(spline) > 1:
                    self.wigl = spline[1:]
                continue
            elif line[:4] == 'WPDB':
                # WPDB n[1]
                self.wpdb = spline[1:]
                continue
            elif line[:4] == 'XNPD':
                # XNPD Umin[-0.001]
                self.xnpd = spline[1:]
                continue
            elif line[:4] == 'BEDE':
                # Later...
                continue
            elif line[:4] == 'LONE':
                # Later...
                continue
            elif line[:4] == 'MOLE':
                # print('*** MOLE is deprecated! Do not use it! ***')
                pass
            elif line[:4] == 'HOPE':
                # print('*** HOPE is deprecated! Do not use it! ***')
                pass
            elif line[:1] == '+':
                pass
            else:
                if DEBUG:
                    print(line)
                    raise ParseUnknownParam

    @time_this_method
    def run_after_parse(self):
        if self.sump:
            for x in self.sump:
                for y in x:
                    self.fvars.set_fvar_usage(int(y[1]))
        for r in self.restraints:
            if r.name == "DFIX" or r.name == "DANG":
                if abs(r.d) > 4:
                    fvar, value = split_fvar_and_parameter(r.d)
                    self.fvars.set_fvar_usage(fvar)
        if self.abin:
            if len(self.abin) > 1:
                self.fvars.set_fvar_usage(self.abin[0])
                self.fvars.set_fvar_usage(self.abin[1])
            else:
                self.fvars.set_fvar_usage(self.abin[0])
        # Check if basf parameters are consistent:
        if self.basf:
            if self.twin:
                basfs = flatten(self.basf)
                if int(self.twin.allowed_N) != len(basfs):
                    if DEBUG:
                        print('*** Invalid TWIN instruction! BASF with wrong number of parameters. ***')

    def restore_acta_card(self, acta: str):
        """
        Place ACTA after UNIT
        """
        self.add_line(self.unit.line_numbers[-1] + 1, acta)

    def orthogonal_matrix(self):
        """
        Converts von fractional to cartesian by .
        Invert the matrix to do the opposite.

        Old tests:
        #>>> import mpmath as mpm
        #>>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
        #>>> coord = (-0.186843,   0.282708,   0.526803)
        #>>> print(mpm.nstr(A*mpm.matrix(coord)))
        [-2.74151]
        [ 5.90959]
        [ 10.7752]
        #>>> cartcoord = mpm.matrix([['-2.74150542399906'], ['5.909586678'], ['10.7752007008937']])
        #>>> print(mpm.nstr(A**-1*cartcoord))
        [-0.186843]
        [ 0.282708]
        [ 0.526803]
        """
        return Matrix([[self.a, self.b * cos(self.gamma), self.c * cos(self.beta)],
                       [0, self.b * sin(self.gamma),
                        (self.c * (cos(self.alpha) - cos(self.beta) * cos(self.gamma)) / sin(self.gamma))],
                       [0, 0, self.V / (self.a * self.b * sin(self.gamma))]])

    @staticmethod
    def vol_unitcell(a, b, c, al, be, ga) -> float:
        """
        calculates the volume of a unit cell
        >>> v = ShelXlFile.vol_unitcell(2, 2, 2, 90, 90, 90)
        >>> print(v)
        8.0
        """
        ca, cb, cg = cos(radians(al)), cos(radians(be)), cos(radians(ga))
        v = a * b * c * sqrt(1 + 2 * ca * cb * cg - ca ** 2 - cb ** 2 - cg ** 2)
        return v

    def __repr__(self):
        """
        Represents the shelxl object.
        """
        resl = []
        for num, line in enumerate(self._reslist):
            if num in self.delete_on_write:
                if DEBUG:
                    pass
                    # print('Deleted line {}'.format(num + 1))
                continue
            if line == '' and self._reslist[num + 1] == '':
                continue
            line = self.wrap_line(line)
            resl.append(line)
        return "\n".join(resl)

    def wrap_line(self, line: str) -> str:
        # Generally, all Shelx opbjects have no line wrap. I do this now:
        line = textwrap.wrap(str(line), 78, subsequent_indent='  ', drop_whitespace=False, replace_whitespace=False)
        if len(line) > 1:
            newline = []
            for n, l in enumerate(line):
                if n < len(line) - 1:
                    l += ' =\n'
                newline.append(l)
            line = ' '.join(newline)
        else:
            line = ''.join(line)
        return line

    # @time_this_method
    def read_file_to_list(self, resfile: str) -> list:
        """
        Read in shelx file and returns a list without line endings. +include files are inserted
        also.
        :param resfile: The path to a SHLEL .res or .ins file.
        """
        reslist = []
        #resnodes = ResList()
        includefiles = []
        try:
            with open(resfile, 'r') as f:
                reslist = f.read().splitlines(keepends=False)
                #for ll in reslist:
                    #resnodes.append(ll)
                for n, line in enumerate(reslist):
                    if line.startswith('+'):
                        try:
                            include_filename = line.split()[1]
                            # Detect recoursive file inclusion:
                            if include_filename in includefiles:
                                raise ValueError('*** Recoursive include files detected! ***')
                            includefiles.append(include_filename)
                            newfile = ShelXlFile.read_nested_file_to_list(os.path.join(os.path.dirname(resfile),
                                                                                       include_filename))
                            if newfile:
                                for num, l in enumerate(newfile):
                                    lnum = n + 1 + num
                                    # '+filename' include files are not copied to res file,
                                    #  so I have to delete these lines on write.
                                    # '++filename' copies them to the .res file where appropriate
                                    if l.startswith('+') and l[:2] != '++':
                                        self.delete_on_write.update([lnum])
                                    reslist.insert(lnum, l)
                                continue
                        except IndexError:
                            if DEBUG:
                                print('*** CANNOT READ INCLUDE FILE {} ***'.format(line))
                            #del reslist[n]
        except (IOError) as e:
            print(e)
            print('*** CANNOT READ FILE {} ***'.format(resfile))
        return reslist

    @staticmethod
    def read_nested_file_to_list(resfile: str) -> list:
        """
        Read in shelx file and returns a list without line endings.
        :param resfile: The path to a SHLEL .res or .ins file.
        """
        reslist = []
        try:
            with open(os.path.abspath(resfile), 'r') as f:
                reslist = f.read().splitlines(keepends=False)
        except (IOError) as e:
            if DEBUG:
                print(e)
                print('*** CANNOT OPEN NESTED INPUT FILE {} ***'.format(resfile))
            return reslist
        return reslist

    def reload(self, resfile):
        """
        Reloads the shelx file and parses it again.
        """
        if resfile:
            if DEBUG:
                print('loading file:', resfile)
            #self.write_shelx_file(resfile)
            self.__init__(self.resfile)
        else:
            if DEBUG:
                print('*** Can not read empty file data! ***')

    def write_shelx_file(self, filename=None, verbose=False):
        if not filename:
            filename = self.resfile
        with open(filename, 'w') as f:
            for num, line in enumerate(self._reslist):
                if num in self.delete_on_write:
                    if DEBUG:
                        pass
                        #print('Deleted line {}'.format(num + 1))
                    continue
                if line == '' and self._reslist[num + 1] == '':
                    continue
                line = self.wrap_line(line)
                f.write(line + '\n')
        if verbose or DEBUG:
            print('File successfully written to {}'.format(os.path.abspath(filename)))
            return True
        return True

    def append_card(self, obj, card, line_num):
        obj.append(card)
        self._reslist[line_num] = card

    def assign_card(self, card, line_num):
        self._reslist[line_num] = card
        return card

    @staticmethod
    def get_resi_definition_dict(resi: list) -> dict:
        """
        Returns the residue number and class of a string like 'RESI TOL 1'
        or 'RESI 1 TOL'

        Residue names may now begin with a digit.
        They must however contain at least one letter

        Allowed residue numbers is now from -999 to 9999 (2017/1)

        TODO: support alias

        :param resi: ['number', 'class']
        :type resi: list or string

        >>> sorted(list(ShelXlFile.get_resi_definition_dict('RESI 1 TOL'.split()[1:])))
        ['ID', 'class', 'number']
        >>> sorted(ShelXlFile.get_resi_definition_dict('RESI 1 TOL'.split()[1:]).items())
        [('ID', None), ('class', 'TOL'), ('number', 1)]
        >>> ShelXlFile.get_resi_definition_dict('RESI 1 TOL'.split()[1:])
        {'class': 'TOL', 'number': 1, 'ID': None}
        >>> ShelXlFile.get_resi_definition_dict('RESI A:100 TOL'.split()[1:])
        {'class': 'TOL', 'number': 100, 'ID': 'A'}
        >>> ShelXlFile.get_resi_definition_dict('RESI -10 TOL'.split()[1:])
        {'class': 'TOL', 'number': -10, 'ID': None}
        >>> ShelXlFile.get_resi_definition_dict('RESI b:-10 TOL'.split()[1:])
        {'class': 'TOL', 'number': -10, 'ID': 'b'}
        """
        resi_dict = {'class': None, 'number': None, 'ID': None}
        for x in resi:
            if re.search('[a-zA-Z]', x):
                if ':' in x:
                    # contains :, must be a chain-id+number
                    resi_dict['ID'], resi_dict['number'] = x.split(':')[0], int(x.split(':')[1])
                else:
                    # contains letters, must be a name
                    resi_dict['class'] = x
            else:
                # everything else can only be a number
                resi_dict['number'] = int(x)
        return resi_dict

    @staticmethod
    def get_partnumber(partstring: str) -> (int, float):
        """
        get the part number from a string like PART 1 oder PART 2 -21

        PART n sof
        partstring: string like 'PART 2 -21'

        >>> ShelXlFile.get_partnumber(partstring='PART 2 -21')
        (2, -21.0)
        """
        part = partstring.upper().split()
        sof = 0
        try:
            partnum = int(part[1])
        except(ValueError, IndexError):
            if DEBUG:
                print('*** Wrong PART definition found! Check your PART instructions ***')
            partnum = 0
        if len(part) > 2:
            sof = float(part[2])
        return partnum, sof

    @staticmethod
    def get_afix_numbers(spline: list, line_num: int) -> ((int, float, float, float), [int, float, float, float]):
        """
        Returns a tuple with the AFIX instructions. afixnum, d, sof, u

        AFIX mn d[#] sof[11] U[10.08]

        >>> ShelXlFile.get_afix_numbers(['AFIX', 137], 1)
        (137, 0.0, 11, 10.08)
        >>> ShelXlFile.get_afix_numbers(['AFIX', '137b'], 1)
        *** Wrong AFIX definition in line 1. Check your AFIX instructions ***
        (0, 0.0, 11, 10.08)
        >>> ShelXlFile.get_afix_numbers(['AFIX', 13, 1.234], 1)
        (13, 1.234, 11, 10.08)
        >>> ShelXlFile.get_afix_numbers(['AFIX', 13, 1.234, -21, 10.05], 1)
        (13, 1.234, -21.0, 10.05)
        """
        d = 0.0
        sof = 11
        u = 10.08
        try:
            afixnum = int(spline[1])
        except(ValueError, IndexError):
            if DEBUG:
                print('*** Wrong AFIX definition in line {}. Check your AFIX instructions ***'.format(line_num))
            afixnum = 0
        if len(spline) > 2:
            d = float(spline[2])
        if len(spline) > 3:
            sof = float(spline[3])
        if len(spline) > 4:
            u = float(spline[4])
        return afixnum, d, sof, u

    @staticmethod
    def is_atom(atomline: str) -> bool:
        """
        Returns True is line contains an atom.
        atomline:  'O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...'
        >>> ShelXlFile.is_atom(atomline = 'O1    3    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        True
        >>> ShelXlFile.is_atom(atomline = 'O1    0.120080   0.336659   0.494426  11.00000   0.01445 ...')
        False
        >>> ShelXlFile.is_atom(atomline = 'O1  4  0.120080    0.494426  11.00000   0.01445 ...')
        True
        >>> ShelXlFile.is_atom("AFIX 123")
        False
        >>> ShelXlFile.is_atom("AFIX")
        False
        >>> ShelXlFile.is_atom('O1    3    0.120080   0.336659   0.494426')
        True
        """
        # no empty line, not in cards and not space at start:
        if atomline[:4].upper() not in SHX_CARDS:  # exclude all non-atom cards
            # Too few parameter for an atom:
            if len(atomline.split()) < 5:
                return False
            # means sfac number is missing:
            if '.' in atomline.split()[1]:
                return False
            return True
        else:
            return False

    def elem2sfac(self, atom_type: str) -> int:
        """
        returns an sfac-number for the element given in "atom_type"
        >>> shx = ShelXlFile('p21c.res')
        >>> shx.elem2sfac('O')
        3
        >>> shx.elem2sfac('c')
        1
        >>> shx.elem2sfac('Ar')

        """
        for num, element in enumerate(self.sfac_table, 1):
            if atom_type.upper() == element.upper():
                return num  # return sfac number

    def sfac2elem(self, sfacnum: int) -> str:
        """
        returns an element and needs an sfac-number
        :param sfacnum: string like '2'
        >>> shx = ShelXlFile('./p21c.res')
        >>> shx.sfac2elem(1)
        'C'
        >>> shx.sfac2elem(2)
        'H'
        >>> shx.sfac2elem(3)
        'O'
        >>> shx.sfac2elem(5)
        'Al'
        >>> shx.sfac2elem(8)
        ''
        >>> shx.sfac2elem(0)
        ''
        """
        try:
            elem = self.sfac_table[int(sfacnum)]
        except IndexError:
            return ''
        return elem

    def add_line(self, linenum: int, obj):
        """
        Adds a new SHELX card to the reslist after linenum.
        """
        self._reslist.insert(linenum + 1, obj)

    def replace_line(self, obj, new_line: str):
        """
        Replaces a single line in the res file with new_line.
        """
        self._reslist[self.index_of(obj)] = new_line

    def index_of(self, obj):
        return self._reslist.index(obj)

    def insert_frag_fend_entry(self, dbatoms: list, cell: list):
        """
        Inserts the FRAG ... FEND entry in the res file.
        :param dbatoms:   list of atoms in the database entry
        :param cell:  string with "FRAG 17 cell" from the database entry
        """
        dblist = []
        for line in dbatoms:
            dblist.append("{:4} {:<4} {:>8}  {:>8}  {:>8}".format(*line))
        dblines = ' The following is from DSR:\n'
        dblines = dblines + 'FRAG 17 {} {} {} {} {} {}'.format(*cell) + '\n'
        dblines = dblines + '\n'.join(dblist)
        dblines = dblines + '\nFEND\n'
        # insert the db entry right after FVAR
        print(self.fvars.line_number)
        self.add_line(self.fvars.line_number, dblines)



if __name__ == "__main__":
    def runall():
        """
        >>> file = r'p21c.res'
        >>> try:
        >>>     shx = ShelXlFile(file)
        >>> except Exception:
        >>>    raise
        """
        pass

    """
    def get_commands():
        url = "http://shelx.uni-goettingen.de/shelxl_html.php"
        response = urlopen('{}/version.txt'.format(url))
        html = response.read().decode('UTF-8')
        #res = BeautifulSoup(html, "html5lib")
        tags = res.findAll("p", {"class": 'instr'})
        for l in tags:
            if l:
                print(str(l).split(">")[1].split("<")[0])
    """
    #get_commands()
    #sys.exit()

    file = r'test-data/p21c.res'
    try:
        shx = ShelXlFile(file)
    except Exception:
        raise

    print(shx.atoms.distance('Ga1', 'Al1'))
    print(shx.hklf)
    print(shx.sfac_table.is_exp(shx.sfac_table[1]))
    print(shx.sfac_table)
    shx.sfac_table.add_element('Zn')
    print(shx.unit)
    #print(shx.sfac_table.remove_element('In'))
    print(shx.sfac_table)
    print(shx.unit)
    print(shx.fvars.line_number)
    shx.cycles.set_refine_cycles(33)
    shx.write_shelx_file(r'./test.ins')
    print('\n\n')
    print(shx)
    print('######################')
    sys.exit()
    # for x in shx.atoms:
    #    print(x)
    # shx.reload()
    # for x in shx.restraints:
    #    print(x)
    # for x in shx.rem:
    #    print(x)
    # print(shx.size)
    # for x in shx.sump:
    #    print(x)
    # print(float(shx.temp)+273.15)
    # print(shx.atoms.atoms_in_class('CCF3'))
    #sys.exit()
    files = walkdir(r'/Users/daniel', 'res')
    print('finished')
    for f in files:
        if "dsrsaves" in str(f) or ".olex" in str(f) or 'ED' in str(f) or 'shelXlesaves' in str(f):
            continue
        #path = f.parent
        #file = f.name
        #print(path.joinpath(file))
        #id = id_generator(size=4)
        #copy(str(f), Path(r"d:/Github/testresfiles/").joinpath(id+file))
        #print('copied', str(f.name))
        print(f)
        shx = ShelXlFile(f, debug=False)
        # print(len(shx.atoms), f)


"""
SHELXL cards:

ABIN n1 n2 
ACTA 2θfull[#]
AFIX mn d[#] sof[11] U[10.08]
ANIS n 
ANIS names
ANSC six coefficients
ANSR anres[0.001]
BASF scale factors
BIND atom1 atom2
BIND m n
BLOC n1 n2 atomnames 
BOND atomnames
BUMP s [0.02]
CELL λ a b c α β γ
CGLS nls[0] nrf[0] nextra[0]
CHIV V[0] s[0.1] atomnames
CONF atomnames max_d[1.9] max_a[170]
CONN bmax[12] r[#] atomnames or CONN bmax[12]
DAMP damp[0.7] limse[15]
DANG d s[0.04] atom pairs
DEFS sd[0.02] sf[0.1] su[0.01] ss[0.04] maxsof[1]
DELU s1[0.01] s2[0.01] atomnames
DFIX d s[0.02] atom pairs
DISP E f' f"[#] mu[#]
EADP atomnames
END
EQIV $n symmetry operation
EXTI x[0] 
EXYZ atomnames
FEND
FLAT s[0.1] four or more atoms
FMAP code[2] axis[#] nl[53]
FRAG code[17] a[1] b[1] c[1] α[90] β[90] γ[90]
FREE atom1 atom2
FVAR osf[1] free variables
GRID sl[#] sa[#] sd[#] dl[#] da[#] dd[#]
HFIX mn U[#] d[#] atomnames
HKLF N[0] S[1] r11...r33[1 0 0 0 1 0 0 0 1] sm[1] m[0]
HTAB dh[2.0]
HTAB donor-atom acceptor-atom
ISOR s[0.1] st[0.2] atomnames
LATT N[1]
LAUE E
LIST m[#] mult[1]
L.S. nls[0] nrf[0] nextra[0]
MERG n[2]
MORE m[1]
MOVE dx[0] dy[0] dz[0] sign[1]
MPLA na atomnames
NCSY DN sd[0.1] su[0.05] atoms
NEUT
OMIT atomnames
OMIT s[-2] 2θ(lim)[180]
OMIT h k l
PART n sof
PLAN npeaks[20] d1[#] d2[#]
PRIG p[#]
REM
RESI class[ ] number[0] alias
RIGU s1[0.004] s2[0.004] atomnames
RTAB codename atomnames
SADI s[0.02] pairs of atoms
SAME s1[0.02] s2[0.04] atomnames
SFAC elements
SFAC E a1 b1 a2 b2 a3 b3 a4 b4 c f' f" mu r wt
SHEL lowres[infinite] highres[0]
SIMU s[0.04] st[0.08] dmax[2.0] atomnames
SIZE dx dy dz
SPEC del[0.2]
STIR sres step[0.01]
SUMP c sigma c1 m1 c2 m2 ... 
SWAT g[0] U[2] 
SYMM symmetry operation
TEMP T[20]
TITL [ ]
TWIN 3x3 matrix [-1 0 0 0 -1 0 0 0 -1] N[2]
TWST N[0]
UNIT n1 n2 ...
WGHT a[0.1] b[0] c[0] d[0] e[0] f[.33333]
WIGL del[0.2] dU[0.2]
WPDB n[1]
XNPD Umin[-0.001]
ZERR Z esd(a) esd(b) esd(c) esd(α) esd(β) esd(γ)
"""