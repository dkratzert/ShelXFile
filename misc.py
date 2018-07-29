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
import os
import re
import textwrap
import time
from shutil import get_terminal_size

from dsrmath import frac_to_cart, subtract_vect, determinante

DEBUG = True
PROFILE = False

dsr_regex = re.compile(r'^rem\s+DSR\s+(PUT|REPLACE).*', re.IGNORECASE)


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
    def __init__(self, arg=None):
        if DEBUG:
            if arg:
                print(arg)
            print("*** UNKNOWN PARAMETER ***")


class ParseSyntaxError(Exception):
    def __init__(self):
        print("*** Syntax Error ***")


try:
    width, height = get_terminal_size()  # @UnusedVariable
except():
    width = 80
sep_line = (width-1)*'-'


def remove_file(filename, exit_dsr=False):
    """
    removes the file "filename" from disk
    program exits when exit is true
    platon gets terminated if terminate is true

    >>> remove_file('foobar')
    """
    if os.path.isfile(filename):
        try:
            os.remove(filename)
        except(IOError, OSError):
            return False
        return True


def find_line(inputlist: list, regex: str, start: int = None):
    """
    returns the index number of the line where regex is found in the inputlist
    if stop is true, stop searching with first line found
    :param inputlist: list of strings
    :type inputlist: list
    :param regex: regular expression to search
    :type regex: string
    :param start: line number where to start the search
    :param start: start searching at line start
    :type start: string or int
    >>> input = ['Hallo blub', 'foo bar blub', '123', '1 blub 2 3 4']
    >>> find_line(input, '.*blub.*')
    0
    >>> input = [['foo'],['bar']]
    >>> find_line(input, '.*blub.*') #doctest: +REPORT_NDIFF +NORMALIZE_WHITESPACE +ELLIPSIS
    Traceback (most recent call last):
        ...
    TypeError: expected string or ...
    """
    if start:
        start = int(start)
        inputlist_slice = inputlist[start:]
        for i, string in enumerate(inputlist_slice, start):
            if re.match(regex, string, re.IGNORECASE):
                return i  # returns the index number if regex found
    else:
        for i, string in enumerate(inputlist):
            if re.match(regex, string, re.IGNORECASE):
                return i  # returns the index number if regex found
    return False  # returns False if no regex found (xt solution has no fvar)


def which(name: str, flags=os.X_OK, exts=None) -> list:
    """
    Search PATH for executable files with the given name.

    On MS-Windows the only flag that has any meaning is os.F_OK. Any other
    flags will be ignored.
    """
    if exts is None:
        exts = ['.exe', '.EXE', '.bat']
    result = []
    # exts = filter(None, os.environ.get('PATHEXT', '').split(os.pathsep))
    path = os.getenv('PATH', None)
    if path is None:
        return []
    for p in os.getenv('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if os.access(p, flags):
            result.append(p)
        for e in exts:
            pext = p + e
            if os.access(pext, flags):
                result.append(pext)
    return result


def split_fvar_and_parameter(parameter: float) -> tuple:
    """
    Returns the free variable and value of a given parameter e.g. 30.5 for the occupancy.
    :return (fvar: int, value: float)

    >>> split_fvar_and_parameter(30.5)
    (3, 0.5)
    >>> split_fvar_and_parameter(31.0)
    (3, 1.0)
    >>> split_fvar_and_parameter(-30.5)
    (-3, -0.5)
    >>> split_fvar_and_parameter(11.0)
    (1, 1.0)
    >>> split_fvar_and_parameter(-11.0)
    (-1, -1.0)
    >>> split_fvar_and_parameter(-10.33333333)
    (-1, -0.33333333)
    """
    fvar = abs(int(str(parameter).split('.')[0])) // 10  # The free variable number e.g. 2
    value = abs(float(parameter)) % 10  # The value with which the free variable was multiplied e.g. 0.5
    if parameter < 0:
        value *= -1
        fvar *= -1
    return fvar, round(value, 8)


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
        if line.startswith("REM") and not dsr_regex.match(line):
            return False
        return True
    else:
        return False


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
    Contains the lines of the res file as unordered linked list.
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


def wrap_line(line: str) -> str:
    """
    Wraps long lines according to SHELXL syntax with = at end and space characters before the next line.
    The wrapping will only be at whitespace, not inside words.

    >>> wrap_line("This is a really long line with over 79 characters. Shelxl wants it to be wrapped.")
    'This is a really long line with over 79 characters. Shelxl wants it to be  =\\n   wrapped.'
    """
    line = textwrap.wrap(line, 79, subsequent_indent='  ', drop_whitespace=False, replace_whitespace=False)
    if len(line) > 1:
        newline = []
        for n, ln in enumerate(line):
            if n < len(line) - 1:
                ln += ' =\n'
            newline.append(ln)
        line = ' '.join(newline)
    else:
        line = ''.join(line)
    return line


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


def walkdir(rootdir, include="", exclude=""):
    """
    Returns a list of files in all subdirectories with full path.
    :param rootdir: base path from which walk should start
    :param filter: list of file endings to include only e.g. ['.py', '.res']
    :return: list of files

    #>>> walkdir("../docs") #doctest: +REPORT_NDIFF +NORMALIZE_WHITESPACE +ELLIPSIS
    #['../docs/test.txt']
    #>>> walkdir("../setup/modpath.iss")
    #['../setup/modpath.iss']
    #>>> walkdir("../setup/modpath.iss", exclude=['.iss'])
    #[]
    #>>> walkdir("../docs", exclude=['.txt']) #doctest: +REPORT_NDIFF +NORMALIZE_WHITESPACE +ELLIPSIS
    #[]
    """
    results = []
    if not os.path.isdir(rootdir):
        if os.path.splitext(rootdir)[1] in exclude:
            return []
        return [rootdir]
    for root, subFolders, files in os.walk(rootdir):
        for file in files:
            fullfilepath = os.path.join(root, file)
            if exclude:
                if os.path.splitext(fullfilepath)[1] in exclude:
                    continue
            if include:
                if os.path.splitext(fullfilepath)[1] in include:
                    results.append(os.path.normpath(fullfilepath).replace('\\', '/'))
            else:
                results.append(os.path.normpath(fullfilepath).replace('\\', '/'))
    return results