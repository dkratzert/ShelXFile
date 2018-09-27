# /usr/bin/env python
# -*- encoding: utf-8 -*-
# möp
# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42):
# <daniel.kratzert@ac.uni-freiburg.de> wrote this file. As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some day,
# and you think this stuff is worth it, you can buy me a beer in return.
# Daniel Kratzert
# ----------------------------------------------------------------------------
#

import re

num2element = {
    1:   'H',
    2:   'He',
    3:   'Li',
    4:   'Be',
    5:   'B',
    6:   'C',
    7:   'N',
    8:   'O',
    9:   'F',
    10:  'Ne',
    11:  'Na',
    12:  'Mg',
    13:  'Al',
    14:  'Si',
    15:  'P',
    16:  'S',
    17:  'Cl',
    18:  'Ar',
    19:  'K',
    20:  'Ca',
    21:  'Sc',
    22:  'Ti',
    23:  'V',
    24:  'Cr',
    25:  'Mn',
    26:  'Fe',
    27:  'Co',
    28:  'Ni',
    29:  'Cu',
    30:  'Zn',
    31:  'Ga',
    32:  'Ge',
    33:  'As',
    34:  'Se',
    35:  'Br',
    36:  'Kr',
    37:  'Rb',
    38:  'Sr',
    39:  'Y',
    40:  'Zr',
    41:  'Nb',
    42:  'Mo',
    43:  'Tc',
    44:  'Ru',
    45:  'Rh',
    46:  'Pd',
    47:  'Ag',
    48:  'Cd',
    49:  'In',
    50:  'Sn',
    51:  'Sb',
    52:  'Te',
    53:  'I',
    54:  'Xe',
    55:  'Cs',
    56:  'Ba',
    57:  'La',
    58:  'Ce',
    59:  'Pr',
    60:  'Nd',
    61:  'Pm',
    62:  'Sm',
    63:  'Eu',
    64:  'Gd',
    65:  'Tb',
    66:  'Dy',
    67:  'Ho',
    68:  'Er',
    69:  'Tm',
    70:  'Yb',
    71:  'Lu',
    72:  'Hf',
    73:  'Ta',
    74:  'W',
    75:  'Re',
    76:  'Os',
    77:  'Ir',
    78:  'Pt',
    79:  'Au',
    80:  'Hg',
    81:  'Tl',
    82:  'Pb',
    83:  'Bi',
    84:  'Po',
    85:  'At',
    86:  'Rn',
    87:  'Fr',
    88:  'Ra',
    89:  'Ac',
    90:  'Th',
    91:  'Pa',
    92:  'U',
    93:  'Np',
    94:  'Pu',
    95:  'Am',
    96:  'Cm',
    97:  'Bk',
    98:  'Cf',
    99:  'Es',
    100: 'Fm',
    101: 'Md',
    102: 'No',
    103: 'Lr',
    104: 'Rf',
    105: 'Db',
    106: 'Sg',
    107: 'Bh',
    108: 'Hs',
    109: 'Mt',
    110: 'Ds',
    111: 'Rg',
    112: 'Cn',
    114: 'Uuq',
    116: 'Uuh'
}

element2num = {
    'H':  1,
    'He': 2,
    'Li': 3,
    'Be': 4,
    'B':  5,
    'C':  6,
    'N':  7,
    'O':  8,
    'F':  9,
    'Ne': 10,
    'Na': 11,
    'Mg': 12,
    'Al': 13,
    'Si': 14,
    'P':  15,
    'S':  16,
    'Cl': 17,
    'Ar': 18,
    'K':  19,
    'Ca': 20,
    'Sc': 21,
    'Ti': 22,
    'V':  23,
    'Cr': 24,
    'Mn': 25,
    'Fe': 26,
    'Co': 27,
    'Ni': 28,
    'Cu': 29,
    'Zn': 30,
    'Ga': 31,
    'Ge': 32,
    'As': 33,
    'Se': 34,
    'Br': 35,
    'Kr': 36,
    'Rb': 37,
    'Sr': 38,
    'Y':  39,
    'Zr': 40,
    'Nb': 41,
    'Mo': 42,
    'Tc': 43,
    'Ru': 44,
    'Rh': 45,
    'Pd': 46,
    'Ag': 47,
    'Cd': 48,
    'In': 49,
    'Sn': 50,
    'Sb': 51,
    'Te': 52,
    'I':  53,
    'Xe': 54,
    'Cs': 55,
    'Ba': 56,
    'La': 57,
    'Ce': 58,
    'Pr': 59,
    'Nd': 60,
    'Pm': 61,
    'Sm': 62,
    'Eu': 63,
    'Gd': 64,
    'Tb': 65,
    'Dy': 66,
    'Ho': 67,
    'Er': 68,
    'Tm': 69,
    'Yb': 70,
    'Lu': 71,
    'Hf': 72,
    'Ta': 73,
    'W':  74,
    'Re': 75,
    'Os': 76,
    'Ir': 77,
    'Pt': 78,
    'Au': 79,
    'Hg': 80,
    'Tl': 81,
    'Pb': 82,
    'Bi': 83,
    'Po': 84,
    'At': 85,
    'Rn': 86,
    'Fr': 87,
    'Ra': 88,
    'Ac': 89,
    'Th': 90,
    'Pa': 91,
    'U':  92,
    'Np': 93,
    'Pu': 94,
    'Am': 95,
    'Cm': 96,
    'Bk': 97,
    'Cf': 98,
}

atoms = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
         'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
         'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',
         'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
         'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
         'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',
         'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
         'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf']

num2covradius = {
    1:  0.55,
    2:  1.00,
    3:  1.23,
    4:  0.90,
    5:  0.80,
    6:  0.77,
    7:  0.74,
    8:  0.71,
    9:  0.72,
    10: 1.00,
    11: 1.54,
    12: 1.49,
    13: 1.18,
    14: 1.11,
    15: 1.06,
    16: 1.02,
    17: 0.99,
    18: 1.00,
    19: 2.03,
    20: 1.74,
    21: 1.44,
    22: 1.32,
    23: 1.22,
    24: 1.18,
    25: 1.17,
    26: 1.17,
    27: 1.16,
    28: 1.15,
    29: 1.17,
    30: 1.25,
    31: 1.26,
    32: 1.22,
    33: 1.20,
    34: 1.16,
    35: 1.14,
    36: 1.12,
    37: 2.16,
    38: 1.91,
    39: 1.62,
    40: 1.45,
    41: 1.34,
    42: 1.30,
    43: 1.27,
    44: 1.25,
    45: 1.25,
    46: 1.28,
    47: 1.34,
    48: 1.48,
    49: 1.44,
    50: 1.41,
    51: 1.40,
    52: 1.36,
    53: 1.33,
    54: 1.31,
    55: 2.35,
    56: 1.98,
    57: 1.69,
    58: 1.65,
    59: 1.65,
    60: 1.64,
    61: 1.63,
    62: 1.62,
    63: 1.85,
    64: 1.61,
    65: 1.59,
    66: 1.59,
    67: 1.58,
    68: 1.57,
    69: 1.56,
    70: 1.74,
    71: 1.56,
    72: 1.44,
    73: 1.34,
    74: 1.30,
    75: 1.28,
    76: 1.26,
    77: 1.27,
    78: 1.30,
    79: 1.34,
    80: 1.49,
    81: 1.48,
    82: 1.47,
    83: 1.46,
    84: 1.46,
    85: 1.45,
    86: 1.00,
    87: 1.00,
    88: 1.00,
    89: 1.88,
    90: 1.65,
    91: 1.61,
    92: 1.42,
    93: 1.30,
    94: 1.51,
    95: 1.82,
    96: 1.20,
    97: 1.20,
    98: 1.20,
    99: 1.20
}


def get_radius(atomic_number: int) -> float:
    """
    Get the covalent radius in pm for the element.

    >>> get_radius(6)
    0.77
    """
    return num2covradius[atomic_number]


def get_radius_from_element(element: str) -> float:
    """
    Returns the radius of an atom by its element name.

    >>> get_radius_from_element('F')
    0.72
    """
    return get_radius(element2num[element])


def get_atomic_number(element: str) -> int:
    """
    returns the atomic number from the element symbol

    >>> get_atomic_number('F')
    9
    """
    return element2num[element]


def get_element(atomic_number: int) -> str:
    """
    returns the element symbol from the atomic number

    >>> get_element(7)
    'N'
    """
    return num2element[atomic_number]


def get_atomlabel(input_atom: str) -> str:
    """
    converts an atom name like C12 to the element symbol C.

    >>> get_atomlabel('C12')
    'C'
    """
    atom = ''
    for x in input_atom:  # iterate over characters in i
        if re.match(r'^[A-Za-z#]', x):  # Alphabet and "#" as allowed characters in names
            atom = atom + x.upper()  # add characters to atoms until numbers occur
        else:  # now we have atoms like C, Ca, but also Caaa
            break
    try:
        if atom[0:2].capitalize() in atoms:  # fixes names like Caaa to be just Ca
            return atom[0:2].capitalize()  # atoms first, search for all two-letter atoms
        elif atom[0].upper() in atoms:
            return atom[0]  # then for all one-letter atoms
        else:
            print('*** {} is not a valid atom!! ***'.format(atom))
            raise KeyError
    except(IndexError):
        print('*** {} is not a valid atom! ***'.format(atom))
        raise KeyError


if __name__ == '__main__':
    import doctest

    failed, attempted = doctest.testmod()  # verbose=True)
    if failed == 0:
        print('passed all {} tests!'.format(attempted))
    else:
        print('{} of {} tests failed'.format(failed, attempted))
