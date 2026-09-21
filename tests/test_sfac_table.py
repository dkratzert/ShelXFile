"""The ``SFAC`` table and its positional relationship to ``UNIT``.

Two forms exist: *"SFAC elements"* and the explicit
*"SFAC E a1 b1 a2 b2 a3 b3 a4 b4 c f' f" mu r wt"*.

Both are **positional**: an atom's second field *"refers to the list
defined by the SFAC instruction(s)"*, and ``UNIT`` gives one number per
entry in that same order.  Anything that reorders, merges or drops an
entry renumbers every atom in the file.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile

BODY = (
    'C1    1    0.100000    0.200000    0.300000    11.00000    0.0500\n'
    'HKLF 4\n'
    'END\n'
)


def _shx(sfac: str, unit: str) -> Shelxfile:
    shx = Shelxfile()
    shx.read_string(
        'TITL t\n'
        'CELL 0.71073 10.0 11.0 12.0 90 90 90\n'
        'ZERR 4 0.001 0.001 0.001 0 0 0\n'
        'LATT 1\n'
        f'{sfac}\n{unit}\n' + BODY
    )
    return shx


# ------------------------------------------------------- repeated elements

def test_the_same_element_twice_is_two_slots() -> None:
    """``Cl`` and ``CL`` are one element but two scattering-factor slots."""
    shx = _shx('SFAC  C  H  O  F  AL  SI  BR  AG  Cl  CL', 'UNIT  1 2 3 4 5 6 7 8 9 10')
    assert len(shx.sfac_table.sfac_table) == 10
    assert len(shx.unit.values) == 10


def test_repeated_elements_survive_a_write() -> None:
    """Merging them would shift every later sfac number by one."""
    shx = _shx('SFAC  C  H  O  F  AL  SI  BR  AG  Cl  CL', 'UNIT  1 2 3 4 5 6 7 8 9 10')
    line = [x for x in shx.dumps().splitlines() if x.startswith('SFAC')][0]
    assert line.split()[1:] == ['C', 'H', 'O', 'F', 'Al', 'Si', 'Br', 'Ag', 'Cl', 'Cl']


def test_repeated_elements_round_trip() -> None:
    shx = _shx('SFAC  C  H  O  F  AL  SI  BR  AG  Cl  CL', 'UNIT  1 2 3 4 5 6 7 8 9 10')
    once = shx.dumps()
    again = Shelxfile()
    again.read_string(once)
    assert again.dumps() == once
    assert len(again.sfac_table.sfac_table) == 10


def test_sfac_numbering_still_matches_the_atoms() -> None:
    shx = _shx('SFAC  C  H  O  F  AL  SI  BR  AG  Cl  CL', 'UNIT  1 2 3 4 5 6 7 8 9 10')
    assert shx.sfac2elem(9) == 'Cl'
    assert shx.sfac2elem(10) == 'Cl'
    assert shx.sfac2elem(1) == 'C'


# ---------------------------------------------------------- explicit form

EXPLICIT = ('SFAC C 2.31 20.8439 1.02 10.2075 1.5886 0.5687 0.865 51.6512 '
            '0.2156 0.0033 0.0016 0.00 0.77 12.01')


def _sfac_lines(shx: Shelxfile) -> str:
    """The SFAC instruction(s) as one logical string, continuations joined."""
    joined = shx.dumps().replace('=\n', '')
    return '\n'.join(x for x in joined.splitlines() if x.startswith('SFAC'))


def test_explicit_scattering_factors_are_kept() -> None:
    """The full form carries data that cannot be regenerated."""
    shx = _shx(EXPLICIT, 'UNIT 4')
    line = _sfac_lines(shx)
    for value in ('2.31', '20.8439', '51.6512', '12.01'):
        assert value in line, f'{value} was dropped from {line!r}'


def test_explicit_form_round_trips() -> None:
    shx = _shx(EXPLICIT, 'UNIT 4')
    once = shx.dumps()
    again = Shelxfile()
    again.read_string(once)
    assert again.dumps() == once


def test_explicit_and_plain_forms_can_be_mixed() -> None:
    """*"There may be more than one SFAC ... instruction."*"""
    shx = _shx(EXPLICIT + '\nSFAC H O', 'UNIT 4 8 2')
    assert len(shx.sfac_table.sfac_table) == 3
    once = shx.dumps()
    again = Shelxfile()
    again.read_string(once)
    assert again.dumps() == once
    assert len(again.sfac_table.sfac_table) == 3


@pytest.mark.parametrize('element', ['C', 'H', 'O'])
def test_mixed_forms_keep_every_element(element) -> None:
    shx = _shx(EXPLICIT + '\nSFAC H O', 'UNIT 4 8 2')
    assert shx.sfac_table.has_element(element)
