"""Atom parameters written in SHELXL's ``10*m + p`` form.

*"To fix any atom parameter, add 10."* -- so a coordinate held at 0.6666
is written ``10.666600``.  Such a line is an ordinary atom, but the value
sits far outside the plain fractional range that
:meth:`Shelxfile.is_atom` uses to tell atoms from other instructions.

Getting this wrong is expensive twice over: the atom is not parsed, *and*
its ``=`` continuation line is then swallowed as if it belonged to the
line above, so a second atom disappears with it.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit import ShelxDocument

HEADER = (
    'TITL fixed coordinates\n'
    'CELL 0.71073 10.0 11.0 12.0 90 90 90\n'
    'ZERR 4 0.001 0.001 0.001 0 0 0\n'
    'LATT -1\n'
    'SFAC C N\n'
    'UNIT 8 4\n'
    'FVAR 0.40559 0.84667\n'
)
FOOTER = 'HKLF 4\nEND\n'

#: An atom on a special position: all three coordinates fixed, and an
#: anisotropic U list long enough to need a continuation line.
FIXED_ANISO = (
    'MG1   1   10.666600   10.333330   10.833330    10.16667    0.00743    0.00743 =\n'
    '         0.00400  -10.00000  -10.00000    0.00371\n'
)
PLAIN_ANISO = (
    'N1    2    0.557056    0.407220    0.776816    11.00000    0.00863    0.00963 =\n'
    '         0.00476   -0.00022    0.00010    0.00555\n'
)


def _shx(body: str) -> Shelxfile:
    shx = Shelxfile()
    shx.read_string(HEADER + body + FOOTER)
    return shx


# ------------------------------------------------------------ recognition

@pytest.mark.parametrize('line', [
    'MG1   1   10.666600   10.333330   10.833330    10.16667    0.05',
    'C1    1    0.100000   10.250000    0.300000    11.00000    0.05',
    'C2    1  -10.000000    0.200000    0.300000    11.00000    0.05',
    'C3    1   20.500000    0.200000    0.300000    11.00000    0.05',
    # Coordinates tied to free variables 2, 3 and 4 with a coefficient of
    # 1.0 -- there is no free variable 1, so no ambiguity with a sof.
    'C4    1   21.000000   31.000000   41.000000    10.25000    0.05',
])
def test_encoded_coordinates_are_recognised_as_atoms(line) -> None:
    assert Shelxfile.is_atom(line) is True


@pytest.mark.parametrize('line', [
    # A coordinate is missing, so the sof slid into the z column.
    # ``11.00000`` decodes to exactly 1.0, which is an occupancy far more
    # often than a coordinate fixed at the cell edge.
    'O1  4  0.120080    0.494426  11.00000   0.01445 ...',
    # Decodes to 5.5/6.6/7.7 -- no fractional coordinate looks like that.
    'XX  1  55.5  66.6  77.7  11.00000  0.05',
    # Below the encoding range and far outside the plain one.
    'XX  1  8.0  8.0  8.0  11.00000  0.05',
])
def test_implausible_lines_are_still_rejected(line) -> None:
    assert Shelxfile.is_atom(line) is False


def test_both_atom_checks_agree(  ) -> None:
    """``is_atom`` and the fast ``is_atom_spline`` must not diverge."""
    lines = [
        'MG1   1   10.666600   10.333330   10.833330    10.16667    0.05',
        'C1    1    0.100000    0.200000    0.300000    11.00000    0.05',
        'O1  4  0.120080    0.494426  11.00000   0.01445 ...',
    ]
    for line in lines:
        spline = line.split()
        assert Shelxfile.is_atom(line) == Shelxfile.is_atom_spline(
            line[:4].upper(), spline
        ), line


# --------------------------------------------------------------- parsing

def test_fixed_coordinates_are_decoded() -> None:
    atom = _shx(FIXED_ANISO).atoms.get_atom_by_name('MG1')
    assert atom is not None
    assert atom.frac_coords == pytest.approx((0.6666, 0.33333, 0.83333))


def test_fixed_coordinates_are_written_back_encoded() -> None:
    """Dropping the code would silently release the constraint."""
    shx = _shx(FIXED_ANISO)
    written = str(shx.atoms.get_atom_by_name('MG1'))
    assert '10.666600' in written
    assert '10.333330' in written
    assert '10.833330' in written


def test_a_moved_atom_keeps_its_fixing_code() -> None:
    """The value moves; the ``m`` that fixes it does not."""
    shx = _shx(FIXED_ANISO)
    atom = shx.atoms.get_atom_by_name('MG1')
    atom.x = 0.25
    assert '10.250000' in str(atom)


def test_free_variable_coordinates_round_trip() -> None:
    """``21.000000`` means ``1.0 * fv2`` and must be written back as such."""
    body = (
        'C2    1   21.000000   31.000000   41.000000    10.25000    0.01848    0.02578 =\n'
        '         0.02696    0.00000    0.02039    0.00000\n'
    )
    shx = _shx(body)
    atom = shx.atoms.get_atom_by_name('C2')
    assert atom is not None
    written = str(atom)
    assert '21.000000' in written
    assert '31.000000' in written
    assert '41.000000' in written


def test_an_unparsed_atom_would_swallow_the_next_instruction() -> None:
    """The continuation ``=`` eats whatever follows if the atom is missed."""
    body = (
        'C2    1   21.000000   31.000000   41.000000    10.25000    0.01848    0.02578 =\n'
        '         0.02696    0.00000    0.02039    0.00000\n'
        'AFIX 13\n'
        'H2    2    0.446417    0.341452    0.043044    10.25000   -1.20000\n'
        'AFIX 0\n'
    )
    shx = _shx(body)
    assert 'AFIX 13' in shx.dumps()
    assert [a.name for a in shx.atoms] == ['C2', 'H2']


def test_plain_coordinates_are_written_plainly() -> None:
    written = str(_shx(PLAIN_ANISO).atoms.get_atom_by_name('N1'))
    assert '0.557056' in written
    assert '10.557056' not in written


# ---------------------------------------------- the knock-on consequence

def test_the_following_atom_is_not_swallowed() -> None:
    """An unparsed atom's ``=`` would eat the next line on re-read."""
    shx = _shx(FIXED_ANISO + PLAIN_ANISO)
    assert [a.name for a in shx.atoms] == ['MG1', 'N1']


def test_round_trip_keeps_every_atom() -> None:
    first = _shx(FIXED_ANISO + PLAIN_ANISO)
    once = first.dumps()
    second = Shelxfile()
    second.read_string(once)
    assert [a.name for a in second.atoms] == ['MG1', 'N1']
    assert second.dumps() == once


def test_document_sees_the_fixed_atom() -> None:
    doc = ShelxDocument.from_string(HEADER + FIXED_ANISO + PLAIN_ANISO + FOOTER)
    atom = doc.shelxfile.atoms.get_atom_by_name('MG1')
    assert doc.line_of(atom) is not None
