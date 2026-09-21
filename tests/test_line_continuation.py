"""Line continuation with a trailing ``=``.

*"Continuation lines are flagged by '=' at the end of a line, the
instruction being continued on the next line which must start with one or
more spaces."*

The qualifier matters: *"all characters following '!' or '=' in an
instruction line are ignored"*, so a ``=`` in the middle of a line is not
a continuation marker.  Treating it as one makes a line such as
``TITL foo  R = 0.40`` glue the following line onto itself -- and since
the next line is very often the ``created by SHELXL`` comment, that
comment is lost on every round-trip.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.misc.misc import multiline_test

HEADER_TAIL = (
    'CELL 0.71073 10.0 11.0 12.0 90 90 90\n'
    'ZERR 4 0.001 0.001 0.001 0 0 0\n'
    'LATT 1\n'
    'SFAC C\n'
    'UNIT 4\n'
)
BODY = (
    'C1    1    0.100000    0.200000    0.300000    11.00000    0.0500\n'
    'HKLF 4\n'
    'END\n'
)


@pytest.mark.parametrize('line', [
    'O1 3 -0.01453 1.66590 0.10966 11.00 0.05 =',
    'SADI 0.02 C1 C2 C3 C4 =',
    'C1 1 0.1 0.2 0.3 11.0 0.05 0.05 =   ',
])
def test_trailing_equals_continues(line) -> None:
    assert multiline_test(line) is True


@pytest.mark.parametrize('line', [
    'TITL foo      P n m a      R = 0.40 New:Pnma',
    'CELL 0.71073 10 11 12 90 90 90',
    'SADI 0.02 C1 C2',
])
def test_an_equals_in_the_middle_does_not_continue(line) -> None:
    assert multiline_test(line) is False


def test_rem_is_never_a_continuation() -> None:
    """*"A '=' character in a rem line is not a line break."*"""
    assert multiline_test('REM wR2 = 0.12 for all data =') is False


def test_a_title_containing_equals_keeps_the_next_line() -> None:
    """The regression this fix exists for."""
    source = (
        'TITL converted      P n m a      R = 0.40 New:Pnma\n'
        '    created by SHELXL-2019/3 at 10:36:25 on 26-Nov-2025\n'
        + HEADER_TAIL + BODY
    )
    shx = Shelxfile()
    shx.read_string(source)
    assert 'created by SHELXL-2019/3' in shx.dumps()


def test_comment_lines_survive_repeated_round_trips() -> None:
    """One comment was lost per round-trip, so two round-trips lost two."""
    source = (
        'TITL converted      P n m a      R = 0.40 New:Pnma\n'
        '    converted_a_pl.res\n'
        '    created by SHELXL-2019/3 at 10:36:25 on 26-Nov-2025\n'
        + HEADER_TAIL + BODY
    )
    text = source
    for _ in range(3):
        shx = Shelxfile()
        shx.read_string(text)
        text = shx.dumps()
    assert 'converted_a_pl.res' in text
    assert 'created by SHELXL-2019/3' in text


def test_a_real_continuation_is_still_joined() -> None:
    """An anisotropic atom spans two physical lines but is one atom."""
    source = (
        'TITL t\n' + HEADER_TAIL
        + 'C1    1    0.100000    0.200000    0.300000    11.00000    0.05 =\n'
          '         0.05    0.05    0.00    0.00    0.00\n'
          'HKLF 4\nEND\n'
    )
    shx = Shelxfile()
    shx.read_string(source)
    atom = shx.atoms.get_atom_by_name('C1')
    assert atom is not None
    assert atom.is_anisotropic
    assert len(shx.atoms) == 1
