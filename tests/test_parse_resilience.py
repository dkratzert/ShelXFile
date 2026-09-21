"""One bad instruction must not cost the whole file.

``parse_cards`` used to let any exception out of ``_parse_cards``, catch
it, and return -- handing the caller a model silently missing every card
and atom below the failure.  Nothing indicated that anything was wrong.

Measured on the reference corpus before the fix: **456 of 5970 files
(7.6 %)** aborted before reaching their own ``HKLF`` line.  One
``TEMP -100.0C`` discarded six ``MPLA`` cards and 71 atoms.

SHELXL keeps reading past an instruction it cannot use, so ShelXFile now
records the line, quarantines it, and carries on.  The line stays in the
file as text so it still round-trips.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile

HEADER = (
    'TITL t\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C\n'
    'UNIT 4\n'
)
ATOMS = (
    'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
    'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
)
FOOTER = 'HKLF 4\nEND\n'
BAD = 'TEMP totally-not-a-number\n'


def _read(body: str, **kwargs) -> Shelxfile:
    shx = Shelxfile(**kwargs)
    shx.read_string(HEADER + body + FOOTER)
    return shx


# ------------------------------------------------- parsing continues

def test_atoms_after_a_bad_line_are_still_parsed() -> None:
    shx = _read(BAD + ATOMS)
    assert [a.name for a in shx.atoms] == ['C1', 'C2']


def test_cards_after_a_bad_line_are_still_parsed() -> None:
    shx = _read(BAD + 'SADI C1 C2\n' + ATOMS)
    assert len(list(shx.restraints)) == 1


def test_the_end_of_the_file_is_reached() -> None:
    shx = _read(BAD + ATOMS)
    assert shx.hklf is not None


def test_a_bad_line_early_on_does_not_hide_the_cell() -> None:
    shx = _read(BAD + ATOMS)
    assert shx.cell is not None


# --------------------------------------------------- the error is reported

def test_the_failure_is_recorded() -> None:
    shx = _read(BAD + ATOMS)
    assert len(shx.parse_errors) == 1
    assert 'TEMP' in shx.parse_errors[0]


def test_the_report_names_the_line_number() -> None:
    shx = _read(BAD + ATOMS)
    assert 'line 7' in shx.parse_errors[0]


def test_the_report_names_the_underlying_error() -> None:
    shx = _read(BAD + ATOMS)
    assert 'ValueError' in shx.parse_errors[0]


def test_a_clean_file_reports_nothing() -> None:
    shx = _read(ATOMS)
    assert shx.parse_errors == []


def test_verbose_mode_prints_the_problem(capsys: pytest.CaptureFixture) -> None:
    _read(BAD + ATOMS, verbose=True)
    assert 'Could not parse' in capsys.readouterr().out


def test_debug_mode_still_raises() -> None:
    """Debug exists to surface problems, not to paper over them."""
    shx = Shelxfile(debug=True)
    with pytest.raises(Exception):
        shx.read_string(HEADER + BAD + ATOMS + FOOTER)


# ------------------------------------------------- the line still round-trips

def test_the_bad_line_is_preserved_in_the_output() -> None:
    """Quarantined, not discarded: dropping it would lose user data."""
    shx = _read(BAD + ATOMS)
    assert 'TEMP totally-not-a-number' in shx.dumps()


def test_round_trip_keeps_the_bad_line(tmp_path) -> None:
    shx = _read(BAD + ATOMS)
    out = tmp_path / 'rt.res'
    shx.write_shelx_file(out)
    reread = Shelxfile()
    reread.read_file(out)
    assert 'TEMP totally-not-a-number' in reread.dumps()
    assert [a.name for a in reread.atoms] == ['C1', 'C2']


# ------------------------------------------------------- several bad lines

def test_multiple_bad_lines_are_each_reported() -> None:
    shx = _read('TEMP nonsense\n' + 'ZERR also nonsense here\n' + ATOMS)
    assert len(shx.parse_errors) >= 1
    assert [a.name for a in shx.atoms] == ['C1', 'C2']


def test_a_file_of_nothing_but_junk_gives_up() -> None:
    """Past a point, this is not a SHELXL file and should not be retried forever."""
    junk = ''.join(f'TEMP bad{i}\n' for i in range(Shelxfile.MAX_UNPARSABLE_LINES + 5))
    shx = _read(junk + ATOMS)
    assert len(shx.parse_errors) <= Shelxfile.MAX_UNPARSABLE_LINES + 1


# --------------------------------------------------------- TEMP specifics

@pytest.mark.parametrize('value, expected', [
    ('20', 20.0),
    ('-100.0', -100.0),
    ('-100.0C', -100.0),      # stray unit, seen in real files
    ('-173.15c', -173.15),
    ('20(2)', 20.0),          # esd in parentheses
])
def test_temp_tolerates_how_it_is_written(value, expected) -> None:
    shx = _read(f'TEMP {value}\n' + ATOMS)
    assert shx.temp == pytest.approx(expected)
    assert shx.parse_errors == []


def test_temp_in_kelvin_follows() -> None:
    shx = _read('TEMP -100.0C\n' + ATOMS)
    assert shx.temp_in_kelvin == pytest.approx(173.15)


# ------------------------------------------------ HKLF after an open AFIX

def test_hklf_is_parsed_even_with_an_afix_still_open() -> None:
    """Closing a dangling ``AFIX`` must not swallow the dispatch.

    Reaching ``HKLF`` with an ``AFIX`` still open used to match the
    "close the group" guard, and because that guard headed an ``elif``
    chain the ``HKLF`` branch below never ran.  ``shx.hklf`` stayed
    ``None`` in 411 corpus files, so callers saw a structure with no
    reflection-data instruction at all.
    """
    shx = Shelxfile()
    shx.read_string(
        HEADER
        + 'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
          'AFIX 43\n'
          'H1    1     0.15  0.25  0.35  11.00000  0.05\n'
        + FOOTER          # HKLF with the AFIX never closed
    )
    assert shx.hklf is not None


def test_the_dangling_afix_is_still_closed() -> None:
    """The original purpose of the guard must survive the fix."""
    shx = Shelxfile()
    shx.read_string(
        HEADER
        + 'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
          'AFIX 43\n'
          'H1    1     0.15  0.25  0.35  11.00000  0.05\n'
        + FOOTER
    )
    assert shx.afix.mn == 0


def test_hklf_after_a_properly_closed_afix_still_works() -> None:
    shx = Shelxfile()
    shx.read_string(
        HEADER
        + 'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
          'AFIX 43\n'
          'H1    1     0.15  0.25  0.35  11.00000  0.05\n'
          'AFIX 0\n'
        + FOOTER
    )
    assert shx.hklf is not None


def test_afix_cards_are_still_parsed_normally() -> None:
    shx = Shelxfile()
    shx.read_string(
        HEADER
        + 'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
          'AFIX 137\n'
          'H1    1     0.15  0.25  0.35  11.00000  0.05\n'
          'AFIX 0\n'
        + FOOTER
    )
    riding = [a for a in shx.atoms if a.afix and a.afix.mn == 137]
    assert riding, 'the AFIX branch must still run for ordinary AFIX lines'
