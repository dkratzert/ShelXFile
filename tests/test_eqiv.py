"""``EQIV`` is a real card, not a raw string (plan item D-6a).

``EQIV $n <symop>`` defines a symmetry operation that other instructions
reference as ``C1_$2``.  ShelXFile used to stash these as bare
``spline[1:]`` lists and never call ``_append_card``, so the ``_reslist``
entry stayed a plain ``str``.  That made ``EQIV`` invisible to
identity-based removal and to the atom-linking graph.

Two manual rules are encoded here:

* *"Such a symmetry operation must be defined before it is used"*
* *"The same $n may not appear on two separate EQIV instructions"*
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.shelx.cards import EQIV

HEADER = (
    'TITL eqiv test\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C\n'
    'UNIT 4\n'
)
FOOTER = 'HKLF 4\nEND\n'


def _read(body: str) -> Shelxfile:
    shx = Shelxfile()
    shx.read_string(HEADER + body + FOOTER)
    return shx


@pytest.fixture
def shx() -> Shelxfile:
    return _read(
        'EQIV $1 1-x, y, 1-z\n'
        'EQIV $2 x+1/2, -y+1/2, z+1/2\n'
        'DFIX 1.5 C1 C2_$2\n'
        'C1    1     0.1  0.2  0.3  11.00000  0.05\n'
        'C2    1     0.4  0.5  0.6  11.00000  0.05\n'
    )


# ---------------------------------------------------------------- parsing

def test_eqiv_cards_are_objects(shx: Shelxfile) -> None:
    assert len(shx.eqiv) == 2
    assert all(isinstance(card, EQIV) for card in shx.eqiv)


def test_eqiv_is_present_in_the_reslist(shx: Shelxfile) -> None:
    """Required for identity-based removal to reach it at all."""
    for card in shx.eqiv:
        assert any(item is card for item in shx._reslist)


def test_eqiv_exposes_id_number_and_symmop(shx: Shelxfile) -> None:
    first = shx.eqiv[0]
    assert first.id == '$1'
    assert first.number == 1
    assert first.symmop == '1-x, y, 1-z'


def test_eqiv_renders_unchanged(shx: Shelxfile) -> None:
    assert str(shx.eqiv[1]) == 'EQIV $2 x+1/2, -y+1/2, z+1/2'


def test_eqiv_survives_a_roundtrip(shx: Shelxfile, tmp_path) -> None:
    out = tmp_path / 'eqiv.res'
    shx.write_shelx_file(out)
    reread = Shelxfile()
    reread.read_file(out)
    assert [c.id for c in reread.eqiv] == ['$1', '$2']
    assert [c.symmop for c in reread.eqiv] == [c.symmop for c in shx.eqiv]


# ------------------------------------------------------------- back-compat

def test_entries_still_index_like_the_old_token_lists(shx: Shelxfile) -> None:
    """``shx.eqiv`` held ``spline[1:]`` lists; ``entry[0]`` must still work."""
    assert shx.eqiv[0][0] == '$1'
    assert shx.eqiv[0][1] == '1-x,'
    assert len(shx.eqiv[0]) == 4


# ---------------------------------------------------------------- lookup

def test_lookup_accepts_both_token_and_number(shx: Shelxfile) -> None:
    assert shx.eqiv_by_id('$2') is shx.eqiv[1]
    assert shx.eqiv_by_id(2) is shx.eqiv[1]


def test_lookup_returns_none_for_unknown_ids(shx: Shelxfile) -> None:
    assert shx.eqiv_by_id('$9') is None
    assert shx.eqiv_by_id('$nonsense') is None


# ------------------------------------------------------------- validation

def test_reference_to_a_defined_eqiv_is_not_reported(shx: Shelxfile) -> None:
    assert not [w for w in shx.restraint_errors if 'EQIV' in w]


def test_reference_to_an_undefined_eqiv_is_reported() -> None:
    shx = _read(
        'EQIV $1 1-x, y, 1-z\n'
        'DFIX 1.5 C1 C2_$7\n'
        'C1    1     0.1  0.2  0.3  11.00000  0.05\n'
        'C2    1     0.4  0.5  0.6  11.00000  0.05\n'
    )
    assert any('EQIV' in w for w in shx.restraint_errors)


def test_duplicate_ids_are_reported(capsys: pytest.CaptureFixture) -> None:
    """The manual forbids the same $n on two EQIV instructions."""
    shx = Shelxfile(verbose=True)
    shx.read_string(
        HEADER
        + 'EQIV $1 1-x, y, 1-z\n'
          'EQIV $1 x, 1-y, z\n'
          'C1    1     0.1  0.2  0.3  11.00000  0.05\n'
        + FOOTER
    )
    assert 'Duplicate EQIV' in capsys.readouterr().out


def test_malformed_eqiv_is_ignored_not_crashing() -> None:
    shx = _read(
        'EQIV nonsense\n'
        'C1    1     0.1  0.2  0.3  11.00000  0.05\n'
    )
    assert shx.eqiv == []


# ---------------------------------------------------------------- removal

def test_eqiv_can_be_removed_by_identity(shx: Shelxfile) -> None:
    card = shx.eqiv[0]
    shx.remove_from_reslist(card)
    shx.eqiv.remove(card)
    assert not any(item is card for item in shx._reslist)
    assert 'EQIV $1' not in shx.dumps()
    assert 'EQIV $2' in shx.dumps(), 'removed the wrong EQIV'
