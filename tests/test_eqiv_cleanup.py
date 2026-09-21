"""``EQIV`` housekeeping and count-coupled cards (plan items D-6c/d, B3).

Two small but easy-to-get-wrong pieces of consistency:

* a ``C1_$3`` token depends on ``EQIV $3`` as well as on ``C1``, so an
  ``EQIV`` can be left behind with nothing pointing at it;
* ``MPLA na`` counts how many of the named atoms define the plane, so
  ``na`` is coupled to the atom list rather than independent of it.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit import ShelxDocument
from shelxfile.edit.eqiv_cleanup import (
    EqivCleaner,
    symmetry_suffix,
    validate_symmetry_arity,
)
from shelxfile.edit.graph import AtomRestraintGraph
from shelxfile.shelx.cards import MPLA

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
    'C3    1     0.20  0.30  0.40  11.00000  0.05\n'
    'C4    1     0.25  0.35  0.45  11.00000  0.05\n'
)
FOOTER = 'HKLF 4\nEND\n'


def _doc(body: str) -> ShelxDocument:
    return ShelxDocument.from_string(HEADER + body + ATOMS + FOOTER)


def _cleaner(body: str) -> tuple[Shelxfile, EqivCleaner]:
    shx = Shelxfile()
    shx.read_string(HEADER + body + ATOMS + FOOTER)
    graph = AtomRestraintGraph(shx)
    graph.rebuild()
    return shx, EqivCleaner(shx, graph)


# -------------------------------------------------------- suffix detection

@pytest.mark.parametrize('token, expected', [
    ('C1_$3', '$3'),
    ('C1_2_$11', '$11'),
    ('C1_2', None),
    ('C1', None),
    ('C1_*', None),
])
def test_symmetry_suffix_detection(token, expected) -> None:
    assert symmetry_suffix(token) == expected


# ---------------------------------------------------------- orphan removal

def test_unreferenced_eqiv_is_found() -> None:
    shx, cleaner = _cleaner(
        'EQIV $1 1-x, y, 1-z\n'
        'EQIV $2 x, 1-y, z\n'
        'DFIX 1.5 C1 C2_$1\n'
    )
    assert [e.id for e in cleaner.orphans()] == ['$2']


def test_removing_orphans_reports_the_reason() -> None:
    shx, cleaner = _cleaner(
        'EQIV $1 1-x, y, 1-z\n'
        'EQIV $2 x, 1-y, z\n'
        'DFIX 1.5 C1 C2_$1\n'
    )
    report = cleaner.remove_orphans()
    assert len(report.cards) == 1
    assert report.cards[0].reason.name == 'EQIV_UNREFERENCED'
    assert 'EQIV $2' not in shx.dumps()
    assert 'EQIV $1' in shx.dumps()


def test_deleting_the_last_user_collects_the_eqiv() -> None:
    doc = _doc('EQIV $1 1-x, y, 1-z\nDFIX 1.5 C1 C2_$1\n')
    assert 'EQIV $1' in doc.text

    restraint = next(iter(doc.shelxfile.restraints))
    report = doc.remove_restraint(restraint)

    assert 'EQIV $1' not in doc.text
    assert any(c.reason.name == 'EQIV_UNREFERENCED' for c in report.cards)


def test_an_eqiv_with_other_users_survives() -> None:
    doc = _doc(
        'EQIV $1 1-x, y, 1-z\n'
        'DFIX 1.5 C1 C2_$1\n'
        'SADI C3 C4_$1\n'
    )
    doc.remove_restraint(next(iter(doc.shelxfile.restraints)))
    assert 'EQIV $1' in doc.text, 'another card still references $1'


def test_ids_are_never_renumbered() -> None:
    """A surviving reference must keep meaning the same operation."""
    doc = _doc(
        'EQIV $1 1-x, y, 1-z\n'
        'EQIV $2 x, 1-y, z\n'
        'EQIV $3 -x, -y, -z\n'
        'DFIX 1.5 C1 C2_$3\n'
    )
    doc.delete_atom(doc.shelxfile.atoms.get_atom_by_name('C4_0'))
    remaining = [e.id for e in doc.shelxfile.eqiv]
    assert '$3' in remaining, 'the referenced id must be untouched'
    assert 'EQIV $3 -x, -y, -z' in doc.text


def test_post_end_eqivs_are_left_alone() -> None:
    shx = Shelxfile()
    shx.read_string(
        HEADER + ATOMS + FOOTER + 'EQIV $9 1-x, y, 1-z\nHTAB C1 C2_$9\n'
    )
    graph = AtomRestraintGraph(shx)
    graph.rebuild()
    report = EqivCleaner(shx, graph).remove_orphans()
    assert report.is_empty
    assert 'EQIV $9' in shx.dumps()


# ------------------------------------------------------- defined before use

def test_eqiv_defined_before_its_user() -> None:
    shx, cleaner = _cleaner('EQIV $1 1-x, y, 1-z\nDFIX 1.5 C1 C2_$1\n')
    assert cleaner.is_defined_before_use(shx.eqiv[0])


def test_first_use_index_points_at_the_referencing_card() -> None:
    shx, cleaner = _cleaner('EQIV $1 1-x, y, 1-z\nDFIX 1.5 C1 C2_$1\n')
    restraint = next(iter(shx.restraints))
    assert cleaner.first_use_index(shx.eqiv[0]) == shx.index_of(restraint)


def test_unused_eqiv_has_no_first_use() -> None:
    shx, cleaner = _cleaner('EQIV $1 1-x, y, 1-z\n')
    assert cleaner.first_use_index(shx.eqiv[0]) is None
    assert cleaner.is_defined_before_use(shx.eqiv[0])


# ---------------------------------------------------- D-6d: arity limits

def test_free_rejects_two_symmetry_operands() -> None:
    """"Only one of the two atoms may be an equivalent atom"."""
    shx = Shelxfile()
    shx.read_string(
        HEADER + 'EQIV $1 1-x, y, 1-z\nFREE C1_$1 C2_$1\n' + ATOMS + FOOTER
    )
    assert validate_symmetry_arity(shx.free[0]) is not None


def test_free_accepts_one_symmetry_operand() -> None:
    shx = Shelxfile()
    shx.read_string(
        HEADER + 'EQIV $1 1-x, y, 1-z\nFREE C1 C2_$1\n' + ATOMS + FOOTER
    )
    assert validate_symmetry_arity(shx.free[0]) is None


def test_htab_allows_symmetry_only_on_the_acceptor() -> None:
    """"Only the acceptor atom may specify a symmetry operation"."""
    shx = Shelxfile()
    shx.read_string(
        HEADER + 'EQIV $1 1-x, y, 1-z\nHTAB C1_$1 C2\n' + ATOMS + FOOTER
    )
    assert validate_symmetry_arity(shx.htab) is not None


def test_htab_acceptor_symmetry_is_fine() -> None:
    shx = Shelxfile()
    shx.read_string(
        HEADER + 'EQIV $1 1-x, y, 1-z\nHTAB C1 C2_$1\n' + ATOMS + FOOTER
    )
    assert validate_symmetry_arity(shx.htab) is None


def test_cards_without_a_limit_are_not_flagged() -> None:
    shx = Shelxfile()
    shx.read_string(
        HEADER + 'EQIV $1 1-x, y, 1-z\nSADI C1_$1 C2_$1\n' + ATOMS + FOOTER
    )
    assert validate_symmetry_arity(next(iter(shx.restraints))) is None


# ------------------------------------------------ B3: MPLA count coupling

def _mpla(instruction: str) -> MPLA:
    shx = Shelxfile()
    shx.read_string(HEADER + instruction + ATOMS + FOOTER)
    return next(i for i in shx._reslist if isinstance(i, MPLA))


def test_na_within_the_atom_list_is_consistent() -> None:
    assert _mpla('MPLA 3 C1 C2 C3 C4\n').na_is_consistent


def test_na_larger_than_the_atom_list_is_not() -> None:
    card = _mpla('MPLA 4 C1 C2 C3 C4\n')
    card.atoms.remove('C4')
    assert card.na_is_consistent is False


def test_na_below_three_is_not_consistent() -> None:
    """"na must be at least 3"."""
    assert _mpla('MPLA 2 C1 C2 C3\n').na_is_consistent is False


def test_clamping_brings_na_back_into_range() -> None:
    card = _mpla('MPLA 4 C1 C2 C3 C4\n')
    card.atoms.remove('C4')
    assert card.clamp_na() is True
    assert card.na == 3
    assert card.na_is_consistent


def test_clamping_a_consistent_card_is_a_no_op() -> None:
    card = _mpla('MPLA 3 C1 C2 C3 C4\n')
    assert card.clamp_na() is False
    assert card.na == 3


def test_clamping_fails_when_too_few_atoms_remain() -> None:
    """Below three atoms there is no plane to fit; the card must go."""
    card = _mpla('MPLA 4 C1 C2 C3 C4\n')
    card.atoms[:] = ['C1', 'C2']
    assert card.clamp_na() is False
    assert card.na_is_consistent is False


def test_clamping_rewrites_the_rendered_line() -> None:
    card = _mpla('MPLA 4 C1 C2 C3 C4\n')
    card.atoms.remove('C4')
    card.clamp_na()
    assert str(card).split() == ['MPLA', '3', 'C1', 'C2', 'C3']


def test_missing_na_is_always_consistent() -> None:
    """"If na is omitted the plane is fitted to all the atoms specified"."""
    assert _mpla('MPLA C1 C2 C3\n').na_is_consistent
