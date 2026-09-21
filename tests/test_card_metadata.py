"""Card classification metadata (plan items D-5, B2, B3).

Replaces the magic strings the original edit plan proposed with enums, so
the type checker can catch typos and a ``match`` can be checked for
exhaustiveness.

The important behaviour here is **B2**: an empty atom list is a valid,
meaningful state, and what it means depends on the card.  A minimum atom
count may only ever be applied to a card that actually named atoms.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit.card_meta import (
    AfixDependency,
    AtomGrouping,
    AtomListSemantics,
    CardLifetime,
)

HEADER = (
    'TITL t\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C O\n'
    'UNIT 4 4\n'
)
ATOMS = (
    'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
    'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
    'C3    1     0.20  0.30  0.40  11.00000  0.05\n'
    'C4    1     0.25  0.35  0.45  11.00000  0.05\n'
)
FOOTER = 'HKLF 4\nEND\n'


def _card(instruction: str):
    shx = Shelxfile()
    shx.read_string(HEADER + instruction + '\n' + ATOMS + FOOTER)
    return next(iter(shx.restraints))


# -------------------------------------------------------------- enum shape

def test_members_are_not_strings() -> None:
    """The point of the enums is to stop string comparison."""
    for member in (AtomGrouping.FLAT, AtomListSemantics.EXPLICIT,
                   AfixDependency.RIDING, CardLifetime.INPUT):
        assert not isinstance(member, str)


def test_members_compare_by_identity() -> None:
    assert AtomGrouping.PAIRS is AtomGrouping.PAIRS
    assert AtomGrouping.PAIRS != AtomGrouping.FLAT


# ------------------------------------------------------------- grouping

@pytest.mark.parametrize('instruction, expected', [
    ('DFIX 1.5 C1 C2', AtomGrouping.PAIRS),
    ('DANG 2.4 C1 C3', AtomGrouping.PAIRS),
    ('SADI C1 C2 C1 C3', AtomGrouping.PAIRS),
    ('SIMU C1 C2', AtomGrouping.FLAT),
    ('ISOR C1 C2', AtomGrouping.FLAT),
    ('EADP C1 C2', AtomGrouping.FLAT),
    ('FLAT C1 C2 C3 C4', AtomGrouping.FLAT),
])
def test_atom_grouping_per_card(instruction, expected) -> None:
    assert _card(instruction).atom_grouping is expected


def test_plain_same_is_matched_against_following_atoms() -> None:
    """"...compared with the same number of atoms which follow"."""
    assert _card('SAME C1 C2').atom_grouping is AtomGrouping.SAME_FOLLOWING_ATOMS


def test_same_with_residue_number_still_uses_following_atoms() -> None:
    assert _card('SAME_1 C1 C2').atom_grouping is AtomGrouping.SAME_FOLLOWING_ATOMS


def test_same_with_residue_class_uses_the_other_mode() -> None:
    """``SAME_XYZ`` "no longer uses the following atoms"."""
    assert _card('SAME_TOL C1 C2').atom_grouping is AtomGrouping.SAME_RESIDUE_CLASS


# ------------------------------------------------- B2: empty-list semantics

@pytest.mark.parametrize('card_name', ['SIMU', 'ISOR', 'DELU', 'RIGU'])
def test_bare_card_means_all_non_hydrogen_atoms(card_name) -> None:
    """1217 corpus files rely on this form."""
    card = _card(card_name)
    assert card.atom_semantics is AtomListSemantics.GLOBAL_WHEN_EMPTY


def test_bare_sadi_is_a_directive_not_a_restraint() -> None:
    """A bare SADI asks SHELXL to emit SAME-derived restraints after END."""
    assert _card('SADI').atom_semantics is AtomListSemantics.DIRECTIVE_WHEN_EMPTY


@pytest.mark.parametrize('instruction', [
    'SIMU C1 C2', 'ISOR C1 C2', 'DELU C1 C2', 'RIGU C1 C2', 'SADI C1 C2',
])
def test_naming_atoms_always_wins(instruction) -> None:
    """A card that named atoms is EXPLICIT regardless of its class."""
    assert _card(instruction).atom_semantics is AtomListSemantics.EXPLICIT


@pytest.mark.parametrize('card_name', ['SIMU', 'ISOR', 'DELU', 'RIGU', 'SADI'])
def test_bare_cards_are_not_atom_linked(card_name) -> None:
    """They reference no specific atoms, so deletion must leave them alone."""
    assert _card(card_name).is_atom_linked is False


def test_explicit_cards_are_atom_linked() -> None:
    assert _card('SIMU C1 C2').is_atom_linked is True


def test_post_end_cards_are_never_atom_linked() -> None:
    shx = Shelxfile()
    shx.read_string(
        HEADER + ATOMS + FOOTER + 'SADI 0.02 C1 C2\n'
    )
    post = [r for r in shx.restraints
            if r.lifetime is CardLifetime.POST_END_OUTPUT]
    assert post and all(not r.is_atom_linked for r in post)


# ----------------------------------------------------- B3: minimum atoms

@pytest.mark.parametrize('instruction, expected', [
    ('CHIV C1', 1),           # one chiral volume per named atom
    ('FLAT C1 C2 C3 C4', 4),  # "four or more atoms"
    ('EADP C1 C2', 2),
    ('EXYZ C1 C2', 2),
    ('DFIX 1.5 C1 C2', 2),
    ('SADI C1 C2', 2),
    ('NCSY 1 C1 C2', 1),
])
def test_minimum_atom_counts(instruction, expected) -> None:
    assert _card(instruction).MIN_ATOMS == expected


def test_chiv_minimum_is_one_not_four() -> None:
    """The original plan put CHIV at 4 alongside FLAT; only FLAT needs four."""
    assert _card('CHIV C1').MIN_ATOMS == 1
    assert _card('FLAT C1 C2 C3 C4').MIN_ATOMS == 4


# --------------------------------------------------------- AFIX dependency

@pytest.mark.parametrize('code, expected', [
    (0, AfixDependency.NONE),
    (1, AfixDependency.NONE),
    (2, AfixDependency.NONE),
    (13, AfixDependency.RIDING),     # AFIX 13: tertiary C-H, riding
    (43, AfixDependency.RIDING),     # aromatic C-H, riding
    (4, AfixDependency.RIDING),
    (5, AfixDependency.RIGID_DEPENDENT),
    (6, AfixDependency.RIGID_PIVOT),
    (66, AfixDependency.RIGID_PIVOT),  # rigid hexagon
    (9, AfixDependency.RIGID_PIVOT),
    (137, AfixDependency.ROTATING),  # rotating methyl
    (8, AfixDependency.ROTATING),
    (None, AfixDependency.NONE),
])
def test_afix_dependency_from_code(code, expected) -> None:
    assert AfixDependency.from_code(code) is expected


def test_only_the_last_digit_selects_the_constraint() -> None:
    """``m`` picks the geometry, ``n`` the constraint."""
    assert AfixDependency.from_code(3) is AfixDependency.from_code(123)


@pytest.mark.parametrize('dependency, outside', [
    (AfixDependency.RIDING, True),
    (AfixDependency.ROTATING, True),
    (AfixDependency.RIGID_PIVOT, False),
    (AfixDependency.RIGID_DEPENDENT, False),
    (AfixDependency.NONE, False),
])
def test_pivot_location(dependency, outside) -> None:
    """Riding groups depend on an atom before the card; rigid ones don't."""
    assert dependency.pivot_is_outside_bracket is outside


def test_afix_card_exposes_its_dependency() -> None:
    shx = Shelxfile()
    shx.read_string(
        HEADER
        + 'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
          'AFIX 43\n'
          'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
          'AFIX 0\n'
        + FOOTER
    )
    riding = [a for a in shx.atoms if a.afix and a.afix.mn == 43]
    assert riding, 'precondition: an AFIX 43 atom was parsed'
    assert riding[0].afix.dependency is AfixDependency.RIDING
