"""Atom references outside restraints (plan decision D-1).

Restraints are not the only instructions that name atoms.  ``HFIX``,
``ANIS``, ``BOND``, ``MPLA``, ``CONN``, ``CONF``, ``BLOC``, ``BIND``,
``FREE`` and ``HTAB`` do too, and would be left pointing at deleted atoms
if only restraints were maintained.

``AtomReferencingCard`` gives them all the same classification surface.
The shapes differ -- a flat ``atoms`` list, or named fields such as
``FREE``'s ``atom1``/``atom2`` -- so ``referenced_atoms`` normalises
access.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit.card_meta import (
    AtomGrouping,
    AtomListSemantics,
    AtomReferencingCard,
)
from shelxfile.shelx import cards

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
)
FOOTER = 'HKLF 4\nEND\n'

COVERED = ['BOND', 'ANIS', 'CONF', 'CONN', 'BLOC', 'MPLA',
           'HFIX', 'FREE', 'HTAB', 'BIND']


def _shx(instructions: str) -> Shelxfile:
    shx = Shelxfile()
    shx.read_string(HEADER + instructions + ATOMS + FOOTER)
    return shx


# ------------------------------------------------------------- coverage

@pytest.mark.parametrize('name', COVERED)
def test_card_implements_the_protocol(name) -> None:
    assert issubclass(getattr(cards, name), AtomReferencingCard)


@pytest.mark.parametrize('name', COVERED + ['Restraint'])
def test_protocol_surface_is_complete(name) -> None:
    card_class = getattr(cards, name)
    for attribute in ('ATOM_GROUPING', 'MIN_ATOMS', 'EMPTY_MEANS',
                      'referenced_atoms', 'atom_grouping',
                      'atom_semantics', 'is_atom_linked'):
        assert hasattr(card_class, attribute), f'{name} lacks {attribute}'


def test_restraints_use_the_same_protocol() -> None:
    assert issubclass(cards.Restraint, AtomReferencingCard)


# ------------------------------------------------- normalised atom access

def test_flat_list_cards_expose_their_atoms() -> None:
    shx = _shx('HFIX 43 C1 C2\n')
    assert shx.hfixes[0].referenced_atoms == ['C1', 'C2']


def test_free_exposes_its_named_operands() -> None:
    """``FREE`` stores ``atom1``/``atom2``, not a list."""
    shx = _shx('FREE C1 C2\n')
    assert shx.free[0].referenced_atoms == ['C1', 'C2']


def test_htab_exposes_donor_then_acceptor() -> None:
    """Order matters: only the acceptor may carry a symmetry operation."""
    shx = _shx('HTAB C1 C2\n')
    assert shx.htab.referenced_atoms == ['C1', 'C2']


def test_symmetry_suffixed_operands_are_preserved() -> None:
    shx = _shx('EQIV $1 1-x, y, 1-z\nFREE C1 C2_$1\n')
    assert shx.free[0].referenced_atoms == ['C1', 'C2_$1']


# --------------------------------------------- B2 semantics for these cards

@pytest.mark.parametrize('instruction, card_name', [
    ('ANIS\n', 'ANIS'),
    ('BOND\n', 'BOND'),
])
def test_bare_card_means_all_atoms(instruction, card_name) -> None:
    shx = _shx(instruction)
    card = getattr(cards, card_name)
    found = [i for i in shx._reslist if isinstance(i, card)]
    assert found[0].atom_semantics is AtomListSemantics.GLOBAL_WHEN_EMPTY
    assert found[0].is_atom_linked is False


def test_bind_with_part_numbers_is_a_different_instruction() -> None:
    """``BIND m n`` links two PARTs and names no atoms."""
    shx = _shx('BIND 1 2\n')
    assert shx.bind[0].referenced_atoms == []
    assert shx.bind[0].atom_semantics is AtomListSemantics.DIRECTIVE_WHEN_EMPTY
    assert shx.bind[0].is_atom_linked is False


def test_bind_with_atoms_is_explicit() -> None:
    shx = _shx('BIND C1 C2\n')
    assert shx.bind[0].referenced_atoms == ['C1', 'C2']
    assert shx.bind[0].is_atom_linked is True


def test_bare_htab_sets_a_distance_and_names_nothing() -> None:
    shx = _shx('HTAB 2.2\n')
    assert shx.htab.referenced_atoms == []
    assert shx.htab.is_atom_linked is False


def test_bare_conn_changes_a_default() -> None:
    shx = _shx('CONN 0\n')
    found = [i for i in shx._reslist if isinstance(i, cards.CONN)]
    assert found[0].atom_semantics is AtomListSemantics.DIRECTIVE_WHEN_EMPTY


def test_named_cards_are_atom_linked() -> None:
    shx = _shx('HFIX 43 C1\n')
    assert shx.hfixes[0].is_atom_linked is True


# ------------------------------------------------------------- grouping

@pytest.mark.parametrize('name, expected', [
    ('FREE', AtomGrouping.PAIRS),
    ('BIND', AtomGrouping.PAIRS),
    ('HTAB', AtomGrouping.PAIRS),
    ('HFIX', AtomGrouping.FLAT),
    ('ANIS', AtomGrouping.FLAT),
    ('MPLA', AtomGrouping.FLAT),
])
def test_grouping_per_card(name, expected) -> None:
    assert getattr(cards, name).ATOM_GROUPING is expected


def test_mpla_needs_at_least_three_atoms() -> None:
    """"na must be at least 3" for a least-squares plane."""
    assert cards.MPLA.MIN_ATOMS == 3


def test_hfix_needs_at_least_one_atom() -> None:
    assert cards.HFIX.MIN_ATOMS == 1


def test_bond_pairs_are_not_pairs() -> None:
    """``BOND`` lists atoms of interest, not bonded couples."""
    assert cards.BOND.ATOM_GROUPING is AtomGrouping.FLAT


# ------------------------------------------------------------ regression

def test_adding_the_protocol_did_not_break_parsing() -> None:
    shx = Shelxfile()
    shx.read_file('tests/resources/p21c.res')
    assert len(shx.atoms) > 0
    assert shx.restraint_errors == [] or all(
        isinstance(w, str) for w in shx.restraint_errors
    )
