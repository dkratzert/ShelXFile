"""Atom/card link graph (plan items D-1, D-6c, C3).

The graph is *derived* data.  Rather than hooking every mutation -- which
risks a stale graph, and a stale graph is worse than none -- it records
the model version it was built at and rebuilds when that moves on.

Two link directions matter beyond the obvious one:

* a ``SAME`` reaches atoms it never names (see B1), and
* a ``C1_$3`` reference depends on ``EQIV $3`` as well as on ``C1``,
  so an ``EQIV`` can become orphaned (see D-6c).
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit import ShelxDocument
from shelxfile.edit.graph import AtomRestraintGraph

HEADER = (
    'TITL t\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C H\n'
    'UNIT 8 8\n'
)
ATOMS = (
    'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
    'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
    'C3    1     0.20  0.30  0.40  11.00000  0.05\n'
)
FOOTER = 'HKLF 4\nEND\n'


def _graph(body: str) -> tuple[Shelxfile, AtomRestraintGraph]:
    shx = Shelxfile()
    shx.read_string(HEADER + body + FOOTER)
    return shx, AtomRestraintGraph(shx)


@pytest.fixture
def doc() -> ShelxDocument:
    return ShelxDocument.from_file('tests/resources/p21c.res')


# ------------------------------------------------------------ basic links

def test_atom_maps_to_the_cards_that_name_it() -> None:
    shx, graph = _graph('DFIX 1.5 C1 C2\n' + ATOMS)
    cards = graph.cards_for_atom('C1_0')
    assert len(cards) == 1
    assert type(cards[0]).__name__ == 'DFIX'


def test_card_maps_to_its_atoms() -> None:
    shx, graph = _graph('DFIX 1.5 C1 C2\n' + ATOMS)
    card = next(iter(shx.restraints))
    assert graph.atoms_for_card(card) == ['C1_0', 'C2_0']


def test_unreferenced_atom_has_no_cards() -> None:
    shx, graph = _graph('DFIX 1.5 C1 C2\n' + ATOMS)
    assert graph.cards_for_atom('C3_0') == []


def test_ranges_are_expanded_into_links() -> None:
    shx, graph = _graph('SIMU C1 > C3\n' + ATOMS)
    assert graph.cards_for_atom('C2_0'), 'the middle of a range must be linked'


def test_accepts_atom_objects_and_names() -> None:
    shx, graph = _graph('DFIX 1.5 C1 C2\n' + ATOMS)
    atom = shx.atoms.get_atom_by_name('C1_0')
    assert graph.cards_for_atom(atom) == graph.cards_for_atom('C1_0')


# --------------------------------------------------- B2/B9 exclusions

def test_bare_cards_are_not_indexed() -> None:
    """A card meaning "all non-hydrogen atoms" links to nothing specific."""
    shx, graph = _graph('SIMU\n' + ATOMS)
    assert graph.cards == []


def test_post_end_cards_are_not_indexed() -> None:
    shx = Shelxfile()
    shx.read_string(HEADER + ATOMS + FOOTER + 'SADI 0.02 C1 C2\n')
    graph = AtomRestraintGraph(shx)
    assert graph.cards == []


# -------------------------------------------------- B1: SAME reaches out

def test_same_links_atoms_it_never_names() -> None:
    """The B1 case: the counterpart fragment is not on the card."""
    shx, graph = _graph(
        'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
        'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
        'SAME C1 C2\n'
        'C1B   1     0.40  0.50  0.60  11.00000  0.05\n'
        'C2B   1     0.45  0.55  0.65  11.00000  0.05\n'
    )
    card = next(iter(shx.restraints))
    assert 'C1B' not in card.atoms
    assert 'C1B_0' in graph.atoms_for_card(card)
    assert graph.cards_for_atom('C1B_0') == [card]


def test_deleting_an_unnamed_same_partner_is_visible() -> None:
    """Which is the whole reason the graph tracks them."""
    shx, graph = _graph(
        'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
        'SAME C1\n'
        'C1B   1     0.40  0.50  0.60  11.00000  0.05\n'
    )
    assert graph.cards_for_atom('C1B_0'), (
        'an atom following SAME must be reported as affected'
    )


# ----------------------------------------------- D-6c: EQIV back-references

def _eqiv_body() -> str:
    return (
        'EQIV $1 1-x, y, 1-z\n'
        'EQIV $2 x, 1-y, z\n'
        'DFIX 1.5 C1 C2_$1\n'
        + ATOMS
    )


def test_eqiv_reference_is_indexed() -> None:
    shx, graph = _graph(_eqiv_body())
    cards = graph.cards_for_eqiv('$1')
    assert len(cards) == 1
    assert type(cards[0]).__name__ == 'DFIX'


def test_unused_eqiv_is_reported_as_orphaned() -> None:
    shx, graph = _graph(_eqiv_body())
    orphans = [e.id for e in graph.unreferenced_eqivs()]
    assert orphans == ['$2']


def test_eqiv_becomes_orphaned_when_its_last_user_goes() -> None:
    shx, graph = _graph(_eqiv_body())
    assert '$1' not in [e.id for e in graph.unreferenced_eqivs()]

    next(iter(shx.restraints)).delete()

    assert '$1' in [e.id for e in graph.unreferenced_eqivs()]


def test_symmetry_reference_still_links_the_base_atom() -> None:
    shx, graph = _graph(_eqiv_body())
    assert graph.cards_for_atom('C2_0'), 'C2_$1 depends on C2 as well as $1'


def test_post_end_eqivs_are_never_reported_as_orphaned() -> None:
    """SHELXL writes generated ``HTAB``/``EQIV`` pairs after ``END``.

    Both are inert output, so the ``EQIV`` only *looks* unreferenced --
    the ``HTAB`` using it is excluded from the graph too.  Reporting it
    would invite an edit to a block that must stay verbatim (D-10).
    """
    shx = Shelxfile()
    shx.read_string(
        HEADER + ATOMS + FOOTER
        + 'EQIV $9 1-x, y, 1-z\n'
          'HTAB C1 C2_$9\n'
    )
    graph = AtomRestraintGraph(shx)
    assert any(e.id == '$9' for e in shx.eqiv), 'precondition: the EQIV parsed'
    assert [e.id for e in graph.unreferenced_eqivs()] == []


def test_atom_lookup_is_case_insensitive() -> None:
    """Names keep the file's casing, so both map directions must agree."""
    shx, graph = _graph(
        'DFIX 1.5 C1 H9a\n'
        'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
        'H9a   2     0.15  0.25  0.35  11.00000  0.05\n'
    )
    card = next(iter(shx.restraints))
    for name in graph.atoms_for_card(card):
        assert card in graph.cards_for_atom(name)
        assert card in graph.cards_for_atom(name.lower())


# ------------------------------------------------------ C3: staleness

def test_graph_is_fresh_after_building(doc: ShelxDocument) -> None:
    assert doc.graph.is_stale is False


def test_deleting_an_atom_marks_the_graph_stale(doc: ShelxDocument) -> None:
    graph = doc.graph
    doc.shelxfile.atoms.all_atoms[0].delete()
    assert graph.is_stale


def test_graph_rebuilds_itself_on_access(doc: ShelxDocument) -> None:
    before = len(doc.graph.referencing_atoms())
    doc.shelxfile.atoms.get_atom_by_name('C1_4').delete()
    after = len(doc.graph.referencing_atoms())
    assert after < before, 'the rebuilt graph must reflect the deletion'


def test_adding_an_atom_marks_the_graph_stale(doc: ShelxDocument) -> None:
    graph = doc.graph
    assert not graph.is_stale
    doc.shelxfile.add_atom(name='C99', coordinates=[0.1, 0.1, 0.1])
    assert graph.is_stale


def test_removing_a_card_marks_the_graph_stale(doc: ShelxDocument) -> None:
    graph = doc.graph
    assert not graph.is_stale
    next(iter(doc.shelxfile.restraints)).delete()
    assert graph.is_stale


# -------------------------------------------------------------- real file

def test_real_file_builds_a_populated_graph(doc: ShelxDocument) -> None:
    graph = doc.graph
    assert len(graph) > 0
    assert len(graph.referencing_atoms()) > 0


def test_every_indexed_card_is_atom_linked(doc: ShelxDocument) -> None:
    assert all(card.is_atom_linked for card in doc.graph.cards)


def test_links_are_symmetric(doc: ShelxDocument) -> None:
    """Whatever a card claims must claim it back."""
    graph = doc.graph
    for card in graph.cards:
        for name in graph.atoms_for_card(card):
            assert card in graph.cards_for_atom(name)
