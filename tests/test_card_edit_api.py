"""Adding and removing atom-referencing cards (plan item D-7c/d/e).

The surface a viewer needs: the user clicks two atoms -- one possibly a
symmetry image -- and asks for a bond.  ShelXFile works out how to name
them, creates the ``EQIV`` if that operation has not been named yet, and
writes a valid instruction.

Two arity limits from the manual are enforced rather than assumed:
``FREE``/``BIND`` allow only one symmetry operand, and ``HTAB`` allows one
only on the acceptor.
"""

from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.edit import ShelxDocument

RES = 'tests/resources/p21c.res'


@pytest.fixture
def doc() -> ShelxDocument:
    return ShelxDocument.from_file(RES)


def _mate(doc: ShelxDocument):
    return next(a for a in doc.shelxfile.pack() if a.symmgen)


# --------------------------------------------------------- plain atoms

def test_bind_between_two_plain_atoms(doc: ShelxDocument) -> None:
    report = doc.add_bind(
        doc.shelxfile.atoms.get_atom_by_name('C1_4'),
        doc.shelxfile.atoms.get_atom_by_name('C2_4'),
    )
    assert str(report.added[0]) == 'BIND C1_4 C2_4'
    assert 'BIND C1_4 C2_4' in doc.text


def test_free_between_two_plain_atoms(doc: ShelxDocument) -> None:
    doc.add_free(
        doc.shelxfile.atoms.get_atom_by_name('C1_4'),
        doc.shelxfile.atoms.get_atom_by_name('C2_4'),
    )
    assert 'FREE C1_4 C2_4' in doc.text


def test_atoms_may_be_named_as_strings(doc: ShelxDocument) -> None:
    """A viewer may already hold the name rather than the object."""
    doc.add_bind('C1_4', 'C2_4')
    assert 'BIND C1_4 C2_4' in doc.text


def test_the_card_reparses(doc: ShelxDocument) -> None:
    doc.add_bind('C1_4', 'C2_4')
    reread = Shelxfile()
    reread.read_string(doc.text)
    assert any(str(c) == 'BIND C1_4 C2_4' for c in reread.bind)


# ----------------------------------------------------- symmetry images

def test_a_symmetry_image_gets_an_eqiv_reference(doc: ShelxDocument) -> None:
    report = doc.add_bind(doc.shelxfile.atoms.get_atom_by_name('C1_4'), _mate(doc))
    rendered = str(report.added[0])
    assert '_$' in rendered
    assert [ln for ln in doc.text.splitlines() if ln.startswith('EQIV')]


def test_the_eqiv_is_defined_before_the_card(doc: ShelxDocument) -> None:
    """"Such a symmetry operation must be defined before it is used"."""
    doc.add_bind(doc.shelxfile.atoms.get_atom_by_name('C1_4'), _mate(doc))
    lines = doc.text.splitlines()
    eqiv_at = next(i for i, ln in enumerate(lines) if ln.startswith('EQIV'))
    bind_at = next(i for i, ln in enumerate(lines) if ln.startswith('BIND'))
    assert eqiv_at < bind_at


def test_one_eqiv_serves_several_cards(doc: ShelxDocument) -> None:
    mates = [a for a in doc.shelxfile.pack() if a.symmgen]
    same_op = [m for m in mates
               if m.symm_mate.symm_number == mates[0].symm_mate.symm_number][:3]
    assert len(same_op) > 1

    for mate in same_op:
        doc.add_bind(doc.shelxfile.atoms.get_atom_by_name('C1_4'), mate)

    eqivs = [ln for ln in doc.text.splitlines() if ln.startswith('EQIV')]
    assert len(eqivs) == 1


def test_the_symmetry_card_reparses_without_warnings(doc: ShelxDocument) -> None:
    doc.add_bind(doc.shelxfile.atoms.get_atom_by_name('C1_4'), _mate(doc))
    reread = Shelxfile()
    reread.read_string(doc.text)
    assert not [w for w in reread.restraint_errors if 'EQIV' in w]


# ------------------------------------------------------- D-7d: arity

def test_bind_refuses_two_symmetry_operands(doc: ShelxDocument) -> None:
    """"Only one of the two atoms may be an equivalent atom"."""
    mate = _mate(doc)
    with pytest.raises(ValueError, match='only one of its two atoms'):
        doc.add_bind(mate, mate)


def test_free_refuses_two_symmetry_operands(doc: ShelxDocument) -> None:
    mate = _mate(doc)
    with pytest.raises(ValueError, match='only one of its two atoms'):
        doc.add_free(mate, mate)


def test_htab_refuses_a_symmetry_donor(doc: ShelxDocument) -> None:
    """"Only the acceptor atom may specify a symmetry operation"."""
    with pytest.raises(ValueError, match='acceptor'):
        doc.add_htab(_mate(doc), doc.shelxfile.atoms.get_atom_by_name('C1_4'))


def test_htab_accepts_a_symmetry_acceptor(doc: ShelxDocument) -> None:
    report = doc.add_htab(
        doc.shelxfile.atoms.get_atom_by_name('C1_4'), _mate(doc)
    )
    assert '_$' in str(report.added[0])


def test_a_refused_card_leaves_no_trace(doc: ShelxDocument) -> None:
    mate = _mate(doc)
    before = doc.text
    with pytest.raises(ValueError):
        doc.add_bind(mate, mate)
    assert 'BIND' not in doc.text
    assert doc.text.count('EQIV') <= before.count('EQIV') + 1, (
        'at most the EQIV itself may remain, never a half-written card'
    )


# ------------------------------------------------------------- removal

def test_removing_a_card_keeps_the_atoms(doc: ShelxDocument) -> None:
    """D-9b holds here too."""
    report = doc.add_bind('C1_4', 'C2_4')
    before = len(doc.shelxfile.atoms)

    doc.remove_card(report.added[0])

    assert len(doc.shelxfile.atoms) == before
    assert 'BIND' not in doc.text


def test_removing_the_last_user_collects_its_eqiv(doc: ShelxDocument) -> None:
    report = doc.add_bind(doc.shelxfile.atoms.get_atom_by_name('C1_4'), _mate(doc))
    assert 'EQIV' in doc.text

    removal = doc.remove_card(report.added[0])

    assert 'EQIV' not in doc.text
    assert any(c.reason.name == 'EQIV_UNREFERENCED' for c in removal.cards)


def test_an_eqiv_with_other_users_survives(doc: ShelxDocument) -> None:
    mate = _mate(doc)
    first = doc.add_bind(doc.shelxfile.atoms.get_atom_by_name('C1_4'), mate)
    doc.add_bind(doc.shelxfile.atoms.get_atom_by_name('C2_4'), mate)

    doc.remove_card(first.added[0])

    assert 'EQIV' in doc.text


def test_removing_an_absent_card_is_harmless(doc: ShelxDocument) -> None:
    report = doc.add_bind('C1_4', 'C2_4')
    card = report.added[0]
    doc.remove_card(card)
    assert doc.remove_card(card).is_empty


# ------------------------------------------------------------ plumbing

def test_observers_are_notified(doc: ShelxDocument) -> None:
    seen: list[object] = []
    doc.subscribe(seen.append)
    doc.add_bind('C1_4', 'C2_4')
    assert seen == [doc]


def test_an_unresolvable_operand_is_rejected(doc: ShelxDocument) -> None:
    class Stranger:
        symm_mate = None
        fullname_short = None

    with pytest.raises(ValueError, match='atom reference'):
        doc.add_bind('C1_4', Stranger())


def test_restraints_still_go_through_their_own_path(doc: ShelxDocument) -> None:
    """``add_restraint`` stays the route for restraints."""
    report = doc.add_restraint('SADI C1_4 C2_4')
    assert type(report.added[0]).__name__ == 'SADI'
