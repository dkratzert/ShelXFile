"""Cascading deletion (plan items D-3, D-9, B1, B2).

This is where the original goal finally happens: deleting an atom cleans
up the instructions that referenced it, instead of leaving dangling
names behind.

Two properties matter as much as the cleanup itself:

* the consequences are computed **before** anything is touched, so a
  cascade cannot be half-applied;
* cascades run one way only -- an atom deletion may remove cards, but
  removing a card never deletes atoms (**D-9**).
"""

from __future__ import annotations

from shelxfile.edit import ShelxDocument
from shelxfile.edit.reports import RemovalReason

HEADER = (
    'TITL t\n'
    'CELL 0.71073 10 10 10 90 90 90\n'
    'ZERR 4 0 0 0 0 0 0\n'
    'LATT 1\n'
    'SFAC C H\n'
    'UNIT 8 8\n'
)
FOOTER = 'HKLF 4\nEND\n'
ATOMS = ''.join(
    f'C{i:<4} 1     0.{i:02d}0  0.{i:02d}5  0.{i:02d}9  11.00000  0.05\n'
    for i in range(1, 7)
)


def _doc(cards: str, atoms: str = ATOMS) -> ShelxDocument:
    return ShelxDocument.from_string(HEADER + cards + atoms + FOOTER)


def _delete(doc: ShelxDocument, *names: str):
    return doc.delete_atoms(
        [doc.shelxfile.atoms.get_atom_by_name(n) for n in names]
    )


# -------------------------------------------------- trimming what survives

def test_flat_card_loses_only_the_deleted_atom() -> None:
    doc = _doc('SIMU C1 C2 C3\n')
    report = _delete(doc, 'C2_0')
    assert len(report.edited) == 1
    assert 'SIMU C1 C3' in doc.text


def test_pair_card_loses_the_whole_pair() -> None:
    """"the first and second named atom, the third and fourth ..."."""
    doc = _doc('SADI C1 C2 C3 C4\n')
    _delete(doc, 'C2_0')
    assert 'SADI C3 C4' in doc.text
    assert 'C1' not in doc.text.split('SADI')[1].split('\n')[0]


def test_pair_removal_works_from_either_side() -> None:
    doc = _doc('SADI C1 C2 C3 C4\n')
    _delete(doc, 'C3_0')
    assert 'SADI C1 C2' in doc.text


def test_trimmed_card_reports_what_it_looked_like_before() -> None:
    doc = _doc('SIMU C1 C2 C3\n')
    report = _delete(doc, 'C2_0')
    assert report.edited[0].before == 'SIMU C1 C2 C3'
    assert 'trimmed' in report.summary()


# ------------------------------------------------- removing what cannot stay

def test_card_below_its_minimum_is_removed() -> None:
    doc = _doc('DFIX 1.5 C1 C2\n')
    report = _delete(doc, 'C2_0')
    assert 'DFIX' not in doc.text
    assert report.cards[0].reason is RemovalReason.BELOW_MIN_ATOMS


def test_flat_card_emptied_would_become_global() -> None:
    """B2: a bare ``SIMU`` restrains every non-hydrogen atom."""
    doc = _doc('SIMU C1 C2\n')
    report = _delete(doc, 'C1_0', 'C2_0')
    assert 'SIMU' not in doc.text
    assert report.cards[0].reason is RemovalReason.WOULD_BECOME_GLOBAL


def test_a_batch_is_judged_as_a_batch() -> None:
    """Deleting both atoms at once is emptying, not two near-misses."""
    doc = _doc('SIMU C1 C2\n')
    report = _delete(doc, 'C1_0', 'C2_0')
    assert {c.reason for c in report.cards} == {RemovalReason.WOULD_BECOME_GLOBAL}


def test_chiv_survives_on_a_single_atom() -> None:
    """B3: each named atom gets its own chiral volume restraint."""
    doc = _doc('CHIV C1 C2\n')
    _delete(doc, 'C2_0')
    assert 'CHIV C1' in doc.text


def test_flat_needs_four_atoms() -> None:
    doc = _doc('FLAT C1 C2 C3 C4\n')
    report = _delete(doc, 'C4_0')
    assert 'FLAT' not in doc.text
    assert report.cards[0].reason is RemovalReason.BELOW_MIN_ATOMS


# ------------------------------------- B8: residue-scoped cards are tolerant

def test_a_residue_scoped_card_survives_losing_one_residue() -> None:
    """"the instruction is simply ignored for that residue"."""
    doc = _doc(
        'SADI_TOL C1 C2\n',
        atoms=(
            'RESI 1 TOL\n'
            'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
            'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
            'RESI 2 TOL\n'
            'C1    1     0.40  0.50  0.60  11.00000  0.05\n'
            'C2    1     0.45  0.55  0.65  11.00000  0.05\n'
            'RESI 0\n'
        ),
    )
    report = _delete(doc, 'C1_2')
    assert 'SADI_TOL C1 C2' in doc.text, 'residue 1 still satisfies the card'
    assert report.edited == []


def test_a_token_with_no_survivors_is_dropped() -> None:
    doc = _doc(
        'SADI_TOL C1 C2\n',
        atoms=(
            'RESI 1 TOL\n'
            'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
            'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
            'RESI 0\n'
        ),
    )
    _delete(doc, 'C1_1')
    assert 'SADI_TOL' not in doc.text


def test_losing_one_residue_is_not_reported_as_an_unknown_atom() -> None:
    """The validator must match SHELXL's tolerance, not exceed it.

    A class-scoped instruction is applied once per residue, and one that
    cannot be satisfied for a given residue is "simply ignored for that
    residue".  Flagging it made 72 of 388 real deletions look like they
    had left dangling references behind, when they had not.
    """
    doc = _doc(
        'SADI_TOL C1 C2\n',
        atoms=(
            'RESI 1 TOL\n'
            'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
            'C2    1     0.15  0.25  0.35  11.00000  0.05\n'
            'RESI 2 TOL\n'
            'C1    1     0.40  0.50  0.60  11.00000  0.05\n'
            'C2    1     0.45  0.55  0.65  11.00000  0.05\n'
            'RESI 0\n'
        ),
    )
    _delete(doc, 'C1_2')
    reparsed = ShelxDocument.from_string(doc.text)
    unknown = [w for w in reparsed.shelxfile.restraint_errors
               if 'Unknown atom' in w]
    assert unknown == []


def test_a_name_missing_from_every_residue_is_still_reported() -> None:
    """Relaxing the check must not hide a genuine typo."""
    doc = _doc(
        'SADI_TOL C1 CZZ\n',
        atoms=(
            'RESI 1 TOL\n'
            'C1    1     0.10  0.20  0.30  11.00000  0.05\n'
            'RESI 2 TOL\n'
            'C1    1     0.40  0.50  0.60  11.00000  0.05\n'
            'RESI 0\n'
        ),
    )
    unknown = [w for w in doc.shelxfile.restraint_errors if 'CZZ' in w]
    assert unknown, 'a name that resolves nowhere is a real problem'


# ------------------------------------------------------ ranges stay coherent

def test_deleting_inside_a_range_leaves_it_alone() -> None:
    doc = _doc('SIMU C1 > C4\n')
    report = _delete(doc, 'C2_0')
    assert report.edited == [], 'the endpoints still resolve'
    assert 'SIMU C1 > C4' in doc.text


def test_losing_a_range_endpoint_drops_the_marker() -> None:
    doc = _doc('SIMU C1 > C4\n')
    _delete(doc, 'C4_0')
    rendered = [line for line in doc.text.splitlines() if line.startswith('SIMU')]
    assert rendered == [] or '>' not in rendered[0], (
        'a range marker must not be left without an endpoint'
    )


# ----------------------------------------------------- D-9: direction

def test_removing_a_card_never_deletes_atoms() -> None:
    doc = _doc('SADI C1 C2 C3 C4\n')
    before = len(doc.shelxfile.atoms)
    doc.remove_restraint(next(iter(doc.shelxfile.restraints)))
    assert len(doc.shelxfile.atoms) == before


def test_a_card_removed_by_a_cascade_does_not_delete_further_atoms() -> None:
    """The indirect path must be as safe as the direct one."""
    doc = _doc('DFIX 1.5 C1 C2\nSIMU C3 C4\n')
    report = _delete(doc, 'C2_0')
    removed = {a.item.name for a in report.atoms}
    assert removed == {'C2'}, 'removing the DFIX must not take C1 with it'
    assert doc.shelxfile.atoms.get_atom_by_name('C1_0') is not None


# ------------------------------------------------------- plan before apply

def test_planning_changes_nothing() -> None:
    doc = _doc('DFIX 1.5 C1 C2\n')
    before = doc.text
    plan = doc.plan_deletion([doc.shelxfile.atoms.get_atom_by_name('C2_0')])
    assert not plan.is_empty
    assert doc.text == before, 'planning must not touch the model'


def test_a_plan_describes_the_whole_closure() -> None:
    doc = _doc('DFIX 1.5 C1 C2\n')
    plan = doc.plan_deletion([doc.shelxfile.atoms.get_atom_by_name('C2_0')])
    assert len(plan.atoms) == 1
    assert len(plan.cards) == 1


def test_planning_a_symmetry_mate_plans_nothing() -> None:
    doc = ShelxDocument.from_file('tests/resources/p21c.res')
    packed = [a for a in doc.shelxfile.pack() if a.symmgen][:2]
    assert packed
    assert doc.plan_deletion(packed).is_empty


# ------------------------------------------------------------- termination

def test_cascade_reaches_a_fixpoint() -> None:
    doc = _doc(
        'SIMU C1 C2 C3\n'
        'SADI C3 C4 C5 C6\n'
        'DFIX 1.5 C1 C6\n'
    )
    report = _delete(doc, 'C1_0')
    assert len(report) > 0
    assert doc.text  # completed rather than looping


def test_no_atom_is_removed_twice() -> None:
    doc = _doc('SIMU C1 C2 C3\n')
    report = _delete(doc, 'C1_0', 'C2_0')
    names = [a.item.name for a in report.atoms]
    assert len(names) == len(set(names))


def test_deleting_everything_leaves_a_valid_file() -> None:
    doc = _doc('SIMU C1 C2 C3\nSADI C1 C2 C3 C4\n')
    _delete(doc, *[f'C{i}_0' for i in range(1, 7)])
    assert doc.shelxfile.atoms.all_atoms == []
    assert 'HKLF' in doc.text


# ---------------------------------------------------------------- reporting

def test_every_removal_carries_a_reason() -> None:
    doc = _doc('DFIX 1.5 C1 C2\n')
    report = _delete(doc, 'C2_0')
    assert all(item.reason is not None for item in report.atoms + report.cards)


def test_report_length_counts_everything() -> None:
    doc = _doc('SIMU C1 C2 C3\n')
    report = _delete(doc, 'C2_0')
    assert len(report) == len(report.atoms) + len(report.cards) + len(report.edited)


def test_edited_cards_round_trip() -> None:
    doc = _doc('SIMU C1 C2 C3\n')
    _delete(doc, 'C2_0')
    reparsed = ShelxDocument.from_string(doc.text)
    simu = next(r for r in reparsed.shelxfile.restraints
                if type(r).__name__ == 'SIMU')
    assert simu.atoms == ['C1', 'C3']
