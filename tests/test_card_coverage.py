"""Registry-driven coverage of every atom-referencing SHELXL card.

The plan's **coverage gate**: rather than testing whichever cards happen
to be common, these tests enumerate the card *registries* and demand a
sample, a round-trip test and a deletion test for each.  Adding a card
class to ``shelxfile`` fails this file until it is covered, so coverage
cannot drift toward whatever the reference corpus happens to contain.

The samples and their expected outcomes live in ``tests/card_catalog.py``,
each justified by a quoted sentence from the SHELXL manual.
"""

from __future__ import annotations

import inspect

import pytest

from shelxfile import Shelxfile
from shelxfile.edit import ShelxDocument
from shelxfile.edit.card_meta import AtomListSemantics, AtomReferencingCard
from shelxfile.shelx import cards as cards_module
from tests.card_catalog import CARD_SAMPLES, COVERED_KEYWORDS, CardSample, Expect

#: Base classes, not instructions in their own right.
ABSTRACT = frozenset({'AtomReferencingCard', 'Restraint'})

SAMPLE_IDS = sorted(CARD_SAMPLES)


def _registry_keywords() -> set[str]:
    """Every concrete card class that declares it references atoms."""
    found = set()
    for name, obj in vars(cards_module).items():
        if (inspect.isclass(obj)
                and issubclass(obj, AtomReferencingCard)
                and name not in ABSTRACT):
            found.add(name)
    return found


def _card_in(shx: Shelxfile, keyword: str):
    for item in shx._reslist:
        if type(item).__name__ == keyword:
            return item
    return None


# ------------------------------------------------------------ the gate

def test_every_registry_card_has_a_sample() -> None:
    """A new atom-referencing card class must arrive with a fixture."""
    missing = sorted(_registry_keywords() - COVERED_KEYWORDS)
    assert not missing, (
        'These card classes reference atoms but have no sample in '
        f'tests/card_catalog.py, so nothing tests what deleting one of '
        f'their atoms does: {missing}'
    )


def test_every_restraint_registry_card_has_a_sample() -> None:
    """``add_restraint`` offers these keywords, so they must be covered."""
    missing = sorted(set(Shelxfile.RESTRAINT_CARD_CLASSES) - COVERED_KEYWORDS)
    assert not missing, (
        f'RESTRAINT_CARD_CLASSES entries without a sample: {missing}'
    )


def test_every_atom_card_registry_entry_has_a_sample() -> None:
    """The symmetry-aware editing API's registry, likewise."""
    missing = sorted(set(ShelxDocument.ATOM_CARD_CLASSES) - COVERED_KEYWORDS)
    assert not missing, (
        f'ATOM_CARD_CLASSES entries without a sample: {missing}'
    )


def test_catalog_names_only_real_cards() -> None:
    """Guard the gate itself: a typo must not count as coverage."""
    unknown = sorted(COVERED_KEYWORDS - _registry_keywords())
    assert not unknown, (
        f'card_catalog.py names cards that do not exist: {unknown}'
    )


# ------------------------------------------------------- per-sample use

@pytest.mark.parametrize('key', SAMPLE_IDS)
def test_sample_parses(key) -> None:
    """Each hand-written sample must be valid SHELXL we can read."""
    sample: CardSample = CARD_SAMPLES[key]
    shx = Shelxfile()
    shx.read_string(sample.source())
    assert shx.cell, f'{sample.instruction!r} broke parsing'
    assert _card_in(shx, sample.keyword) is not None, (
        f'{sample.instruction!r} did not produce a {sample.keyword} card'
    )


@pytest.mark.parametrize('key', SAMPLE_IDS)
def test_sample_round_trips(key) -> None:
    """Writing and re-reading a sample must not change it (I2)."""
    sample: CardSample = CARD_SAMPLES[key]
    first = Shelxfile()
    first.read_string(sample.source())
    once = first.dumps()

    second = Shelxfile()
    second.read_string(once)
    assert second.dumps() == once


@pytest.mark.parametrize('key', SAMPLE_IDS)
def test_sample_is_echoed_verbatim_when_untouched(key) -> None:
    """An unedited card is echoed exactly, never reflowed (A1)."""
    sample: CardSample = CARD_SAMPLES[key]
    shx = Shelxfile()
    shx.read_string(sample.source())
    assert sample.instruction in shx.dumps()


@pytest.mark.parametrize('key', SAMPLE_IDS)
def test_sample_deletion_matches_the_manual(key) -> None:
    """Deleting an atom must do what the manual implies for this card."""
    sample: CardSample = CARD_SAMPLES[key]
    doc = ShelxDocument.from_string(sample.source())
    atom = doc.shelxfile.atoms.get_atom_by_name(sample.delete)
    assert atom is not None, f'fixture has no atom {sample.delete}'

    doc.delete_atoms([atom])
    lines = [line for line in doc.text.splitlines()
             if line.split() and line.split()[0].split('_')[0] == sample.keyword]

    if sample.expect is Expect.REMOVED:
        assert not lines, (
            f'{sample.instruction!r} should be removed when {sample.delete} '
            f'goes. {sample.why} Got: {lines}'
        )
        return

    assert lines, (
        f'{sample.instruction!r} should survive losing {sample.delete}. '
        f'{sample.why}'
    )
    if sample.expect is Expect.UNTOUCHED:
        assert lines == [sample.instruction], (
            f'{sample.instruction!r} should be left alone. {sample.why}'
        )
    else:
        assert lines != [sample.instruction], (
            f'{sample.instruction!r} should have been trimmed. {sample.why}'
        )


@pytest.mark.parametrize('key', SAMPLE_IDS)
def test_deletion_leaves_no_dangling_reference(key) -> None:
    """No surviving card may name an atom that is gone (the whole point)."""
    sample: CardSample = CARD_SAMPLES[key]
    doc = ShelxDocument.from_string(sample.source())
    atom = doc.shelxfile.atoms.get_atom_by_name(sample.delete)
    doc.delete_atoms([atom])

    reread = Shelxfile()
    reread.read_string(doc.text)
    card = _card_in(reread, sample.keyword)
    if card is None:
        return
    if card.atom_semantics is not AtomListSemantics.EXPLICIT:
        return
    named = [token.split('_')[0].upper() for token in card.referenced_atoms]
    assert sample.delete.upper() not in named, (
        f'{sample.instruction!r} still names the deleted {sample.delete} '
        f'after the edit: {card}'
    )


@pytest.mark.parametrize('key', SAMPLE_IDS)
def test_edit_stays_local(key) -> None:
    """Only the sampled card and the deleted atom may change (I3)."""
    sample: CardSample = CARD_SAMPLES[key]
    doc = ShelxDocument.from_string(sample.source())
    before = doc.text.splitlines()
    atom = doc.shelxfile.atoms.get_atom_by_name(sample.delete)
    doc.delete_atoms([atom])
    after = doc.text.splitlines()

    vanished = [line for line in before if line not in after]
    allowed = {sample.instruction}
    for line in vanished:
        first = line.split()[0] if line.split() else ''
        assert line in allowed or first == sample.delete, (
            f'deleting {sample.delete} also removed an unrelated line: '
            f'{line!r}'
        )


@pytest.mark.parametrize('key', SAMPLE_IDS)
def test_rename_reaches_every_card(key) -> None:
    """Renaming an atom must reach the card text, whatever its shape.

    Cards keep their operands in different places -- a flat ``atoms``
    list, or named fields such as ``HTAB``'s donor/acceptor -- and a card
    that fails to regenerate would be left naming an atom that no longer
    exists.
    """
    sample: CardSample = CARD_SAMPLES[key]
    doc = ShelxDocument.from_string(sample.source())
    atom = doc.shelxfile.atoms.get_atom_by_name(sample.delete)
    doc.rename_atom(atom, 'Z9')

    lines = [line for line in doc.text.splitlines()
             if line.split() and line.split()[0].split('_')[0] == sample.keyword]
    assert lines, f'{sample.instruction!r} vanished during a rename'

    if sample.expect is Expect.UNTOUCHED:
        assert lines == [sample.instruction], (
            f'{sample.instruction!r} names no atoms, so a rename must not '
            f'touch it. {sample.why}'
        )
        return

    text = ' '.join(lines)
    assert 'Z9' in text, (
        f'{sample.instruction!r} did not follow {sample.delete} -> Z9; '
        f'the card still reads {text!r}'
    )
    named = [token.split('_')[0].upper()
             for line in lines for token in line.split()[1:]]
    assert sample.delete.upper() not in named, (
        f'{sample.instruction!r} still names the old {sample.delete}: {text!r}'
    )


@pytest.mark.parametrize('key', SAMPLE_IDS)
def test_result_is_still_readable(key) -> None:
    """The edited file must reparse without losing the atom list (I1)."""
    sample: CardSample = CARD_SAMPLES[key]
    doc = ShelxDocument.from_string(sample.source())
    atom = doc.shelxfile.atoms.get_atom_by_name(sample.delete)
    doc.delete_atoms([atom])

    reread = Shelxfile()
    reread.read_string(doc.text)
    assert reread.cell
    assert len(reread.atoms) == len(doc.shelxfile.atoms)
