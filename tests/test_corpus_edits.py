"""Destructive edit invariants over the real-structure corpus (I3-I5).

Companion to ``test_corpus.py``, which covers parsing and round-tripping
(I1/I2).  Here every sampled structure loses one atom, and the result is
held to three rules:

* **I3 edit locality.** Only lines belonging to the deleted atom, or to a
  card the report names, may change.  This is the real form of the
  guarantee the original plan stated as "byte-identical output" -- an
  invariant that cannot hold, because reading and writing normalises
  formatting.
* **I4 report completeness.** Nothing disappears silently.  Everything
  that vanished between the two dumps must appear in the
  :class:`DeletionReport`.
* **I5 no semantic escalation.** No card may be left having named atoms
  but naming none.  An emptied ``ISOR`` does not mean "no atoms"; it
  means *"all non-hydrogen atoms"*, so trimming one down to nothing would
  silently widen it to the whole structure.

Editing costs more than parsing, so the whole sweep is held to
``--corpus-time-budget`` (900 s by default).  ``--corpus-edit-sample``
bounds this sweep alone; ``--corpus-sample`` bounds every corpus test.
"""

from __future__ import annotations

import random
import time
from pathlib import Path

import pytest

from shelxfile import Shelxfile
from shelxfile.edit import ShelxDocument
from shelxfile.edit.card_meta import AtomReferencingCard, CardLifetime

pytestmark = pytest.mark.corpus

#: Same seed for every run, so a failure is reproducible.
_PICK_SEED = 20260921


def _normalise(line: str) -> str:
    """Compare lines by content, not by how they happen to be wrapped."""
    stripped = line.strip()
    if stripped.endswith('='):
        stripped = stripped[:-1]
    return ' '.join(stripped.split())


def _attributable(line: str, allowed: str) -> bool:
    """Whether *line* is part of one of the texts the report accounts for.

    Substring rather than equality because a long instruction is written
    across several physical lines; each fragment belongs to the same
    logical card.
    """
    norm = _normalise(line)
    return bool(norm) and norm in allowed


def _pick_atom(doc: ShelxDocument):
    """An atom that at least one live card refers to, or ``None``.

    Deleting an atom nothing references exercises almost nothing, so the
    sweep looks for one that will actually make the cascade work.
    """
    shx = doc.shelxfile
    referenced: list[str] = []
    for item in shx._reslist:
        if not isinstance(item, AtomReferencingCard):
            continue
        if not getattr(item, 'is_atom_linked', False):
            continue
        referenced.extend(
            token for token in item.referenced_atoms
            if token not in ('>', '<', '=') and not token.startswith('$')
        )
    rng = random.Random(_PICK_SEED)
    rng.shuffle(referenced)
    for token in referenced:
        atom = shx.atoms.get_atom_by_name(token.split('_')[0])
        if atom is not None and not atom.qpeak:
            return atom
    for atom in shx.atoms:
        if not atom.qpeak:
            return atom
    return None


def _escalated_cards(shx: Shelxfile) -> list[str]:
    """Cards left having named atoms but naming none (**I5**).

    The test is *"had atoms, has none"*, not merely *"has none"*: a card
    authored without atom names is perfectly valid -- ``HFIX 0`` and bare
    ``ISOR`` both occur in real files -- and must not be mistaken for one
    that was emptied by an edit.
    """
    offenders = []
    for item in shx._reslist:
        if not isinstance(item, AtomReferencingCard):
            continue
        if not getattr(item, 'lifetime', CardLifetime.INPUT):
            continue  # post-END output is inert
        if item.referenced_atoms:
            continue
        if not getattr(item, 'atoms_were_edited', False):
            continue  # authored empty, not emptied
        offenders.append(str(item))
    return offenders


@pytest.fixture(scope='session')
def edit_sweep(corpus_edit_files: list[Path]) -> dict:
    """Delete one atom from every sampled structure and record the fallout."""
    result: dict[str, list] = {
        'edited': [], 'skipped': [], 'errors': [],
        'unreported': [], 'unexpected_additions': [], 'escalated': [],
        'reparse_errors': [], 'atom_count_drift': [],
    }
    started = time.perf_counter()
    for path in corpus_edit_files:
        key = str(path)
        try:
            doc = ShelxDocument.from_file(path)
        except Exception as exc:  # noqa: BLE001 - measuring, not handling
            result['errors'].append(f'{key}: {type(exc).__name__}: {exc}')
            continue
        atom = _pick_atom(doc)
        if atom is None:
            result['skipped'].append(key)
            continue

        atom_line = str(atom)
        before_text = doc.text
        atoms_before = len(doc.shelxfile.atoms)
        try:
            report = doc.delete_atoms([atom])
        except Exception as exc:  # noqa: BLE001
            result['errors'].append(f'{key}: delete raised {type(exc).__name__}: {exc}')
            continue
        after_text = doc.text
        result['edited'].append(key)

        accounted = ' || '.join(
            [_normalise(atom_line)]
            + [_normalise(str(item.item)) for item in report.atoms]
            + [_normalise(item.text) for item in report.atoms if item.text]
            + [_normalise(str(item.item)) for item in report.cards]
            + [_normalise(item.text) for item in report.cards if item.text]
            + [_normalise(edit.before) for edit in report.edited]
            + [_normalise(str(edit.card)) for edit in report.edited]
        )

        before_lines = before_text.splitlines()
        after_lines = after_text.splitlines()
        remaining = list(after_lines)
        for line in before_lines:
            if line in remaining:
                remaining.remove(line)
                continue
            if not _attributable(line, accounted):
                result['unreported'].append(f'{key}: vanished line {line!r}')
                break
        leftover = list(before_lines)
        for line in after_lines:
            if line in leftover:
                leftover.remove(line)
                continue
            if not _attributable(line, accounted):
                result['unexpected_additions'].append(f'{key}: new line {line!r}')
                break

        escalated = _escalated_cards(doc.shelxfile)
        if escalated:
            result['escalated'].append(f'{key}: {escalated[:3]}')

        try:
            reread = Shelxfile()
            reread.read_string(after_text)
        except Exception as exc:  # noqa: BLE001
            result['reparse_errors'].append(f'{key}: {type(exc).__name__}: {exc}')
            continue
        if len(reread.atoms) != atoms_before - 1 and not report.atoms[1:]:
            result['atom_count_drift'].append(
                f'{key}: {atoms_before} -> {len(reread.atoms)} atoms but the '
                f'report names only {len(report.atoms)}'
            )
    result['elapsed'] = time.perf_counter() - started
    return result


def test_the_sweep_actually_edited_something(edit_sweep: dict) -> None:
    """Guard the guard: a sweep that edits nothing proves nothing."""
    assert edit_sweep['edited'], (
        'no corpus structure could be edited, so I3-I5 were never exercised. '
        f'Skipped {len(edit_sweep["skipped"])}, errored {len(edit_sweep["errors"])}.'
    )


def test_deleting_an_atom_never_raises(edit_sweep: dict) -> None:
    """A deletion that throws leaves the model in an unknown state."""
    assert not edit_sweep['errors'], (
        f'{len(edit_sweep["errors"])} structures failed to edit: '
        f'{edit_sweep["errors"][:5]}'
    )


def test_nothing_vanishes_unreported(edit_sweep: dict) -> None:
    """**I4** - everything removed is named in the DeletionReport."""
    assert not edit_sweep['unreported'], (
        f'{len(edit_sweep["unreported"])} structures lost a line that the '
        f'report does not mention: {edit_sweep["unreported"][:5]}'
    )


def test_edits_stay_local(edit_sweep: dict) -> None:
    """**I3** - no line appears that the report does not explain."""
    assert not edit_sweep['unexpected_additions'], (
        f'{len(edit_sweep["unexpected_additions"])} structures gained an '
        f'unexplained line: {edit_sweep["unexpected_additions"][:5]}'
    )


def test_no_card_is_emptied_instead_of_removed(edit_sweep: dict) -> None:
    """**I5** - a card that named atoms must never be left naming none."""
    assert not edit_sweep['escalated'], (
        f'{len(edit_sweep["escalated"])} structures kept a card that lost '
        f'every atom it named, silently widening its scope: '
        f'{edit_sweep["escalated"][:5]}'
    )


def test_edited_output_reparses(edit_sweep: dict) -> None:
    """An edit that produces unreadable output is data loss."""
    assert not edit_sweep['reparse_errors'], (
        f'{len(edit_sweep["reparse_errors"])} edited structures cannot be '
        f'read back: {edit_sweep["reparse_errors"][:5]}'
    )


def test_atom_count_matches_the_report(edit_sweep: dict) -> None:
    """Atoms may only disappear via the report (**I4** for atoms)."""
    assert not edit_sweep['atom_count_drift'], (
        f'{len(edit_sweep["atom_count_drift"])} structures lost a different '
        f'number of atoms than reported: {edit_sweep["atom_count_drift"][:5]}'
    )


def test_edit_sweep_stays_within_its_time_budget(edit_sweep: dict,
                                                 time_budget: float) -> None:
    """The opt-in suite is only useful if people are willing to run it."""
    if not time_budget:
        pytest.skip('time budget disabled (--corpus-time-budget 0)')
    assert edit_sweep['elapsed'] <= time_budget, (
        f'the edit sweep took {edit_sweep["elapsed"]:.0f}s over '
        f'{len(edit_sweep["edited"])} structures, above the budget of '
        f'{time_budget:.0f}s. Lower --corpus-edit-sample or make the '
        f'cascade cheaper.'
    )
