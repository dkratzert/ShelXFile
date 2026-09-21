"""Opt-in sweep over a corpus of real SHELXL structures.

Run with::

    pytest --corpus D:/frames/Workordner -m corpus

Skipped entirely when no corpus is configured, so CI and contributors
without the (private) data are unaffected.  See ``tests/conftest.py``.

Covers invariants **I1** (parse stability) and **I2** (``dumps()``
idempotence) from the plan.  Both are enforced as *ratchets*: the
aggregate budgets in ``tests/resources/corpus_expectations.json`` may only
ever be lowered.  The corpus itself is never committed, and neither are
the file names -- only the aggregate counts.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from shelxfile import Shelxfile

pytestmark = pytest.mark.corpus

EXPECTATIONS_FILE = Path(__file__).resolve().parent / 'resources' / 'corpus_expectations.json'


@pytest.fixture(scope='session')
def expectations() -> dict:
    return json.loads(EXPECTATIONS_FILE.read_text(encoding='utf-8'))


@pytest.fixture(scope='session')
def is_sampled(request: pytest.FixtureRequest) -> bool:
    """True when ``--corpus-sample`` limited the sweep to a subset."""
    return bool(request.config.getoption('--corpus-sample'))


@pytest.fixture(scope='session')
def sweep(corpus_files: list[Path]) -> dict:
    """Parse and round-trip every corpus structure once, for all tests here."""
    result: dict[str, list[str]] = {
        'parsed': [], 'errors': [], 'non_idempotent': [], 'reparse_errors': [],
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        scratch = Path(tmpdir)
        for path in corpus_files:
            key = path.name
            try:
                first = Shelxfile()
                first.read_file(str(path))
                dumped = first.dumps()
            except Exception as exc:  # noqa: BLE001 - measuring, not handling
                result['errors'].append(f'{key}: {type(exc).__name__}: {exc}')
                continue
            result['parsed'].append(key)

            roundtrip = scratch / path.name
            roundtrip.write_text(dumped, encoding='utf-8', errors='replace')
            try:
                second = Shelxfile()
                second.read_file(str(roundtrip))
            except Exception as exc:  # noqa: BLE001
                result['reparse_errors'].append(f'{key}: {type(exc).__name__}: {exc}')
                continue
            if second.dumps() != dumped:
                result['non_idempotent'].append(key)
    return result


def test_corpus_is_not_empty(corpus_files: list[Path]) -> None:
    assert corpus_files, 'corpus directory contains no .res/.ins files'


def test_parse_errors_within_budget(sweep: dict, corpus_files: list[Path],
                                    expectations: dict, is_sampled: bool) -> None:
    """I1 - parsing real structures must not regress."""
    if is_sampled:
        pytest.skip('ratchet budgets are only meaningful over the full sweep; '
                    'drop --corpus-sample to enforce them')
    budget = expectations['max_parse_error_ratio'] * len(corpus_files)
    errors = sweep['errors']
    assert len(errors) <= budget, (
        f'{len(errors)} parse errors exceeds the budget of {budget:.1f} '
        f'({expectations["max_parse_error_ratio"]:.2%} of {len(corpus_files)} files). '
        f'First few: {errors[:5]}'
    )


def test_dumps_is_idempotent_within_budget(sweep: dict, expectations: dict,
                                           is_sampled: bool) -> None:
    """I2 - ``dumps(parse(dumps(parse(f)))) == dumps(parse(f))``.

    Known remaining causes are catalogued in
    ``corpus_expectations.json``; the budget may only be lowered.

    Enforced only on a full sweep: on a small sample a single failure
    easily exceeds a sub-percent ratio, which would be noise rather than
    a regression.
    """
    if is_sampled:
        pytest.skip('ratchet budgets are only meaningful over the full sweep; '
                    'drop --corpus-sample to enforce them')
    parsed = len(sweep['parsed'])
    if not parsed:
        pytest.skip('no structures parsed')
    failures = sweep['non_idempotent']
    ratio = len(failures) / parsed
    assert ratio <= expectations['max_non_idempotent_ratio'], (
        f'{len(failures)}/{parsed} structures ({ratio:.2%}) are not '
        f'dumps()-idempotent, above the ratchet of '
        f'{expectations["max_non_idempotent_ratio"]:.2%}. '
        f'First few: {failures[:5]}'
    )


def test_roundtrip_output_always_reparses(sweep: dict) -> None:
    """Whatever ``dumps()`` writes must itself be parsable.

    Weaker than idempotence but absolute: a file we wrote and cannot read
    back is data loss, not a formatting quirk.
    """
    assert not sweep['reparse_errors'], (
        'ShelXFile produced output it cannot parse back: '
        f'{sweep["reparse_errors"][:5]}'
    )
