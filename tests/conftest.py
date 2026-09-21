"""Shared pytest configuration for the ShelXFile test suite.

Adds the opt-in *corpus* facility: a large collection of real
``.res``/``.ins`` structures living outside the repository, used to check
the parser and the editing layer against production data.

The corpus is private and is never committed.  Tests that need it are
marked ``@pytest.mark.corpus`` and are skipped automatically unless a
path is supplied::

    pytest --corpus D:/frames/Workordner -m corpus
    SHELXFILE_CORPUS=D:/frames/Workordner pytest -m corpus

Use ``--corpus-sample N`` to run against a deterministic random subset
while iterating; the full sweep is the default.
"""

from __future__ import annotations

import os
import random
from pathlib import Path

import pytest

CORPUS_ENV_VAR = 'SHELXFILE_CORPUS'

#: Deterministic seed so ``--corpus-sample`` selects the same files every run.
_SAMPLE_SEED = 20260921


def pytest_addoption(parser: pytest.Parser) -> None:
    group = parser.getgroup('shelxfile')
    group.addoption(
        '--corpus',
        action='store',
        default=None,
        metavar='PATH',
        help='Directory tree of real .res/.ins structures to test against. '
             f'Falls back to the {CORPUS_ENV_VAR} environment variable.',
    )
    group.addoption(
        '--corpus-sample',
        action='store',
        type=int,
        default=0,
        metavar='N',
        help='Use a deterministic random sample of N corpus files instead of '
             'the full sweep. 0 (the default) means all files.',
    )


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        'markers',
        'corpus: test requires the external structure corpus '
        '(--corpus PATH or $' + CORPUS_ENV_VAR + '); skipped otherwise.',
    )


def _corpus_root(config: pytest.Config) -> Path | None:
    """Resolve the corpus directory from the CLI option or the environment."""
    raw = config.getoption('--corpus') or os.environ.get(CORPUS_ENV_VAR)
    if not raw:
        return None
    path = Path(raw).expanduser()
    return path if path.is_dir() else None


def pytest_collection_modifyitems(config: pytest.Config,
                                  items: list[pytest.Item]) -> None:
    if _corpus_root(config) is not None:
        return
    skip = pytest.mark.skip(
        reason=f'no structure corpus configured (--corpus PATH or ${CORPUS_ENV_VAR})'
    )
    for item in items:
        if 'corpus' in item.keywords:
            item.add_marker(skip)


@pytest.fixture(scope='session')
def corpus_root(request: pytest.FixtureRequest) -> Path:
    """Root directory of the external structure corpus."""
    root = _corpus_root(request.config)
    if root is None:
        pytest.skip(f'no structure corpus configured (--corpus PATH or ${CORPUS_ENV_VAR})')
    return root


@pytest.fixture(scope='session')
def corpus_files(request: pytest.FixtureRequest, corpus_root: Path) -> list[Path]:
    """All corpus structures, or a deterministic sample of ``--corpus-sample``.

    Sorted so that the full sweep is reproducible, and seeded so that a
    sample is the same set on every run.
    """
    files = sorted(
        [p for p in corpus_root.rglob('*.res') if p.is_file()]
        + [p for p in corpus_root.rglob('*.ins') if p.is_file()]
    )
    limit = int(request.config.getoption('--corpus-sample') or 0)
    if 0 < limit < len(files):
        rng = random.Random(_SAMPLE_SEED)
        files = sorted(rng.sample(files, limit))
    return files
