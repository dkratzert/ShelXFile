"""``ShelxDocument.refine()``: SHELXL through the edit layer's "new instance" pattern.

Mirrors :meth:`ShelxDocument.apply`: a successful refinement replaces the
model bound to the document with a freshly parsed one (read from the
``.res`` file SHELXL produced) instead of mutating the old one in place. A
failed refinement restores the document to its pre-refine state -- including
the ``ACTA`` card, which has to be removed before SHELXL runs -- and raises
:class:`RuntimeError` instead of the underlying code's ``sys.exit()``.

Tests that actually invoke SHELXL are skipped when no ``shelxl``/``xl``
executable is on ``PATH``, matching the guard already used by
``tests/test_refine.py``.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from shelxfile.edit import ShelxDocument
from shelxfile.refine.refine import find_shelxl_exe

MODEL_FINISHED = Path('tests/resources/model_finished')


def _clean_refine_files(basename: str = 'p21c') -> None:
    for suffix in ('.res', '.ins', '.lst', '.fcf', '.fcf6', '.cif', '.hkl', '.shx-bak'):
        Path(f'{basename}{suffix}').unlink(missing_ok=True)
    shxsaves = Path('shxsaves')
    if shxsaves.exists():
        shutil.rmtree(shxsaves, ignore_errors=True)


@pytest.fixture
def doc() -> ShelxDocument:
    if not find_shelxl_exe():
        pytest.skip('SHELXL not found')
    shutil.copy(MODEL_FINISHED / 'p21c.res', '.')
    shutil.copy(MODEL_FINISHED / 'p21c.hkl', '.')
    document = ShelxDocument.from_file('p21c.res', debug=True)
    yield document
    _clean_refine_files('p21c')


def test_refine_replaces_the_model_with_the_refined_result(doc: ShelxDocument) -> None:
    old_shx = doc.shelxfile
    assert doc.refine(cycles=2) is True
    assert doc.shelxfile is not old_shx
    assert 'L.S. 2' in doc.text


def test_refine_notifies_observers(doc: ShelxDocument) -> None:
    seen = []
    doc.subscribe(lambda d: seen.append(d))
    doc.refine(cycles=1)
    assert seen and seen[-1] is doc


def test_refine_without_resfile_raises() -> None:
    document = ShelxDocument.from_string(Path('tests/resources/p21c.res').read_text())
    with pytest.raises(RuntimeError, match='no file path'):
        document.refine()


def test_refine_missing_hkl_restores_previous_state(doc: ShelxDocument) -> None:
    Path('p21c.hkl').unlink()
    before_text = doc.text
    with pytest.raises(RuntimeError):
        doc.refine(cycles=1)
    assert doc.text == before_text
