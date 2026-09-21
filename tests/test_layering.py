"""Invariant I6 — architectural layering guard (plan decision D-8).

The SHELXL model must stay independent of Qt, and the Qt widget must stay
a pure view:

* ``shelxfile/edit/``  — Qt-free editing layer.  May not import Qt.
* ``shelxfile/gui/``   — display and user interaction only.  May not touch
  the private ``_reslist`` and may not make SHELXL parsing decisions.

Dependency direction is ``gui -> edit -> shelx/atoms``, never the reverse.

These are cheap textual checks on purpose: they run everywhere, cost
nothing, and stop the boundary eroding again.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parent.parent / 'src' / 'shelxfile'
EDIT_DIR = SRC / 'edit'
GUI_DIR = SRC / 'gui'

QT_ROOTS = frozenset({'qtpy', 'PyQt5', 'PyQt6', 'PySide2', 'PySide6', 'Qt'})

#: Modules under ``gui/`` that are allowed to keep SHELXL knowledge.
#: Syntax highlighting legitimately needs the instruction vocabulary.
GUI_PARSING_EXEMPT = frozenset({'syntax_highlighter.py'})


def _python_files(directory: Path) -> list[Path]:
    if not directory.is_dir():
        return []
    return sorted(p for p in directory.rglob('*.py') if '__pycache__' not in p.parts)


def _imported_roots(path: Path) -> set[str]:
    """Top-level package names imported by *path*, via AST (not text)."""
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    roots: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                roots.add(alias.name.split('.')[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level == 0 and node.module:
                roots.add(node.module.split('.')[0])
    return roots


def test_edit_layer_exists_or_is_not_yet_created() -> None:
    """Documents the expected layout; informative rather than enforcing."""
    if not EDIT_DIR.is_dir():
        pytest.skip('shelxfile/edit/ not created yet (plan task: edit-layer)')
    assert (EDIT_DIR / '__init__.py').is_file(), 'edit layer must be a package'


def test_edit_layer_is_qt_free() -> None:
    """D-8c: nothing under ``shelxfile/edit/`` may import Qt."""
    files = _python_files(EDIT_DIR)
    if not files:
        pytest.skip('shelxfile/edit/ not created yet (plan task: edit-layer)')
    offenders = {
        str(path.relative_to(SRC)): sorted(_imported_roots(path) & QT_ROOTS)
        for path in files
        if _imported_roots(path) & QT_ROOTS
    }
    assert not offenders, (
        'The edit layer must stay Qt-free (D-8c). Qt imports found in: '
        f'{offenders}'
    )


def test_edit_layer_does_not_import_gui() -> None:
    """D-8c: the dependency direction is gui -> edit, never the reverse."""
    files = _python_files(EDIT_DIR)
    if not files:
        pytest.skip('shelxfile/edit/ not created yet (plan task: edit-layer)')
    offenders = [
        str(path.relative_to(SRC))
        for path in files
        if 'shelxfile.gui' in path.read_text(encoding='utf-8')
    ]
    assert not offenders, (
        f'The edit layer must not depend on the GUI layer (D-8c): {offenders}'
    )


@pytest.mark.xfail(
    reason='D-8b not implemented yet: editor_widget.py still reaches into '
           '_reslist. Remove this xfail as part of the widget-slimming task.',
    strict=False,
)
def test_gui_does_not_touch_reslist() -> None:
    """D-8b: the widget is a pure view and must not index ``_reslist``.

    ``AGENTS.md``: "Never access ``_reslist`` directly in user/test code."
    """
    offenders: dict[str, list[int]] = {}
    for path in _python_files(GUI_DIR):
        lines = path.read_text(encoding='utf-8').splitlines()
        hits = [
            num for num, line in enumerate(lines, start=1)
            if '_reslist' in line and not line.lstrip().startswith(('#', '*'))
            and '"""' not in line
        ]
        if hits:
            offenders[str(path.relative_to(SRC))] = hits
    assert not offenders, (
        'The GUI layer must route all model access through ShelxDocument '
        f'(D-8b). Direct _reslist access found at: {offenders}'
    )


@pytest.mark.xfail(
    reason='D-8b not implemented yet: editor_widget.py still parses and '
           'dispatches SHELXL instructions. Remove with widget-slimming.',
    strict=False,
)
def test_gui_does_not_make_parsing_decisions() -> None:
    """D-8b: card semantics belong to the edit layer, not the widget."""
    markers = ('RESTRAINT_CARD_CLASSES', '.split()', 'add_restraint(')
    offenders: dict[str, list[str]] = {}
    for path in _python_files(GUI_DIR):
        if path.name in GUI_PARSING_EXEMPT:
            continue
        text = path.read_text(encoding='utf-8')
        hits = [m for m in markers if m in text]
        if hits:
            offenders[str(path.relative_to(SRC))] = hits
    assert not offenders, (
        'The GUI layer must not parse or dispatch SHELXL instructions '
        f'(D-8b): {offenders}'
    )
