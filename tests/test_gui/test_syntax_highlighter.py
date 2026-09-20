from __future__ import annotations

import pytest

from shelxfile.gui.syntax_highlighter import SHELX_KEYWORDS, ShelxSyntaxHighlighter

pytestmark = pytest.mark.usefixtures('qtbot')


@pytest.fixture
def editor(qtbot):
    from qtpy.QtWidgets import QPlainTextEdit
    widget = QPlainTextEdit()
    qtbot.addWidget(widget)
    widget.highlighter = ShelxSyntaxHighlighter(widget.document())
    return widget


def _format_at(editor, block_no: int, column: int):
    from qtpy.QtWidgets import QApplication
    QApplication.processEvents()
    block = editor.document().findBlockByNumber(block_no)
    for fmt_range in block.layout().formats():
        if fmt_range.start <= column < fmt_range.start + fmt_range.length:
            return fmt_range.format
    return block.charFormat()


def test_shelx_keywords_deduplicated_and_stripped():
    assert 'TITL' in SHELX_KEYWORDS
    assert 'REM' in SHELX_KEYWORDS
    # No entry should carry the internal fixed-width padding.
    assert all(kw == kw.strip() for kw in SHELX_KEYWORDS)


def test_keyword_line_is_highlighted(editor):
    editor.setPlainText('CELL 0.71073 10 10 10 90 90 90')
    fmt = _format_at(editor, 0, 1)
    assert fmt.foreground().color().name() == '#0057b7'
    assert fmt.fontWeight() > 50  # bold


def test_rem_line_is_comment(editor):
    editor.setPlainText('REM this is a comment')
    fmt = _format_at(editor, 0, 5)
    assert fmt.fontItalic()


def test_restraint_with_residue_class_suffix_is_highlighted(editor):
    editor.setPlainText('SADI_CCF3 0.02 C1 C2 C1 C3')
    fmt = _format_at(editor, 0, 1)
    assert fmt.foreground().color().name() == '#0057b7'


def test_indented_non_continuation_line_is_comment(editor):
    editor.setPlainText('CELL 0.71073 10 10 10 90 90 90\n   this looks like a comment')
    fmt = _format_at(editor, 1, 3)
    assert fmt.fontItalic()


def test_continuation_marker_is_highlighted(editor):
    editor.setPlainText('SADI 0.02 C1 C2 =\n C1 C3')
    fmt = _format_at(editor, 0, 16)  # the '=' character
    assert fmt.foreground().color().name() == '#800080'


def test_inline_bang_comment_is_highlighted(editor):
    editor.setPlainText('CELL 0.71073 10 10 10 90 90 90 ! a comment')
    bang_pos = len('CELL 0.71073 10 10 10 90 90 90 ')
    fmt = _format_at(editor, 0, bang_pos)
    assert fmt.fontItalic()
