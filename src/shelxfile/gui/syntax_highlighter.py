"""SHELXL syntax highlighting for the optional GUI text editor.

Highlights instruction/restraint keywords (e.g. ``TITL``, ``CELL``, ``SADI``),
``REM`` comment lines, ``!`` inline comments (valid anywhere on a line per the
SHELXL manual), indented lines that are comments rather than ``=``-
continuations, and the ``=`` line-continuation marker.

Adapted from the ``ShelxSyntaxHighlighter`` originally written for FinalCif
(https://github.com/dkratzert/FinalCif) by the same author, reusing this
project's own :data:`~shelxfile.shelx.shelx.SHX_CARDS` keyword list instead of
duplicating it.
"""

from __future__ import annotations

from qtpy.QtGui import QColor, QFont, QSyntaxHighlighter, QTextCharFormat

from shelxfile.shelx.shelx import SHX_CARDS

#: SHELXL instruction keywords, deduplicated and stripped of the padding
#: spaces ``SHX_CARDS`` uses internally for fixed-width comparisons.
SHELX_KEYWORDS = frozenset(keyword.strip() for keyword in SHX_CARDS)


def _make_format(color: str | None = None, bold: bool = False) -> QTextCharFormat:
    """Build a QTextCharFormat."""
    fmt = QTextCharFormat()
    if color is not None:
        fmt.setForeground(QColor(color))
    if bold:
        fmt.setFontWeight(QFont.Weight.Bold)
    return fmt


class ShelxSyntaxHighlighter(QSyntaxHighlighter):
    """Syntax highlighter for SHELXL .ins/.res files."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)

        self.keyword_format = _make_format("#0057B7", bold=True)
        self.comment_format = _make_format("#808080")
        self.comment_format.setFontItalic(True)
        self.continuation_format = _make_format("#800080", bold=True)
        self.value_format = _make_format("#763127")

    @staticmethod
    def _strip_inline_comment(text: str) -> str:
        """Return *text* up to (excluding) a '!' comment, if any.

        Per the SHELXL manual, "All characters following '!' ... in an
        instruction line are ignored", so a trailing '=' continuation marker
        must be detected *before* any '!' comment on the same line.
        """
        bang_pos = text.find('!')
        return text if bang_pos == -1 else text[:bang_pos]

    def _ends_with_continuation_marker(self, text: str) -> bool:
        """Whether *text* is a '='-terminated continuation line, ignoring
        any trailing '!' comment."""
        return self._strip_inline_comment(text).rstrip().endswith('=')

    def _is_continuation_line(self) -> bool:
        """Whether the current block continues a preceding '='-terminated line.

        Reads the previous block's own text directly (rather than relying on
        ``previousBlockState()``) since that stays reliable across the
        multiple internal reformat passes Qt may perform on a block.
        """
        prev_block = self.currentBlock().previous()
        return prev_block.isValid() and self._ends_with_continuation_marker(prev_block.text())

    def _continuation_root_keyword(self) -> str:
        """Return the upper-cased first word of the line that started the
        current '='-continuation chain (the "root" instruction line), or ''
        if that root line is empty/unavailable.

        A continuation chain can itself span several '='-terminated lines,
        so this walks back to the earliest block in the unbroken chain.
        """
        block = self.currentBlock().previous()
        while True:
            earlier = block.previous()
            if earlier.isValid() and self._ends_with_continuation_marker(earlier.text()):
                block = earlier
            else:
                break
        root_stripped = block.text().strip()
        return root_stripped.split(None, 1)[0].upper() if root_stripped else ''

    def highlightBlock(self, text: str) -> None:
        stripped = text.strip()
        if not stripped:
            return

        leading_ws = len(text) - len(text.lstrip())
        is_continuation = self._is_continuation_line()

        # ---------- Indented lines that are not continuations are comments
        # (per the SHELXL manual: "Other lines beginning with spaces are
        # treated as comments"). A line following a trailing '=' is a real
        # continuation of the previous instruction, not a comment - it is
        # colored the same as the line it continues, as if there was no
        # line break at all. ----------

        if leading_ws > 0 and not is_continuation:
            self.setFormat(0, len(text), self.comment_format)
            return

        if is_continuation:
            root_word = self._continuation_root_keyword()
            if root_word == 'REM':
                self.setFormat(0, len(text), self.comment_format)
            elif root_word.split('_', 1)[0] in SHELX_KEYWORDS:
                self.setFormat(0, len(text), self.value_format)
        else:
            first_word = stripped.split(None, 1)[0]
            upper_word = first_word.upper()

            # ---------- Comments ----------

            if upper_word == 'REM':
                self.setFormat(0, len(text), self.comment_format)
                return

            # ---------- Keywords (restraint names may carry a residue-class
            # suffix, e.g. SADI_CCF3, so only the part before '_' is matched). ----------

            base_keyword = upper_word.split('_', 1)[0]
            if base_keyword in SHELX_KEYWORDS:
                self.setFormat(leading_ws, len(base_keyword), self.keyword_format)

                # ---------- Rest of the line following the keyword ----------

                rest_start = leading_ws + len(upper_word)
                if rest_start < len(text):
                    self.setFormat(rest_start, len(text) - rest_start, self.value_format)

        # ---------- Line continuation (ignoring any trailing '!' comment) ----------

        stripped_before_comment = self._strip_inline_comment(text).rstrip()
        if stripped_before_comment.endswith('='):
            self.setFormat(len(stripped_before_comment) - 1, 1, self.continuation_format)

        # ---------- Trailing '!' comment (anywhere on the line) ----------

        bang_pos = text.find('!')
        if bang_pos != -1:
            self.setFormat(bang_pos, len(text) - bang_pos, self.comment_format)

        self.setCurrentBlockState(0)
