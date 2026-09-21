"""A reusable Qt6/qtpy widget for viewing and editing SHELX ``.res``/``.ins``
files, backed by a :class:`~shelxfile.edit.document.ShelxDocument`.

The widget is a **pure view** (plan decision D-8b).  It renders what the
document gives it, reports what the user pointed at, and asks the document
to make changes; it never indexes the model's private line list and never
decides what a SHELXL instruction means.  All such knowledge lives in
``shelxfile.edit``, which is Qt-free and testable without a display.

Every mutation made through the widget's toolbar actions (add/delete atom,
add/delete restraint) goes through the document, which then regenerates the
text — so the model stays the single source of truth.

Free-hand text edits are only applied to the model when the user (or host
application) explicitly calls :meth:`ShelxEditorWidget.apply`.

This widget is designed to be embedded by a host application (e.g.
`Fastmolwidget <https://github.com/dkratzert/Fastmolwidget>`_) that renders
the same structure in 3D. Two hooks make that integration a single
``connect()`` call in each direction:

* :meth:`jump_to_atom` — a slot accepting an atom's
  ``fullname_short`` (e.g. ``"C1_1"``), which moves the cursor to that atom's
  line and highlights it. Fastmolwidget's ``MoleculeWidget.atomClicked``
  signal already emits this same string, so a host can simply do
  ``mol_widget.atomClicked.connect(editor.jump_to_atom)``.
* :attr:`atom_selected` — a signal emitting an atom's ``fullname_short``
  whenever the cursor moves onto that atom's line, so a host can highlight
  the same atom in its 3D view.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from qtpy.QtCore import Signal
from qtpy.QtGui import QColor, QFont, QTextCharFormat, QTextCursor
from qtpy.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from shelxfile.atoms.atom import Atom
from shelxfile.edit import ShelxDocument, restraint_keywords
from shelxfile.gui.syntax_highlighter import ShelxSyntaxHighlighter
from shelxfile.shelx.cards import Restraint

if TYPE_CHECKING:
    from shelxfile.shelx.shelx import Shelxfile

__all__ = ['ShelxEditorWidget', 'AddAtomDialog', 'AddRestraintDialog', 'main']


class ShelxEditorWidget(QWidget):
    """
    A syntax-highlighted SHELX text editor bound to a :class:`ShelxDocument`.
    """

    #: Emitted with the (possibly new) :class:`Shelxfile` instance after a
    #: successful :meth:`apply`, or after any add/delete action.
    model_changed = Signal(object)

    #: Emitted with an error message and a 0-based line number (or -1 if
    #: unknown) whenever :meth:`apply` fails to parse the edited text.
    parse_error = Signal(str, int)

    #: Emitted with an atom's ``fullname_short`` (e.g. ``"C1_1"``) whenever
    #: the cursor moves onto that atom's line. Not emitted while
    #: :meth:`jump_to_atom` is programmatically moving the cursor, to avoid
    #: a feedback loop with a host's 3D viewer.
    atom_selected = Signal(str)

    def __init__(self, shelxfile: Shelxfile | ShelxDocument | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._document: ShelxDocument | None = None
        self._dirty_since_apply: bool = False
        self._updating_text: bool = False
        self._jumping: bool = False
        self._build_ui()
        if shelxfile is not None:
            self.set_shelxfile(shelxfile)

    # ------------------------------------------------------------------ UI

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)

        toolbar = QHBoxLayout()
        self.apply_button = QPushButton('Apply')
        self.add_atom_button = QPushButton('Add atom…')
        self.delete_atom_button = QPushButton('Delete selected atom(s)')
        self.add_restraint_button = QPushButton('Add restraint…')
        self.delete_restraint_button = QPushButton('Delete selected restraint')
        for button in (
                self.apply_button, self.add_atom_button, self.delete_atom_button,
                self.add_restraint_button, self.delete_restraint_button,
        ):
            toolbar.addWidget(button)
        toolbar.addStretch(1)
        layout.addLayout(toolbar)

        self.error_label = QLabel('')
        self.error_label.setStyleSheet('color: #a00000;')
        self.error_label.setWordWrap(True)
        self.error_label.hide()
        layout.addWidget(self.error_label)

        self.editor = QPlainTextEdit()
        self.editor.setFont(QFont('Courier New', 10))
        self.editor.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.highlighter = ShelxSyntaxHighlighter(self.editor.document())
        layout.addWidget(self.editor)

        self.apply_button.clicked.connect(self.apply)
        self.add_atom_button.clicked.connect(self._on_add_atom_clicked)
        self.delete_atom_button.clicked.connect(self.delete_selected_atoms)
        self.add_restraint_button.clicked.connect(self._on_add_restraint_clicked)
        self.delete_restraint_button.clicked.connect(self.delete_selected_restraint)
        self.editor.textChanged.connect(self._on_text_changed)
        self.editor.cursorPositionChanged.connect(self._on_cursor_position_changed)

    # --------------------------------------------------------------- model

    @property
    def document(self) -> ShelxDocument | None:
        """The :class:`ShelxDocument` currently bound to this widget."""
        return self._document

    @property
    def shelxfile(self) -> Shelxfile | None:
        """The :class:`Shelxfile` model behind the bound document."""
        return None if self._document is None else self._document.shelxfile

    def set_document(self, document: ShelxDocument) -> None:
        """Bind *document* to this widget and refresh the text view from it."""
        if self._document is not None:
            self._document.unsubscribe(self._on_document_changed)
        self._document = document
        document.subscribe(self._on_document_changed)
        self._hide_error()
        self._refresh_text_from_model(preserve_cursor=False)

    def set_shelxfile(self, shx: Shelxfile | ShelxDocument) -> None:
        """Bind a model (or a document) to this widget and refresh the view."""
        self.set_document(shx if isinstance(shx, ShelxDocument) else ShelxDocument(shx))

    def _on_document_changed(self, _document: ShelxDocument) -> None:
        self._refresh_text_from_model()

    def text(self) -> str:
        """The editor's current (possibly unapplied) text content."""
        return self.editor.toPlainText()

    @property
    def is_dirty(self) -> bool:
        """Whether the text was hand-edited since the last successful :meth:`apply`."""
        return self._dirty_since_apply

    # ------------------------------------------------------ text <-> model

    def _refresh_text_from_model(self, *, preserve_cursor: bool = True) -> None:
        if self._document is None:
            return
        cursor = self.editor.textCursor()
        block_no = cursor.blockNumber()
        column = cursor.positionInBlock()
        vscroll = self.editor.verticalScrollBar().value()

        self._updating_text = True
        try:
            self.editor.setPlainText(self._document.text)
        finally:
            self._updating_text = False
        self._dirty_since_apply = False

        if preserve_cursor:
            doc = self.editor.document()
            block = doc.findBlockByNumber(min(block_no, doc.blockCount() - 1))
            new_cursor = QTextCursor(block)
            new_cursor.movePosition(
                QTextCursor.MoveOperation.Right,
                QTextCursor.MoveMode.MoveAnchor,
                min(column, len(block.text())),
            )
            self.editor.setTextCursor(new_cursor)
            self.editor.verticalScrollBar().setValue(vscroll)

    def _selected_blocks(self) -> range:
        cursor = self.editor.textCursor()
        doc = self.editor.document()
        start = doc.findBlock(cursor.selectionStart()).blockNumber()
        end = doc.findBlock(cursor.selectionEnd()).blockNumber()
        return range(start, end + 1)

    def _selected_items(self) -> list[object]:
        """Whatever the document says sits on the selected lines."""
        if self._document is None:
            return []
        blocks = self._selected_blocks()
        return self._document.items_in_lines(blocks.start, blocks.stop - 1)

    # ------------------------------------------------------------- signals

    def _on_text_changed(self) -> None:
        if not self._updating_text:
            self._dirty_since_apply = True

    def _on_cursor_position_changed(self) -> None:
        if self._jumping or self._document is None:
            return
        item = self._document.item_at_line(self.editor.textCursor().blockNumber())
        if isinstance(item, Atom):
            self.atom_selected.emit(item.fullname_short)

    # --------------------------------------------------------------- apply

    def apply(self) -> bool:
        """
        Re-parse the current editor text into a fresh document.

        On success, the new document replaces the previously bound one, the
        text is regenerated from it (normalising formatting), and
        :attr:`model_changed` is emitted. On failure, the text and the old
        document are left untouched and :attr:`parse_error` is emitted
        together with an inline error message.
        """
        if self._document is None:
            return False
        attempt = ShelxDocument.try_from_string(self.editor.toPlainText())
        if attempt.document is None:
            self._show_error(attempt.error or 'Could not parse SHELX file.',
                             attempt.error_line)
            return False
        self.set_document(attempt.document)
        self.model_changed.emit(attempt.document.shelxfile)
        return True

    def _show_error(self, message: str, line_num: int) -> None:
        if line_num is not None and line_num >= 0:
            message = f'Line {line_num + 1}: {message}'
        self.error_label.setText(message)
        self.error_label.show()
        self.parse_error.emit(message, line_num if line_num is not None else -1)

    def _hide_error(self) -> None:
        self.error_label.hide()
        self.error_label.setText('')

    def _require_applied_model(self) -> bool:
        """Guard against mutating a model while the text has unapplied hand-edits."""
        if self._dirty_since_apply:
            self._show_error('Apply your text edits before adding or deleting atoms/restraints.', -1)
            return False
        return True

    # ---------------------------------------------------------------- atoms

    def _on_add_atom_clicked(self) -> None:
        if self._document is None or not self._require_applied_model():
            return
        dialog = AddAtomDialog(self._document, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self._document.add_atom(**dialog.values())
            self.model_changed.emit(self._document.shelxfile)

    def delete_selected_atoms(self) -> bool:
        """Delete every :class:`Atom` touched by the current text selection.

        The document decides what else cannot survive the deletion; the
        widget only says which atoms the user pointed at.
        """
        if self._document is None or not self._require_applied_model():
            return False
        atoms = [item for item in self._selected_items() if isinstance(item, Atom)]
        if not atoms:
            return False
        self._document.delete_atoms(atoms)
        self.model_changed.emit(self._document.shelxfile)
        return True

    # ----------------------------------------------------------- restraints

    def _on_add_restraint_clicked(self) -> None:
        if self._document is None or not self._require_applied_model():
            return
        dialog = AddRestraintDialog(self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            line = dialog.instruction_line()
            if not line:
                return
            try:
                self._document.add_restraint(line)
            except ValueError as exc:
                self._show_error(str(exc), -1)
                return
            self.model_changed.emit(self._document.shelxfile)

    def delete_selected_restraint(self) -> bool:
        """Remove every :class:`Restraint` touched by the current selection.

        Removing a restraint never deletes atoms (plan decision D-9).
        """
        if self._document is None or not self._require_applied_model():
            return False
        restraints = [item for item in self._selected_items() if isinstance(item, Restraint)]
        if not restraints:
            return False
        for restraint in restraints:
            self._document.remove_card(restraint)
        self.model_changed.emit(self._document.shelxfile)
        return True

    # -------------------------------------------- Fastmolwidget integration

    def jump_to_atom(self, fullname_short: str) -> bool:
        """
        Move the cursor to, and highlight, the line of the atom identified by
        *fullname_short* (e.g. ``"C1_1"``, as produced by
        :attr:`shelxfile.atoms.atom.Atom.fullname_short` and emitted by
        Fastmolwidget's ``MoleculeWidget.atomClicked`` signal).

        Intended to be connected directly to that signal, e.g.::

            mol_widget.atomClicked.connect(editor.jump_to_atom)

        :returns: ``True`` if the atom was found and the cursor moved,
                  ``False`` otherwise (unknown atom, or unapplied hand-edits
                  make the line map unreliable).
        """
        if self._document is None or self._dirty_since_apply:
            return False
        block_no = self._document.line_of_atom(fullname_short)
        if block_no is None:
            return False

        self._jumping = True
        try:
            doc = self.editor.document()
            block = doc.findBlockByNumber(block_no)
            cursor = QTextCursor(block)
            cursor.movePosition(QTextCursor.MoveOperation.EndOfBlock, QTextCursor.MoveMode.KeepAnchor)
            self.editor.setTextCursor(cursor)
            self.editor.centerCursor()
            self._highlight_line(block_no)
        finally:
            self._jumping = False
        return True

    def _highlight_line(self, block_no: int) -> None:
        """Give the line at *block_no* a temporary background highlight."""
        selection = QTextEdit.ExtraSelection()
        selection.format = QTextCharFormat()
        selection.format.setBackground(QColor('#fff2a8'))
        selection.format.setProperty(QTextCharFormat.Property.FullWidthSelection, True)
        doc = self.editor.document()
        cursor = QTextCursor(doc.findBlockByNumber(block_no))
        selection.cursor = cursor
        self.editor.setExtraSelections([selection])


class AddAtomDialog(QDialog):
    """Small dialog collecting the parameters for :meth:`ShelxDocument.add_atom`."""

    def __init__(self, document: ShelxDocument, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle('Add atom')
        self._document = document

        from qtpy.QtWidgets import QDoubleSpinBox, QLineEdit, QSpinBox

        self.name_edit = QLineEdit(document.unused_atom_name('C'))
        self.element_edit = QLineEdit('C')
        self.x_edit = QDoubleSpinBox()
        self.y_edit = QDoubleSpinBox()
        self.z_edit = QDoubleSpinBox()
        for spin in (self.x_edit, self.y_edit, self.z_edit):
            spin.setDecimals(6)
            spin.setRange(-10.0, 10.0)
        self.occupancy_edit = QDoubleSpinBox()
        self.occupancy_edit.setRange(0.0, 1.0)
        self.occupancy_edit.setValue(1.0)
        self.resi_edit = QSpinBox()
        self.resi_edit.setRange(0, 9999)
        self.part_edit = QSpinBox()
        self.part_edit.setRange(-9, 9)

        form = QFormLayout()
        form.addRow('Name', self.name_edit)
        form.addRow('Element', self.element_edit)
        form.addRow('x', self.x_edit)
        form.addRow('y', self.y_edit)
        form.addRow('z', self.z_edit)
        form.addRow('Occupancy', self.occupancy_edit)
        form.addRow('Residue number', self.resi_edit)
        form.addRow('Part', self.part_edit)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def values(self) -> dict:
        """The dialog's contents as ``document.add_atom(**values())`` keyword arguments."""
        return {
            'name': self.name_edit.text().strip(),
            'coordinates': [self.x_edit.value(), self.y_edit.value(), self.z_edit.value()],
            'element': self.element_edit.text().strip() or 'C',
            'occupancy': self.occupancy_edit.value(),
            'resi': self.resi_edit.value(),
            'part': self.part_edit.value(),
        }


class AddRestraintDialog(QDialog):
    """Small dialog assembling a single restraint instruction line."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle('Add restraint')

        from qtpy.QtWidgets import QComboBox, QLineEdit

        self.keyword_combo = QComboBox()
        self.keyword_combo.addItems(restraint_keywords())
        self.params_edit = QLineEdit()
        self.params_edit.setPlaceholderText('e.g. 0.02 C1 C2 C1 C3')

        form = QFormLayout()
        form.addRow('Restraint', self.keyword_combo)
        form.addRow('Parameters / atoms', self.params_edit)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def instruction_line(self) -> str:
        """The assembled instruction line, e.g. ``'SADI 0.02 C1 C2 C1 C3'``."""
        params = self.params_edit.text().strip()
        keyword = self.keyword_combo.currentText()
        return f'{keyword} {params}'.strip()


def _default_demo_file() -> Path | None:
    """Fall back to the bundled test resource, if this is a checkout of the repo."""
    candidate = Path(__file__).resolve().parents[3] / 'tests' / 'resources' / 'p21c.res'
    return candidate if candidate.exists() else None


def main(argv: list[str] | None = None) -> int:
    """
    Run :class:`ShelxEditorWidget` as a standalone window.

    Usage::

        python -m shelxfile.gui.editor_widget [path/to/file.res]

    If no path is given, the bundled ``tests/resources/p21c.res`` demo file is
    used when available. This entry point exists for manual smoke-testing of
    the widget during development; it is not part of the public API used by
    host applications (which should embed :class:`ShelxEditorWidget` directly).
    """
    import sys

    from qtpy.QtWidgets import QApplication, QMainWindow

    argv = sys.argv[1:] if argv is None else argv
    resfile = Path(argv[0]) if argv else _default_demo_file()
    if resfile is None or not resfile.exists():
        print('Usage: python -m shelxfile.gui.editor_widget <path/to/file.res>')
        return 1

    app = QApplication(sys.argv[:1])

    document = ShelxDocument.from_file(resfile, verbose=True)

    editor = ShelxEditorWidget(document)
    editor.model_changed.connect(lambda model: print(f'*** Model updated ({len(model.atoms)} atoms) ***'))
    editor.parse_error.connect(lambda message, line: print(f'*** Parse error: {message} ***'))
    editor.atom_selected.connect(lambda name: print(f'Selected atom: {name}'))

    window = QMainWindow()
    window.setWindowTitle(f'ShelxEditorWidget \u2014 {resfile}')
    window.setCentralWidget(editor)
    window.resize(900, 700)
    window.show()

    return app.exec()


if __name__ == '__main__':
    raise SystemExit(main())
