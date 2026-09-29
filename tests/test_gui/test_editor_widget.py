"""The Qt editor widget, exercised as the pure view it is meant to be.

Every assertion here goes through the public API — :class:`ShelxDocument`
for the model side, the widget's signals and slots for the view side.
Nothing touches the model's private line list, which is precisely the
boundary ``tests/test_layering.py`` guards (plan decision D-8b).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from qtpy.QtWidgets import QVBoxLayout, QWidget

from shelxfile import Shelxfile
from shelxfile.edit import ShelxDocument
from shelxfile.gui.editor_widget import (
    AddAtomDialog,
    AddRestraintDialog,
    ShelxEditorToolbar,
    ShelxEditorWidget,
)

RESOURCE = 'tests/resources/p21c.res'


@pytest.fixture
def document() -> ShelxDocument:
    return ShelxDocument.from_file(RESOURCE)


@pytest.fixture
def shx(document) -> Shelxfile:
    return document.shelxfile


@pytest.fixture
def widget(qtbot, document) -> ShelxEditorWidget:
    w = ShelxEditorWidget(document)
    qtbot.addWidget(w)
    return w


def _select_line(widget, block_no) -> None:
    from qtpy.QtGui import QTextCursor

    cursor = QTextCursor(widget.editor.document().findBlockByNumber(block_no))
    widget.editor.setTextCursor(cursor)


# ----------------------------------------------------------------- binding


def test_set_shelxfile_populates_editor_text(widget, document):
    assert widget.editor.toPlainText() == document.text
    assert widget.editor.document().blockCount() == document.line_count


def test_accepts_a_bare_shelxfile_for_backwards_compatibility(qtbot, shx):
    w = ShelxEditorWidget(shx)
    qtbot.addWidget(w)
    assert w.shelxfile is shx
    assert isinstance(w.document, ShelxDocument)
    assert w.editor.toPlainText() == shx.dumps()


def test_document_and_shelxfile_agree(widget, document):
    assert widget.document is document
    assert widget.shelxfile is document.shelxfile


def test_clear_unbinds_the_document(widget, document):
    widget.clear()
    assert widget.document is None
    assert widget.editor.toPlainText() == ''
    assert widget.editor.isReadOnly()
    assert not widget.is_dirty
    assert all(not button.isEnabled() for button in widget.toolbar.buttons)

    # The view must no longer follow the document it was unbound from.
    document.add_restraint('SADI 0.02 C1 C2')
    assert widget.editor.toPlainText() == ''


def test_binding_again_after_clear_works(widget, document):
    widget.clear()
    widget.set_document(document)
    assert widget.editor.toPlainText() == document.text
    assert not widget.editor.isReadOnly()
    assert all(button.isEnabled() for button in widget.toolbar.buttons)


def test_an_unbound_editor_is_disabled(qtbot):
    w = ShelxEditorWidget()
    qtbot.addWidget(w)
    assert w.document is None
    assert w.editor.isReadOnly()
    assert all(not button.isEnabled() for button in w.toolbar.buttons)


# --------------------------------------------------------------- toolbar


def test_toolbar_is_a_widget_of_its_own(widget):
    """So a host can place it outside the 80-column-wide text area."""
    assert isinstance(widget.toolbar, ShelxEditorToolbar)
    assert widget.toolbar.parent() is widget
    assert widget.apply_button is widget.toolbar.apply_button


def test_toolbar_buttons_still_work_after_reparenting(qtbot, widget):
    host = QWidget()
    qtbot.addWidget(host)
    QVBoxLayout(host).addWidget(widget.toolbar)
    assert widget.toolbar.parent() is host

    widget.editor.setPlainText(widget.editor.toPlainText() + '\nREM from the toolbar\n')
    widget.apply_button.click()
    assert 'REM from the toolbar' in widget.document.text


def test_character_width_covers_84_columns(widget):
    from qtpy.QtGui import QFontMetricsF

    width = widget.character_width()
    text = QFontMetricsF(widget.editor.font()).horizontalAdvance('0' * 84)
    # The frame, document margins and scrollbar come on top of the glyphs,
    # but the toolbar's much larger width must not.
    assert width > text
    assert width < text + 120
    assert widget.character_width(40) < width


# ------------------------------------------------------------- navigation


def test_jump_to_atom_moves_cursor_and_highlights(widget, document, shx):
    atom = shx.atoms.all_atoms[0]
    assert widget.jump_to_atom(atom.fullname_short) is True
    block_no = widget.editor.textCursor().blockNumber()
    assert block_no == document.line_of(atom)
    assert len(widget.editor.extraSelections()) == 1


def test_jump_to_atom_unknown_name_returns_false(widget):
    assert widget.jump_to_atom('DOES_NOT_EXIST_1') is False


def test_atom_selected_emitted_on_cursor_move(qtbot, widget, document, shx):
    atom = shx.atoms.all_atoms[2]
    with qtbot.waitSignal(widget.atom_selected, timeout=1000) as blocker:
        _select_line(widget, document.line_of(atom))
    assert blocker.args == [atom.fullname_short]


def test_atom_selected_not_emitted_during_jump(widget, shx):
    atom = shx.atoms.all_atoms[3]
    received = []
    widget.atom_selected.connect(received.append)
    widget.jump_to_atom(atom.fullname_short)
    assert received == []


# ---------------------------------------------------------------- editing


def test_delete_selected_atoms_updates_model(widget, document, shx):
    n_before = len(shx.atoms)
    atom = shx.atoms.all_atoms[0]
    _select_line(widget, document.line_of(atom))

    assert widget.delete_selected_atoms() is True
    assert len(shx.atoms) == n_before - 1
    assert document.line_of(atom) is None


def test_delete_selected_atoms_refreshes_the_view(widget, document, shx):
    atom = shx.atoms.all_atoms[0]
    _select_line(widget, document.line_of(atom))
    widget.delete_selected_atoms()
    assert widget.editor.toPlainText() == document.text


def test_delete_selected_atoms_blocked_while_dirty(widget, shx):
    widget.editor.insertPlainText(' ')  # hand-edit -> dirty
    assert widget.is_dirty is True
    n_before = len(shx.atoms)
    assert widget.delete_selected_atoms() is False
    assert len(shx.atoms) == n_before
    assert not widget.error_label.isHidden()


def test_add_and_delete_restraint_updates_model(widget, document, shx):
    n_before = len(shx.restraints._restraints)
    report = document.add_restraint('SADI 0.02 C1 C2')
    restraint = report.added[0]
    assert len(shx.restraints._restraints) == n_before + 1
    assert widget.editor.toPlainText() == document.text  # observer refreshed it

    _select_line(widget, document.line_of(restraint))
    assert widget.delete_selected_restraint() is True
    assert len(shx.restraints._restraints) == n_before


def test_deleting_a_restraint_keeps_its_atoms(widget, document, shx):
    """D-9: card removal is never a statement about the atoms it named."""
    report = document.add_restraint('SADI 0.02 C1 C2')
    restraint = report.added[0]
    n_atoms = len(shx.atoms)

    _select_line(widget, document.line_of(restraint))
    widget.delete_selected_restraint()
    assert len(shx.atoms) == n_atoms


def test_model_changed_emitted_on_delete(qtbot, widget, document, shx):
    atom = shx.atoms.all_atoms[0]
    _select_line(widget, document.line_of(atom))
    with qtbot.waitSignal(widget.model_changed, timeout=1000) as blocker:
        widget.delete_selected_atoms()
    assert blocker.args == [shx]


# ------------------------------------------------------------------ apply


def test_apply_keeps_the_document_and_its_history(widget):
    """``apply`` swaps the *model*, in place, and keeps everything around it."""
    old_document = widget.document
    assert widget.apply() is True
    assert widget.document is old_document
    assert isinstance(widget.document, ShelxDocument)
    assert isinstance(widget.shelxfile, Shelxfile)


def test_apply_does_not_bury_an_earlier_edit(widget):
    """A step made before ``apply`` must still be reachable by undo."""
    widget.document.add_restraint('SADI 0.02 C1 C2')
    assert widget.document.history.undo_label == 'Add restraint'

    widget.editor.setPlainText(widget.editor.toPlainText() + '\nREM applied\n')
    assert widget.apply() is True
    assert widget.document.history.undo_label == 'Edit text'

    assert widget.document.undo().label == 'Edit text'
    assert widget.document.undo().label == 'Add restraint'
    assert 'SADI 0.02 C1 C2' not in widget.document.text


def test_apply_keeps_the_observer_bound(widget):
    """The view keeps following the document it was bound to."""
    document = widget.document
    widget.apply()
    document.add_restraint('SADI 0.02 C1 C2')
    assert widget.editor.toPlainText() == document.text


def test_apply_normalises_the_text_by_default(widget):
    widget.editor.setPlainText(widget.editor.toPlainText() + '\nREM normalise me\n')
    assert widget.apply() is True
    assert widget.editor.toPlainText() == widget.document.text


def test_apply_can_leave_the_text_alone(widget):
    """A live editor must not have its own input reformatted under the cursor."""
    typed = widget.editor.toPlainText() + '\nREM   kept    as-is\n'
    widget.editor.setPlainText(typed)
    assert widget.apply(normalise=False) is True
    assert widget.editor.toPlainText() == typed
    assert 'kept    as-is' in widget.editor.toPlainText()
    assert not widget.is_dirty


def test_coalescing_applies_fold_into_one_undo_step(widget):
    base = widget.editor.toPlainText()
    for n in range(3):
        widget.editor.setPlainText(f'{base}\nREM burst {n}\n')
        assert widget.apply(normalise=False, coalesce=True) is True

    assert widget.document.history.undo_labels.count('Edit text') == 1
    widget.document.undo()
    assert 'REM burst' not in widget.document.text


def test_apply_on_broken_text_keeps_old_model(widget):
    old_document = widget.document
    old_text = old_document.text
    widget.editor.setPlainText('this is not a shelx file at all')
    assert widget.apply() is False
    assert widget.document is old_document
    assert widget.document.text == old_text
    assert not widget.error_label.isHidden()


def test_parse_error_signal_carries_a_line_number(qtbot, widget):
    widget.editor.setPlainText('this is not a shelx file at all')
    with qtbot.waitSignal(widget.parse_error, timeout=1000) as blocker:
        widget.apply()
    message, line = blocker.args
    assert message
    assert isinstance(line, int)


# --------------------------------------------------------------- refine


def test_refine_without_a_file_path_shows_an_error(qtbot):
    """No ``resfile`` -> a clear inline error, no attempt to run SHELXL."""
    document = ShelxDocument.from_string(Path(RESOURCE).read_text())
    w = ShelxEditorWidget(document)
    qtbot.addWidget(w)
    assert w.refine() is False
    assert not w.error_label.isHidden()


def test_refine_applies_pending_text_edits_first(qtbot, widget):
    """A dirty editor should be committed, exactly like clicking Apply."""
    widget.editor.setPlainText('this is not a shelx file at all')
    assert widget.refine() is False
    assert not widget.error_label.isHidden()


@pytest.fixture
def refine_widget(qtbot):
    from shelxfile.refine.refine import find_shelxl_exe

    if not find_shelxl_exe():
        pytest.skip('SHELXL not found')
    resource_dir = Path('tests/resources/model_finished')
    shutil.copy(resource_dir / 'p21c.res', '.')
    shutil.copy(resource_dir / 'p21c.hkl', '.')
    w = ShelxEditorWidget(ShelxDocument.from_file('p21c.res'))
    qtbot.addWidget(w)
    yield w
    for suffix in ('.res', '.ins', '.lst', '.fcf', '.fcf6', '.cif', '.hkl', '.shx-bak'):
        Path(f'p21c{suffix}').unlink(missing_ok=True)
    shutil.rmtree('shxsaves', ignore_errors=True)


def test_refine_button_loads_the_refined_result(refine_widget):
    old_shx = refine_widget.shelxfile
    assert refine_widget.refine(cycles=1) is True
    assert refine_widget.shelxfile is not old_shx
    assert refine_widget.editor.toPlainText() == refine_widget.document.text
    assert 'L.S. 1' in refine_widget.editor.toPlainText()


def test_refine_emits_model_changed(qtbot, refine_widget):
    with qtbot.waitSignal(refine_widget.model_changed, timeout=30000):
        refine_widget.refine(cycles=1)


# ---------------------------------------------------------------- dialogs


def test_add_atom_dialog_values(qtbot, document):
    dialog = AddAtomDialog(document)
    qtbot.addWidget(dialog)
    dialog.name_edit.setText('C999')
    dialog.x_edit.setValue(0.1)
    dialog.y_edit.setValue(0.2)
    dialog.z_edit.setValue(0.3)
    values = dialog.values()
    assert values['name'] == 'C999'
    assert values['coordinates'] == [0.1, 0.2, 0.3]
    n_before = len(document.shelxfile.atoms)
    document.add_atom(**values)
    assert len(document.shelxfile.atoms) == n_before + 1


def test_add_atom_dialog_proposes_a_free_name(qtbot, document):
    dialog = AddAtomDialog(document)
    qtbot.addWidget(dialog)
    proposed = dialog.name_edit.text()
    assert proposed
    assert document.shelxfile.atoms.get_atom_by_name(proposed) is None


def test_add_restraint_dialog_instruction_line(qtbot):
    dialog = AddRestraintDialog()
    qtbot.addWidget(dialog)
    idx = dialog.keyword_combo.findText('SADI')
    dialog.keyword_combo.setCurrentIndex(idx)
    dialog.params_edit.setText('0.02 C1 C2')
    assert dialog.instruction_line() == 'SADI 0.02 C1 C2'


def test_add_restraint_dialog_offers_the_documented_vocabulary(qtbot):
    from shelxfile.edit import restraint_keywords

    dialog = AddRestraintDialog()
    qtbot.addWidget(dialog)
    offered = [dialog.keyword_combo.itemText(i) for i in range(dialog.keyword_combo.count())]
    assert offered == restraint_keywords()
