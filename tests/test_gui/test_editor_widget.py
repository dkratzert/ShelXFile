from __future__ import annotations

import pytest

from shelxfile import Shelxfile
from shelxfile.gui.editor_widget import AddAtomDialog, AddRestraintDialog, ShelxEditorWidget

RESOURCE = 'tests/resources/p21c.res'


@pytest.fixture
def shx() -> Shelxfile:
    s = Shelxfile()
    s.read_file(RESOURCE)
    return s


@pytest.fixture
def widget(qtbot, shx) -> ShelxEditorWidget:
    w = ShelxEditorWidget(shx)
    qtbot.addWidget(w)
    return w


def test_set_shelxfile_populates_editor_text(widget, shx):
    assert widget.editor.toPlainText() == shx.dumps()
    assert len(widget._line_index_map) == widget.editor.document().blockCount()


def test_jump_to_atom_moves_cursor_and_highlights(widget, shx):
    atom = shx.atoms.all_atoms[0]
    assert widget.jump_to_atom(atom.fullname_short) is True
    block_no = widget.editor.textCursor().blockNumber()
    assert widget._line_index_map[block_no] == shx._reslist.index(atom)
    assert len(widget.editor.extraSelections()) == 1


def test_jump_to_atom_unknown_name_returns_false(widget):
    assert widget.jump_to_atom('DOES_NOT_EXIST_1') is False


def test_atom_selected_emitted_on_cursor_move(qtbot, widget, shx):
    from qtpy.QtGui import QTextCursor

    atom = shx.atoms.all_atoms[2]
    block_no = widget._line_index_map.index(shx._reslist.index(atom))

    with qtbot.waitSignal(widget.atom_selected, timeout=1000) as blocker:
        cursor = QTextCursor(widget.editor.document().findBlockByNumber(block_no))
        widget.editor.setTextCursor(cursor)
    assert blocker.args == [atom.fullname_short]


def test_atom_selected_not_emitted_during_jump(qtbot, widget, shx):
    atom = shx.atoms.all_atoms[3]
    received = []
    widget.atom_selected.connect(received.append)
    widget.jump_to_atom(atom.fullname_short)
    assert received == []


def test_delete_selected_atoms_updates_model(widget, shx):
    from qtpy.QtGui import QTextCursor

    n_before = len(shx.atoms)
    atom = shx.atoms.all_atoms[0]
    block_no = widget._line_index_map.index(shx._reslist.index(atom))
    cursor = QTextCursor(widget.editor.document().findBlockByNumber(block_no))
    widget.editor.setTextCursor(cursor)

    assert widget.delete_selected_atoms() is True
    assert len(shx.atoms) == n_before - 1
    assert atom not in shx._reslist


def test_delete_selected_atoms_blocked_while_dirty(widget, shx):
    widget.editor.insertPlainText(' ')  # hand-edit -> dirty
    assert widget.is_dirty is True
    n_before = len(shx.atoms)
    assert widget.delete_selected_atoms() is False
    assert len(shx.atoms) == n_before
    assert not widget.error_label.isHidden()


def test_add_and_delete_restraint_updates_model(widget, shx):
    from qtpy.QtGui import QTextCursor

    n_before = len(shx.restraints._restraints)
    restraint = shx.add_restraint('SADI 0.02 C1 C2')
    widget._refresh_text_from_model()
    assert len(shx.restraints._restraints) == n_before + 1

    block_no = widget._line_index_map.index(shx._reslist.index(restraint))
    cursor = QTextCursor(widget.editor.document().findBlockByNumber(block_no))
    widget.editor.setTextCursor(cursor)

    assert widget.delete_selected_restraint() is True
    assert len(shx.restraints._restraints) == n_before


def test_apply_reparses_and_swaps_model(widget, shx):
    old_model = widget.shelxfile
    assert widget.apply() is True
    assert widget.shelxfile is not old_model
    assert isinstance(widget.shelxfile, Shelxfile)


def test_apply_on_broken_text_keeps_old_model(widget):
    old_model = widget.shelxfile
    widget.editor.setPlainText('this is not a shelx file at all')
    assert widget.apply() is False
    assert widget.shelxfile is old_model
    assert not widget.error_label.isHidden()


def test_add_atom_dialog_values(qtbot, shx):
    dialog = AddAtomDialog(shx)
    qtbot.addWidget(dialog)
    dialog.name_edit.setText('C999')
    dialog.x_edit.setValue(0.1)
    dialog.y_edit.setValue(0.2)
    dialog.z_edit.setValue(0.3)
    values = dialog.values()
    assert values['name'] == 'C999'
    assert values['coordinates'] == [0.1, 0.2, 0.3]
    n_before = len(shx.atoms)
    shx.add_atom(**values)
    assert len(shx.atoms) == n_before + 1


def test_add_restraint_dialog_instruction_line(qtbot):
    dialog = AddRestraintDialog()
    qtbot.addWidget(dialog)
    idx = dialog.keyword_combo.findText('SADI')
    dialog.keyword_combo.setCurrentIndex(idx)
    dialog.params_edit.setText('0.02 C1 C2')
    assert dialog.instruction_line() == 'SADI 0.02 C1 C2'
