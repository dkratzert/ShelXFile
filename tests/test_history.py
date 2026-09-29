"""Snapshot-based undo and redo on :class:`ShelxDocument`."""

from __future__ import annotations

import pytest

from shelxfile.edit import ShelxDocument

RES = 'tests/resources/jkd77.res'


@pytest.fixture
def doc() -> ShelxDocument:
    return ShelxDocument.from_file(RES)


def names(doc: ShelxDocument) -> list[str]:
    return [a.name for a in doc.shelxfile.atoms]


def test_fresh_document_has_no_history(doc):
    assert not doc.can_undo
    assert not doc.can_redo
    assert not doc.is_modified
    assert doc.undo() is None
    assert doc.redo() is None


def test_undo_and_redo_restore_the_text(doc):
    original = doc.text
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1'), 'C99')
    edited = doc.text
    assert doc.is_modified

    result = doc.undo()
    assert result.label == 'Rename atom'
    assert doc.text == original
    assert not doc.is_modified
    assert doc.can_redo

    doc.redo()
    assert doc.text == edited
    assert doc.is_modified


def test_failed_edit_leaves_no_step(doc):
    report = doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1'), 'C2')
    assert not report.ok
    assert not doc.can_undo


def test_batch_is_one_step_with_one_notification(doc):
    calls = []
    doc.subscribe(calls.append)
    with doc.batch('Two renames'):
        doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1'), 'C91')
        doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C2'), 'C92')
    assert len(calls) == 1
    assert doc.history.undo_labels == ['Two renames']
    doc.undo()
    assert 'C1' in names(doc) and 'C2' in names(doc)


def test_failing_batch_is_rolled_back(doc):
    original = doc.text
    with pytest.raises(RuntimeError):
        with doc.batch('broken'):
            doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1'), 'C91')
            raise RuntimeError('boom')
    assert doc.text == original
    assert not doc.can_undo


def test_new_edit_clears_redo(doc):
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1'), 'C91')
    doc.undo()
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C2'), 'C92')
    assert not doc.can_redo


def test_payload_follows_the_restored_state(doc):
    state = {'value': 'before'}
    doc.set_state_provider(lambda: dict(state))
    with doc.batch('edit'):
        doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1'), 'C91')
        state['value'] = 'after'
    assert doc.undo().payload == {'value': 'before'}
    assert doc.redo().payload == {'value': 'after'}


def test_saving_marks_the_state(doc, tmp_path):
    doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1'), 'C91')
    doc.write(tmp_path / 'out.ins')
    assert not doc.is_modified
    doc.undo()
    assert doc.is_modified
    doc.redo()
    assert not doc.is_modified


def test_undo_replaces_the_model_object(doc):
    before = doc.shelxfile
    doc.rename_atom(before.atoms.get_atom_by_name('C1'), 'C91')
    doc.undo()
    assert doc.shelxfile is not before
    assert doc.shelxfile.resfile == before.resfile


def test_undo_inside_batch_is_refused(doc):
    with pytest.raises(RuntimeError):
        with doc.batch('x'):
            doc.undo()


def test_edit_that_mutates_without_touching_is_rolled_back(doc):
    """A primitive that raises after mutating must leave nothing behind.

    ``set_uvals`` used to mutate its atoms one by one and only ``touch()``
    the model afterwards, so a mid-loop failure sidestepped the rollback and
    left the model flattened while ``is_modified`` claimed otherwise.
    """
    original = doc.text
    first = doc.shelxfile.atoms.get_atom_by_name('C1')
    second = doc.shelxfile.atoms.get_atom_by_name('C2')
    before = list(first.uvals)

    with pytest.raises(ValueError):
        doc.set_uvals([(first, [0.035]), (second, [0.1, 0.2])])

    assert doc.text == original
    assert doc.shelxfile.dumps() == original
    assert doc.shelxfile.atoms.get_atom_by_name('C1').uvals == before
    assert not doc.is_modified
    assert not doc.can_undo


def test_a_batch_that_cancels_itself_out_leaves_no_step(doc):
    with doc.batch('Rename and back'):
        doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C1'), 'C99')
        doc.rename_atom(doc.shelxfile.atoms.get_atom_by_name('C99'), 'C1')
    assert not doc.can_undo
    assert not doc.is_modified
