"""Optional Qt6/qtpy GUI components for editing SHELX files.

This subpackage is **not** imported by ``shelxfile/__init__.py`` — it requires
the ``gui`` extra (``pip install shelxfile[gui]``) plus a Qt binding
(PyQt6/PySide6/PyQt5/PySide2) supported by ``qtpy``. Import the pieces you
need explicitly, e.g.::

    from shelxfile.gui.editor_widget import ShelxEditorWidget
"""
