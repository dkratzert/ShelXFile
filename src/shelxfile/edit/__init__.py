"""Qt-free editing layer for ShelXFile.

Everything that understands SHELXL semantics and mutates the model lives
here.  Nothing under this package may import Qt: the GUI depends on the
edit layer, never the other way round.

See ``tests/test_layering.py`` for the guards that enforce this.
"""

from shelxfile.edit.card_meta import CardLifetime

__all__ = ['CardLifetime']
