"""Make the GUI test suite run headless (no real display needed)."""

import os

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
