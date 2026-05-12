"""Test package: inserts the in-tree ``src/`` directory at the front of
sys.path so ``from PyConSolv...`` resolves to the working copy rather
than the version installed into site-packages.
"""
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
