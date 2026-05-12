"""Shared helpers for the PyConSolv test suite."""
import os
import shutil
import tempfile


TEST_ROOT = os.path.dirname(os.path.abspath(__file__))
TESTFILES = os.path.join(TEST_ROOT, 'Testfiles')
REFERENCE = os.path.join(TEST_ROOT, 'TestReference')


def src_root() -> str:
    return os.path.abspath(os.path.join(TEST_ROOT, '..', 'src'))


def db_dir() -> str:
    return os.path.join(src_root(), 'PyConSolv', 'db')


def radii_files() -> tuple:
    return (os.path.join(db_dir(), 'atom-radius.txt'),
            os.path.join(db_dir(), 'metal-radius.txt'))


class TempDir:
    """Context-manager temporary directory. Also restores cwd on exit, so
    interfaces that chdir into the created directory don't leave the test
    process stranded after cleanup.
    """

    def __init__(self):
        self.path = None
        self._orig_cwd = None

    def __enter__(self):
        try:
            self._orig_cwd = os.getcwd()
        except FileNotFoundError:
            self._orig_cwd = TEST_ROOT
        self.path = tempfile.mkdtemp(prefix='pyconsolv-test-')
        return self.path

    def __exit__(self, exc_type, exc, tb):
        try:
            os.chdir(self._orig_cwd or TEST_ROOT)
        except (FileNotFoundError, OSError):
            os.chdir(TEST_ROOT)
        if self.path and os.path.isdir(self.path):
            shutil.rmtree(self.path, ignore_errors=True)


def write(path: str, text: str) -> str:
    with open(path, 'w') as f:
        f.write(text)
    return path


SAMPLE_XYZ = """3
test water
O  0.000  0.000  0.000
H  0.957  0.000  0.000
H -0.240  0.927  0.000
"""


SAMPLE_METAL_XYZ = """5
Fe tetrahedral dummy
Fe  0.000  0.000  0.000
N   2.000  0.000  0.000
N  -2.000  0.000  0.000
N   0.000  2.000  0.000
N   0.000 -2.000  0.000
"""
