import os
import unittest

import tests  # noqa: F401
from PyConSolv.utils.copier import Copier
from tests.helpers import TempDir, write


class TestCopier(unittest.TestCase):
    def test_copy_preserves_content(self):
        with TempDir() as d:
            src = write(os.path.join(d, 'a.txt'), 'hello')
            dst = os.path.join(d, 'b.txt')
            Copier(src, dst).copy()
            self.assertTrue(os.path.isfile(dst))
            with open(dst) as f:
                self.assertEqual(f.read(), 'hello')

    def test_copy_overwrites_existing(self):
        with TempDir() as d:
            src = write(os.path.join(d, 'a.txt'), 'new')
            dst = write(os.path.join(d, 'b.txt'), 'old')
            Copier(src, dst).copy()
            with open(dst) as f:
                self.assertEqual(f.read(), 'new')


if __name__ == '__main__':
    unittest.main()
