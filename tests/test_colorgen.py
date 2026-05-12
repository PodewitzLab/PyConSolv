import unittest

import tests  # noqa: F401 - sys.path fix
from PyConSolv.utils.colorgen import Color


class TestColor(unittest.TestCase):
    def test_has_expected_constants(self):
        for attr in ('RED', 'GREEN', 'YELLOW', 'END'):
            self.assertTrue(hasattr(Color, attr))
            self.assertIsInstance(getattr(Color, attr), str)

    def test_end_resets(self):
        self.assertIn('0', Color.END)

    def test_codes_are_escapes(self):
        for attr in ('RED', 'GREEN', 'YELLOW'):
            self.assertTrue(getattr(Color, attr).startswith('\033['))


if __name__ == '__main__':
    unittest.main()
