import os
import unittest

import tests  # noqa: F401
from PyConSolv.misc.restart import RestartFile
from tests.helpers import TempDir, write


class TestRestartFile(unittest.TestCase):
    def test_default_state(self):
        with TempDir() as d:
            write(os.path.join(d, 'pyconsolv.restart'), 'unknown_keyword\n')
            r = RestartFile(d)
            self.assertEqual(r.getstate(), 0)

    def test_state_progression(self):
        expected = {
            'setup': 1, 'orca': 2, 'antechamber': 3, 'frcmod': 4,
            'multiwfn': 5, 'mcpb': 6, 'tleap': 7, 'equilibration': 8,
            'DONE': 9,
        }
        for keyword, state in expected.items():
            with TempDir() as d:
                write(os.path.join(d, 'pyconsolv.restart'), keyword + '\n')
                self.assertEqual(RestartFile(d).getstate(), state,
                                 'keyword {} should map to {}'.format(keyword, state))

    def test_write_creates_file(self):
        with TempDir() as d:
            r = RestartFile(d)
            r.write('orca')
            with open(os.path.join(d, 'pyconsolv.restart')) as f:
                self.assertEqual(f.read(), 'orca')

    def test_write_overwrites(self):
        with TempDir() as d:
            r = RestartFile(d)
            r.write('setup')
            r.write('tleap')
            self.assertEqual(r.getstate(), 7)

    def test_parseInput_is_noop(self):
        with TempDir() as d:
            self.assertIsNone(RestartFile(d).parseInput())
            self.assertIsNone(RestartFile(d).writeInput())


if __name__ == '__main__':
    unittest.main()
