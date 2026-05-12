import os
import unittest
from unittest import mock

import tests  # noqa: F401
from PyConSolv.interfaces.cgenff import CGenFFInterface
from tests.helpers import TempDir, write


SAMPLE_STR = """* CGenFF stream file
read rtf card append
27 1
MASS -1 CG331 12.011
RESI LIG 0.00 ! penalty= 2.500
ATOM C1 CG331 -0.27 ! penalty 12.000
BOND C1 H1
END

read para card flex append
BONDS
CG331 HGA1 322.0 1.111 ! penalty= 8.500
END
"""


class TestCGenFF(unittest.TestCase):
    def test_checkpath_missing(self):
        iface = CGenFFInterface(cgenff_path='/nonexistent/cgenff-binary')
        self.assertFalse(iface.checkpath())
        self.assertEqual(iface.status, 0)

    @mock.patch('shutil.which', return_value='/usr/local/bin/cgenff')
    def test_checkpath_found(self, _which):
        iface = CGenFFInterface()
        self.assertTrue(iface.checkpath())

    def test_parse_penalties(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'a.str'), SAMPLE_STR)
            iface = CGenFFInterface()
            pens = iface.parsePenalties(path)
            self.assertGreater(len(pens), 0)
            # all values must be numeric
            for v in pens.values():
                self.assertIsInstance(v, float)

    def test_split_str(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'a.str'), SAMPLE_STR)
            iface = CGenFFInterface()
            rtf, prm = iface.splitSTR(path)
            self.assertIn('RESI LIG', rtf)
            self.assertIn('MASS', rtf)
            self.assertIn('BONDS', prm)
            self.assertIn('CG331', prm)

    def test_parametrize_runs_subprocess(self):
        with TempDir() as d:
            mol2 = write(os.path.join(d, 'x.mol2'), 'dummy')
            iface = CGenFFInterface(cgenff_path='cgenff')
            with mock.patch('subprocess.run') as run:
                run.return_value = mock.Mock(returncode=0)
                # fabricate the expected output
                with open(os.path.join(d, 'x.str'), 'w') as f:
                    f.write(SAMPLE_STR)
                out = iface.parametrize(mol2, output_dir=d)
            self.assertTrue(out.endswith('x.str'))
            self.assertEqual(iface.status, 1)
            # command contained both input and output paths
            called = run.call_args[0][0][0]
            self.assertIn('x.mol2', called)
            self.assertIn('x.str', called)

    def test_parametrize_failure(self):
        with TempDir() as d:
            mol2 = write(os.path.join(d, 'x.mol2'), 'dummy')
            iface = CGenFFInterface(cgenff_path='cgenff')
            with mock.patch('subprocess.run') as run:
                run.return_value = mock.Mock(returncode=1)
                out = iface.parametrize(mol2, output_dir=d)
            self.assertEqual(out, '')
            self.assertEqual(iface.status, 0)


if __name__ == '__main__':
    unittest.main()
