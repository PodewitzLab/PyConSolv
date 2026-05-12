import os
import unittest

import tests  # noqa: F401

try:
    from PyConSolv.misc.fragmenting import Fragmentor
    HAS_FRAGMENTOR = True
except Exception:
    HAS_FRAGMENTOR = False

from tests.helpers import TempDir, write


SAMPLE_PT = """11
Pt complex
Pt -2.87200000 0.14700000 0.17500000
N  -0.89000000 0.63600000 0.53600000
H  -0.77000000 1.54500000 1.01000000
H  -0.31200000 0.65800000 -0.31800000
H  -0.49200000 -0.09000000 1.15800000
N  -3.19200000 1.80700000 -1.02600000
H  -3.05300000 2.70200000 -0.53200000
H  -4.17500000 1.78100000 -1.34800000
H  -2.59300000 1.82500000 -1.86600000
Cl -2.48100000 -1.71900000 1.54000000
Cl -5.11600000 -0.38200000 -0.25300000
"""


@unittest.skipUnless(HAS_FRAGMENTOR, 'RDKit not available')
class TestFragmentor(unittest.TestCase):
    def setUp(self):
        self.tmp = TempDir()
        self.d = self.tmp.__enter__()
        self.xyz_path = write(os.path.join(self.d, 'input.xyz'), SAMPLE_PT)
        self.frag = Fragmentor(path=self.xyz_path, radius=4.0)

    def tearDown(self):
        self.tmp.__exit__(None, None, None)

    def test_findMetalIndices_returns_platinum(self):
        self.assertEqual(self.frag.metal_indices, [0])

    def test_checkRadius_keeps_all_within_radius(self):
        self.frag.checkRadius()
        self.assertIn(0, self.frag.keep)
        # Pt-N bonds < 4 Å so the Ns should be included
        self.assertIn(1, self.frag.keep)
        self.assertIn(5, self.frag.keep)

    def test_prepareXYZ_populates_metal_links(self):
        # After __init__ prepareXYZ has been called
        # Metal atom's linkList should contain at least one ligand index
        self.assertTrue(len(self.frag.xyz.linkList[0]) > 0)

    def test_writeXYZ_creates_file(self):
        coords = '2\ntest\nH 0 0 0\nH 1 0 0\n'
        self.frag.writeXYZ(coords, filename='frag.xyz')
        out = os.path.join(self.d, 'frag.xyz')
        self.assertTrue(os.path.isfile(out))
        with open(out) as f:
            self.assertEqual(f.read(), coords)


if __name__ == '__main__':
    unittest.main()
