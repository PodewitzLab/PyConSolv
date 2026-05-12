import os
import unittest

import tests  # noqa: F401
from PyConSolv.utils.charge import ChargeChanger
from tests.helpers import TempDir, write


SAMPLE_MOL2 = """@<TRIPOS>MOLECULE
OLD
 2 1 1 0 1
SMALL
bcc

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 c3        1 TEMP       0.000000
      2 H1          1.0900    0.0000    0.0000 h1        1 TEMP       0.000000
@<TRIPOS>BOND
     1    1    2 1
@<TRIPOS>SUBSTRUCTURE
     1 TEMP              0 ****  ****    0 ROOT
"""


class TestChargeChanger(unittest.TestCase):
    def test_change_replaces_resname_and_charges(self):
        with TempDir() as d:
            fin = write(os.path.join(d, 'in.mol2'), SAMPLE_MOL2)
            fout = os.path.join(d, 'out.mol2')
            # charges is list-of-lists (iterator picks [-1])
            charges = [[0.123], [-0.123]]
            changer = ChargeChanger()
            changer.change(fin, fout, 'LIG', charges)
            with open(fout) as f:
                text = f.read()
            self.assertIn('LIG', text)
            self.assertIn('0.123', text)
            self.assertIn('-0.123', text)
            # bcc gets replaced by RESP Charge
            self.assertIn('RESP Charge', text)
            # molecule name switched to resname
            self.assertIn('\nLIG\n', text)

    def test_change_increments_iterator(self):
        changer = ChargeChanger()
        self.assertEqual(changer.iterator, 0)
        with TempDir() as d:
            fin = write(os.path.join(d, 'in.mol2'), SAMPLE_MOL2)
            fout = os.path.join(d, 'out.mol2')
            changer.change(fin, fout, 'LIG', [[0.0], [0.0]])
            self.assertEqual(changer.iterator, 2)


if __name__ == '__main__':
    unittest.main()
