import unittest

import tests  # noqa: F401
from PyConSolv.misc.ions import ionlib


class TestIonLib(unittest.TestCase):
    def setUp(self):
        self.ions = ionlib()

    def test_common_cations(self):
        for key in ('Na+', 'K+', 'Li+', 'Cs+'):
            self.assertIn(key, self.ions.ionsinAmber)
            self.assertEqual(self.ions.ionsinAmber[key][1], 1)

    def test_common_anions(self):
        for key in ('F-', 'Cl-', 'Br-'):
            self.assertEqual(self.ions.ionsinAmber[key][1], -1)

    def test_transition_metal_charges(self):
        self.assertEqual(self.ions.ionsinAmber['Zn2+'][1], 2)
        self.assertEqual(self.ions.ionsinAmber['Fe3+'][1], 3)
        self.assertEqual(self.ions.ionsinAmber['Ce4+'][1], 4)

    def test_resnames_are_strings(self):
        for name, (resname, charge) in self.ions.ionsinAmber.items():
            self.assertIsInstance(resname, str)
            self.assertIsInstance(charge, int)


if __name__ == '__main__':
    unittest.main()
