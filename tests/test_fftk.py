import os
import unittest

import numpy as np

import tests  # noqa: F401
from PyConSolv.interfaces.fftk import (
    parseOrcaHessian, seminarioBond, seminarioAngle, metalAtomType,
    loadMetalNonbonded, FFTKBondedInterface,
    HA_BOHR2_TO_KCAL_ANG2, BOHR_TO_ANG,
)
from tests.helpers import TempDir, write


def diatomicHessian(unit_bond_vec, k_hartree_bohr2):
    """Exact Hessian for two atoms connected by a harmonic bond.

    E = k (r - r0)^2, so H_AB = -2 k u u^T (cross-block), and
    H_AA = H_BB = +2 k u u^T.
    """
    u = np.asarray(unit_bond_vec, dtype=float)
    u = u / np.linalg.norm(u)
    outer = 2.0 * k_hartree_bohr2 * np.outer(u, u)
    H = np.zeros((6, 6))
    H[0:3, 0:3] = outer
    H[3:6, 3:6] = outer
    H[0:3, 3:6] = -outer
    H[3:6, 0:3] = -outer
    return H


SAMPLE_HESS_FILE = """$orca_hessian_file

$hessian
3
                  0          1          2
      0    1.000000     0.200000     0.300000
      1    0.200000     2.000000     0.400000
      2    0.300000     0.400000     3.000000

$end
"""


class TestOrcaHessianParser(unittest.TestCase):
    def test_parseOrcaHessian_reads_square_matrix(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'orca.hess'), SAMPLE_HESS_FILE)
            H = parseOrcaHessian(path)
            self.assertEqual(H.shape, (3, 3))
            self.assertAlmostEqual(H[0, 0], 1.0)
            self.assertAlmostEqual(H[1, 2], 0.4)

    def test_parseOrcaHessian_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            parseOrcaHessian('/nonexistent/path/orca.hess')


class TestSeminarioBond(unittest.TestCase):
    def test_diatomic_recovers_force_constant(self):
        """For E = Kb (r-r0)^2, Seminario should return exactly Kb."""
        k_true_hartree_bohr2 = 0.5  # arbitrary
        expected_Kb = 0.5 * k_true_hartree_bohr2 * HA_BOHR2_TO_KCAL_ANG2 * 2
        # Note: CHARMM convention has no 1/2; our seminarioBond divides by 2.
        # Expected Kb(CHARMM) = k_true * HA_BOHR2_TO_KCAL_ANG2. But we built
        # the Hessian with H = 2k outer(u,u), matching E = k(r-r0)^2 with no 1/2,
        # so the Seminario projection returns 2k; divide by 2 -> k.
        expected_Kb = k_true_hartree_bohr2 * HA_BOHR2_TO_KCAL_ANG2

        u = np.array([1.0, 0.0, 0.0])
        H = diatomicHessian(u, k_true_hartree_bohr2)
        coords_ang = np.array([[0.0, 0.0, 0.0],
                                [1.5, 0.0, 0.0]])
        Kb, b0 = seminarioBond(H, coords_ang, 0, 1)
        self.assertAlmostEqual(b0, 1.5)
        self.assertAlmostEqual(Kb, expected_Kb, places=3)

    def test_bond_along_arbitrary_axis(self):
        k_true = 0.8
        u = np.array([1.0, 1.0, 1.0])
        u = u / np.linalg.norm(u)
        H = diatomicHessian(u, k_true)
        r0 = 2.0
        coords = np.array([[0.0, 0.0, 0.0], 2.0 * u])
        Kb, b0 = seminarioBond(H, coords, 0, 1)
        self.assertAlmostEqual(b0, r0, places=6)
        self.assertAlmostEqual(Kb, k_true * HA_BOHR2_TO_KCAL_ANG2, places=3)


class TestSeminarioAngle(unittest.TestCase):
    def test_angle_geometry(self):
        """Triangle with 90-degree apex; check theta0 is 90 and Ktheta >= 0."""
        coords = np.array([[1.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0],   # apex
                            [0.0, 1.0, 0.0]])
        # Hessian of zero yields Ktheta=0 but correct theta0.
        H = np.zeros((9, 9))
        Kth, th0 = seminarioAngle(H, coords, 0, 1, 2)
        self.assertAlmostEqual(th0, 90.0, places=4)
        self.assertEqual(Kth, 0.0)

    def test_linear_angle_returns_zero(self):
        coords = np.array([[0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                            [2.0, 0.0, 0.0]])
        H = np.zeros((9, 9))
        Kth, th0 = seminarioAngle(H, coords, 0, 1, 2)
        self.assertEqual(Kth, 0.0)
        self.assertAlmostEqual(th0, 180.0, places=4)


class TestMetalAtomType(unittest.TestCase):
    def test_two_letter_element(self):
        self.assertEqual(metalAtomType('Fe'), 'MFE')
        self.assertEqual(metalAtomType('cu'), 'MCU')

    def test_single_letter_element(self):
        self.assertEqual(metalAtomType('W'), 'MW')

    def test_truncates_to_four(self):
        self.assertEqual(metalAtomType('Abcd'), 'MABC')


class TestMetalNonbondedDB(unittest.TestCase):
    def test_parses_db_file(self):
        db_path = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                                'src', 'PyConSolv', 'db',
                                'charmm_metal_nonbonded.txt')
        table = loadMetalNonbonded(db_path)
        self.assertIn('FE', table)
        self.assertIn('ZN', table)
        eps, rmin2 = table['FE']
        self.assertLess(eps, 0.0)
        self.assertGreater(rmin2, 0.0)

    def test_missing_file_returns_empty(self):
        self.assertEqual(loadMetalNonbonded('/nonexistent'), {})


class TestFFTKBondedInterface(unittest.TestCase):
    def test_scope_selection_metal_forces_refit(self):
        with TempDir() as d:
            iface = FFTKBondedInterface(d)
            self.assertTrue(iface._shouldRefitBond(0, 1, {0}, {}))
            self.assertFalse(iface._shouldRefitBond(2, 3, {0}, {}))

    def test_scope_selection_penalty(self):
        with TempDir() as d:
            iface = FFTKBondedInterface(d, penalty_threshold=50.0)
            self.assertTrue(iface._shouldRefitBond(
                2, 3, set(), {2: 75.0, 3: 5.0}))
            self.assertFalse(iface._shouldRefitBond(
                2, 3, set(), {2: 5.0, 3: 5.0}))

    def test_deriveParameters_emits_rtf_and_prm(self):
        with TempDir() as d:
            iface = FFTKBondedInterface(d)
            coords = np.array([[0.0, 0.0, 0.0],
                                [2.0, 0.0, 0.0]])
            # Synthetic Hessian from diatomic with k=0.5.
            u = np.array([1.0, 0.0, 0.0])
            H = diatomicHessian(u, 0.5)
            rtf, prm = iface.deriveParameters(
                xyz_coords_ang=coords,
                elements=['Fe', 'N'],
                hessian=H,
                bonds=[(0, 1)],
                angles=[],
                metal_indices=[0],
                atom_types=['MFE', 'NG311'],
                penalties={},
                out_basename='testrun',
            )
            self.assertTrue(os.path.isfile(rtf))
            self.assertTrue(os.path.isfile(prm))
            with open(rtf) as f:
                rtf_txt = f.read()
            self.assertIn('MASS', rtf_txt)
            self.assertIn('MFE', rtf_txt)
            with open(prm) as f:
                prm_txt = f.read()
            self.assertIn('BONDS', prm_txt)
            # The one bond we asked about should appear with its types.
            self.assertIn('MFE', prm_txt)
            self.assertIn('NG311', prm_txt)
            self.assertIn('NONBONDED', prm_txt)

    def test_deriveParameters_skips_bonds_outside_scope(self):
        with TempDir() as d:
            iface = FFTKBondedInterface(d, penalty_threshold=50.0)
            coords = np.array([[0.0, 0.0, 0.0],
                                [2.0, 0.0, 0.0],
                                [4.0, 0.0, 0.0]])
            H = np.zeros((9, 9))
            rtf, prm = iface.deriveParameters(
                xyz_coords_ang=coords,
                elements=['C', 'C', 'C'],
                hessian=H,
                bonds=[(0, 1), (1, 2)],
                angles=[],
                metal_indices=[],
                atom_types=['CG331', 'CG331', 'CG331'],
                penalties={},
                out_basename='nometal',
            )
            with open(prm) as f:
                prm_txt = f.read()
            # No metals, no penalties -> BONDS section omitted.
            self.assertNotIn('BONDS', prm_txt)


if __name__ == '__main__':
    unittest.main()
