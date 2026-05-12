import unittest

import numpy as np

import tests  # noqa: F401
from PyConSolv.misc.charge_fitter import (
    fitCharges, mmCrossEnergy, FitResult,
    TIP3P_CHARGES, TIP3P_LJ, COULOMB_K,
    _buildExpansion, _initialReduced,
)
from PyConSolv.interfaces.water_interaction import InteractionResult


def _fakeResult(idx, ligand_coords, water_xyz, e_scaled):
    geometry = np.vstack([ligand_coords, water_xyz])
    return InteractionResult(
        site_label='synth_{}'.format(idx),
        probe_idx=idx,
        geometry=geometry,
        n_ligand_atoms=len(ligand_coords),
        e_complex=-100.0,
        e_ligand=-50.0,
        e_water=-50.0,
        e_int_kcal=e_scaled / 1.16,
        e_int_scaled=e_scaled,
        distance=2.0,
    )


METHANOL_COORDS = np.array([
    [0.000,  0.000,  0.000],   # C
    [1.430,  0.000,  0.000],   # O
    [1.760,  0.945,  0.000],   # H on O
    [-0.355,  1.010, 0.000],   # H on C
    [-0.355, -0.505, 0.875],   # H on C
    [-0.355, -0.505, -0.875],  # H on C
])
METHANOL_LJ = np.array([
    [-0.0780, 2.0500],   # C  (CG331-ish)
    [-0.1921, 1.7650],   # O  (OG311)
    [-0.0460, 0.2245],   # H (polar, HGP1)
    [-0.0240, 1.3400],   # H (aliphatic)
    [-0.0240, 1.3400],
    [-0.0240, 1.3400],
])


def _waterAt(O_xyz):
    """Build a roughly TIP3P-shaped water with O at the given position."""
    return np.array([O_xyz,
                     O_xyz + np.array([0.9572, 0.0, 0.0]),
                     O_xyz + np.array([-0.2400, 0.9266, 0.0])])


class TestMMCrossEnergy(unittest.TestCase):
    def test_pure_coulomb_two_charges(self):
        ca = np.array([[0.0, 0.0, 0.0]])
        cb = np.array([[1.0, 0.0, 0.0]])
        qa = np.array([1.0])
        qb = np.array([-1.0])
        zero_lj = np.array([[0.0, 0.0]])
        e = mmCrossEnergy(ca, qa, zero_lj, cb, qb, zero_lj)
        self.assertAlmostEqual(e, -COULOMB_K, places=4)

    def test_lj_at_minimum_equals_neg_eps(self):
        eps_a, rmin_a = 0.1, 1.0
        eps_b, rmin_b = 0.1, 1.0
        ca = np.array([[0.0, 0.0, 0.0]])
        cb = np.array([[rmin_a + rmin_b, 0.0, 0.0]])
        zero_q = np.array([0.0])
        lja = np.array([[-eps_a, rmin_a]])
        ljb = np.array([[-eps_b, rmin_b]])
        e = mmCrossEnergy(ca, zero_q, lja, cb, zero_q, ljb)
        self.assertAlmostEqual(e, -np.sqrt(eps_a * eps_b), places=6)

    def test_lj_uses_absolute_epsilon(self):
        """Repulsive distances should give positive LJ regardless of sign."""
        ca = np.array([[0.0, 0.0, 0.0]])
        cb = np.array([[0.5, 0.0, 0.0]])  # well inside Rmin
        zero_q = np.array([0.0])
        lja = np.array([[-0.1, 1.0]])
        ljb = np.array([[-0.1, 1.0]])
        e = mmCrossEnergy(ca, zero_q, lja, cb, zero_q, ljb)
        self.assertGreater(e, 0.0)


class TestBuildExpansion(unittest.TestCase):
    def test_no_groups_returns_identity(self):
        independent, sizes, expand = _buildExpansion(4, None)
        self.assertEqual(independent, [0, 1, 2, 3])
        self.assertEqual(len(sizes), 0)
        q = expand(np.array([0.1, 0.2, 0.3, 0.4]))
        np.testing.assert_allclose(q, [0.1, 0.2, 0.3, 0.4])

    def test_groups_share_charges(self):
        independent, sizes, expand = _buildExpansion(5, [[1, 2], [3, 4]])
        self.assertEqual(independent, [0])
        np.testing.assert_allclose(sizes, [2, 2])
        q = expand(np.array([0.5, -0.3, 0.1]))   # g0=0.5, g1=-0.3, indep=0.1
        np.testing.assert_allclose(q, [0.1, 0.5, 0.5, -0.3, -0.3])

    def test_overlapping_groups_raise(self):
        with self.assertRaises(ValueError):
            _buildExpansion(5, [[0, 1], [1, 2]])


class TestFitChargesRoundTrip(unittest.TestCase):
    def _make_results(self, q_true, lig_coords, lig_lj):
        # Probe waters dropped at three different positions around methanol O.
        water_positions = [
            np.array([3.6, 0.5, 0.0]),
            np.array([1.5, 2.5, 0.5]),
            np.array([0.0, -2.5, 0.5]),
        ]
        results = []
        for i, O_pos in enumerate(water_positions):
            wxyz = _waterAt(O_pos)
            e_mm = mmCrossEnergy(lig_coords, q_true, lig_lj,
                                  wxyz, TIP3P_CHARGES, TIP3P_LJ)
            results.append(_fakeResult(i, lig_coords, wxyz, e_mm))
        return results

    def test_recovers_known_charges_no_symmetry(self):
        q_true = np.array([0.10, -0.65, 0.43, 0.04, 0.04, 0.04])
        results = self._make_results(q_true, METHANOL_COORDS, METHANOL_LJ)
        q_init = np.array([0.0, -0.50, 0.30, 0.05, 0.05, 0.05])
        res = fitCharges(
            qm_results=results, ligand_coords=METHANOL_COORDS,
            ligand_elements=['C', 'O', 'H', 'H', 'H', 'H'],
            initial_charges=q_init, ligand_lj=METHANOL_LJ,
            total_charge=0.0, regularizer=0.0,
            symmetry_groups=[[3, 4, 5]],
        )
        self.assertTrue(res.success, res.message)
        np.testing.assert_allclose(res.charges, q_true, atol=5e-3)
        self.assertLess(res.rmsd, 1e-3)

    def test_total_charge_constraint(self):
        q_true = np.array([0.10, -0.65, 0.43, 0.04, 0.04, 0.04])
        results = self._make_results(q_true, METHANOL_COORDS, METHANOL_LJ)
        q_init = np.array([0.2, -0.5, 0.4, 0.1, 0.1, 0.1])
        res = fitCharges(
            qm_results=results, ligand_coords=METHANOL_COORDS,
            ligand_elements=['C', 'O', 'H', 'H', 'H', 'H'],
            initial_charges=q_init, ligand_lj=METHANOL_LJ,
            total_charge=-1.0, regularizer=0.0,
            symmetry_groups=[[3, 4, 5]],
        )
        self.assertAlmostEqual(float(np.sum(res.charges)), -1.0, places=4)

    def test_symmetry_groups_enforced(self):
        q_true = np.array([0.10, -0.65, 0.43, 0.04, 0.04, 0.04])
        results = self._make_results(q_true, METHANOL_COORDS, METHANOL_LJ)
        res = fitCharges(
            qm_results=results, ligand_coords=METHANOL_COORDS,
            ligand_elements=['C', 'O', 'H', 'H', 'H', 'H'],
            initial_charges=np.zeros(6), ligand_lj=METHANOL_LJ,
            total_charge=0.0, regularizer=0.0,
            symmetry_groups=[[3, 4, 5]],
        )
        self.assertAlmostEqual(res.charges[3], res.charges[4], places=8)
        self.assertAlmostEqual(res.charges[4], res.charges[5], places=8)


class TestRegulariser(unittest.TestCase):
    def test_pulls_toward_initial_when_data_is_weak(self):
        # Single tiny-weight probe -> regulariser dominates.
        q_init = np.array([0.0, -0.5, 0.3, 0.05, 0.05, 0.1])
        # Force the total charge to match q_init so the constraint is happy.
        target = float(np.sum(q_init))
        wxyz = _waterAt(np.array([3.6, 0.5, 0.0]))
        # Set a "QM target" that's wildly inconsistent with q_init.
        results = [_fakeResult(0, METHANOL_COORDS, wxyz, e_scaled=-50.0)]
        res = fitCharges(
            qm_results=results, ligand_coords=METHANOL_COORDS,
            ligand_elements=['C', 'O', 'H', 'H', 'H', 'H'],
            initial_charges=q_init, ligand_lj=METHANOL_LJ,
            total_charge=target,
            weights=[1e-6],         # almost-zero data weight
            regularizer=1.0,        # large regulariser
        )
        np.testing.assert_allclose(res.charges, q_init, atol=5e-3)


class TestBounds(unittest.TestCase):
    def test_per_atom_bounds_respected(self):
        q_true = np.array([0.10, -0.65, 0.43, 0.04, 0.04, 0.04])
        wxyz = _waterAt(np.array([3.6, 0.5, 0.0]))
        results = [_fakeResult(0, METHANOL_COORDS, wxyz,
                                e_scaled=mmCrossEnergy(METHANOL_COORDS,
                                                        q_true, METHANOL_LJ,
                                                        wxyz, TIP3P_CHARGES,
                                                        TIP3P_LJ))]
        res = fitCharges(
            qm_results=results, ligand_coords=METHANOL_COORDS,
            ligand_elements=['C', 'O', 'H', 'H', 'H', 'H'],
            initial_charges=np.zeros(6), ligand_lj=METHANOL_LJ,
            total_charge=0.0,
            charge_bounds=(-0.3, 0.3),     # tighter than q_true[1] = -0.65
            regularizer=0.0,
        )
        self.assertTrue(np.all(np.abs(res.charges) <= 0.3 + 1e-6))


if __name__ == '__main__':
    unittest.main()
