import os
import unittest

import numpy as np

import tests  # noqa: F401
from PyConSolv.misc.polar_sites import (
    PolarSite, Probe,
    detectPolarSites, placeWater, placeWaters,
    writeProbePDB, writeAllProbes,
    TIP3P_OH, TIP3P_HOH_DEG, DONOR_O_DIST, ACCEPTOR_H_DIST,
    _unit, _perpendicular, _rotate,
)
from tests.helpers import TempDir


def methanol():
    coords = np.array([
        [0.000,  0.000,  0.000],   # C
        [1.430,  0.000,  0.000],   # O
        [1.760,  0.945,  0.000],   # H on O
        [-0.355,  1.010, 0.000],   # H on C
        [-0.355, -0.505, 0.875],   # H on C
        [-0.355, -0.505, -0.875],  # H on C
    ])
    elements = ['C', 'O', 'H', 'H', 'H', 'H']
    bonds = [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5)]
    types = ['CG331', 'OG311', 'HGP1', 'HGA3', 'HGA3', 'HGA3']
    return coords, elements, bonds, types


def methylamine():
    coords = np.array([
        [0.000,  0.000,  0.000],   # C
        [1.470,  0.000,  0.000],   # N
        [1.825,  0.935,  0.000],   # H on N
        [1.825, -0.470,  0.810],   # H on N
        [-0.355,  1.010, 0.000],   # H on C
        [-0.355, -0.505, 0.875],   # H on C
        [-0.355, -0.505, -0.875],  # H on C
    ])
    elements = ['C', 'N', 'H', 'H', 'H', 'H', 'H']
    bonds = [(0, 1), (1, 2), (1, 3), (0, 4), (0, 5), (0, 6)]
    types = ['CG331', 'NG321', 'HGPAM2', 'HGPAM2', 'HGA3', 'HGA3', 'HGA3']
    return coords, elements, bonds, types


def acetone():
    coords = np.array([
        [0.000,  1.220,  0.000],   # O
        [0.000,  0.000,  0.000],   # C (carbonyl)
        [1.290, -0.745,  0.000],   # C methyl 1
        [-1.290, -0.745, 0.000],   # C methyl 2
        [1.290, -1.395,  0.890],
        [1.290, -1.395, -0.890],
        [2.180, -0.110,  0.000],
        [-1.290, -1.395, 0.890],
        [-1.290, -1.395, -0.890],
        [-2.180, -0.110, 0.000],
    ])
    elements = ['O', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
    bonds = [(0, 1), (1, 2), (1, 3),
             (2, 4), (2, 5), (2, 6),
             (3, 7), (3, 8), (3, 9)]
    types = ['OG2D3', 'CG2O5', 'CG331', 'CG331',
             'HGA3', 'HGA3', 'HGA3', 'HGA3', 'HGA3', 'HGA3']
    return coords, elements, bonds, types


def pyridine():
    angles = np.linspace(0, 2 * np.pi, 7)[:6]
    radius = 1.39
    ring = np.array([[radius * np.cos(a), radius * np.sin(a), 0.0]
                     for a in angles])
    elements = ['N', 'C', 'C', 'C', 'C', 'C']
    h_coords = []
    for i in range(1, 6):
        out_dir = ring[i] / np.linalg.norm(ring[i])
        h_coords.append(ring[i] + 1.08 * out_dir)
    coords = np.vstack([ring, np.array(h_coords)])
    elements += ['H'] * 5
    bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
             (1, 6), (2, 7), (3, 8), (4, 9), (5, 10)]
    types = ['NG2R60'] + ['CG2R61'] * 5 + ['HGR61'] * 5
    return coords, elements, bonds, types


def formamide():
    """HC(=O)NH2 - a sp2 N donor *and* a sp2 O acceptor."""
    coords = np.array([
        [0.000,  0.000,  0.000],   # C
        [0.000,  1.220,  0.000],   # O
        [1.190, -0.690,  0.000],   # N
        [-0.940, -0.490, 0.000],   # H on C
        [1.190, -1.695,  0.000],   # H on N
        [2.060, -0.190,  0.000],   # H on N
    ])
    elements = ['C', 'O', 'N', 'H', 'H', 'H']
    bonds = [(0, 1), (0, 2), (0, 3), (2, 4), (2, 5)]
    types = ['CG2O1', 'OG2D1', 'NG2S1', 'HGR52', 'HGP1', 'HGP1']
    return coords, elements, bonds, types


class TestPolarSiteDetection(unittest.TestCase):
    def test_methanol(self):
        coords, elements, bonds, types = methanol()
        sites = detectPolarSites(coords, elements, bonds, types)
        donors = [s for s in sites if s.role == 'donor']
        acceptors = [s for s in sites if s.role == 'acceptor']
        self.assertEqual(len(donors), 1)
        self.assertEqual(donors[0].heavy_idx, 1)
        self.assertEqual(donors[0].h_idx, 2)
        self.assertEqual(len(acceptors), 2)
        self.assertTrue(all(a.heavy_idx == 1 for a in acceptors))

    def test_methylamine(self):
        coords, elements, bonds, types = methylamine()
        sites = detectPolarSites(coords, elements, bonds, types)
        donors = [s for s in sites if s.role == 'donor']
        acceptors = [s for s in sites if s.role == 'acceptor']
        self.assertEqual(len(donors), 2)
        self.assertEqual(len(acceptors), 1)
        self.assertEqual(acceptors[0].heavy_idx, 1)

    def test_acetone(self):
        coords, elements, bonds, types = acetone()
        sites = detectPolarSites(coords, elements, bonds, types)
        donors = [s for s in sites if s.role == 'donor']
        acceptors = [s for s in sites if s.role == 'acceptor']
        self.assertEqual(len(donors), 0)
        self.assertEqual(len(acceptors), 2)
        for acc in acceptors:
            self.assertAlmostEqual(abs(acc.direction[2]), 0.0, places=4)

    def test_pyridine(self):
        coords, elements, bonds, types = pyridine()
        sites = detectPolarSites(coords, elements, bonds, types)
        donors = [s for s in sites if s.role == 'donor']
        acceptors = [s for s in sites if s.role == 'acceptor']
        self.assertEqual(len(donors), 0)
        self.assertEqual(len(acceptors), 1)
        self.assertAlmostEqual(abs(acceptors[0].direction[2]), 0.0, places=4)
        # outward from ring centre
        self.assertGreater(acceptors[0].direction[0], 0.99)

    def test_formamide(self):
        coords, elements, bonds, types = formamide()
        sites = detectPolarSites(coords, elements, bonds, types)
        donors = [s for s in sites if s.role == 'donor']
        acceptors_O = [s for s in sites
                        if s.role == 'acceptor' and s.heavy_idx == 1]
        # Two N-H donors on the amide nitrogen
        self.assertEqual(len(donors), 2)
        # sp2 carbonyl O -> two in-plane LPs
        self.assertEqual(len(acceptors_O), 2)

    def test_fallback_without_types(self):
        coords, elements, bonds, _ = methanol()
        sites = detectPolarSites(coords, elements, bonds)
        donors = [s for s in sites if s.role == 'donor']
        acceptors = [s for s in sites if s.role == 'acceptor']
        self.assertEqual(len(donors), 1)
        self.assertEqual(len(acceptors), 2)

    def test_skips_nonpolar_atoms(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        sites = detectPolarSites(coords, ['C', 'C'], [(0, 1)], ['CG331', 'CG331'])
        self.assertEqual(sites, [])


class TestWaterPlacement(unittest.TestCase):
    def test_donor_distance_and_geometry(self):
        coords, elements, bonds, types = methanol()
        sites = detectPolarSites(coords, elements, bonds, types)
        donor = next(s for s in sites if s.role == 'donor')
        water = placeWater(donor, coords)
        O, H1, H2 = water[0], water[1], water[2]
        H_pos = coords[donor.h_idx]
        self.assertAlmostEqual(np.linalg.norm(O - H_pos), DONOR_O_DIST,
                                places=4)
        self.assertAlmostEqual(np.linalg.norm(O - H1), TIP3P_OH, places=4)
        self.assertAlmostEqual(np.linalg.norm(O - H2), TIP3P_OH, places=4)
        v1 = _unit(H1 - O)
        v2 = _unit(H2 - O)
        angle = np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1, 1)))
        self.assertAlmostEqual(angle, TIP3P_HOH_DEG, places=3)

    def test_donor_water_O_lies_on_XH_extension(self):
        coords, elements, bonds, types = methanol()
        sites = detectPolarSites(coords, elements, bonds, types)
        donor = next(s for s in sites if s.role == 'donor')
        water = placeWater(donor, coords)
        O_w = water[0]
        H_pos = coords[donor.h_idx]
        X_pos = coords[donor.heavy_idx]
        cos_a = (np.dot(O_w - H_pos, H_pos - X_pos)
                 / (np.linalg.norm(O_w - H_pos)
                    * np.linalg.norm(H_pos - X_pos)))
        self.assertAlmostEqual(cos_a, 1.0, places=4)

    def test_acceptor_distance_and_geometry(self):
        coords, elements, bonds, types = acetone()
        sites = detectPolarSites(coords, elements, bonds, types)
        acc = next(s for s in sites if s.role == 'acceptor')
        water = placeWater(acc, coords)
        O_w, H1, H2 = water[0], water[1], water[2]
        O_acc = coords[acc.heavy_idx]
        self.assertAlmostEqual(np.linalg.norm(H1 - O_acc), ACCEPTOR_H_DIST,
                                places=4)
        self.assertAlmostEqual(np.linalg.norm(O_w - H1), TIP3P_OH, places=4)
        self.assertAlmostEqual(np.linalg.norm(O_w - H2), TIP3P_OH, places=4)
        v1 = _unit(H1 - O_w)
        v2 = _unit(H2 - O_w)
        angle = np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1, 1)))
        self.assertAlmostEqual(angle, TIP3P_HOH_DEG, places=3)

    def test_acceptor_water_H_lies_along_LP(self):
        coords, elements, bonds, types = acetone()
        sites = detectPolarSites(coords, elements, bonds, types)
        acc = next(s for s in sites if s.role == 'acceptor')
        water = placeWater(acc, coords)
        H1 = water[1]
        cos_a = (np.dot(H1 - coords[acc.heavy_idx], acc.direction)
                 / np.linalg.norm(H1 - coords[acc.heavy_idx]))
        self.assertAlmostEqual(cos_a, 1.0, places=4)

    def test_placeWaters_returns_one_probe_per_site(self):
        coords, elements, bonds, types = methanol()
        sites = detectPolarSites(coords, elements, bonds, types)
        probes = placeWaters(sites, coords)
        self.assertEqual(len(probes), len(sites))
        self.assertTrue(all(isinstance(p, Probe) for p in probes))
        self.assertTrue(all(p.water_coords.shape == (3, 3) for p in probes))


class TestPDBWriter(unittest.TestCase):
    def test_atom_count_matches_ligand_plus_water(self):
        with TempDir() as d:
            coords, elements, bonds, types = methanol()
            sites = detectPolarSites(coords, elements, bonds, types)
            probes = placeWaters(sites, coords)
            paths = writeAllProbes(probes, coords, elements, d)
            self.assertEqual(len(paths), len(sites))
            for path in paths:
                with open(path) as f:
                    text = f.read()
                atoms = [l for l in text.splitlines() if l.startswith('ATOM')]
                self.assertEqual(len(atoms), len(elements) + 3)
                self.assertIn('TIP3', text)
                self.assertIn('LIG', text)
                self.assertTrue(text.rstrip().endswith('END'))

    def test_creates_directory(self):
        with TempDir() as d:
            sub = os.path.join(d, 'probes')
            coords, elements, bonds, types = methanol()
            sites = detectPolarSites(coords, elements, bonds, types)
            probes = placeWaters(sites, coords)
            writeAllProbes(probes, coords, elements, sub)
            self.assertTrue(os.path.isdir(sub))


class TestGeometryHelpers(unittest.TestCase):
    def test_unit_normalises(self):
        v = _unit(np.array([3.0, 4.0, 0.0]))
        self.assertAlmostEqual(np.linalg.norm(v), 1.0)
        self.assertAlmostEqual(v[0], 0.6)
        self.assertAlmostEqual(v[1], 0.8)

    def test_unit_rejects_zero(self):
        with self.assertRaises(ValueError):
            _unit(np.zeros(3))

    def test_perpendicular_is_orthogonal(self):
        for v in (np.array([1.0, 2.0, 3.0]),
                  np.array([1.0, 0.0, 0.0]),
                  np.array([0.0, 1.0, 0.0]),
                  np.array([0.0, 0.0, 1.0])):
            p = _perpendicular(v)
            self.assertAlmostEqual(np.dot(p, _unit(v)), 0.0, places=6)
            self.assertAlmostEqual(np.linalg.norm(p), 1.0, places=6)

    def test_rotate_90_around_z(self):
        v = np.array([1.0, 0.0, 0.0])
        out = _rotate(v, np.array([0.0, 0.0, 1.0]), np.pi / 2)
        self.assertAlmostEqual(out[0], 0.0, places=6)
        self.assertAlmostEqual(out[1], 1.0, places=6)
        self.assertAlmostEqual(out[2], 0.0, places=6)

    def test_rotate_preserves_norm(self):
        v = np.array([1.3, -0.7, 0.5])
        out = _rotate(v, np.array([1.0, 1.0, 1.0]), 1.234)
        self.assertAlmostEqual(np.linalg.norm(out), np.linalg.norm(v),
                                places=6)


if __name__ == '__main__':
    unittest.main()