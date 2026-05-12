import os
import unittest

import numpy as np

import tests  # noqa: F401
from PyConSolv.misc.inputparser import XYZ
from tests.helpers import TempDir, radii_files, write, SAMPLE_XYZ, SAMPLE_METAL_XYZ


class TestXYZIsMetal(unittest.TestCase):
    def setUp(self):
        db, metal = radii_files()
        self.xyz = XYZ(db, metal)

    def test_metal_detection_common(self):
        for m in ('Fe', 'FE', 'Pt', 'Cu', 'Zn', 'ru'):
            self.assertTrue(self.xyz.isMetal(m), '{} should be metal'.format(m))

    def test_non_metals(self):
        for e in ('C', 'H', 'N', 'O', 'S', 'P', 'Cl', 'F'):
            self.assertFalse(self.xyz.isMetal(e), '{} should not be metal'.format(e))


class TestXYZReadAndBonds(unittest.TestCase):
    def setUp(self):
        db, metal = radii_files()
        self.tmp = TempDir()
        self.d = self.tmp.__enter__()
        self.water = write(os.path.join(self.d, 'water.xyz'), SAMPLE_XYZ)
        self.metal = write(os.path.join(self.d, 'fe.xyz'), SAMPLE_METAL_XYZ)
        self.xyz = XYZ(db, metal)

    def tearDown(self):
        self.tmp.__exit__(None, None, None)

    def test_readXYZ_shapes(self):
        self.xyz.readXYZ(self.water)
        self.assertEqual(self.xyz.atoms.shape, (3,))
        self.assertEqual(self.xyz.coords.shape, (3, 3))
        self.assertEqual(self.xyz.atoms[0], 'O')

    def test_distance_matrix_symmetric(self):
        self.xyz.readXYZ(self.water)
        self.xyz.calculateDistanceMatrix()
        self.assertTrue(np.allclose(self.xyz.Dmat, self.xyz.Dmat.T))
        self.assertAlmostEqual(float(self.xyz.Dmat[0][0]), 0.0)

    def test_adjacency_diagonal_zero(self):
        self.xyz.readXYZ(self.water)
        self.xyz.calculateDistanceMatrix()
        self.xyz.generateAdjacencyMatrix()
        for i in range(self.xyz.Adjmat.shape[0]):
            self.assertEqual(self.xyz.Adjmat[i][i], 0.0)

    def test_link_list_length_matches_atoms(self):
        self.xyz.readXYZ(self.water)
        self.xyz.calculateDistanceMatrix()
        self.xyz.generateAdjacencyMatrix()
        self.xyz.generateLinkList()
        self.assertEqual(len(self.xyz.linkList), len(self.xyz.atoms))

    def test_metal_bonds_populated(self):
        self.xyz.readXYZ(self.metal)
        self.xyz.calculateDistanceMatrix()
        self.xyz.generateAdjacencyMatrix()
        # Fe-N bonds at 2 Å should be detected as metal bonds.
        # hasMetal is only set by writePDBFiles; check metalBonds directly.
        self.assertGreater(len(self.xyz.metalBonds), 0)

    def test_connected_components_single_water(self):
        self.xyz.readXYZ(self.water)
        self.xyz.calculateDistanceMatrix()
        self.xyz.generateAdjacencyMatrix()
        self.xyz.generateLinkList()
        self.xyz.connectedCompponents()
        # 3 atoms, all connected via O
        flat = [a for comp in self.xyz.connected for a in comp]
        self.assertEqual(sorted(flat), [0, 1, 2])


if __name__ == '__main__':
    unittest.main()
