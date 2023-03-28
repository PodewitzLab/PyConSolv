import os
import unittest
import filecmp
import numpy as np
from ..src.PyConSolv.misc.inputparser import XYZ


class TestSetup(unittest.TestCase):
    def setUp(self):
        self.testRefPath = os.path.split(__file__)[0] + '/TestReference/'
        self.testFilePath = os.path.split(__file__)[0] + '/Testfiles/'
        db_file = os.path.split(__file__)[0] + '/../src/PyConSolv/db/atom-radius.txt'
        db_metal_file = os.path.split(__file__)[0] + '/../src/PyConSolv/db/metal-radius.txt'
        #run functions of inputparser
        self.xyz = XYZ(db_file,db_metal_file)
        self.xyz.readXYZ(os.path.split(__file__)[0] + '/Testfiles/input.xyz')
        self.xyz.calculateDistanceMatrix()
        self.xyz.generateAdjacencyMatrix()
        self.xyz.generateLinkList()
        self.xyz.connectedCompponents()
        self.xyz.assignChain()
        self.xyz.createPDB()
        # checked reference values
        self.coords_ref = np.array([[-2.872, 0.148, 0.176],
                                    [-0.889, 0.636, 0.537],
                                    [-0.768, 1.543, 1.015],
                                    [-0.31, 0.66, -0.317],
                                    [-0.491, -0.092, 1.155],
                                    [-3.193, 1.807, -1.027],
                                    [-3.056, 2.703, -0.534],
                                    [-4.176, 1.781, -1.35],
                                    [-2.594, 1.826, -1.868],
                                    [-2.482, -1.719, 1.541],
                                    [-5.116, -0.382, -0.252]])
        self.atoms_ref = np.array(['Pt', 'N', 'H', 'H', 'H', 'N', 'H', 'H', 'H', 'Cl', 'Cl'], dtype='<U2')
        self.Amat_ref = np.array([['Pt-Pt', 'Pt-N', 'Pt-H', 'Pt-H', 'Pt-H', 'Pt-N', 'Pt-H', 'Pt-H',
        'Pt-H', 'Pt-Cl', 'Pt-Cl'],
       ['N-Pt', 'N-N', 'N-H', 'N-H', 'N-H', 'N-N', 'N-H', 'N-H', 'N-H',
        'N-Cl', 'N-Cl'],
       ['H-Pt', 'H-N', 'H-H', 'H-H', 'H-H', 'H-N', 'H-H', 'H-H', 'H-H',
        'H-Cl', 'H-Cl'],
       ['H-Pt', 'H-N', 'H-H', 'H-H', 'H-H', 'H-N', 'H-H', 'H-H', 'H-H',
        'H-Cl', 'H-Cl'],
       ['H-Pt', 'H-N', 'H-H', 'H-H', 'H-H', 'H-N', 'H-H', 'H-H', 'H-H',
        'H-Cl', 'H-Cl'],
       ['N-Pt', 'N-N', 'N-H', 'N-H', 'N-H', 'N-N', 'N-H', 'N-H', 'N-H',
        'N-Cl', 'N-Cl'],
       ['H-Pt', 'H-N', 'H-H', 'H-H', 'H-H', 'H-N', 'H-H', 'H-H', 'H-H',
        'H-Cl', 'H-Cl'],
       ['H-Pt', 'H-N', 'H-H', 'H-H', 'H-H', 'H-N', 'H-H', 'H-H', 'H-H',
        'H-Cl', 'H-Cl'],
       ['H-Pt', 'H-N', 'H-H', 'H-H', 'H-H', 'H-N', 'H-H', 'H-H', 'H-H',
        'H-Cl', 'H-Cl'],
       ['Cl-Pt', 'Cl-N', 'Cl-H', 'Cl-H', 'Cl-H', 'Cl-N', 'Cl-H', 'Cl-H',
        'Cl-H', 'Cl-Cl', 'Cl-Cl'],
       ['Cl-Pt', 'Cl-N', 'Cl-H', 'Cl-H', 'Cl-H', 'Cl-N', 'Cl-H', 'Cl-H',
        'Cl-H', 'Cl-Cl', 'Cl-Cl']], dtype='<U5')
        self.Dmat_ref = np.array([[0.        , 2.07382593, 2.66021841, 2.65876607, 2.58557576,
        2.07425432, 2.6581913 , 2.58762072, 2.65911715, 2.34542406,
        2.34512686],
       [2.07382593, 0.        , 1.03236331, 1.03205281, 1.03455884,
        3.02088613, 3.18047465, 3.95931345, 3.17917442, 3.01524294,
        4.41886569],
       [2.66021841, 1.03236331, 0.        , 1.66243105, 1.6641977 ,
        3.1812081 , 2.99665564, 4.15503706, 3.42433264, 3.7222461 ,
        4.92097734],
       [2.65876607, 1.03205281, 1.66243105, 0.        , 1.66284365,
        3.18298571, 3.42949763, 4.15568117, 2.9969673 , 3.71878865,
        4.9180916 ],
       [2.58557576, 1.03455884, 1.6641977 , 1.66284365, 0.        ,
        3.95829875, 4.15258606, 4.83346449, 4.15209128, 2.60003962,
        4.84297161],
       [2.07425432, 3.02088613, 3.1812081 , 3.18298571, 3.95829875,
        0.        , 1.03181103, 1.03503333, 1.03268727, 4.41959512,
        3.01500829],
       [2.6581913 , 3.18047465, 2.99665564, 3.42949763, 4.15258606,
        1.03181103, 0.        , 1.66443384, 1.6619654 , 4.9182502 ,
        3.72026195],
       [2.58762072, 3.95931345, 4.15503706, 4.15568117, 4.83346449,
        1.03503333, 1.66443384, 0.        , 1.66525464, 4.84536036,
        2.60149438],
       [2.65911715, 3.17917442, 3.42433264, 2.9969673 , 4.15209128,
        1.03268727, 1.6619654 , 1.66525464, 0.        , 4.91943594,
        3.72118315],
       [2.34542406, 3.01524294, 3.7222461 , 3.71878865, 2.60003962,
        4.41959512, 4.9182502 , 4.84536036, 4.91943594, 0.        ,
        3.45548463],
       [2.34512686, 4.41886569, 4.92097734, 4.9180916 , 4.84297161,
        3.01500829, 3.72026195, 2.60149438, 3.72118315, 3.45548463,
        0.        ]])

        self.Adjmat_ref = np.array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 1., 1., 1., 0., 0., 0., 0., 0., 0.],
           [0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 1., 1., 1., 0., 0.],
           [0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]])

        self.linkList_ref = [[], [2, 3, 4], [1], [1], [1], [6, 7, 8], [5], [5], [5], [], []]
        self.connected_ref = [[0], [1, 2, 3, 4], [5, 6, 7, 8], [9], [10]]
        self.files_ref = [['HETATM    1 Pt     A 1   1      -2.872   0.148   0.176  1.00  0.00          Pt\n', 'TER\n'],
                          ['HETATM    1 N      B 2   2      -0.889   0.636   0.537  1.00  0.00           N\n',
                           'HETATM    2 H      B 2   2      -0.768   1.543   1.015  1.00  0.00           H\n',
                           'HETATM    3 H      B 2   2       -0.31    0.66  -0.317  1.00  0.00           H\n',
                           'HETATM    4 H      B 2   2      -0.491  -0.092   1.155  1.00  0.00           H\n',
                           'TER\n'],
                          ['HETATM    1 N      C 3   3      -3.193   1.807  -1.027  1.00  0.00           N\n',
                           'HETATM    2 H      C 3   3      -3.056   2.703  -0.534  1.00  0.00           H\n',
                           'HETATM    3 H      C 3   3      -4.176   1.781   -1.35  1.00  0.00           H\n',
                           'HETATM    4 H      C 3   3      -2.594   1.826  -1.868  1.00  0.00           H\n',
                           'TER\n'],
                          ['HETATM    1 Cl     D 4   4      -2.482  -1.719   1.541  1.00  0.00          Cl\n',
                           'TER\n'],
                          ['HETATM    1 Cl     E 5   5      -5.116  -0.382  -0.252  1.00  0.00          Cl\n',
                           'TER\n']]
        self.chain_ref = np.array([[ 0,  0],
       [ 1,  1],
       [ 2,  1],
       [ 3,  1],
       [ 4,  1],
       [ 5,  2],
       [ 6,  2],
       [ 7,  2],
       [ 8,  2],
       [ 9,  3],
       [10,  4]])

    def test_readXYZ_coords(self):
        self.assertTrue(np.array_equal(self.coords_ref,self.xyz.coords),'Coords arrays should be equal')

    def test_readXYZ_atoms(self):
        self.assertTrue(np.array_equal(self.atoms_ref,self.xyz.atoms),'Atoms arrays should be equal')

    def test_calculateDistanceMatrix_Amat(self):
        self.assertTrue(np.array_equal(self.xyz.Amat,self.Amat_ref),'Amats should be equal')

    def test_calculateDistanceMatrix_Dmat(self):
        self.assertTrue(np.allclose(self.xyz.Dmat, self.Dmat_ref), 'Dmats should be within numerical error')

    def test_generateAdjacencyMatrix(self):
        self.assertTrue(np.array_equal(self.xyz.Adjmat, self.Adjmat_ref), 'Adjmat should be identical')

    def test_generateLinkList(self):
        self.assertListEqual(self.xyz.linkList,self.linkList_ref,'Linklists should be identical')

    def test_connectedComponents(self):
        self.assertListEqual(self.xyz.connected,self.connected_ref,'Connected should be identical')

    def test_assignChain(self):
            self.assertTrue(np.array_equal(self.xyz.chain, self.chain_ref), 'Chains should be identical')
    def test_createPDB(self):
        self.assertListEqual(self.xyz.files,self.files_ref,'files should be identical')

    def test_writePDBFiles(self):
        self.xyz.writePDBFiles(self.testFilePath)
        self.assertTrue(True)

    def test_writeMol2Files(self):
        self.xyz.path = self.testRefPath + '/MCPB_setup/'
        self.xyz.writeMol2Files()
        sum = 0
        if filecmp.cmp(self.testFilePath + '/PT.mol2',
                    self.testRefPath + '/MCPB_setup/PT.mol2'):
            sum += 1
        if filecmp.cmp(self.testFilePath + '/B.mol2',
                    self.testRefPath + '/MCPB_setup/B.mol2'):
            sum += 1
        if filecmp.cmp(self.testFilePath + '/C.mol2',
                    self.testRefPath + '/MCPB_setup/C.mol2'):
            sum += 1
        if filecmp.cmp(self.testFilePath + '/D.mol2',
                    self.testRefPath + '/MCPB_setup/D.mol2'):
            sum += 1
        if filecmp.cmp(self.testFilePath + '/E.mol2',
                    self.testRefPath + '/MCPB_setup/E.mol2'):
            sum += 1
        self.assertEqual(sum, 5, 'Mol2 files should be identical to references')

if __name__ == '__main__':
    unittest.main()