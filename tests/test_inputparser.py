import os
import shutil
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
        self.coords_ref = np.array([[-2.87230799,  0.14721348,  0.17546035],
       [-0.88989958,  0.63642335,  0.53637083],
       [-0.77026521,  1.5453448 ,  1.0103898 ],
       [-0.31174366,  0.65752973, -0.31832161],
       [-0.49249057, -0.08959617,  1.1578487 ],
       [-3.19190674,  1.80661442, -1.02584679],
       [-3.05341969,  2.70203475, -0.5316748 ],
       [-4.17512018,  1.78134192, -1.34812185],
       [-2.59336275,  1.82503979, -1.86635777],
       [-2.48097143, -1.71899766,  1.5396982 ],
       [-5.11551222, -0.38194842, -0.25344504]])
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
        self.Dmat_ref = np.array([[0.        , 2.07352979, 2.65903408, 2.65720398, 2.58547803,
        2.07337736, 2.65705727, 2.58630974, 2.65742974, 2.34457528,
        2.34434156],
       [2.07352979, 0.        , 1.03205841, 1.03208959, 1.0350232 ,
        3.01812989, 3.17620705, 3.95661795, 3.17611406, 3.01432971,
        4.41777007],
       [2.65903408, 1.03205841, 0.        , 1.66250764, 1.66491287,
        3.17472345, 2.98809124, 4.14865149, 3.41724897, 3.72325326,
        4.9186288 ],
       [2.65720398, 1.03208959, 1.66250764, 0.        , 1.66431528,
        3.18061737, 3.4267051 , 4.15320605, 2.99420797, 3.71559809,
        4.91537549],
       [2.58547803, 1.0350232 , 1.66491287, 1.66431528, 0.        ,
        3.956133  , 4.1480177 , 4.83136198, 4.15033969, 2.59900257,
        4.84247344],
       [2.07337736, 3.01812989, 3.17472345, 3.18061737, 3.956133  ,
        0.        , 1.03206695, 1.03499206, 1.0320141 , 4.4178491 ,
        3.01441036],
       [2.65705727, 3.17620705, 2.98809124, 3.4267051 , 4.1480177 ,
        1.03206695, 0.        , 1.66507449, 1.66197208, 4.91566992,
        3.72029429],
       [2.58630974, 3.95661795, 4.14865149, 4.15320605, 4.83136198,
        1.03499206, 1.66507449, 0.        , 1.66506292, 4.84376112,
        2.6004768 ],
       [2.65742974, 3.17611406, 3.41724897, 2.99420797, 4.15033969,
        1.0320141 , 1.66197208, 1.66506292, 0.        , 4.91671135,
        3.71934434],
       [2.34457528, 3.01432971, 3.72325326, 3.71559809, 2.59900257,
        4.4178491 , 4.91566992, 4.84376112, 4.91671135, 0.        ,
        3.45599024],
       [2.34434156, 4.41777007, 4.9186288 , 4.91537549, 4.84247344,
        3.01441036, 3.72029429, 2.6004768 , 3.71934434, 3.45599024,
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
        self.files_ref = [['HETATM    1 Pt     A 1   1      -2.872   0.147   0.175  1.00  0.00          Pt\n', 'TER\n'], ['HETATM    1 N      B 2   2       -0.89   0.636   0.536  1.00  0.00           N\n', 'HETATM    2 H      B 2   2       -0.77   1.545    1.01  1.00  0.00           H\n', 'HETATM    3 H      B 2   2      -0.312   0.658  -0.318  1.00  0.00           H\n', 'HETATM    4 H      B 2   2      -0.492   -0.09   1.158  1.00  0.00           H\n', 'TER\n'], ['HETATM    1 N      C 3   3      -3.192   1.807  -1.026  1.00  0.00           N\n', 'HETATM    2 H      C 3   3      -3.053   2.702  -0.532  1.00  0.00           H\n', 'HETATM    3 H      C 3   3      -4.175   1.781  -1.348  1.00  0.00           H\n', 'HETATM    4 H      C 3   3      -2.593   1.825  -1.866  1.00  0.00           H\n', 'TER\n'], ['HETATM    1 Cl     D 4   4      -2.481  -1.719    1.54  1.00  0.00          Cl\n', 'TER\n'], ['HETATM    1 Cl     E 5   5      -5.116  -0.382  -0.253  1.00  0.00          Cl\n', 'TER\n']]
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
        self.assertTrue(np.allclose(self.coords_ref,self.xyz.coords),'Coords arrays should be equal')

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
        sum = 0
        if filecmp.cmp(self.testFilePath + 'B.pdb',
                    self.testRefPath + '/MCPB_setup/B.pdb'):
            sum += 1
        if filecmp.cmp(self.testFilePath + 'C.pdb',
                    self.testRefPath + '/MCPB_setup/C.pdb'):

            sum += 1
        if filecmp.cmp(self.testFilePath + 'D.pdb',
                    self.testRefPath + '/MCPB_setup/D.pdb'):

            sum += 1
        if filecmp.cmp(self.testFilePath + 'E.pdb',
                    self.testRefPath + '/MCPB_setup/E.pdb'):

            sum += 1
        if filecmp.cmp(self.testFilePath + 'PT.pdb',
                    self.testRefPath + '/MCPB_setup/PT.pdb'):

            sum += 1
        if filecmp.cmp(self.testFilePath + 'Full_PDB.pdb',
                    self.testRefPath + '/MCPB_setup/Full_PDB.pdb'):

            sum += 1
        os.remove(self.testFilePath + 'B.pdb')
        os.remove(self.testFilePath + 'C.pdb')
        os.remove(self.testFilePath + 'D.pdb')
        os.remove(self.testFilePath + 'E.pdb')
        os.remove(self.testFilePath + 'PT.pdb')
        os.remove(self.testFilePath + 'Full_PDB.pdb')
        os.remove(self.testFilePath + 'filenames.restart')
        self.assertEqual(sum, 6, 'PDB files should be identical to references')

    def test_writeMol2Files(self):
        shutil.copytree(self.testRefPath + '/MCPB_setup/', self.testFilePath + '/MCPB_setup/')
        self.xyz.path = self.testFilePath + '/MCPB_setup/'
        self.xyz.writeMol2Files()
        sum = 0
        if filecmp.cmp(self.testFilePath + '/MCPB_setup//PT.mol2',
                    self.testRefPath + '/MCPB_setup/PT.mol2'):
            sum += 1
        if filecmp.cmp(self.testFilePath + '/MCPB_setup//B.mol2',
                    self.testRefPath + '/MCPB_setup/B.mol2'):
            sum += 1
        if filecmp.cmp(self.testFilePath + '/MCPB_setup//C.mol2',
                    self.testRefPath + '/MCPB_setup/C.mol2'):
            sum += 1
        if filecmp.cmp(self.testFilePath + '/MCPB_setup//D.mol2',
                    self.testRefPath + '/MCPB_setup/D.mol2'):
            sum += 1
        if filecmp.cmp(self.testFilePath + '/MCPB_setup//E.mol2',
                    self.testRefPath + '/MCPB_setup/E.mol2'):
            sum += 1
        shutil.rmtree(self.testFilePath + '/MCPB_setup/')

        self.assertEqual(sum, 5, 'Mol2 files should be identical to references')

if __name__ == '__main__':
    unittest.main()