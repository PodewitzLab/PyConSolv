import os
import unittest
import filecmp
from ..src.PyConSolv.misc.filestructure import Setup


class TestSetup(unittest.TestCase):
    def setUp(self):
        self.testRefPath = os.path.split(__file__)[0] + '/TestReference/'
        self.testFilePath = os.path.split(__file__)[0] + '/Testfiles/'
        self.setup = Setup(self.testFilePath + '/input.xyz', charge=0)

    def test_createFolders(self):
        self.setup.createFolders()
        sum = 0
        if os.path.exists(self.testFilePath + '/equilibration'):
            sum += 1
            os.rmdir(self.testFilePath + '/equilibration')
        if os.path.exists(self.testFilePath + '/simulation'):
            sum += 1
            os.rmdir(self.testFilePath + '/simulation')
        if os.path.exists(self.testFilePath + '/MCPB_setup'):
            sum += 1
            os.rmdir(self.testFilePath + '/MCPB_setup')
        if os.path.exists(self.testFilePath + '/orca_calculations/opt'):
            sum += 1
            os.rmdir(self.testFilePath + '/orca_calculations/opt')
        if os.path.exists(self.testFilePath + '/orca_calculations/freq'):
            sum += 1
            os.rmdir(self.testFilePath + '/orca_calculations/freq')
        if os.path.exists(self.testFilePath + '/orca_calculations'):
            sum += 1
            os.rmdir(self.testFilePath + '/orca_calculations')
        self.assertEqual(sum,6,'There should be 6 checks completed for 5 directories')
    def test_createFiles(self):
        method = 'BP86'
        basis = 'def2-TZVP'
        dsp = 'D3'
        solvent = 'CH2Cl2'
        cpu = 8

        self.setup.Method(method = method, basis = basis, DSP = dsp, CPCM=solvent, CPU = cpu)
        os.mkdir(self.testFilePath + '/orca_calculations/')
        os.mkdir(self.testFilePath + '/orca_calculations/opt')
        os.mkdir(self.testFilePath + '/orca_calculations/freq')
        self.setup.createFiles()
        sum = 0
        if filecmp.cmp(self.testFilePath + '/orca_calculations/opt/orca_opt.inp',
                       self.testRefPath + '/orca_calculations/opt/orca_opt.inp'):
            sum += 1
            os.remove(self.testFilePath + '/orca_calculations/opt/orca_opt.inp')

        if filecmp.cmp(self.testFilePath + '/orca_calculations/freq/orca_freq.inp',
                       self.testRefPath + '/orca_calculations/freq/orca_freq.inp'):
            sum += 1
            os.remove(self.testFilePath + '/orca_calculations/freq/orca_freq.inp')

        if os.path.exists(self.testFilePath + '/orca_calculations/opt/input.xyz'):
            if filecmp.cmp(self.testFilePath + '/orca_calculations/opt/input.xyz',
                           self.testRefPath + '/orca_calculations/opt/input.xyz'):
                sum += 1
            os.remove(self.testFilePath + '/orca_calculations/opt/input.xyz')

        os.rmdir(self.testFilePath + '/orca_calculations/opt')
        os.rmdir(self.testFilePath + '/orca_calculations/freq')
        os.rmdir(self.testFilePath + '/orca_calculations/')

        self.assertEqual(sum, 3, 'Not all orca input files are identical')

if __name__ == '__main__':
    unittest.main()