import os
import shutil
import unittest
import filecmp
from ..src.PyConSolv.interfaces.amber import amberInterface


class TestSetup(unittest.TestCase):
    def setUp(self):
        self.testRefPath = os.path.split(__file__)[0] + '/TestReference/'
        self.testFilePath = os.path.split(__file__)[0] + '/Testfiles/'
        db_file = os.path.split(__file__)[0] + '/../src/PyConSolv/db/atom-radius.txt'
        db_metal_file = os.path.split(__file__)[0] + '/../src/PyConSolv/db/metal-radius.txt'
        # run functions of inputparser
        self.amber = amberInterface(self.testFilePath)
        self.metalbonds_ref = ['0-1 @N1', '0-5 @N5', '0-9 @Cl9', '0-10 @Cl10']
        self.connections_ref = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

    def test_inputFileGenerator(self):
        metals = 'PT'
        ligands = ['B','C','D','E']
        self.amber.inputFileGenerator(metals,ligands)
        result = filecmp.cmp(self.testFilePath + 'input.in',
                    self.testRefPath + '/MCPB_setup/input_reduced.in')
        os.remove(self.testFilePath + 'input.in')
        self.assertTrue(result)

    def test_readMetalBonds(self):
        metalbonds = self.amber.readMetalBonds(self.testRefPath + '/MCPB_setup/')
        self.assertListEqual(self.metalbonds_ref, metalbonds, 'Metalbonds should be the same')

    def test_readConnections(self):
        connections = self.amber.readConnections(self.testRefPath + '/MCPB_setup/')
        self.assertListEqual(self.connections_ref, connections, 'Connections should be the same')

    def test_checkMCPBBonds(self):
        shutil.copytree(self.testRefPath + '/MCPB_setup/', self.testFilePath + '/MCPB_setup/')
        shutil.copyfile(self.testRefPath + '/MCPB_setup/input_reduced.in', self.testFilePath + '/MCPB_setup/input.in')
        result = self.amber.checkMCPBBonds(self.testFilePath + '/MCPB_setup/')
        result =  filecmp.cmp(self.testFilePath + '/MCPB_setup/input.in',self.testRefPath + '/MCPB_setup/input.in')
        shutil.rmtree(self.testFilePath + '/MCPB_setup/')
        self.assertTrue(result)

if __name__ == '__main__':
    unittest.main()