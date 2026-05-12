import os
import unittest
from unittest import mock

import tests  # noqa: F401
from PyConSolv.interfaces.amber import amberInterface
from tests.helpers import TempDir


class TestAmberInterface(unittest.TestCase):
    def test_input_file_generator_multiple_metals(self):
        with TempDir() as d:
            iface = amberInterface(d)
            iface.inputFileGenerator(['FE', 'CU'], ['B', 'C'])
            with open(os.path.join(d, 'input.in')) as f:
                text = f.read()
            self.assertIn('ion_ids 1 2', text)
            self.assertIn('FE.mol2 CU.mol2', text)
            self.assertIn('B.mol2', text)
            self.assertIn('C.frcmod', text)
            self.assertIn('group_name LIG', text)

    def test_input_file_generator_single_metal_string(self):
        with TempDir() as d:
            iface = amberInterface(d)
            iface.inputFileGenerator('PT', ['B'])
            with open(os.path.join(d, 'input.in')) as f:
                text = f.read()
            self.assertIn('ion_ids 1', text)
            self.assertIn('PT.mol2', text)

    def test_antechamber_success(self):
        with TempDir() as d:
            iface = amberInterface(d)
            with mock.patch('subprocess.run') as run:
                run.return_value = mock.Mock(returncode=0)
                status = iface.antechamber('LIG', 0)
            self.assertEqual(status, 1)
            cmd = run.call_args[0][0][0]
            self.assertIn('antechamber', cmd)
            self.assertIn('LIG.pdb', cmd)
            self.assertIn('LIG.mol2', cmd)

    def test_antechamber_failure(self):
        with TempDir() as d:
            iface = amberInterface(d)
            with mock.patch('subprocess.run') as run:
                run.return_value = mock.Mock(returncode=2)
                self.assertEqual(iface.antechamber('LIG', 0), 0)

    def test_runMCPB_builds_command(self):
        with TempDir() as d:
            iface = amberInterface(d)
            with mock.patch('subprocess.run') as run:
                run.return_value = mock.Mock(returncode=0)
                iface.runMCPB('1')
            cmd = run.call_args[0][0][0]
            self.assertIn('MCPB.py', cmd)
            self.assertIn('-s 1', cmd)

    def test_runParmchk2_builds_command(self):
        with TempDir() as d:
            iface = amberInterface(d)
            with mock.patch('subprocess.run') as run:
                run.return_value = mock.Mock(returncode=0)
                iface.runParmchk2(['ignored', 'LIG'])
            cmd = run.call_args[0][0][0]
            self.assertIn('parmchk2', cmd)
            self.assertIn('LIG.mol2', cmd)

    def test_tleap_no_metal_writes_file(self):
        with TempDir() as d:
            iface = amberInterface(d)
            iface.tleapNoMetal(d, name='LIG')
            with open(os.path.join(d, 'LIG_tleap.in')) as f:
                text = f.read()
            self.assertIn('source leaprc.gaff', text)
            self.assertIn('loadmol2', text)

    def test_tleap_no_metal_solv_writes_box_size(self):
        with TempDir() as d:
            iface = amberInterface(d)
            iface.defaultbox = 15
            iface.tleapNoMetalSolv(d, name='LIG')
            with open(os.path.join(d, 'LIG_tleap.in')) as f:
                text = f.read()
            self.assertIn('solvatebox LIG TIP3PBOX 15', text)

    def test_tleapChecker_removes_duplicates(self):
        with TempDir() as d:
            path = os.path.join(d, 'LIG_tleap.in')
            with open(path, 'w') as f:
                f.write('source leaprc.gaff\n')
                f.write('source leaprc.gaff\n')
                f.write('quit\n')
            iface = amberInterface(d)
            iface.tleapChecker(d)
            with open(path) as f:
                lines = f.readlines()
            self.assertEqual(lines.count('source leaprc.gaff\n'), 1)

    def test_readMetalBonds_and_connections(self):
        with TempDir() as d:
            with open(os.path.join(d, 'metalConnections'), 'w') as f:
                f.write('0 1 N\n0 5 N\n')
            with open(os.path.join(d, 'Connections'), 'w') as f:
                f.write('0 1 2 3\n4 5 6\n')
            iface = amberInterface(d)
            self.assertEqual(iface.readMetalBonds(d),
                             ['0-N 1', '0-N 5'])
            self.assertEqual(iface.readConnections(d),
                             ['0', '1', '2', '3', '4', '5', '6'])

    def test_equil_writes_21_files(self):
        with TempDir() as d:
            os.makedirs(os.path.join(d, 'equilibration'))
            iface = amberInterface(d)
            iface.equil(d)
            files = sorted(os.listdir(os.path.join(d, 'equilibration')))
            self.assertEqual(len(files), 21)


if __name__ == '__main__':
    unittest.main()
