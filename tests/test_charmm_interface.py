import os
import unittest
from unittest import mock

import numpy as np

import tests  # noqa: F401
from PyConSolv.interfaces.charmm import CharmmInterface
from tests.helpers import TempDir, write


SAMPLE_STR = """* stream
read rtf card append
MASS -1 CG331 12.011
RESI LIG 0.0
ATOM C1 CG331 -0.27
END

read para card flex append
BONDS
CG331 HGA1 322.0 1.111
END
"""


class TestCharmmInterface(unittest.TestCase):
    def test_init_creates_workdir(self):
        with TempDir() as d:
            sub = os.path.join(d, 'charmm_dir')
            iface = CharmmInterface(sub)
            self.assertTrue(os.path.isdir(sub))

    def test_generateLigandParams_writes_files(self):
        with TempDir() as d:
            mol2 = write(os.path.join(d, 'LIG.mol2'), 'dummy')
            iface = CharmmInterface(d)
            with mock.patch.object(iface.cgenff, 'parametrize') as para:
                str_path = os.path.join(d, 'LIG.str')
                write(str_path, SAMPLE_STR)
                para.return_value = str_path
                rtf, prm = iface.generateLigandParams(mol2)
            self.assertTrue(os.path.isfile(rtf))
            self.assertTrue(os.path.isfile(prm))
            self.assertIn(rtf, iface.rtf_files)
            self.assertIn(prm, iface.prm_files)

    def test_generateMetalParams_delegates_to_fftk(self):
        with TempDir() as d:
            iface = CharmmInterface(d)
            coords = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
            with mock.patch('PyConSolv.interfaces.charmm.parseOrcaHessian',
                            return_value=np.zeros((6, 6))), \
                 mock.patch.object(iface.fftk, 'deriveParameters',
                                   return_value=('m.rtf', 'm.prm')) as derive:
                rtf, prm = iface.generateMetalParams(
                    xyz_coords_ang=coords,
                    elements=['Fe', 'N'],
                    hessian_file='ignored',
                    bonds=[(0, 1)],
                    angles=[],
                    metal_indices=[0],
                    atom_types=['MFE', 'NG311'],
                )
            derive.assert_called_once()
            self.assertEqual((rtf, prm), ('m.rtf', 'm.prm'))

    def test_mergeParameters(self):
        with TempDir() as d:
            iface = CharmmInterface(d)
            rtf = write(os.path.join(d, 'a.rtf'),
                        '* header\nMASS -1 C 12.0\nRESI L 0.0\nATOM X C 0\nEND\n')
            prm = write(os.path.join(d, 'a.prm'),
                        '* header\nBONDS\nC C 300 1.5\nEND\n')
            iface.rtf_files = [rtf]
            iface.prm_files = [prm]
            rtf_out, prm_out = iface.mergeParameters(out_base='sys')
            self.assertTrue(os.path.isfile(rtf_out))
            self.assertTrue(os.path.isfile(prm_out))

    def test_solvate_without_psf_returns_empty(self):
        with TempDir() as d:
            iface = CharmmInterface(d)
            self.assertEqual(iface.solvate(), ('', ''))

    def test_checkDependencies_short_circuits(self):
        with TempDir() as d:
            iface = CharmmInterface(d)
            with mock.patch.object(iface.cgenff, 'checkpath', return_value=False), \
                 mock.patch.object(iface.packmol, 'checkpath', return_value=True):
                self.assertFalse(iface.checkDependencies())


if __name__ == '__main__':
    unittest.main()
