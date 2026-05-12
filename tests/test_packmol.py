import os
import unittest
from unittest import mock

import tests  # noqa: F401
from PyConSolv.interfaces.packmol import (PackmolInterface, SOLVENT_PROPERTIES,
                                          AVOGADRO)
from tests.helpers import TempDir, write


SOLUTE_PDB = """\
ATOM      1  N   LIG     1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  C   LIG     1       1.500   0.000   0.000  1.00  0.00           C
ATOM      3  O   LIG     1       0.000   1.500   0.000  1.00  0.00           O
END
"""


class TestPackmolStatics(unittest.TestCase):
    def test_solvent_count_water(self):
        # 1 nm^3 = 1000 A^3; ~33 water molecules at 1 g/mL
        n = PackmolInterface.solventCount(1000.0, 18.015, 1.0)
        self.assertGreaterEqual(n, 30)
        self.assertLessEqual(n, 36)

    def test_solvent_count_scales_linearly(self):
        a = PackmolInterface.solventCount(1000.0, 18.015, 1.0)
        b = PackmolInterface.solventCount(2000.0, 18.015, 1.0)
        self.assertAlmostEqual(b / a, 2.0, delta=0.05)

    def test_ion_count(self):
        # 1 L = 1e27 A^3 at 0.15 M => 0.15 * N_A ions
        n = PackmolInterface.ionCount(1e27, 0.15)
        expected = round(0.15 * AVOGADRO)
        self.assertEqual(n, expected)

    def test_known_solvents(self):
        for key in ('water', 'acetonitrile', 'methanol', 'dmso'):
            self.assertIn(key, SOLVENT_PROPERTIES)
            self.assertGreater(SOLVENT_PROPERTIES[key]['mw'], 0)
            self.assertGreater(SOLVENT_PROPERTIES[key]['density'], 0)


class TestPackmolFunctions(unittest.TestCase):
    def test_checkpath_missing(self):
        iface = PackmolInterface(packmol_cmd='nonexistent-packmol-x9')
        self.assertFalse(iface.checkpath())

    @mock.patch('shutil.which', return_value='/usr/bin/packmol')
    def test_checkpath_found(self, _which):
        self.assertTrue(PackmolInterface().checkpath())

    def test_boxDimensions(self):
        with TempDir() as d:
            pdb = write(os.path.join(d, 's.pdb'), SOLUTE_PDB)
            iface = PackmolInterface()
            box = iface.boxDimensions(pdb, padding=5.0)
            self.assertEqual(len(box), 6)
            # mins are below 0 after padding, maxes above solute extents
            self.assertLess(box[0], 0)
            self.assertGreater(box[3], 1.5)

    def test_boxDimensions_empty_pdb_raises(self):
        with TempDir() as d:
            pdb = write(os.path.join(d, 'e.pdb'), 'REMARK empty\nEND\n')
            iface = PackmolInterface()
            with self.assertRaises(ValueError):
                iface.boxDimensions(pdb, padding=5.0)

    def test_generateScript_contains_key_directives(self):
        iface = PackmolInterface()
        script = iface.generateScript(
            solute_pdb='solute.pdb',
            components=[{'pdb': 'water.pdb', 'count': 100}],
            box=(-5., -5., -5., 5., 5., 5.),
            output_pdb='out.pdb')
        self.assertIn('tolerance', script)
        self.assertIn('output out.pdb', script)
        self.assertIn('structure solute.pdb', script)
        self.assertIn('structure water.pdb', script)
        self.assertIn('number 100', script)
        self.assertIn('inside box', script)

    def test_execute_success(self):
        with TempDir() as d:
            iface = PackmolInterface()
            with mock.patch('subprocess.run') as run:
                run.return_value = mock.Mock(returncode=0)
                ok = iface.execute('tolerance 2.0\n', d)
            self.assertEqual(ok, 1)
            self.assertEqual(iface.status, 1)
            self.assertTrue(os.path.isfile(os.path.join(d, 'packmol.inp')))

    def test_execute_failure(self):
        with TempDir() as d:
            iface = PackmolInterface()
            with mock.patch('subprocess.run') as run:
                run.return_value = mock.Mock(returncode=1)
                self.assertEqual(iface.execute('x', d), 0)
            self.assertEqual(iface.status, 0)

    def test_buildBox_end_to_end(self):
        with TempDir() as d:
            solute = write(os.path.join(d, 's.pdb'), SOLUTE_PDB)
            water = write(os.path.join(d, 'w.pdb'), SOLUTE_PDB)  # stand-in
            iface = PackmolInterface()
            with mock.patch('subprocess.run') as run:
                run.return_value = mock.Mock(returncode=0)
                # fabricate expected output so isfile() check passes
                out = os.path.join(d, 'solvated.pdb')
                write(out, 'ATOM\n')
                result = iface.buildBox(solute, [water], d, padding=5.0,
                                        solvent_name='water')
            self.assertEqual(result, out)


if __name__ == '__main__':
    unittest.main()
