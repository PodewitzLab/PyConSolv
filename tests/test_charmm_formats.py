import os
import unittest

import tests  # noqa: F401
from PyConSolv.misc import charmm_formats
from tests.helpers import TempDir, write, SAMPLE_XYZ


SAMPLE_RTF = """* header
* line two
27 1

MASS  -1  CG331  12.011 ! alkane carbon
MASS  -1  HGA1    1.008 ! alkane hydrogen

RESI LIG  0.00
GROUP
ATOM  C1  CG331 -0.27
ATOM  H1  HGA1   0.09
BOND  C1 H1

PRES DEPROT  -1.00
ATOM  O1  OG311 -0.50

END
"""

SAMPLE_PRM = """* PRM header

BONDS
CG331 HGA1   322.0   1.111

ANGLES
HGA1 CG331 HGA1   35.5   109.0

DIHEDRALS
HGA1 CG331 CG331 HGA1   0.16  3  0.0

IMPROPERS

NONBONDED nbxmod 5
CG331   0.0  -0.0780  2.0500
HGA1    0.0  -0.0240  1.3400
END
"""


class TestRTF(unittest.TestCase):
    def test_parse_collects_mass_and_resi(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'sample.rtf'), SAMPLE_RTF)
            data = charmm_formats.parseRTF(path)
            self.assertEqual(len(data['mass']), 2)
            self.assertIn('LIG', data['resi'])
            self.assertIn('DEPROT', data['pres'])
            self.assertTrue(any('ATOM  C1' in line for line in data['resi']['LIG']))

    def test_write_roundtrip_contains_resi(self):
        with TempDir() as d:
            src = write(os.path.join(d, 'a.rtf'), SAMPLE_RTF)
            dst = os.path.join(d, 'b.rtf')
            data = charmm_formats.parseRTF(src)
            charmm_formats.writeRTF(data, dst)
            with open(dst) as f:
                text = f.read()
            self.assertIn('RESI LIG', text)
            self.assertIn('PRES DEPROT', text)
            self.assertIn('MASS', text)
            self.assertTrue(text.rstrip().endswith('END'))

    def test_merge_deduplicates_mass(self):
        with TempDir() as d:
            a = write(os.path.join(d, 'a.rtf'), SAMPLE_RTF)
            b = write(os.path.join(d, 'b.rtf'), SAMPLE_RTF.replace('LIG', 'LG2'))
            out = os.path.join(d, 'merged.rtf')
            charmm_formats.mergeRTF([a, b], out)
            data = charmm_formats.parseRTF(out)
            self.assertEqual(len(data['mass']), 2)  # deduped
            self.assertIn('LIG', data['resi'])
            self.assertIn('LG2', data['resi'])


class TestPRM(unittest.TestCase):
    def test_parse_sections(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'a.prm'), SAMPLE_PRM)
            data = charmm_formats.parsePRM(path)
            self.assertTrue(any('CG331 HGA1' in l for l in data['BONDS']))
            self.assertTrue(any('HGA1 CG331 HGA1' in l for l in data['ANGLES']))
            self.assertTrue(any('CG331' in l for l in data['NONBONDED']))

    def test_write_has_section_headers(self):
        with TempDir() as d:
            src = write(os.path.join(d, 'a.prm'), SAMPLE_PRM)
            dst = os.path.join(d, 'b.prm')
            charmm_formats.writePRM(charmm_formats.parsePRM(src), dst)
            with open(dst) as f:
                text = f.read()
            for section in ('BONDS', 'ANGLES', 'DIHEDRALS', 'NONBONDED', 'END'):
                self.assertIn(section, text)

    def test_merge_preserves_unique(self):
        with TempDir() as d:
            a = write(os.path.join(d, 'a.prm'), SAMPLE_PRM)
            b = write(os.path.join(d, 'b.prm'), SAMPLE_PRM.replace('HGA1', 'HXX1'))
            out = os.path.join(d, 'merged.prm')
            charmm_formats.mergePRM([a, b], out)
            data = charmm_formats.parsePRM(out)
            joined = ''.join(data['BONDS'])
            self.assertIn('HGA1', joined)
            self.assertIn('HXX1', joined)


class TestXYZtoPDB(unittest.TestCase):
    def test_pdb_atom_count_matches_xyz(self):
        with TempDir() as d:
            xyz_path = write(os.path.join(d, 'w.xyz'), SAMPLE_XYZ)
            pdb_path = os.path.join(d, 'w.pdb')
            charmm_formats.xyzToCharmmPDB(xyz_path, pdb_path,
                                          resname='WAT', segid='W')
            with open(pdb_path) as f:
                lines = f.readlines()
            atoms = [l for l in lines if l.startswith('ATOM')]
            self.assertEqual(len(atoms), 3)
            self.assertIn('WAT', atoms[0])
            # segid should land in columns 73-76 region
            self.assertIn('W', atoms[0][70:])

    def test_names_truncated_to_4_chars(self):
        with TempDir() as d:
            xyz_path = write(os.path.join(d, 'w.xyz'), SAMPLE_XYZ)
            pdb_path = os.path.join(d, 'w.pdb')
            charmm_formats.xyzToCharmmPDB(xyz_path, pdb_path)
            with open(pdb_path) as f:
                for line in f:
                    if line.startswith('ATOM'):
                        name = line[12:16].strip()
                        self.assertLessEqual(len(name), 4)


class TestRewriteRTFCharges(unittest.TestCase):
    def test_overwrite_by_name(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'lig.rtf'), SAMPLE_RTF)
            charmm_formats.rewriteRTFCharges(
                path, charges=[0.1234, -0.1234, 0.5],
                atom_order=['C1', 'H1', 'O1'])
            with open(path) as f:
                text = f.read()
            self.assertIn('0.1234', text)
            self.assertIn('-0.1234', text)
            self.assertIn('0.5000', text)
            self.assertNotIn('-0.27', text)
            self.assertNotIn(' 0.09', text)

    def test_overwrite_by_position(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'lig.rtf'), SAMPLE_RTF)
            charmm_formats.rewriteRTFCharges(
                path, charges=[0.01, 0.02, 0.03])
            with open(path) as f:
                lines = f.readlines()
            atoms = [l for l in lines if l.strip().startswith('ATOM')]
            self.assertEqual(atoms[0].split()[3], '0.0100')
            self.assertEqual(atoms[1].split()[3], '0.0200')
            self.assertEqual(atoms[2].split()[3], '0.0300')

    def test_target_total_absorbs_drift(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'lig.rtf'), SAMPLE_RTF)
            # charges sum to 0.00005; rounding to 4dp drops it, target_total
            # should force the sum to exactly 0.0000 by adjusting the
            # largest-|q| atom.
            charmm_formats.rewriteRTFCharges(
                path, charges=[0.10001, -0.10002, 0.00006],
                atom_order=['C1', 'H1', 'O1'],
                target_total=0.0)
            with open(path) as f:
                lines = f.readlines()
            qs = [float(l.split()[3]) for l in lines
                  if l.strip().startswith('ATOM')]
            self.assertAlmostEqual(sum(qs), 0.0, places=4)

    def test_unknown_name_ignored(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'lig.rtf'), SAMPLE_RTF)
            charmm_formats.rewriteRTFCharges(
                path, charges=[0.5, 0.6],
                atom_order=['XX', 'YY'])
            with open(path) as f:
                text = f.read()
            # No atoms matched; file should be untouched.
            self.assertIn('-0.27', text)


class TestReadRTFAtoms(unittest.TestCase):
    def test_collects_name_type_charge(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'lig.rtf'), SAMPLE_RTF)
            atoms = charmm_formats.readRTFAtoms(path)
            names = [a[0] for a in atoms]
            self.assertIn('C1', names)
            self.assertIn('H1', names)
            self.assertIn('O1', names)
            by_name = {a[0]: a for a in atoms}
            self.assertEqual(by_name['C1'][1], 'CG331')
            self.assertAlmostEqual(by_name['C1'][2], -0.27)
            self.assertAlmostEqual(by_name['O1'][2], -0.50)

    def test_empty_rtf_returns_empty_list(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'empty.rtf'), '* header\nEND\n')
            self.assertEqual(charmm_formats.readRTFAtoms(path), [])


class TestReadNonbondedLJ(unittest.TestCase):
    def test_parses_eps_and_rmin(self):
        with TempDir() as d:
            path = write(os.path.join(d, 'a.prm'), SAMPLE_PRM)
            lj = charmm_formats.readNonbondedLJ(path)
            self.assertIn('CG331', lj)
            eps, rmin2 = lj['CG331']
            self.assertAlmostEqual(eps, -0.0780)
            self.assertAlmostEqual(rmin2, 2.0500)
            self.assertIn('HGA1', lj)

    def test_ignores_1_4_columns(self):
        prm = ('* header\n\nNONBONDED nbxmod 5\n'
               'CG331  0.0  -0.078  2.05   0.0  -0.01  1.9\n'
               'END\n')
        with TempDir() as d:
            path = write(os.path.join(d, 'a.prm'), prm)
            lj = charmm_formats.readNonbondedLJ(path)
            self.assertEqual(len(lj), 1)
            self.assertAlmostEqual(lj['CG331'][0], -0.078)
            self.assertAlmostEqual(lj['CG331'][1], 2.05)


if __name__ == '__main__':
    unittest.main()
