import os
import unittest
from unittest import mock

import numpy as np

import tests  # noqa: F401
from PyConSolv.interfaces.water_interaction import (
    WaterInteractionCharges, InteractionResult,
    ISOLATED_WATER_COORDS, FFTK_SCALE,
)
from PyConSolv.interfaces.fftk import HARTREE_TO_KCAL
from PyConSolv.misc.polar_sites import (
    detectPolarSites, placeWaters,
)
from tests.helpers import TempDir, write
from tests.test_polar_sites import methanol, acetone


SAMPLE_ORCA_OUT_TEMPLATE = """ORCA output stub
... lots of header omitted ...
SCF iterations converged.

FINAL SINGLE POINT ENERGY      {energy:.10f}

... lots of footer omitted ...
"""


def fake_orca(energy_map):
    """Returns a side_effect for subprocess.run that writes a fake ORCA .out
    next to whatever input file the command references.

    energy_map: dict mapping input-file basename -> Hartree energy.
    """
    def _run(args, **kwargs):
        cmd = args[0] if isinstance(args, list) else args
        parts = cmd.split('>')
        out_path = parts[1].strip() if len(parts) > 1 else None
        # The .inp lives in the command between the orca binary and the '>'.
        left = parts[0].strip().split()
        inp_path = left[-1]
        cwd = kwargs.get('cwd', '.')
        if not os.path.isabs(out_path):
            out_path = os.path.join(cwd, out_path)
        if not os.path.isabs(inp_path):
            inp_path = os.path.join(cwd, inp_path)
        base = os.path.splitext(os.path.basename(inp_path))[0]
        # Resolve the energy: exact basename match, or first matching prefix.
        if base in energy_map:
            e = energy_map[base]
        else:
            e = next((v for k, v in energy_map.items() if base.startswith(k)),
                      -100.0)
        with open(out_path, 'w') as f:
            f.write(SAMPLE_ORCA_OUT_TEMPLATE.format(energy=e))
        return mock.Mock(returncode=0)
    return _run


class TestWriteInput(unittest.TestCase):
    def test_input_contains_level_and_coords(self):
        with TempDir() as d:
            wic = WaterInteractionCharges(d, level='HF 6-31G*', cpu=4,
                                            memory=1000)
            coords = np.array([[0.0, 0.0, 0.0], [0.957, 0.0, 0.0]])
            elements = ['O', 'H']
            path = wic.writeInput('test', coords, elements, charge=0,
                                    multiplicity=1)
            self.assertTrue(os.path.isfile(path))
            with open(path) as f:
                text = f.read()
            self.assertIn('! HF 6-31G* SP', text)
            self.assertIn('%PAL NPROCS 4 END', text)
            self.assertIn('%maxcore 1000', text)
            self.assertIn('* xyz 0 1', text)
            self.assertIn('O ', text)
            self.assertIn('0.95700000', text)

    def test_input_respects_charge_and_multiplicity(self):
        with TempDir() as d:
            wic = WaterInteractionCharges(d)
            path = wic.writeInput('q', np.zeros((1, 3)), ['Fe'],
                                    charge=2, multiplicity=5)
            with open(path) as f:
                text = f.read()
            self.assertIn('* xyz 2 5', text)


class TestParseEnergy(unittest.TestCase):
    def test_extracts_final_sp_energy(self):
        with TempDir() as d:
            out = write(os.path.join(d, 'a.out'),
                        SAMPLE_ORCA_OUT_TEMPLATE.format(energy=-76.1234567))
            e = WaterInteractionCharges.parseEnergy(out)
            self.assertAlmostEqual(e, -76.1234567)

    def test_uses_last_occurrence(self):
        with TempDir() as d:
            text = (SAMPLE_ORCA_OUT_TEMPLATE.format(energy=-1.0)
                    + '\n' + SAMPLE_ORCA_OUT_TEMPLATE.format(energy=-2.0))
            out = write(os.path.join(d, 'multi.out'), text)
            e = WaterInteractionCharges.parseEnergy(out)
            self.assertAlmostEqual(e, -2.0)

    def test_missing_file_returns_none(self):
        self.assertIsNone(WaterInteractionCharges.parseEnergy('/nope.out'))

    def test_missing_line_returns_none(self):
        with TempDir() as d:
            out = write(os.path.join(d, 'a.out'), 'no energy here\n')
            self.assertIsNone(WaterInteractionCharges.parseEnergy(out))


class TestRunORCA(unittest.TestCase):
    def test_returns_out_path_on_success(self):
        with TempDir() as d:
            wic = WaterInteractionCharges(d)
            inp = write(os.path.join(d, 'x.inp'), 'dummy')
            with mock.patch('subprocess.run',
                            return_value=mock.Mock(returncode=0)):
                out = wic.runORCA(inp)
            self.assertTrue(out.endswith('x.out'))

    def test_returns_empty_on_failure(self):
        with TempDir() as d:
            wic = WaterInteractionCharges(d)
            inp = write(os.path.join(d, 'x.inp'), 'dummy')
            with mock.patch('subprocess.run',
                            return_value=mock.Mock(returncode=1)):
                self.assertEqual(wic.runORCA(inp), '')


class TestRunBatch(unittest.TestCase):
    def _setup_methanol_batch(self, d):
        coords, elements, bonds, types = methanol()
        sites = detectPolarSites(coords, elements, bonds, types)
        probes = placeWaters(sites, coords)
        wic = WaterInteractionCharges(d, cpu=2, memory=500)
        return wic, probes, coords, elements

    def test_end_to_end_with_mocked_orca(self):
        with TempDir() as d:
            wic, probes, coords, elements = self._setup_methanol_batch(d)
            energies = {
                'isolated_ligand': -115.0,
                'isolated_water': -76.0,
                'probe_000': -191.01,    # E_int = -0.01 Ha = -6.28 kcal/mol
                'probe_001': -191.005,
                'probe_002': -191.005,
            }
            with mock.patch('subprocess.run',
                            side_effect=fake_orca(energies)):
                results = wic.runBatch(probes, coords, elements,
                                        charge=0, multiplicity=1)
            self.assertEqual(len(results), len(probes))
            r0 = results[0]
            self.assertEqual(r0.n_ligand_atoms, len(elements))
            self.assertAlmostEqual(r0.e_ligand, -115.0)
            self.assertAlmostEqual(r0.e_water, -76.0)
            self.assertAlmostEqual(r0.e_complex, -191.01)
            expected = -0.01 * HARTREE_TO_KCAL
            self.assertAlmostEqual(r0.e_int_kcal, expected, places=4)
            self.assertAlmostEqual(r0.e_int_scaled,
                                    expected * FFTK_SCALE, places=4)

    def test_isolated_references_run_once(self):
        with TempDir() as d:
            wic, probes, coords, elements = self._setup_methanol_batch(d)
            energies = {
                'isolated_ligand': -115.0,
                'isolated_water': -76.0,
                'probe_': -191.0,
            }
            calls = []

            def tracking(args, **kwargs):
                calls.append(args[0] if isinstance(args, list) else args)
                return fake_orca(energies)(args, **kwargs)

            with mock.patch('subprocess.run', side_effect=tracking):
                wic.runBatch(probes, coords, elements)
            isolated_calls = [c for c in calls if 'isolated_' in c]
            # Exactly one isolated_ligand + one isolated_water (3 probes share).
            self.assertEqual(len(isolated_calls), 2)

    def test_skips_failed_probe_continues_batch(self):
        with TempDir() as d:
            wic, probes, coords, elements = self._setup_methanol_batch(d)
            energies = {
                'isolated_ligand': -115.0,
                'isolated_water': -76.0,
                'probe_': -191.0,
            }
            call_count = {'n': 0}

            def flaky(args, **kwargs):
                call_count['n'] += 1
                # Fail the second probe SP. First two calls are isolated refs.
                if call_count['n'] == 4:
                    return mock.Mock(returncode=1)
                return fake_orca(energies)(args, **kwargs)

            with mock.patch('subprocess.run', side_effect=flaky):
                results = wic.runBatch(probes, coords, elements)
            self.assertEqual(len(results), len(probes) - 1)

    def test_isolated_failure_aborts_batch(self):
        with TempDir() as d:
            wic, probes, coords, elements = self._setup_methanol_batch(d)
            with mock.patch('subprocess.run',
                            return_value=mock.Mock(returncode=1)):
                results = wic.runBatch(probes, coords, elements)
            self.assertEqual(results, [])


class TestDistanceCalculation(unittest.TestCase):
    def test_donor_distance_recorded(self):
        with TempDir() as d:
            coords, elements, bonds, types = methanol()
            sites = detectPolarSites(coords, elements, bonds, types)
            probes = placeWaters(sites, coords)
            donor_probe = next(p for p in probes if p.site.role == 'donor')
            wic = WaterInteractionCharges(d)
            energies = {'isolated_ligand': -115.0, 'isolated_water': -76.0,
                        'probe_': -191.0}
            with mock.patch('subprocess.run',
                            side_effect=fake_orca(energies)):
                results = wic.runBatch([donor_probe], coords, elements)
            self.assertEqual(len(results), 1)
            # Donor: distance = H_donor to water O (= DONOR_O_DIST = 2.0)
            self.assertAlmostEqual(results[0].distance, 2.0, places=3)

    def test_acceptor_distance_recorded(self):
        with TempDir() as d:
            coords, elements, bonds, types = acetone()
            sites = detectPolarSites(coords, elements, bonds, types)
            probes = placeWaters(sites, coords)
            acc_probe = next(p for p in probes if p.site.role == 'acceptor')
            wic = WaterInteractionCharges(d)
            energies = {'isolated_ligand': -100.0, 'isolated_water': -76.0,
                        'probe_': -176.0}
            with mock.patch('subprocess.run',
                            side_effect=fake_orca(energies)):
                results = wic.runBatch([acc_probe], coords, elements)
            self.assertEqual(len(results), 1)
            # Acceptor: O_acceptor to water H1 (= ACCEPTOR_H_DIST = 1.85)
            self.assertAlmostEqual(results[0].distance, 1.85, places=3)


class TestCheckPath(unittest.TestCase):
    def test_missing_orca(self):
        with TempDir() as d:
            wic = WaterInteractionCharges(d, orca_cmd='/nonexistent/orca-x9')
            self.assertFalse(wic.checkpath())

    @mock.patch('shutil.which', return_value='/usr/local/bin/orca')
    def test_present_orca(self, _which):
        with TempDir() as d:
            wic = WaterInteractionCharges(d)
            self.assertTrue(wic.checkpath())


if __name__ == '__main__':
    unittest.main()
