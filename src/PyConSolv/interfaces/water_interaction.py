"""ORCA single-point batch driver for FFTK-style water-interaction
charge fitting.

Given a list of TIP3P probes from misc.polar_sites, runs:
  - one isolated-ligand SP (cached)
  - one isolated-water SP (cached, standard TIP3P-like geometry)
  - one ligand+water complex SP per probe

For each probe it computes the interaction energy
    E_int = E(complex) - E(ligand) - E(water)
and applies FFTK's empirical scaling (1.16x on E_int) so the result is
on the same energetic footing CHARMM's TIP3P was originally calibrated
against.

Vacuum SPs by design (FFTK convention). Reference level defaults to
``HF 6-31G*`` for the same reason; pass ``level`` to override.
"""
import os
import shutil
import subprocess
from dataclasses import dataclass

import numpy as np

from .fftk import HARTREE_TO_KCAL
from ..utils.colorgen import Color


FFTK_SCALE = 1.16
FFTK_DIST_SHIFT = -0.2

ISOLATED_WATER_COORDS = np.array([
    [0.0000000,  0.0000000, 0.0000000],
    [0.9572000,  0.0000000, 0.0000000],
    [-0.2399837, 0.9266272, 0.0000000],
])


@dataclass
class InteractionResult:
    site_label: str
    probe_idx: int
    geometry: np.ndarray       # (N_lig + 3, 3) coords passed to QM
    n_ligand_atoms: int
    e_complex: float           # Hartree
    e_ligand: float            # Hartree
    e_water: float             # Hartree
    e_int_kcal: float          # raw, kcal/mol
    e_int_scaled: float        # FFTK-scaled, kcal/mol
    distance: float            # placed heavy-atom-to-probe-atom distance (A)


class WaterInteractionCharges:
    """Orchestrate ORCA single-points for a probe batch."""

    def __init__(self, workdir: str, level: str = 'HF 6-31G*',
                 cpu: int = 12, memory: int = 2000,
                 scale_factor: float = FFTK_SCALE,
                 orca_cmd: str = None):
        self.workdir = workdir
        self.level = level
        self.cpu = cpu
        self.memory = memory
        self.scale_factor = scale_factor
        self.orca_cmd = orca_cmd or shutil.which('orca') or 'orca'
        self.status = 0
        self._cache = {}
        if workdir and not os.path.isdir(workdir):
            os.makedirs(workdir, exist_ok=True)

    def checkpath(self) -> bool:
        if shutil.which(self.orca_cmd) is None and not os.path.isfile(self.orca_cmd):
            print(Color.RED + 'ORCA executable not found '
                              '(set orca_cmd or add to PATH)' + Color.END)
            self.status = 0
            return False
        self.status = 1
        return True

    def writeInput(self, name: str, coords: np.ndarray, elements: list,
                    charge: int = 0, multiplicity: int = 1) -> str:
        inp_path = os.path.join(self.workdir, name + '.inp')
        with open(inp_path, 'w') as f:
            f.write('! {} SP\n\n'.format(self.level))
            if self.cpu > 1:
                f.write('%PAL NPROCS {} END\n'.format(self.cpu))
            f.write('%maxcore {}\n'.format(self.memory))
            f.write('%scf maxiter 200 end\n\n')
            f.write('* xyz {} {}\n'.format(charge, multiplicity))
            for el, xyz in zip(elements, coords):
                f.write('{:<2s} {:>14.8f} {:>14.8f} {:>14.8f}\n'.format(
                    el, xyz[0], xyz[1], xyz[2]))
            f.write('*\n')
        return inp_path

    def runORCA(self, inp_path: str) -> str:
        """Run ORCA on inp_path; return the .out path, empty on failure."""
        out_path = inp_path[:-4] + '.out' if inp_path.endswith('.inp') \
                    else inp_path + '.out'
        cmd = '{} {} > {}'.format(self.orca_cmd, inp_path, out_path)
        calc = subprocess.run([cmd], shell=True, cwd=self.workdir)
        if calc.returncode != 0:
            print(Color.RED + 'ORCA failed for {}'.format(inp_path) + Color.END)
            return ''
        return out_path

    @staticmethod
    def parseEnergy(out_path: str) -> float:
        """Return last FINAL SINGLE POINT ENERGY in Hartree, or None."""
        if not os.path.isfile(out_path):
            return None
        energy = None
        with open(out_path, 'r') as f:
            for line in f:
                if 'FINAL SINGLE POINT ENERGY' in line:
                    try:
                        energy = float(line.split()[-1])
                    except (ValueError, IndexError):
                        pass
        return energy

    def _isolatedLigand(self, coords: np.ndarray, elements: list,
                          charge: int, multiplicity: int) -> float:
        if 'ligand' in self._cache:
            return self._cache['ligand']
        inp = self.writeInput('isolated_ligand', coords, elements,
                                charge, multiplicity)
        out = self.runORCA(inp)
        if not out:
            return None
        e = self.parseEnergy(out)
        self._cache['ligand'] = e
        return e

    def _isolatedWater(self) -> float:
        if 'water' in self._cache:
            return self._cache['water']
        inp = self.writeInput('isolated_water', ISOLATED_WATER_COORDS,
                                ['O', 'H', 'H'], 0, 1)
        out = self.runORCA(inp)
        if not out:
            return None
        e = self.parseEnergy(out)
        self._cache['water'] = e
        return e

    def runBatch(self, probes: list, ligand_coords: np.ndarray,
                  ligand_elements: list, charge: int = 0,
                  multiplicity: int = 1) -> list:
        """Run isolated SPs (cached) and one complex SP per probe.

        Returns list of InteractionResult for probes whose SP succeeded.
        Failed probes are skipped (warning printed); the batch continues.
        """
        e_lig = self._isolatedLigand(ligand_coords, ligand_elements,
                                       charge, multiplicity)
        e_wat = self._isolatedWater()
        if e_lig is None or e_wat is None:
            print(Color.RED + 'Isolated reference SP failed; aborting batch'
                  + Color.END)
            return []

        results = []
        for idx, probe in enumerate(probes):
            combined_coords = np.vstack([ligand_coords, probe.water_coords])
            combined_elements = list(ligand_elements) + ['O', 'H', 'H']
            name = 'probe_{:03d}_{}'.format(idx, probe.site.label)
            inp = self.writeInput(name, combined_coords, combined_elements,
                                    charge, multiplicity)
            out = self.runORCA(inp)
            if not out:
                continue
            e_cx = self.parseEnergy(out)
            if e_cx is None:
                continue
            e_int_kcal = (e_cx - e_lig - e_wat) * HARTREE_TO_KCAL
            e_int_scaled = e_int_kcal * self.scale_factor

            if probe.site.role == 'donor':
                target = ligand_coords[probe.site.h_idx]
            else:
                target = ligand_coords[probe.site.heavy_idx]
            r_min = float(np.linalg.norm(probe.water_coords - target,
                                          axis=1).min())

            results.append(InteractionResult(
                site_label=probe.site.label,
                probe_idx=idx,
                geometry=combined_coords,
                n_ligand_atoms=len(ligand_elements),
                e_complex=e_cx,
                e_ligand=e_lig,
                e_water=e_wat,
                e_int_kcal=e_int_kcal,
                e_int_scaled=e_int_scaled,
                distance=r_min,
            ))
            if (idx + 1) % 10 == 0:
                print(Color.PURPLE + '  ... {} / {} probes complete'.format(
                    idx + 1, len(probes)) + Color.END)
        self.status = 1
        return results
