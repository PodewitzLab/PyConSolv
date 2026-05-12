import os
import shutil
import subprocess

from ..utils.colorgen import Color


# Approximate liquid-phase densities (g/mL) and molar masses (g/mol) for
# common CHARMM solvents. Used to compute molecule counts when the caller
# doesn't specify one explicitly.
SOLVENT_PROPERTIES = {
    'water':        {'density': 1.000, 'mw': 18.015},
    'tip3p':        {'density': 1.000, 'mw': 18.015},
    'tip4p':        {'density': 1.000, 'mw': 18.015},
    'spce':         {'density': 1.000, 'mw': 18.015},
    'acetonitrile': {'density': 0.786, 'mw': 41.053},
    'methanol':     {'density': 0.792, 'mw': 32.042},
    'ethanol':      {'density': 0.789, 'mw': 46.068},
    'dmso':         {'density': 1.100, 'mw': 78.133},
    'dmf':          {'density': 0.944, 'mw': 73.094},
    'chloroform':   {'density': 1.489, 'mw': 119.38},
    'ccl4':         {'density': 1.594, 'mw': 153.82},
    'thf':          {'density': 0.889, 'mw': 72.107},
    'toluene':      {'density': 0.867, 'mw': 92.141},
    'benzene':      {'density': 0.879, 'mw': 78.114},
    'cyclohexane':  {'density': 0.779, 'mw': 84.160},
    'hexane':       {'density': 0.655, 'mw': 86.178},
    'ch2cl2':       {'density': 1.327, 'mw': 84.930},
    'octanol':      {'density': 0.824, 'mw': 130.23},
    'ammonia':      {'density': 0.682, 'mw': 17.031},
}

# Atomic masses for common counter-ion flavours (CHARMM atom names).
ION_MW = {'SOD': 22.990, 'CLA': 35.450, 'POT': 39.098, 'CAL': 40.078,
          'MG':  24.305, 'ZN':  65.380, 'LIT': 6.941, 'NA': 22.990, 'CL': 35.450}

AVOGADRO = 6.02214076e23


class PackmolInterface:
    """Wrapper around the Packmol binary for assembling solvated boxes.

    Packmol works for any solvent (water, organics, mixtures) — all we need
    is a single-molecule PDB for each component and a target box size.
    """

    def __init__(self, packmol_cmd: str = 'packmol'):
        self.packmol_cmd = packmol_cmd
        self.status = 0

    def checkpath(self) -> bool:
        if shutil.which(self.packmol_cmd) is None:
            print(Color.RED + 'Packmol not found in PATH' + Color.END)
            self.status = 0
            return False
        self.status = 1
        return True

    @staticmethod
    def solventCount(box_volume_A3: float, mw: float, density_g_ml: float) -> int:
        """Number of solvent molecules needed to reach target density.

        box_volume_A3: volume in cubic angstrom
        density in g/mL, mw in g/mol.
        """
        # box_volume in mL = A^3 * 1e-24
        mass_g = density_g_ml * box_volume_A3 * 1e-24
        moles = mass_g / mw
        return int(round(moles * AVOGADRO))

    @staticmethod
    def ionCount(box_volume_A3: float, concentration_M: float) -> int:
        """Number of ion pairs for a given molar concentration."""
        # 1 L = 1e27 A^3
        volume_L = box_volume_A3 / 1e27
        return int(round(concentration_M * AVOGADRO * volume_L))

    def boxDimensions(self, solute_pdb: str, padding: float) -> tuple:
        """Compute an axis-aligned box large enough to wrap the solute plus
        padding on every face. Returns (xmin, ymin, zmin, xmax, ymax, zmax).
        """
        xs, ys, zs = [], [], []
        with open(solute_pdb, 'r') as f:
            for line in f:
                if not line.startswith(('ATOM', 'HETATM')):
                    continue
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
        if not xs:
            raise ValueError('No atoms found in {}'.format(solute_pdb))
        return (min(xs) - padding, min(ys) - padding, min(zs) - padding,
                max(xs) + padding, max(ys) + padding, max(zs) + padding)

    def generateScript(self, solute_pdb: str, components: list,
                       box: tuple, output_pdb: str,
                       tolerance: float = 2.0) -> str:
        """Build a packmol input file.

        components: list of dicts with keys {'pdb', 'count'}; the solute is
        added separately as a fixed structure at the origin.
        box: (xmin, ymin, zmin, xmax, ymax, zmax)
        """
        lines = [
            'tolerance {}'.format(tolerance),
            'filetype pdb',
            'output {}'.format(output_pdb),
            '',
            'structure {}'.format(solute_pdb),
            '  number 1',
            '  fixed 0. 0. 0. 0. 0. 0.',
            '  centerofmass',
            'end structure',
            '',
        ]
        for comp in components:
            lines += [
                'structure {}'.format(comp['pdb']),
                '  number {}'.format(comp['count']),
                '  inside box {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} {:.3f}'.format(*box),
                'end structure',
                '',
            ]
        return '\n'.join(lines)

    def execute(self, script: str, workdir: str,
                script_name: str = 'packmol.inp') -> int:
        script_path = os.path.join(workdir, script_name)
        with open(script_path, 'w') as f:
            f.write(script)
        cmd = '{} < {} > packmol.log'.format(self.packmol_cmd, script_name)
        calc = subprocess.run([cmd], shell=True, cwd=workdir)
        if calc.returncode != 0:
            print(Color.RED + 'Packmol failed (see packmol.log)' + Color.END)
            self.status = 0
            return 0
        self.status = 1
        return 1

    def buildBox(self, solute_pdb: str, solvent_pdbs: list,
                 workdir: str, padding: float = 10.0,
                 solvent_name: str = 'water',
                 ion_pdbs: list = None,
                 ion_concentration: float = 0.0,
                 output_pdb: str = 'solvated.pdb') -> str:
        """High-level entry point.

        solvent_pdbs: list of single-molecule PDB paths (mixtures allowed)
        ion_pdbs: list of 2 single-atom PDBs (cation, anion) or None
        Returns absolute path to the packed PDB, or '' on failure.
        """
        box = self.boxDimensions(solute_pdb, padding)
        volume = (box[3] - box[0]) * (box[4] - box[1]) * (box[5] - box[2])

        props = SOLVENT_PROPERTIES.get(solvent_name.lower(),
                                       SOLVENT_PROPERTIES['water'])
        n_solvent = self.solventCount(volume, props['mw'], props['density'])
        # Split the count evenly if multiple solvent PDBs (co-solvent mixtures).
        per_component = max(1, n_solvent // max(1, len(solvent_pdbs)))

        components = [{'pdb': pdb, 'count': per_component} for pdb in solvent_pdbs]
        if ion_pdbs and ion_concentration > 0:
            n_ions = self.ionCount(volume, ion_concentration)
            for ion_pdb in ion_pdbs:
                components.append({'pdb': ion_pdb, 'count': n_ions})

        output_abs = os.path.join(workdir, output_pdb)
        script = self.generateScript(solute_pdb, components, box, output_abs)
        if not self.execute(script, workdir):
            return ''
        return output_abs if os.path.isfile(output_abs) else ''
