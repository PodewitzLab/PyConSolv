import os

import numpy as np

from .cgenff import CGenFFInterface
from .fftk import FFTKBondedInterface, parseOrcaHessian
from .packmol import PackmolInterface
from .charmm_builder import CharmmBuilder
from ..misc import charmm_formats
from ..utils.colorgen import Color


DB_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'db')


class CharmmInterface:
    """Top-level CHARMM36 parametrization driver using ParmEd + Packmol.

    Coordinates CGenFF (organic ligands), the native FFTK Seminario
    back-end (metal centers + high-penalty bonded terms), ParmEd (PSF
    assembly), and Packmol (solvation). No VMD required; no easyPARM.
    """

    def __init__(self, path: str, cgenff_path: str = None,
                 packmol_cmd: str = 'packmol',
                 penalty_threshold: float = 50.0):
        self.path = path
        self.status = 0
        self.cgenff = CGenFFInterface(cgenff_path=cgenff_path)
        self.fftk = FFTKBondedInterface(path=path,
                                        penalty_threshold=penalty_threshold)
        self.packmol = PackmolInterface(packmol_cmd=packmol_cmd)
        self.builder = CharmmBuilder(workdir=path)

        self.rtf_files = []
        self.prm_files = []
        self.psf = None
        self.pdb = None
        self.cgenff_penalties = {}
        if path and not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)

    def checkDependencies(self) -> bool:
        ok = self.cgenff.checkpath()
        ok = self.packmol.checkpath() and ok
        self.status = 1 if ok else 0
        return ok

    def generateLigandParams(self, mol2_file: str) -> tuple:
        str_file = self.cgenff.parametrize(mol2_file, output_dir=self.path)
        if not str_file:
            return '', ''
        rtf_text, prm_text = self.cgenff.splitSTR(str_file)
        base = os.path.splitext(os.path.basename(mol2_file))[0]
        rtf_path = os.path.join(self.path, base + '.rtf')
        prm_path = os.path.join(self.path, base + '.prm')
        with open(rtf_path, 'w') as f:
            f.write(rtf_text)
            if not rtf_text.rstrip().upper().endswith('END'):
                f.write('\nEND\n')
        with open(prm_path, 'w') as f:
            f.write(prm_text)
            if not prm_text.rstrip().upper().endswith('END'):
                f.write('\nEND\n')

        penalties = self.cgenff.parsePenalties(str_file)
        self.cgenff_penalties = penalties
        bad = {k: v for k, v in penalties.items() if v > 50}
        if bad:
            print(Color.YELLOW + 'CGenFF penalties > 50 (will be refined by '
                                  'Seminario/FFTK): {}'.format(bad) + Color.END)

        self.rtf_files.append(rtf_path)
        self.prm_files.append(prm_path)
        return rtf_path, prm_path

    def applyRESPCharges(self, rtf_path: str, resp_charges: list,
                          atom_names: list = None,
                          target_total: float = None) -> str:
        """Overwrite CGenFF charges in ``rtf_path`` with RESP values.

        ``resp_charges`` must be in XYZ atom order. ``atom_names`` (also in
        XYZ order) lets the writer match by name rather than position,
        which is more robust when CGenFF reorders atoms inside the RESI
        block. ``target_total`` is the molecular charge used for rounding
        drift correction.
        """
        return charmm_formats.rewriteRTFCharges(
            rtf_path=rtf_path,
            charges=resp_charges,
            atom_order=atom_names,
            target_total=target_total,
        )

    def generateMetalParams(self, xyz_coords_ang, elements: list,
                            hessian_file: str,
                            bonds: list, angles: list,
                            metal_indices: list,
                            atom_types: list,
                            penalties: dict = None) -> tuple:
        """Native CHARMM36 bonded parametrization via Seminario.

        xyz_coords_ang : (N, 3) numpy array in Angstroms.
        elements       : list of element symbols.
        hessian_file   : path to ORCA .hess file.
        bonds, angles  : connectivity lists (0-based indices).
        metal_indices  : list of metal atom indices.
        atom_types     : list of CGenFF (and placeholder) types per atom.
        penalties      : optional {atom_idx: cgenff_penalty}.
        """
        H = parseOrcaHessian(hessian_file)
        penalties = penalties or self.cgenff_penalties
        rtf, prm = self.fftk.deriveParameters(
            xyz_coords_ang=np.asarray(xyz_coords_ang, dtype=float),
            elements=elements,
            hessian=H,
            bonds=bonds,
            angles=angles,
            metal_indices=metal_indices,
            atom_types=atom_types,
            penalties=penalties,
            out_basename='fftk_metal',
        )
        if rtf:
            self.rtf_files.append(rtf)
        if prm:
            self.prm_files.append(prm)
        return rtf, prm

    def mergeParameters(self, out_base: str = 'system') -> tuple:
        rtf_out = os.path.join(self.path, out_base + '.rtf')
        prm_out = os.path.join(self.path, out_base + '.prm')
        charmm_formats.mergeRTF(self.rtf_files, rtf_out)
        charmm_formats.mergePRM(self.prm_files, prm_out)
        return rtf_out, prm_out

    def buildPSF(self, pdb_file: str, metal_bonds: list = None,
                 basename: str = 'system') -> tuple:
        """Assemble the solute PSF via ParmEd."""
        self.builder.loadParameters(self.rtf_files, self.prm_files)
        self.builder.loadStructure(pdb_file, self.rtf_files)
        if metal_bonds:
            self.builder.addMetalBonds(metal_bonds)
        self.builder.applyParameters()
        psf, pdb = self.builder.writePSF(basename=basename)
        self.psf, self.pdb = psf, pdb
        return psf, pdb

    def solvate(self, padding: float = 10.0, solvent_name: str = 'water',
                solvent_pdbs: list = None,
                solvent_rtf_files: list = None,
                solvent_prm_files: list = None,
                concentration: float = 0.15,
                cation: str = 'SOD', anion: str = 'CLA') -> tuple:
        """Pack solvent and ions around the solute with Packmol, then attach
        parameters and emit the final PSF/PDB pair.
        """
        if not (self.psf and self.pdb):
            print(Color.RED + 'No PSF/PDB built yet' + Color.END)
            return '', ''

        if solvent_pdbs is None:
            solvent_pdbs = [os.path.join(DB_DIR, 'water_tip3p.pdb')]
        ion_pdbs = None
        if concentration > 0:
            ion_pdbs = [os.path.join(DB_DIR, 'ion_sod.pdb'),
                        os.path.join(DB_DIR, 'ion_cla.pdb')]

        packed = self.packmol.buildBox(
            solute_pdb=self.pdb,
            solvent_pdbs=solvent_pdbs,
            workdir=self.path,
            padding=padding,
            solvent_name=solvent_name,
            ion_pdbs=ion_pdbs,
            ion_concentration=concentration,
            output_pdb='solvated.pdb',
        )
        if not packed:
            return '', ''

        solvent_rtf_files = solvent_rtf_files or []
        solvent_prm_files = solvent_prm_files or []
        return self.builder.attachSolvatedBox(packed,
                                              solvent_rtf_files,
                                              solvent_prm_files,
                                              basename='solvated')
