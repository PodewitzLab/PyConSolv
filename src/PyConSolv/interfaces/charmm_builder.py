import os

import parmed as pmd

from ..utils.colorgen import Color


class CharmmBuilder:
    """Assemble CHARMM PSF + parameter files via ParmEd (no VMD)."""

    def __init__(self, workdir: str):
        self.workdir = workdir
        self.structure = None
        self.param_set = None
        self.status = 0

    def loadParameters(self, rtf_files: list, prm_files: list):
        """Load one or more CHARMM topology/parameter files into a ParmEd set."""
        args = list(rtf_files) + list(prm_files)
        self.param_set = pmd.charmm.CharmmParameterSet(*args)
        return self.param_set

    def loadStructure(self, pdb_file: str, rtf_files: list) -> pmd.Structure:
        """Build a ParmEd Structure from a PDB, typed against the supplied RTFs."""
        psf_like = pmd.charmm.CharmmPsfFile()  # empty, we'll populate from PDB
        pdb = pmd.load_file(pdb_file)
        # ParmEd doesn't have a direct 'apply RTF to PDB' path, so we merge the
        # PDB geometry into an empty Psf and then apply parameters.
        for atom in pdb.atoms:
            psf_like.add_atom(atom, resname=atom.residue.name,
                              resnum=atom.residue.number,
                              chain=atom.residue.chain or '')
        psf_like.coordinates = pdb.coordinates
        psf_like.box = pdb.box
        self.structure = psf_like
        return self.structure

    def addMetalBonds(self, metal_bonds: list):
        """Register metal-ligand bonds that wouldn't be picked up from RTF.

        metal_bonds: list of (atom_i_idx, atom_j_idx) 0-based indices.
        """
        if self.structure is None:
            return
        atoms = list(self.structure.atoms)
        for i, j in metal_bonds:
            self.structure.bonds.append(pmd.Bond(atoms[i], atoms[j]))

    def applyParameters(self):
        if self.structure is None or self.param_set is None:
            raise RuntimeError('Structure or parameter set not loaded')
        self.structure.load_parameters(self.param_set)

    def writePSF(self, basename: str = 'system') -> tuple:
        if self.structure is None:
            print(Color.RED + 'No structure to write' + Color.END)
            return '', ''
        psf = os.path.join(self.workdir, basename + '.psf')
        pdb = os.path.join(self.workdir, basename + '.pdb')
        self.structure.save(psf, format='psf', overwrite=True)
        self.structure.save(pdb, format='pdb', overwrite=True)
        self.status = 1
        return psf, pdb

    def attachSolvatedBox(self, solvated_pdb: str,
                          solvent_rtf_files: list,
                          solvent_prm_files: list,
                          basename: str = 'solvated') -> tuple:
        """Re-load a Packmol-assembled PDB, attach solvent parameters, and
        emit the final solvated PSF/PDB pair.
        """
        combined = pmd.load_file(solvated_pdb)
        all_params = pmd.charmm.CharmmParameterSet(
            *(list(solvent_rtf_files) + list(solvent_prm_files)))
        if self.param_set is not None:
            all_params.__iadd__(self.param_set)
        combined_psf = pmd.charmm.CharmmPsfFile()
        for atom in combined.atoms:
            combined_psf.add_atom(atom, resname=atom.residue.name,
                                  resnum=atom.residue.number,
                                  chain=atom.residue.chain or '')
        combined_psf.coordinates = combined.coordinates
        combined_psf.box = combined.box
        combined_psf.load_parameters(all_params)

        psf = os.path.join(self.workdir, basename + '.psf')
        pdb = os.path.join(self.workdir, basename + '.pdb')
        combined_psf.save(psf, format='psf', overwrite=True)
        combined_psf.save(pdb, format='pdb', overwrite=True)
        return psf, pdb
