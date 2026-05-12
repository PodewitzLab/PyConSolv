from pathlib import Path

import numpy as np
import rdkit.Chem
from rdkit import Chem

from .inputparser import XYZ


class Fragmentor:
    def __init__(self, path: str = None, radius: float = 4.0):
        '''
        Class to fragment an input structure and perform parametrization on a subunit of the main structure.
        :param string path: input structure as an XMOL structure
        :param float radius: radius to find atoms which should be kept

         Class variables:
            - self.inputpath = path to input XMOL file
            - self.xyz = XYZ object to calculate the bonding of the structure
            - self.radius = radius within which atoms should be kept
            - self.keep = list of atoms which should be kept in the model system
            - self.rings = list of rings identified in structure
            - self.hydrogenate_list = list of added capping hydrogen atoms
            - self.metal_indices = list of indices of metal atoms in the structure

        '''
        self.inputpath = Path(path) if path else None
        self.xyz = self.initializeXYZ(path)
        self.radius = radius
        self.metal_indices = self.findMetalIndices()
        self.prepareXYZ()
        self.keep = []
        self.rings = None
        self.hydrogenate_list = []

    def findMetalIndices(self) -> list:
        """Find the indices of all metal atoms in the structure"""
        indices = []
        for i, atom in enumerate(self.xyz.atoms):
            if self.xyz.isMetal(str(atom)):
                indices.append(i)
        return indices

    def checkBreakPoint(self, maxringsize: int = 10):
        '''
        Check for clean break of structure and make sure not to cut across rings. The atoms connected to those
         within the radius are transformed into hydrogen and used as capping atoms
         :param int maxringsize: maximum size of rings to be identified
        '''
        self.checkRadius()
        self.findRings(maxsize=maxringsize)
        addedRings = []
        for ring in self.rings:
            addedRings.append([ring, False])

        complete = False
        while not complete:
            complete = True
            for ring in addedRings:
                if ring[1]:
                    continue
                else:
                    for element in ring[0]:
                        if element in self.keep:
                            ring[1] = True
                            self.keep = list(set(list(ring[0]) + self.keep))
                            complete = False
                            continue

        self.hydrogenate_list = []
        for atom in self.keep:
            self.hydrogenate_list += self.xyz.linkList[atom]
        self.hydrogenate_list = [x for x in self.hydrogenate_list if x not in self.keep]
        self.keep = list(set(list(self.hydrogenate_list) + self.keep))

    def cutStructure(self) -> str:
        '''
        Cut the coordinates to those of the substructure.
        Capping hydrogens are repositioned to proper bond distance from their parent atoms.
        :return: XMOL structure block
        '''
        # First, find the parent atom for each capping hydrogen
        # Parent is an atom in keep (but not in hydrogenate_list) that is bonded to the capping H
        capping_parents = {}  # maps hydrogenate_idx -> parent_idx
        for h_idx in self.hydrogenate_list:
            for bonded_idx in self.xyz.linkList[h_idx]:
                if bonded_idx in self.keep and bonded_idx not in self.hydrogenate_list:
                    capping_parents[h_idx] = bonded_idx
                    break

        # Typical X-H bond distance (Angstroms)
        H_BOND_DISTANCE = 1.09

        coords = ('''{}
Substructure extracted by pyconsolv with radius {}
''').format(len(self.keep), self.radius)

        with open(self.inputpath, 'r') as f:
            next(f)
            next(f)
            counter = 0
            for line in f:
                if line.split() == []:
                    continue
                if counter in self.keep:
                    tmp = line.split()
                    if counter in self.hydrogenate_list:
                        # This atom becomes a capping hydrogen
                        tmp[0] = 'H'

                        # Reposition the hydrogen to proper bond distance
                        if counter in capping_parents:
                            parent_idx = capping_parents[counter]

                            # Get coordinates
                            h_coords = np.array(self.xyz.coords[counter], dtype=float)
                            parent_coords = np.array(self.xyz.coords[parent_idx], dtype=float)

                            # Calculate direction vector from parent to H
                            direction = h_coords - parent_coords
                            distance = np.linalg.norm(direction)

                            if distance > 0:
                                # Normalize and scale to H bond distance
                                direction = direction / distance
                                new_h_coords = parent_coords + direction * H_BOND_DISTANCE

                                # Update coordinates
                                tmp[1] = f'{new_h_coords[0]:.6f}'
                                tmp[2] = f'{new_h_coords[1]:.6f}'
                                tmp[3] = f'{new_h_coords[2]:.6f}'

                    coords += ' '.join(tmp) + '\n'
                counter += 1
        coords += '\n'
        return coords

    def writeXYZ(self, coords: str = None, filename: str = 'substructure.xyz'):
        '''
        Write out XMOL file for substructure
        :param coords: XMOL string to write
        :param filename: name of file to write into
        :return:
        '''
        output_path = self.inputpath.parent / filename
        with open(output_path, 'w') as f:
            f.write(coords)

    def initializeXYZ(self, path: str) -> XYZ:
        '''
        Initialize XYZ Molecule
        :param path: path to input xyz file
        :return: XYZ molecule object
        '''
        db_path = Path(__file__).parent.parent / 'db'
        xyz = XYZ(db_file=str(db_path / 'atom-radius.txt'),
                  db_metal_file=str(db_path / 'metal-radius.txt'))
        xyz.readXYZ(path)
        xyz.calculateDistanceMatrix()
        xyz.generateAdjacencyMatrix()
        xyz.generateLinkList()
        xyz.connectedCompponents()
        return xyz

    def parametrize(self):
        '''
        Perform parametrization of subunit
        :return:
        '''
        pass

    def prepareXYZ(self):
        '''
        Performs necessary steps for bond detection for all metals in the structure.
        Populates linkList with metal-ligand bonds from metalBonds.
        :return:
        '''
        # Process each metal in the structure
        for metal_idx in self.metal_indices:
            if self.xyz.linkList[metal_idx] == []:
                # Find bonds for this specific metal from metalBonds
                for bond_info in self.xyz.metalBonds:
                    parts = bond_info.split()
                    if len(parts) >= 3:
                        bond_metal_idx = int(parts[0])
                        ligand_idx = int(parts[-1])
                        # Only add if this bond involves the current metal
                        if bond_metal_idx == metal_idx:
                            if ligand_idx not in self.xyz.linkList[metal_idx]:
                                self.xyz.linkList[metal_idx].append(ligand_idx)
                            if metal_idx not in self.xyz.linkList[ligand_idx]:
                                self.xyz.linkList[ligand_idx].append(metal_idx)

    def checkRadius(self):
        '''
        Check which atoms are within the set radius of ANY metal in the structure.
        This ensures all metals and their coordination spheres are included.
        :return:
        '''
        self.keep = []

        if not self.metal_indices:
            # No metals found, keep nothing
            return

        # For each atom, check if it's within radius of ANY metal
        for i in range(self.xyz.Dmat.shape[1]):
            for metal_idx in self.metal_indices:
                if self.xyz.Dmat[metal_idx][i] < self.radius:
                    if i not in self.keep:
                        self.keep.append(i)
                    break  # No need to check other metals once atom is included

    def findRings(self, maxsize: int = 10):
        '''
        Find and mark rings within structure and keep them
        :param int maxsize: maximum size of rings to be identified
        :return:
        '''
        structure = Chem.MolFromXYZFile(str(self.inputpath))
        structure = Chem.RWMol(structure)
        self.addBonds(structure)
        self.rings = structure.GetRingInfo().AtomRings()
        self.rings = [x for x in self.rings if len(x) < maxsize]

    def addBonds(self, mol: Chem.Mol = None, limit=None, bonds=None, metalbonds=None):
        '''
        Add bonds to a RDKit.Mol structure created from xyz.
        Handles multiple metals properly.
        :param mol: rdkit.Chem.Mol structure
        :param limit: add only bonds involving chosen atom ids
        :param bonds: list of bonds
        :param metalbonds: list of bonds to metal
        :return:
        '''
        if not bonds:
            bonds = [list(b) for b in self.xyz.linkList]  # Deep copy to avoid modifying original
        bonded = []

        if not metalbonds:
            metalbonds = self.xyz.metalBonds

        # Add metal bonds to the bonds list for each metal
        for bond_info in metalbonds:
            parts = bond_info.split()
            if len(parts) >= 3:
                metal_idx = int(parts[0])
                ligand_idx = int(parts[-1])
                if ligand_idx not in bonds[metal_idx]:
                    bonds[metal_idx].append(ligand_idx)

        for i in range(len(bonds)):
            if len(bonds[i]) > 0:
                for element in bonds[i]:
                    if limit:
                        if element not in limit or i not in limit:
                            continue
                    if [element, i] in bonded or [i, element] in bonded:
                        continue
                    else:
                        mol.AddBond(i, element, order=rdkit.Chem.rdchem.BondType.SINGLE)
                        bonded.append([i, element])
        try:
            Chem.SanitizeMol(mol)
        except:
            pass  # Sanitization may fail for metal complexes

    def run(self, maxringsize: int = 10, filename: str = 'substructure.xyz') -> bool:
        '''
        Run the fragmentation process
        :param int maxringsize: maximum size of rings to be identified
        :param str filename: name of output file
        :return: True on success, False on failure
        '''
        try:
            if not self.metal_indices:
                print('Warning: No metals found in structure')

            print(f'Found {len(self.metal_indices)} metal(s) at indices: {self.metal_indices}')
            self.checkBreakPoint(maxringsize=maxringsize)
            coords = self.cutStructure()
            self.writeXYZ(coords, filename=filename)
            return True
        except Exception as e:
            print(f'Fragmentation failed: {e}')
            return False
