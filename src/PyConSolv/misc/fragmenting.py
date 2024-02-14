import rdkit.Chem

from PyConSolv.misc.inputparser import XYZ
from rdkit import Chem
import os
from rdkit.Chem import rdDetermineBonds

class Fragmentor:
    def __init__(self, path: str= None, radius: float = 4):
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

        '''
        self.inputpath = path
        self.xyz = self.initializeXYZ(path)
        self.radius = radius
        self.prepareXYZ()
        self.keep = []
        self.rings = None
        self.hydrogenate_list = []


    def checkBreakPoint(self, maxringsize = 10):
        '''
        Check for clean break of structure and make sure not to cut across rings. The atoms connected to those
         within the radius are transformed into hydrogen and used as capping atoms
         :param int maxsize: maximum size of rings to be identified
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
        self.hydrogenate_list = [x for x in self.hydrogenate_list  if x not in self.keep]
        self.keep = list(set(list(self.hydrogenate_list) + self.keep))




    def cutStructure(self) -> str:
        '''
        Cut the coordinates to those of the substructure
        :return: XMOL structure block
        '''

        coords = ('''{}
Substructure extracted by pyconsolv with radius {}
''').format(len(self.keep),self.radius)
        with open(self.inputpath, 'r') as f:
            next(f)
            next(f)
            counter = 0
            for line in f:
                if counter in self.keep:
                    tmp = line.split()
                    if counter in self.hydrogenate_list:
                        tmp[0] = 'H'
                    coords+=' '.join(tmp) + '\n'
                counter += 1
        coords += '\n'
        return coords


    def writeXYZ(self, coords: str = None, filename = 'substructure.xyz'):
        '''
        Write out XMOL file for substructure
        :param coords: XMOL string to write
        :param filename: name of file to write into
        :return:
        '''
        with open('/'.join(self.inputpath.split('/')[:-1] + [filename]), 'w') as f:
            f.write(coords)

    def initializeXYZ(self, path) -> XYZ:
        '''
        Initialize XYZ Molecule
        :param path: path to input xyz file
        :return: XYZ molecule object
        '''
        xyz = XYZ(db_file=os.path.split(__file__)[0] + '/../db/atom-radius.txt',
                  db_metal_file=os.path.split(__file__)[0] + '/../db/metal-radius.txt')
        xyz.readXYZ(path)
        xyz.calculateDistanceMatrix()
        xyz.generateAdjacencyMatrix()
        xyz.generateLinkList()
        xyz.connectedCompponents()
        return xyz


    def capAtom(self, xyz: XYZ, Atom: str = 'H'): #deprecated, will be removed
        '''
        Cap broken bonds with an atom, default is Hydrogen
        '''
        coords = '{}\nSubstructure\n'.format(xyz.atoms.size)
        for i in range(xyz.atoms.size):
            coords += '{} {} {} {}\n'.format(xyz.atoms[i],*xyz.coords[i])

        structure = Chem.MolFromXYZBlock(coords)
        structure = Chem.RWMol(structure)
        self.addBonds(structure, bonds = xyz.linkList, metalbonds=xyz.metalBonds)
        structure_H = Chem.AddHs(structure,	explicitOnly = False, addCoords=True, onlyOnAtoms=[1])
        # structure_H = structure
        coords = ''
        for i, a in enumerate(structure_H.GetAtoms()):
            positions = structure_H.GetConformer().GetAtomPosition(i)
            if (positions.x == 0.0 and positions.y == 0.0 and positions.z == 0.0):
                continue
            coords += '{} {} {} {}\n'.format(a.GetSymbol(),  round(positions.x,3) , round(positions.y,3), round(positions.z,3))
        coords = '{}\n Substructure generated by PyConSolv\n'.format(i) + coords
        self.writeXYZ(coords)

    def parametrize(self):
        '''
        Perform parametrization of subunit
        :return:
        '''
        pass

    def prepareXYZ(self):
        '''
        performs necessary steps for bond detection //TODO make this work for all metals within the molecule
        :return:
        '''

        if self.xyz.linkList[0] == []:
            for id in self.xyz.metalBonds:
                self.xyz.linkList[0].append(int(id.split()[-1]))

            for atomid in self.xyz.linkList[0]:
                self.xyz.linkList[atomid].append(0)

    def checkRadius(self):
        '''
        Check which atoms are within the set radius
        :return:
        '''
        for i in range(self.xyz.Dmat.shape[1]):
            if self.xyz.Dmat[0][i] < self.radius:
                self.keep.append(i)
        inclusion_list = []
        for atom in self.keep:
            inclusion_list += self.xyz.linkList[atom]
        inclusion_list.sort()
        inclusion_list = list(set(inclusion_list))
        tmp = []
        for element in self.keep:
            if element in inclusion_list:
                tmp.append(element)
        self.keep = tmp

    def findRings(self, maxsize: int = 10):
        '''
        Find and mark rings within structure and keep them
        :param int maxsize: maximum size of rings to be identified
        :return:
        '''
        structure = Chem.MolFromXYZFile(self.inputpath)
        structure = Chem.RWMol(structure)
        self.addBonds(structure)
        self.rings = structure.GetRingInfo().AtomRings()
        self.rings = [x for x in self.rings if len(x) < maxsize]


    def addBonds(self, mol: Chem.Mol = None, limit = None, bonds = None, metalbonds = None):
        '''
        Add bonds to a RDKit.Mol structure created from xyz
        :param mol: rdkit.Chem.Mol structure
        :param limit: add only bonds involving chosen atom ids
        :param bonds: list of bonds
        :param bonds: list of bonds to metal
        :return:
        '''
        if not bonds:
            bonds = self.xyz.linkList
        bonded = []

        if not metalbonds:
            metalbonds = self.xyz.metalBonds
        for id in metalbonds:
            bonds[0].append(int(id.split()[-1]))
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
        Chem.SanitizeMol(mol)

    def run(self, maxringsize: int = 10):
        self.checkBreakPoint(maxringsize=maxringsize)
        coords = self.cutStructure()
        self.writeXYZ(coords)
