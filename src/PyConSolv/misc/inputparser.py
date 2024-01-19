import os.path
import shutil

import numpy as np
import pandas as pd

from ..utils.charge import ChargeChanger


class XYZ:
    def __init__(self, db_file: str, db_metal_file: str):
        """
        Initialize database of ionic radii from db_file

        Parameters:
            :param string db_file: file that contains the ionic radii for all elements
            :param string db_metal_file: file that contains the ionic radii for metals

        Class variables:
            - self.map = contains the mapping for atom types and id
            - self.molNotCreated = list of fragments for which antechamber needs to run for the creation of .mol2 files
            - self.charges = charges for each fragment, read from chargeMap.dat
            - self.path - path where pdb files are generated
            - self.hasMetal - True if structure contains a metal, False otherwise
            - self.files - List of lists containing PDB data for each molecule
            - self.connected - list which stores elements of each subgraph
            - self.chain - [x] length array which contains the chain number to which each atom belongs
            - self.linkList - list containing lists of partners for each atom
            - self.Adjmat - [x,x] np array which contains information about bonding; 1 means bonded, 0 means unbonded
            - self.Dmat - [x,x] shaped numpy array which contains the pair distances
            - self.Amat - [x,x] shaped numpy array which contains the pair elements
            - self.atoms = np array containing the atom labels
            - self.coords = np array containing the atom coordinates
            - self.db = dictionary which contains element names and corresponding ionic radius
            - self.metalList = list of all metal (as suppoterd by MCPB.py) atomic symbols, used to cross-reference
                                whether element is metal or not
            - self.metals = list of metals detected in the xyz file and the corresponding filename
            - self.ligands = list of non-metal ligands detected in the xyz file and the corresponding filename
            - self.filenames = contains the names of the PDB files which are written out
            - self.metalBonds = contains the atoms bonded to the metal center
            - self.metalRadius = contains the radius of the metal ion
        """

        self.map = None
        self.molNotCreated = None
        self.charges = None
        self.hasMetal = False
        self.path = None
        self.files = None
        self.chain = None
        self.connected = None
        self.linkList = None
        self.Adjmat = None
        self.Amat = None
        self.Dmat = None
        self.atoms = None
        self.coords = None
        df = pd.read_csv(db_file, delimiter='\t')
        self.db = {}
        for A, B in zip(df.values[:, 1], df.values[:, 5]):
            self.db[A] = B

        df = pd.read_csv(db_metal_file, sep='\t')
        df = df.drop(columns=['Charge', 'Crystal Radius'])
        df['Ion'] = df['Ion'].apply(str.upper)
        df = df.groupby(['Ion']).max()
        self.metalRadius = df.to_dict().get('Ionic Radius')

        self.metalList = ['LI', 'BE', 'NA', 'MG', 'AL', 'SI', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE',
                          'CO', 'NI', 'CU', 'ZN',
                          'GA', 'GE', 'AS', 'SE', 'BR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG',
                          'CD', 'IN', 'SN', 'SB',
                          'TE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER',
                          'TM', 'YB', 'YB', 'LU', 'HF',
                          'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'FR', 'RA', 'AC',
                          'TH', 'PA', 'U', 'NP',
                          'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO',
                          'LR']  # Metals supported by MCPB, F and CL were moved
        self.metals = []
        self.ligands = []
        self.filenames = []
        self.metalBonds = []

    def readXYZ(self, path: str):
        """
        Read XYZ (XMOL format) file from path

        Parameters:
            :param string path: location of XYZ file

        Class variables:
            - self.coords - [x,3] sized array which contains the coordinates of all atoms in the XYZ file
            - self.atoms - [x] length array which contains atom element labels
        """
        self.coords = []
        self.atoms = []
        f = open(path, 'r')
        counter = 0
        for line in f:
            if counter < 2:
                counter = counter + 1
                continue
            else:
                self.atoms.append(line.split()[0])
                self.coords.append(line.split()[1:])
        f.close()
        self.coords = np.asarray(self.coords, dtype=float)
        self.atoms = np.asarray(self.atoms)

    def calculateDistanceMatrix(self):
        """
        Calculate distances between each atom pair

        Parameters:

        Class variables:
        """
        self.Dmat = np.empty((self.coords.shape[0], self.coords.shape[0]))
        self.Amat = np.empty((self.coords.shape[0], self.coords.shape[0]),
                             dtype='<U5')  # <U is better for string formatting due to less compatibility issues
        for i in range(self.coords.shape[0]):
            for j in range(self.coords.shape[0]):
                self.Dmat[i][j] = np.linalg.norm(self.coords[i] - self.coords[j])
                self.Amat[i][j] = self.atoms[i] + '-' + self.atoms[j]

    def generateAdjacencyMatrix(self):
        """
        Check whether atoms are connected to eachother; atoms are considered bonded, based on ionic radius from
        http://crystalmaker.com/support/tutorials/atomic-radii/index.html and
        http://abulafia.mt.ic.ac.uk/shannon/radius.php#R
        Bonds to metal are recorded in self.metalBonds. This is needed to check if all the residues are correctly found
        in MCPB.py, as Metal-C bonds are not autodetected.

        Parameters:

        Class variables:
        """

        self.Adjmat = np.zeros(self.Dmat.shape)
        for i in range(self.Dmat.shape[0]):
            for j in range(self.Dmat.shape[1]):
                met = False
                a = self.Amat[i][j].split('-')
                dist = self.db[a[0]] + self.db[a[1]]
                if i == j:
                    continue
                if self.Dmat[i][j] <= dist * 0.6:
                    self.Adjmat[i][j] = 1
                if self.isMetal(a[0]) and self.isMetal(a[1]):
                    dist = self.metalRadius.get(a[0].upper()) + self.metalRadius.get(a[1].upper())
                    met = True
                elif self.isMetal(a[0]):
                    dist = self.metalRadius.get(a[0].upper()) + self.db[a[1]]
                    met = True
                if met:
                    self.Adjmat[i][j] = 0
                    if self.Dmat[i][j] <= dist * 1.0:
                        self.metalBonds.append('{} @{}{} {}'.format(str(i), a[1], str(j), str(j)))

    def generateLinkList(self):
        """
        Check whether atoms are connected to eachother; atoms are considered bonded, based on ionic radius from
        http://crystalmaker.com/support/tutorials/atomic-radii/index.html

        Parameters:

        Class variables:
        """

        self.linkList = []
        for i in range(self.Dmat.shape[0]):
            temp = []
            for j in range(self.Dmat.shape[1]):
                a = self.Amat[i][j].split('-')
                dist = self.db[a[0]] + self.db[a[1]]
                if i == j:
                    continue
                if self.Dmat[i][j] <= dist * 0.6:
                    if self.isMetal(a[0]) or self.isMetal(a[1]):
                        continue
                    else:
                        temp.append(j)
            self.linkList.append(temp)

    def connectedCompponents(self):
        """
        DFS algorithm for finding components

        Parameters:

        Class variables:
        """
        visited = []
        self.connected = []
        for i in range(len(self.linkList)):
            visited.append(0)
        for i in range(len(self.linkList)):
            if visited[i] == 0:
                temp = []
                self.connected.append(self.DFSaux(i, temp, visited))
        for elem in self.connected:
            elem.sort()

    def DFSaux(self, vertex, temp, visited):
        """
        Auxiliary helper function to perform DFS

        Parameters:

        Class variables:

        """
        visited[vertex] = 1
        temp.append(vertex)
        for i in self.linkList[vertex]:
            if visited[i] == 0:
                temp = self.DFSaux(i, temp, visited)
        return temp

    def assignChain(self):
        """
        Assign each atom to a chain

        Parameters:

        Class variables:
            - self.chain - contains the chain ID of each atom
        """
        temp = []
        for i in range(len(self.connected)):
            for j in range(len(self.connected[i])):
                temp.append([self.connected[i][j], i])
        self.chain = np.array(temp)
        self.chain = self.chain[np.lexsort((self.chain[:, 0], self.chain[:, 1]))]

    def createPDB(self):
        """
        Create a file following PDB format standards

        Parameters:

        Class variables:
        """
        # todo clean up this function and use string.format() method

        self.files = []

        for e in range(len(self.connected)):
            file = []
            atoms = len(self.connected[e])
            metal = False

            for i in range(atoms):

                atom_pos = self.connected[e][i]
                if self.isMetal(self.atoms[atom_pos]):
                    metal = True

                line = ''
                atom = 'HETATM'
                line = line + atom  # Col 1-6 Identifier HETATM

                atomSN = [' ', ' ', ' ', ' ', ' ']  # col 7-11 Atom Serial Number
                SN = str(i + 1)
                if len(SN) > 5:
                    print('Serial Number is too long for atom %d!' % i)
                for j in range(len(SN)):
                    atomSN[4 - j] = SN[-(j + 1)]
                string = ''
                string = string.join(atomSN)
                line = line + string

                line = line + ' '  # Col 12 is empty

                if len(self.atoms[atom_pos]) > 4:
                    print('Atom name is too long for atom %d!' % i)
                atomName = [' ', ' ', ' ', ' ']  # col 13-16 Atom Name
                for j in range(len(self.atoms[atom_pos])):
                    atomName[j] = self.atoms[atom_pos][j]
                string = ''
                string = string.join(atomName)
                line = line + string

                line = line + ' '  # Col 17 is alternate location indicator

                line = line + ' ' + ' ' + chr(65 + e)  # Col 18-20 Residue Name

                line = line + ' '  # Col 21 is empty

                line = line + str(e + 1)  # Col 22 Chain ID

                line = line + ' ' + ' ' + ' ' + str(e + 1)  # col 23-26 Residue sequence Number - Integer

                line = line + ' '  # Col 27 for insertions of residues

                line = line + ' ' + ' ' + ' '  # Col 28-30 are empty

                x = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']  # col 31-38 Residue sequence Number - real (8.3)
                x_data = str(self.coords[atom_pos][0].round(3))
                if len(x_data) == 0:
                    print('x data is missing for atom %d!' % i)
                for j in range(len(x_data)):
                    x[7 - j] = x_data[-(j + 1)]
                string = ''
                string = string.join(x)
                line = line + string

                y = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']  # col 39-46 Residue sequence Number - real (8.3)
                y_data = str(self.coords[atom_pos][1].round(3))
                if len(y_data) == 0:
                    print('y data is missing for atom %d!' % i)
                for j in range(len(y_data)):
                    y[7 - j] = y_data[-(j + 1)]
                string = ''
                string = string.join(y)
                line = line + string

                z = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']  # col 47-54 Residue sequence Number - real (8.3)
                z_data = str(self.coords[atom_pos][2].round(3))
                if len(z_data) == 0:
                    print('z data is missing for atom %d!' % i)
                for j in range(len(z_data)):
                    z[7 - j] = z_data[-(j + 1)]
                string = ''
                string = string.join(z)
                line = line + string

                oc = [' ', ' ', '1', '.', '0', '0']  # col 55-60 Occupancy - real (6.2)
                string = ''
                string = string.join(oc)
                line = line + string

                tf = [' ', ' ', '0', '.', '0', '0']  # col 61-66 Temperature Factor - real (6.2)
                string = ''
                string = string.join(tf)
                line = line + string

                line = line + ' ' + ' ' + ' ' + ' ' + ' ' + ' '  # Col 67-72 are empty

                line = line + ' ' + ' ' + ' ' + ' '  # Col 73-76 Segment identifier

                if len(self.atoms[atom_pos]) > 2:
                    print('Element Symbol is too long for atom %d!' % i)
                es = [' ', ' ']  # col 77-78 Element Symbol
                for j in range(len(self.atoms[atom_pos])):
                    es[1 - j] = self.atoms[atom_pos][-(j + 1)]
                string = ''
                string = string.join(es)
                line = line + string

                line = line + '\n'

                file.append(line)
            file.append('TER\n')
            if metal:
                self.files = [file] + self.files
            else:
                self.files.append(file)

    def writePDBFiles(self, path: str):
        """
        Write out PDB files

        Parameters:
            :param string path: file output path

        Class variables:
        """

        self.path = path + '/'
        self.hasMetal = False
        chrstart = 65
        atomid = 0
        for i in range(len(self.files)):
            if len(self.files[i]) == 2:  # single atom ligands need to be treated differently to work with MCPB.py
                name = self.files[i][0].split()[-1].upper()

                if self.isMetal(name):
                    if name == 'B':
                        chrstart = 67
                    s = self.files[i][0][:7] + '{:>4}'.format(1) + ' ' + '{:<5}'.format(name) + '{:>3}'.format(name) + \
                        self.files[i][0][20:]
                    if [name + '.pdb', name] in self.metals:
                        self.filenames.append(name + str(len(self.metals)) +'.pdb')
                    else:
                        self.filenames.append(name + '.pdb')
                    self.hasMetal = True
                    s = s.upper()
                    self.metals.append([name + '.pdb', name])
                else:
                    s = self.files[i][0][:7] + '{:>4}'.format(1) + ' ' + '{:<5}'.format(name + str(atomid)) + '{:>3}'.format(chr(chrstart + i)) + \
                        self.files[i][0][20:]
                    atomid += 1
                    self.filenames.append(chr(chrstart + i) + '.pdb')
                    name = chr(chrstart + i)
                    self.ligands.append([chr(chrstart + i) + '.pdb', chr(chrstart + i)])
                f = open(self.path + '/' + name + '.pdb', 'w')
                f.write(s)
                f.write('TER')
                f.close()

            else:
                self.filenames.append(chr(chrstart + i) + '.pdb')
                f = open(path + '/' + chr(chrstart + i) + '.pdb', 'w')
                for line in self.files[i]:
                    if 'TER' in line:
                        f.write('END')
                    else:
                        atomid += 1
                        f.write(line[:7] + '{:>4}'.format(atomid) + ' ' + '{:<5}'.format(
                            line[11:16].strip() + str(atomid)) + line[17:])
                f.close()
                self.ligands.append([chr(chrstart + i) + '.pdb', chr(chrstart + i)])
        tmp = []
        atomid = 0
        for i in range(len(self.filenames)):
            f = open(path + '/' + self.filenames[i], 'r')
            for line in f:
                if 'TER' in line:
                    continue
                elif 'END' in line:
                    continue
                else:
                    atomid += 1
                    tmp.append(line[:7] + '{:>4}'.format(atomid) + line[11:])
            f.close()
        f = open(path + '/' + 'Full_PDB.pdb', 'w')
        for i in range(len(tmp)):
            f.write(tmp[i])
        f.write('END')
        f.close()
        self.writeFilenames(path)

    def writeMol2Files(self):
        """
        Write out mol2 files

        Parameters:

        Class variables:
        """

        print('Reading charges from {}\n'.format(self.path + 'chargeMap.dat'))

        self.charges = []
        self.molNotCreated = []
        atomindex = 0
        f = open(self.path + '/chargeMap.dat', 'r')
        for line in f:
            self.charges.append(line.split())
        f.close()
        for i in range(len(self.files)):

            print(self.charges[i])
            if len(self.files[i]) == 2:  # single atom ligands need to be treated differently to work with MCPB.py
                name = self.files[i][0].split()[-1]
                atomlabel = self.files[i][0].split()[2].upper()
                if not self.isMetal(name):
                    atomtype = name
                    name = self.charges[i][0]
                    atomlabel = atomlabel + str(atomindex-1)
                    metal = False
                else:
                    metal = True
                    name = name.upper()
                    atomtype = name.upper()
                    print('Metal atom found: {}, assigning charge {}'.format(name, self.charges[i][1]))
                coords = self.files[i][0].split()[6:9]
                self.createMol2_new(name, coords=coords, charge=float(self.charges[i][1]),
                                    index='', atomtype = atomtype, atomlabel = atomlabel)  # this index should allow for the parsing multiple metals and fix a bug with the ID not being found for mcpb
            else:
                self.molNotCreated.append(self.charges[i])
            atomindex = atomindex + len(self.files[i]) - 1

    def createMol2_new(self, name: str, coords: list, charge: float = -1, index: str = '', atomtype: str = '', atomlabel: str = ''):
        """
        Write out mol2 file for metal and single atom ligands

        Parameters:
            :param string name: name of the mol2 file
            :param list coords: coordinates of atoms
            :param float charge: charge of atom
            :param string index: ID of residue
            :param string atomtype: Alternate label for atom type

        Class variables:
        """

        layout = '''@<TRIPOS>MOLECULE
{:<3}
1 0 1 0 0
SMALL
USER_CHARGES


@<TRIPOS>ATOM
1 {:<2}        {:>.4f}    {:>.4f}    {:>.4f}   {:<2}    1  {:>2} {:>.4f}
@<TRIPOS>BOND
@<TRIPOS>SUBSTRUCTURE
  1 {:<2}    1 TEMP 0 **** **** 0 ROOT
'''
        print('Writing out {} ... \n'.format(self.path + name + index + '.mol2'))
        f = open(self.path + name + index + '.mol2', 'w')
        f.write(
            layout.format(name, atomlabel, float(coords[0]), float(coords[1]), float(coords[2]),
                          atomtype + index, name, charge, name + index))
        f.close()

    def createFinalMol2(self, path: str):
        """
        Create final mol2 files for tleap, using the charges from the RESP calculation

        Parameters:
            :param string path: folder where files are located/created
        Class variables:

        """
        print('Changing charges')
        if self.filenames == []:
            self.readFilenames(path + '/MCPB_setup')
        self.readRESP(path + '/orca_calculations')
        self.readMAP(path + '/MCPB_setup')
        self.chargeChanger = ChargeChanger()
        for name in self.filenames:
            fname = name.split('.')[0]
            fin = self.path + '/' + fname + '.mol2'
            fout = self.path + '/' + fname + '1.mol2'
            self.chargeChanger.change(fin, fout, fname + '1', charges = self.charges, fragmented = True, mapfile = self.map)


    def createFinalMol2_old(self, path: str):
        """
        Create final mol2 files for tleap, using the charges from the RESP calculation

        Parameters:
            :param string path: folder where files are located/created
        Class variables:

        """

        if self.filenames == []:
            self.readFilenames(path + '/MCPB_setup')
        self.readRESP(path + '/orca_calculations')
        self.readMAP(path + '/MCPB_setup')
        iterator = 0  # keeps track of current atom id
        for name in self.filenames:
            switch = 0
            fname = name.split('.')[0]
            fin = open(self.path + '/' + fname + '.mol2', 'r')  # read mol2 files corresponding to the pdb files
            fout = open(self.path + '/' + fname + '1.mol2', 'w')  # create new mol2 file
            for line in fin:
                if 'USER_CHARGES' in line:
                    fout.write('RESP Charge\n')
                elif '@<TRIPOS>ATOM' in line:
                    switch = 1
                    fout.write(line)
                elif '@<TRIPOS>BOND' in line:
                    fout.write(line)
                    switch = 2
                elif '@<TRIPOS>SUBSTRUCTURE' in line:
                    fout.write(line)
                    switch = 3
                elif '@<TRIPOS>MOLECULE' in line:
                    fout.write(line)
                    switch = 4
                else:
                    if switch == 1:

                        tmp_line = line.split()[:-1]  # remove old charge

                        tmp_line.append(float(self.charges[iterator][-1]))  # add resp charge
                        tmp_line[-2] = fname + '1'
                        tmp_line[5] = self.map[iterator][-1]

                        fout.write(
                            '{:>7} {:<7}    {:>7}    {:>7}    {:>7}   {:<2}    {}  {:>2} {:>9.6f}\n'.format(*tmp_line))
                        iterator += 1
                    elif switch == 3:
                        fout.write('     1 {}         1 TEMP              0 ****  ****    0 ROOT\n'.format(fname + '1'))
                        switch = 0
                    elif switch == 4:
                        fout.write(fname + '1\n')
                        switch = 0
                    else:
                        if 'bcc' in line:
                            fout.write(line.replace('bcc', 'RESP Charge'))
                        else:
                            fout.write(line.replace('ar', '1'))
            fout.close()
            fin.close()

    def readRESP(self, path: str):
        """
        Read RESP charges generated by Multiwfn and saved in orca_freq.molden.chg

        Parameters:
            :param string path: path to where charge file is located

        Class variables:
            - self.charges = contains the RESP charges generated by Multiwfn
        """
        if self.hasMetal:
            if not os.path.exists(path + '/freq/orca_freq.molden.chg'):
                shutil.copyfile(path + '/freq/orca_freq.chg',path + '/freq/orca_freq.molden.chg')
            filepath = path + '/freq/orca_freq.molden.chg'
        else:
            if not os.path.exists(path + '/opt/orca_opt.molden.chg'):
                shutil.copyfile(path + '/opt/orca_opt.chg',path + '/opt/orca_opt.molden.chg')
            filepath = path + '/opt/orca_opt.molden.chg'
        f = open(filepath, 'r')
        self.charges = []
        iterator = 0
        for line in f:
            iterator += 1
            self.charges.append([iterator, line.split()[0], line.split()[-1]])
        f.close()

    def readMAP(self, path: str):
        """
        Read Map generated by MCPB for atom types from fingerprint

        Parameters:
            :param string path: path to where map file is located

        Class variables:
        """
        f = open(path + '/LIG_standard.fingerprint', 'r')
        self.map = []
        iterator = 0
        for line in f:
            iterator += 1
            self.map.append([iterator, line.split()[0], line.split()[-1]])
        f.close()

    def isMetal(self, string: str) -> bool:
        """
        checks if given string is a metal element or not

        Parameters:
            :param string string: string to check

        Class variables:

        Returns:
            - Bool - True if string is metal, False if it is not
        """
        if string.upper() in self.metalList:
            return True
        else:
            return False

    def prepareMCPB(self, inputfile: str):
        """
        Runs main PDB generation functionality of the XYZ class.

        Parameters:
            :param string inputfile: full path to xyz inputfile

        Class variables:
        """
        path = '/'.join(inputfile.split('/')[:-1])
        self.readXYZ(inputfile)
        self.calculateDistanceMatrix()
        self.generateAdjacencyMatrix()
        self.generateLinkList()
        self.connectedCompponents()
        self.assignChain()
        self.createPDB()
        self.writePDBFiles(path)

    def remakeXYZ(self, path: str):
        """
        Creates new input.xyz file, with atom ordering suitable for further use with MCPB.py

        Parameters:
            :param string inputfile: full path to xyz inputfile

        Class variables:
        """
        if not self.files:
            print('Could not create file, no connectivity available\n')
            return
        tmp = [atomid for chain in self.connected for atomid in chain]
        try:
            f = open(path + '/input.xyz', 'w')
            f.write(str(len(tmp)) + '\n')
            f.write('Modified XYZ file for MCPB.py - PyConSolv\n')
            for chain in self.files:
                for el in chain:
                    if 'TER' in el:
                        continue
                    elif 'END' in el:
                        continue
                    else:
                        f.write('{:<2} {:>.8f} {:>.8f} {:>.8f}\n'.format(el.split()[2], float(el.split()[6]),
                                                                         float(el.split()[7]), float(el.split()[8])))
            f.close()
            print('Created new optimized input.xyz file')
        except:
            print('Could not create new XYZ file\n')
            return

    def prepareInput(self, inputfile: str):
        """
        Analyze connectivity and create a structure properly organized for MCPB.py
        The atom names will be re-ordered and a map file will be created

        Parameters:
            :param string inputfile: full path to xyz inputfile

        Class variables:
        """
        path = '/'.join(inputfile.split('/')[:-1])
        try:
            shutil.copyfile(inputfile, path + '/input.xyz.original')
        except:
            print('could not back up original file')
        self.readXYZ(inputfile)
        self.calculateDistanceMatrix()
        self.generateAdjacencyMatrix()
        self.generateLinkList()
        self.connectedCompponents()
        self.createPDB()
        self.remakeXYZ(path)

    def readFilenames(self, path: str):
        """
        Reads filenames from file. This is useful for restarting calculations

        Parameters:
            :param string path: full path to the folder containing the filenames.restart file

        Class variables:
        """
        f = open(path + '/filenames.restart', 'r')
        for line in f:
            self.filenames.append(line)
        f.close()

    def writeFilenames(self, path: str):
        """
        writes filenames from file. This is useful for restarting calculations

        Parameters:
            :param string path: full path to the folder containing the filenames.restart file

        Class variables:
        """
        f = open(path + '/filenames.restart', 'w')
        for filename in self.filenames:
            f.write(filename + '\n')
        f.close()

    def writeMetalConnections(self, path: str):
        """
        Writes metal-atom connections to file. This is useful for restarting calculations

        Parameters:
            :param string path: full path to the folder containing the metalConnections file

        Class variables:
        """
        f = open(path + '/metalConnections', 'w')
        for bond in self.metalBonds:
            f.write(bond + '\n')
        f.close()

    def writeConnections(self, path: str):
        """
        Writes initial atom ids for each chain to file. This is useful for restarting calculations

        Parameters:
            :param string path: full path to the folder containing the Connections file

        Class variables:
        """
        bondlen = 0
        f = open(path + '/Connections_old', 'w')
        f2 = open(path + '/Connections', 'w')
        for bond in self.connected:
            f.write(' '.join(str(x) for x in bond) + '\n')
            f2.write(' '.join(str(x) for x in range(0 + bondlen, len(bond) + bondlen)) + '\n')
            bondlen = bondlen + len(bond)
        f.close()
        f2.close()
