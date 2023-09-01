import numpy as np


class ParameterChecker:

    def __init__(self, folder):
        """
        Check the MCPB.py generated parameters and clean the frcmod file

        Parameters:
            - folder = folder that contains the LIG_frcmod, input.xyz and LIG_standard files

        Class variables:
            - self.frcmodfile = path to LIG_frcmod file
            - self.xyzfile = path to input.xyz file
            - self.fingerprintfile = path to fingerprintfile
            - self.xyz = xyz coordinates of system
            - self.bonds = bonds in frcmod file
            - self.angles = angles in frcmod file
            - self.dihedrals = dihedrals in frcmod file
            - self.atoms = dictionary of new atom types added by MCPB.py
            - self.frcmod = frcmodfile template
        """
        self.frcmodfile = folder + '/LIG_mcpbpy.frcmod'
        self.xyzfile = folder + '/input.xyz'
        self.fingerprintfile = folder + '/LIG_standard.fingerprint'
        self.xyz = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.atoms = {}
        self.frcmod = []

    def checkBonds(self):
        """
        Check if bonds are ok

        Parameters:

        Class variables:
        """
        for el in self.bonds:
            if 'M' in el[0] and 'Y' in el[0]:
                atoms = el[0].split('-')
                p = np.zeros((2, 3))
                p[0] = np.array(self.xyz[self.atoms.get(atoms[0])][1:], dtype=float)
                p[1] = np.array(self.xyz[self.atoms.get(atoms[1])][1:], dtype=float)
                value = self.bondCalculator(p)
                self.replaceValue(el[0], value)

    def checkAngles(self):
        """
        Check if angles are ok

        Parameters:

        Class variables:
        """
        for el in self.angles:
            atoms = el[0].split('-')
            if 'M' in atoms[0] or 'Y' in atoms[0]:
                if 'M' in atoms[1] or 'Y' in atoms[1]:
                    if 'M' in atoms[2] or 'Y' in atoms[2]:
                        p = np.zeros((3, 3))
                        p[0] = np.array(self.xyz[self.atoms.get(atoms[0])][1:], dtype=float)
                        p[1] = np.array(self.xyz[self.atoms.get(atoms[1])][1:], dtype=float)
                        p[2] = np.array(self.xyz[self.atoms.get(atoms[2])][1:], dtype=float)
                        value = self.angleCalculator(p)
                        self.replaceValue(el[0], value)

    def checkDihedrals(self):
        """
        Check if dihedrals are ok

        Parameters:

        Class variables:
        """
        for el in self.dihedrals:
            atoms = el[0].split('-')
            p = np.zeros((3, 3))
            p[0] = np.array(self.xyz[self.atoms.get(atoms[0])][1:], dtype=float)
            p[1] = np.array(self.xyz[self.atoms.get(atoms[1])][1:], dtype=float)
            p[2] = np.array(self.xyz[self.atoms.get(atoms[2])][1:], dtype=float)
            p[3] = np.array(self.xyz[self.atoms.get(atoms[3])][1:], dtype=float)
            value = self.bondCalculator(p)
            self.replaceValue(el[0], value)

    def readFrcmod(self):
        """
        Read frcmod file

        Parameters:

        Class variables:
        """
        f = open(self.frcmodfile, 'r')
        switch = 0
        for line in f:
            self.frcmod.append(line)
            if len(line.split()) == 0:
                continue
            if 'BOND' in line:
                switch = 1
                continue
            if 'ANGL' in line:
                switch = 2
                continue
            if 'DIHE' in line:
                switch = 3
                continue
            if switch == 0:
                continue
            elif switch == 1:
                self.bonds.append(line.split())
            elif switch == 2:
                self.angles.append(line.split())
            elif switch == 3:
                self.dihedrals.append(line.split())
        f.close()

    def readXYZ(self):
        """
        Read input.xyz file

        Parameters:

        Class variables:
        """
        f = open(self.xyzfile, 'r')
        next(f)
        next(f)
        for line in f:
            if len(line.split()) == 4:
                self.xyz.append(line.split())
        f.close()

    def readFingerprint(self):
        """
        Read fingerprint file

        Parameters:

        Class variables:
        """
        f = open(self.fingerprintfile, 'r')
        for line in f:
            if '-> M' in line:
                self.atoms[line.split()[-1]] = int(line.split()[1]) - 1
            if '-> Y' in line:
                self.atoms[line.split()[-1]] = int(line.split()[1]) - 1
        f.close()

    def replaceValue(self, key, value):
        """
        Replace value in frcmod file.

        Parameters:
            :param list<String> key - marks the bond, angle or dihedral for which the value should be replaced
            :param float value - value which should be replaced
        Class variables:
        """
        for i in range(len(self.frcmod)):
            if len(self.frcmod[i].split()) > 2:
                if key == self.frcmod[i].split()[0]:
                    l = self.frcmod[i].split()
                    l[2] = str(value)
                    self.frcmod[i] = '{:6<}   {:7>}    {:7>}    Created by Seminario method using MCPB.py\n'.format(
                        l[0], l[1], l[2])

    def removeUseless(self):
        """
        Remove lines which make no sense in the frcmod file

        Parameters:

        Class variables:
        """
        tmp = []
        for i in range(len(self.frcmod)):
            if len(self.frcmod[i].split()) > 2:
                key = self.frcmod[i].split()[0].split('-')
                count = 0
                for j in range(len(key)):
                    if 'M' in key[j] or 'Y' in key[j]:
                        count += 1
                if len(key) == 1:
                    tmp.append(self.frcmod[i])
                elif count != len(key):
                    tmp.append(self.frcmod[i])
                elif len(key) == len(set(key)):
                    tmp.append(self.frcmod[i])
                else:
                    print(key, set(key))
            else:
                tmp.append(self.frcmod[i])
        self.frcmod = tmp

    def writeFrcmod(self):
        """
        Write new frcmod file

        Parameters:

        Class variables:
        """
        f = open(self.frcmodfile, 'w')
        for el in self.frcmod:
            f.write(el)
        f.close()

    def bondCalculator(self, p):
        """
        Calculate distance between 2 atoms

        Parameters:
            :param np.array(2,3) p - array containing the coordinates for the 2 atoms

        Class variables:
        """
        return np.round(abs(np.linalg.norm(p[0] - p[1])), 4)

    def angleCalculator(self, p):
        """
        Calculate angle between 3 atoms

        Parameters:
           :param np.array(3,3) p - array containing the coordinates for the 3 atoms

        Class variables:

        Returns:
            float angle
        """
        ba = p[0] - p[1]
        bc = p[2] - p[1]
        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        return np.round(np.degrees(angle), 2)

    def dihedralCalculator(self, p):
        """
        Calculate dihedral between 4 atoms

        Parameters:
           :param np.array(4,3) p - array containing the coordinates for the 4 atoms

        Class variables:
        """
        b = p[:-1] - p[1:]
        b[0] *= -1
        v = np.array([v - (v.dot(b[1]) / b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])
        # Normalize vectors
        v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1, 1)
        b1 = b[1] / np.linalg.norm(b[1])
        x = np.dot(v[0], v[1])
        m = np.cross(v[0], b1)
        y = np.dot(m, v[1])
        return np.degrees(np.arctan2(y, x))

    def run(self):
        """
        Run check

        Parameters:

        Class variables:
        """
        self.readFingerprint()
        self.readXYZ()
        self.readFrcmod()
        self.removeUseless()
        self.checkBonds()
        self.checkAngles()
        self.writeFrcmod()