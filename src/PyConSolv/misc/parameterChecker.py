import numpy as np


class ParameterChecker:
    def __init__(self, folder):
        self.frcmodfile = folder + '/LIG_mcpbpy.frcmod'
        self.xyzfile = folder + '/input.xyz'
        self.fingerprintfile = folder + '/LIG_standard.fingerprint'
        self.xyz = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.atoms = {}

    def checkBonds(self):
        for el in self.bonds:
            if 'M' in el[0] and 'Y' in el[0]:
                atoms = el[0].split('-')
                p = np.zeros((2, 3))
                p[0] = np.array(self.xyz[self.atoms.get(atoms[0])][1:], dtype=float)
                p[1] = np.array(self.xyz[self.atoms.get(atoms[1])][1:], dtype=float)
                print(float(el[2]), self.bondCalculator(p))

    def checkAngles(self):
        for el in self.angles:
            atoms = el[0].split('-')
            if 'M' in atoms[0] or 'Y' in atoms[0]:
                if 'M' in atoms[1] or 'Y' in atoms[1]:
                    if 'M' in atoms[2] or 'Y' in atoms[2]:
                        print(atoms)
                        p = np.zeros((3, 3))
                        p[0] = np.array(self.xyz[self.atoms.get(atoms[0])][1:], dtype=float)
                        p[1] = np.array(self.xyz[self.atoms.get(atoms[1])][1:], dtype=float)
                        p[2] = np.array(self.xyz[self.atoms.get(atoms[2])][1:], dtype=float)
                        print(float(el[2]), self.angleCalculator(p))

    def checkDihedrals(self):
        for el in self.dihedrals:
            atoms = el[0].split('-')
            p = np.zeros((3, 3))
            p[0] = np.array(self.xyz[self.atoms.get(atoms[0])][1:], dtype=float)
            p[1] = np.array(self.xyz[self.atoms.get(atoms[1])][1:], dtype=float)
            p[2] = np.array(self.xyz[self.atoms.get(atoms[2])][1:], dtype=float)
            p[3] = np.array(self.xyz[self.atoms.get(atoms[3])][1:], dtype=float)
            print(float(el[2]) - self.dihedralCalculator(p))

    def readFrcmod(self):
        f = open(self.frcmodfile, 'r')
        switch = 0
        for line in f:
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
        f = open(self.xyzfile, 'r')
        next(f)
        next(f)
        for line in f:
            if len(line.split()) == 4:
                self.xyz.append(line.split())
        f.close()

    def readFingerprint(self):
        f = open(self.fingerprintfile, 'r')
        for line in f:
            if '-> M' in line:
                self.atoms[line.split()[-1]] = int(line.split()[1]) - 1
            if '-> Y' in line:
                self.atoms[line.split()[-1]] = int(line.split()[1]) - 1
        f.close()

    def writeFrcmod(self):
        pass

    def bondCalculator(self, p):
        return np.round(abs(np.linalg.norm(p[0] - p[1])), 4)

    def angleCalculator(self, p):
        ba = p[0] - p[1]
        bc = p[2] - p[1]
        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        return np.degrees(angle)

    def dihedralCalculator(self, p):
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