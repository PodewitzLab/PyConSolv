import os

import parmed as pmd
class Converter:
    def __init__(self, path, intype = 'amber', outtype = 'gromacs'):
        """
        Class to facilitate conversion between different simulation engines

        Parameters:

        Class variables:
            - self.path - location of input files
            - self.intype - MD engine input type
            - self.outtype - MD engine output type
            - self.engines - dict containing the supported engine types
            - self.status - status of program (0 for error, 1 for success)
            - self.atoms - number of atoms in solute
            - self.system_size - number of atoms in system
            - self.original_wd
        """

        self.path = path
        self.intype = intype
        self.outtype = outtype
        self.engines = ['amber', 'gromacs']
        self.original_wd = os.getcwd()

        self.atoms = None
        self.system_size = None

    def getAtoms(self):
        os.chdir(self.path)
        f = open('input.xyz', 'r')
        for line in f:
            self.atoms = int(line.split()[0])
            break
        f.close()
        os.chdir(self.original_wd)

    def getSystemSize(self):
        os.chdir(self.path)
        f = open('LIG_solv.prmtop', 'r')
        start = 0
        for line in f:
            if '%FLAG POINTERS' in line:
                start = 1
                continue
            if start == 2:
                self.system_size = int(line.split()[0])
                break
            elif start == 1:
                start = 2
                continue
        f.close()
        os.chdir(self.original_wd)
    def convert(self):
        os.chdir(self.path)
        parm = pmd.load_file('LIG_solv.prmtop', 'LIG_solv.inpcrd')

        if self.outtype == 'gromacs':
            if os.path.exists('LIG_solv.top'):
                os.remove('LIG_solv.top')
            if os.path.exists('LIG_solv.gro'):
                os.remove('LIG_solv.gro')
            if os.path.exists('LIG_dry.gro'):
                os.remove('LIG_dry.gro')
            parm.save('LIG_solv.top')
            parm.save('LIG_solv.gro')

            parm = pmd.load_file('LIG_dry.prmtop', 'LIG_dry.inpcrd')
            parm.save('LIG_dry.top')
        os.chdir(self.original_wd)
