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
        """

        self.path = path
        self.intype = intype
        self.outtype = outtype
        self.engines = ['amber', 'gromacs']
        os.chdir(self.path)

    def convert(self):
        amber = pmd.load_file('LIG_solv.prmtop', 'LIG_solv.inpcrd')
        #
        amber.save('LIG_solv.top')
        amber.save('LIG_solv.gro')