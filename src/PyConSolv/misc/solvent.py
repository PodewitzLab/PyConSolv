from .tleapAdderInterface import TleapAdder


class Solvent(TleapAdder):
    def __init__(self):
        """
        Class to modify the tleap file and manage the different solvents

        Parameters:

        """
        self.solventDict = {'Water': 'WAT',
                            'Acetonitrile': 'ACN',
                            'Acetone': 'ACT',
                            'Benzene': 'BNZ',
                            'Cyclohexane': 'CHX',
                            'Chloroform': 'CL3',
                            'CCl4': 'CL4',
                            'CH2Cl2': 'DCM',
                            'DMF': 'DMF',
                            'DMSO': 'DMS',
                            'Ethanol': 'ETL',
                            'Hexane': 'HEX',
                            'Methanol': 'MTL',
                            'Ammonia': 'NH3',
                            'Octanol': 'OCT',
                            'THF': 'THF',
                            'Toluene': 'TOL',
                            'custom': 'SLV'}
        TleapAdder.__init__(self, '../solvents', self.solventDict, 'solvent')