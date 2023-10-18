from .ions import ionlib
from .tleapAdderInterface import TleapAdder


class Counterion(TleapAdder):
    def __init__(self):
        """
        Class to modify the tleap file and manage the different solvents

        Parameters:

        """
        self.counterionDict = {'OTf-': 'OTF',
                               'BF4-': 'BF4',
                               'BARF-': 'BAR',
                               'PFC-': 'PFC',
                               'ScF6-': 'SCF',
                               'ClO4-': 'CLO',
                               'BPh4-': 'BPH',
                               'custom': 'CTI'}
        ions = ionlib()
        self.counterionDict.update(ions.ionsinAmber)
        TleapAdder.__init__(self, '../counterions', self.counterionDict, 'counterion')