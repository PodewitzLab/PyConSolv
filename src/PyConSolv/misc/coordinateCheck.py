import numpy as np


class XYZMapper:

    def __init__(self, filepath):
        '''
        This class will compare 2 XYZ files and return a dictionary which will aid in mapping 2 different structure orders

        :param filepath str - path to the reference XYZ file:

        class variables:

        - self.xyz = reference xyz coordinates, saved as a [N,3] numpy array
        - self.map = dictionary for mapping atomd ids
        '''
        self.xyz = self.readXYZ(filepath)
        self.map  = {}
        for i in range(self.xyz.shape[0]):
            self.map[i+1] = 9999

    def mapXYZ(self, filepath):
        '''

        :param filepath: path to the XYZ file to compare to the reference
        :return:
        '''
        newXYZ = self.readXYZ(filepath)
        for i in range(newXYZ.shape[0]):
            diff = 9999
            match = 9999
            for j in range(self.xyz.shape[0]):
                if abs(np.linalg.norm(self.xyz[i] - newXYZ[j])) <= diff:
                    diff = abs(np.linalg.norm(self.xyz[i] - newXYZ[j]))
                    match = j
            self.map[i+1] = match+1



    def readXYZ(self, file) -> np.array:
        '''

        :param file str: path to the XYZ file to read
        :return: [N,3] numpy array with coordinates
        '''
        return np.loadtxt(file, skiprows=2, dtype=object)[:,1:].astype(float)

    def mapAtom(self, id: int) -> int:
        '''

        :param id int: id of the atom in the target XYZ file
        :return:
        '''
        for key, value in self.map.items():
            if id == value:
                return (key)

    def mapReference(self, id: int) -> int:
        '''

        :param id int: id of the atom in the reference
        :return:
        '''
        return self.map[id]

