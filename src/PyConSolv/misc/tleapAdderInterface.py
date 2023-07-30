import os
import shutil


class TleapAdder:
    def __init__(self, folder, itemDict, stype):
        """
        Class to modify the tleap file and manage the different solvents

        Parameters:

        Class variables:
            - self.solventFilesPath - location of solvent parameter files
            - self.item - Tleap item to be added to be used
            - self.itemDict - dictionary of supported items
            - self.tleap - list to contain tleap file info
            - self.status - status of program (0 for error, 1 for success)
        """
        self.referenceFilesPath = os.path.split(__file__)[0] + '/{}/'.format(folder)
        self.item = None
        self.itemDict = itemDict
        self.tleap = []
        self.status = 1
        self.stype = stype

    def checkItem(self, item: str) -> int:
        """
        Checks whether the solvent desired is available

        Parameters:
            :param str solvent: solvent to be used for MD simulation, check solvent list for available options

        Class variables:
        """
        if item in self.itemDict:
            self.item = self.itemDict[item]
            return 1
        else:
            print('Selected {} is not supported\n'.format(self.stype))
            return 0

    def readLeapFile(self, path: str) -> int:
        """
        read tleap file

        Parameters:
            :param str path: full path to the tleap file

        Class variables:
        """
        try:
            f = open(path, 'r')
            for line in f:
                self.tleap.append(line)
            f.close()
            return 1
        except:
            print('File read error')
            return 0

    def writeLeapFile(self, path: str, solutename = 'LIG') -> int:
        """
        read tleap file

        Parameters:
            :param str path: full path to the tleap file

        Class variables:
        """
        if self.stype == 'solvent':
            if self.item == self.itemDict['Water']:
                return 1
        flag = 0
        try:
            print('writing into ' + path)
            f = open(path,'w')
            for line in self.tleap:
                if 'loadmol2' in line and flag == 0:
                    f.write('{} = loadmol2 {}.mol2\n'.format(self.item,self.item))
                    flag = 1
                if 'loadmol2' in line and flag == 1:
                    f.write('loadamberparams {}.frcmod\n'.format(self.item))
                    flag = 2
                if self.stype == 'solvent':
                    if 'solvatebox' in line and flag == 2:
                        f.write('solvatebox {} {} 20.0 iso\n'.format(solutename, self.item))
                        flag = 3
                        continue
                elif self.stype == 'counterion':
                    if 'addions' in line:
                        f.write('addions {} {} 1\n'.format(solutename, self.item))
                        flag = 3
                f.write(line)
            f.close()
            return 1
        except:
            print('File write error')
            return 0

    def copyfiles(self, path: str):
        """
        copy frcmod and mol2 files for solvent

        Parameters:
            :param str path: full path to the directory where files should be copied to

        Class variables:
        """
        try:
            dest = shutil.copyfile(self.referenceFilesPath + r'/{}.mol2'.format(self.item),
                                   path + r'/{}.mol2'.format(self.item))
            dest = shutil.copyfile(self.referenceFilesPath + r'/{}.frcmod'.format(self.item),
                                   path + r'/{}.frcmod'.format(self.item))
            self.status = 1
        except:
            print('File copy error')
            self.status = 0
    def changeTleapBox(self, boxsizeadd, leapin):
        """
        Change the size of a tleap box

        Parameters:
            :param float boxsize:
            :param str leapin: path to tleap file to be used as input

        Class variables:
        """
        self.readLeapFile(leapin)
        f = open(leapin, 'w')
        for line in self.tleap:
            if 'solvatebox' in line:
                f.write(' '.join(line.split()[:-1]) + ' {:.1f}\n'.format(boxsizeadd + float(line.split()[-1])))
            else:
                f.write(line)
        f.close()

    def applyItem(self, item: str, leapin: str, leapout: str, path: str, solutename = 'LIG') -> int:
        """
        apply changes to solvent in tleapfile

        Parameters:
            :param str item: Tleap item name
            :param str leapin: path to tleap file to be used as input
            :param str leapout: path to tleap file to be written
            :param str path: full path to the modelling folder containing the tleap file

        Class variables:
        """
        self.status = self.checkItem(item)
        if self.status == 0:
            return 0
        self.status = self.readLeapFile(leapin)
        if self.status == 0:
            return 0
        self.status = self.writeLeapFile(leapout, solutename)
        if self.status == 0:
            return 0
        if item != 'custom':
            self.status = self.copyfiles(path)
        if self.status == 0:
            print('Error: Tleap file not modified!\n')
            return 0
        else:
            print('Tleap file modified!\n')
            return 1