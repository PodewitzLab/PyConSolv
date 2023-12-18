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
            - self.stype - type of item to be added e.g. solvent, counterion
            - self.amberNative - flag to determine if a mol2 and frcmod should be loaded for this item
        """
        self.referenceFilesPath = os.path.split(__file__)[0] + '/{}/'.format(folder)
        self.item = None
        self.itemDict = itemDict
        self.tleap = []
        self.status = 1
        self.stype = stype
        self.amberNative = False
        self.watermodel = None

    def checkItem(self, item: str) -> int:
        """
        Checks whether the solvent desired is available

        Parameters:
            :param str item: item to be used for MD simulation, check item list(solvent/counterions) for available options

        Class variables:
        """
        if 'Water' in item:
            self.watermodel = item
            item = 'Water'
        if item in self.itemDict:
            self.item = self.itemDict[item]
            if type(self.item) is list:
                self.item = self.itemDict[item][0]
                self.amberNative = True
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

    def writeLeapFile(self, path: str, solutename = 'LIG', amount: int = 1) -> int:
        """
        read tleap file

        Parameters:
            :param str path: full path to the tleap file
            :param int amount: number of counterions to add

        Class variables:
        """
        if self.stype == 'solvent':
            if self.item == self.itemDict['Water']:
                self.tleap = [line.replace('water.opc', 'water.tip3p') for line in self.tleap]
                self.tleap = [line.replace('OPCBOX', 'TIP3PBOX') for line in self.tleap]
                self.tleap = [line.replace('frcmod.ionslm_126_opc', 'frcmod.ions1lm_126_tip3p') for line in self.tleap]
                if self.watermodel == 'WaterOPC':
                    self.tleap = [line.replace('.water.tip3p', '.water.opc') for line in self.tleap]
                    self.tleap = [line.replace('TIP3PBOX', 'OPCBOX') for line in self.tleap]
                    self.tleap = [line.replace('frcmod.ions1lm_126_tip3p', 'frcmod.ionslm_126_opc') for line in
                                  self.tleap]
                elif self.watermodel == 'WaterTIP4PEW':
                    self.tleap = [line.replace('.water.tip3p', '.water.tip4pew') for line in self.tleap]
                    self.tleap = [line.replace('TIP3PBOX', 'TIP4PEWBOX') for line in self.tleap]
                    self.tleap = [line.replace('frcmod.ions1lm_126_tip3p', 'frcmod.ions1lm_126_tip4pew') for line in
                                  self.tleap]
                try:
                    f = open(path,'w')
                    for line in self.tleap:
                        f.write(line)
                    f.close()
                    return 1
                except:
                    print('File write error')
                    return 0
        flag = 0
        try:
            print('writing into ' + path)
            f = open(path,'w')
            for line in self.tleap:
                if 'loadmol2' in line and flag == 0:
                    if not self.amberNative:
                        f.write('{} = loadmol2 {}.mol2\n'.format(self.item,self.item))
                    flag = 1
                if 'loadmol2' in line and flag == 1:
                    if not self.amberNative:
                        f.write('loadamberparams {}.frcmod\n'.format(self.item))
                    flag = 2
                if self.stype == 'solvent':
                    if 'solvatebox' in line and flag == 2:
                        f.write('solvatebox {} {} 10.0 iso\n'.format(solutename, self.item))
                        flag = 3
                        continue
                elif self.stype == 'counterion':
                    if 'saveamberparm LIG LIG_dry.prmtop LIG_dry.inpcrd' in line:
                        f.write('addions {} {} {}\n'.format(solutename, self.item, amount))
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
                tmp = line.split()
                tmp[3] = str(float(tmp[3]) + boxsizeadd)
                if 'iso' not in line:
                    tmp.append(' iso') #add iso parameter everywhere for the box
                f.write(' '.join(tmp) + '\n' )
            else:
                f.write(line)
        f.close()

    def applyItem(self, item: str, leapin: str, leapout: str, path: str, solutename = 'LIG', amount: int = 1) -> int:
        """
        apply changes to solvent in tleapfile

        Parameters:
            :param str item: Tleap item name
            :param str leapin: path to tleap file to be used as input
            :param str leapout: path to tleap file to be written
            :param str path: full path to the modelling folder containing the tleap file
            :param int amount: amount of counterions to be added

        Class variables:
        """
        self.status = self.checkItem(item)
        if self.status == 0:
            return 0
        self.status = self.readLeapFile(leapin)
        if self.status == 0:
            return 0
        self.status = self.writeLeapFile(leapout, solutename, amount)
        if self.status == 0:
            return 0
        if item != 'custom':
            if 'Water' not in item:
                self.status = self.copyfiles(path)
        if self.status == 0:
            print('Error: Tleap file not modified!\n')
            return 0
        else:
            print('Tleap file modified!\n')
            return 1