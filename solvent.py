import shutil


class Solvent:
    def __init__(self):
        """
        Class to modify the tleap file and manage the different solvents

        Parameters:
            :param str solvent: solvent to be used for MD simulation, check solvent list for available options

        Class variables:
        """
        self.solventFilesPath = ''
        self.solvent = None
        self.solventList = ['ACN','DCM', 'MTL', 'ETL','TOL','CL3', 'THF']
        self.tleap = []
        self.status = 1

    def checkSolvent(self, solvent):
        """
        Checks whether the solvent desired is available

        Parameters:
            :param str solvent: solvent to be used for MD simulation, check solvent list for available options

        Class variables:
        """
        if solvent in self.solventList:
            self.solvent = solvent
            return 1
        else:
            print('Selected solvent is not supported\n')
            return 0

    def readLeapFile(self, path):
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

    def writeLeapFile(self, path):
        """
        read tleap file

        Parameters:
            :param str path: full path to the tleap file

        Class variables:
        """
        flag = 0
        try:
            f = open(path,'w')
            for line in self.tleap:
                if 'loadmol2' in line and flag == 0:
                    f.write('{} = loadmol2 {}.mol2'.format(self.solvent,self.solvent))
                    flag = 1
                if 'loadmol2' in line and flag == 1:
                    f.write('loadamberparams {}.frcmod'.format(self.solvent))
                    flag = 2
                if 'solvatebox' in line and flag == 2:
                    f.write('solvatebox mol {} 20.0 iso'.format(self.solvent))
                    flag = 3
                    continue
                f.write(line)
            f.close()
            return 1
        except:
            print('File write error')
            return 0

    def copyfiles(self, path):#
        """
        copy frcmod and mol2 files for solvent

        Parameters:
            :param str path: full path to the directory where files should be copied to

        Class variables:
        """
        try:
            dest = shutil.copyfile(self.solventFilesPath + r'/{}.mol2.'.format(self.solvent), path)
            dest = shutil.copyfile(self.solventFilesPath + r'/{}.frcmod'.format(self.solvent), path)
            self.status = 1
        except:
            print('File copy error')
            self.status = 0

    def applySolvent(self, solvent, leapin, leapout, path):
        """
        apply changes to solvent in tleapfile

        Parameters:
            :param str solvent: solvent name
            :param str leapin: path to tleap file to be used as input
            :param str leapout: path to tleap file to be written
            :param str path: full path to the modelling folder containing the tleap file

        Class variables:
        """
        self.status = self.checkSolvent(solvent)
        if self.status == 0:
            return
        self.status = self.readLeapFile(leapin)
        if self.status == 0:
            return
        self.status = self.writeLeapFile(leapout)
        if self.status == 0:
            return
        self.status = self.copyfiles(path)
        if self.status == 0:
            print('Tleap file modified!\n')
            return 0
        else:
            print('Error: Tleap file not modified!\n')
            return 1