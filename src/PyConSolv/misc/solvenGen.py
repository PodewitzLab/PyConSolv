import os
import shutil
import subprocess

import numpy as np

from ..utils.charge import ChargeChanger
from ..utils.colorgen import Color
from ..interfaces.amber import amberInterface
from ..interfaces.calculate import Calculation
from ..misc.inputparser import XYZ
from ..interfaces.multiWFN import MultiWfnInterface


class solventParametrizer:
    """
    Class used to parametrize arbitrary solvents using ORCA/Amber/Multiwfn and RESP charges
    """
    def __init__(self, structurePath : str):
        """
        Class used to parametrize arbitrary solvents using ORCA/Amber/Multiwfn and RESP charges

        Parameters:
            :param string structurePath: path to XYZ file containing the solvent geometry

        Class variables:
            - self.charges = contains the RESP charges generated by Multiwfn
        """
        self.db_file = os.path.split(__file__)[0] + '/../db/atom-radius.txt'
        self.db_metal_file = os.path.split(__file__)[0] + '/../db/metal-radius.txt'

        self.inputpath = '/'.join(structurePath.split('/')[:-1]) + '/solv_param/'
        self.input = self.inputpath + structurePath.split('/')[-1]

        try:
            os.mkdir(self.inputpath)
        except:
            print('The folder contains a previously performed solvent parametrization...\n')
        shutil.copyfile(structurePath, self.input)
        self.outpath = self.inputpath
        self.solventPDB = None
        self.solventmol2 = None
        self.orca_inp = '''{}
{}

%PAL NPROCS {} END

%CPCM       EPSILON      {}
            REFRAC       {}
END

%scf
maxiter 350
end

* xyzfile 0 1 input.xyz
'''

    def parseXYZ(self):
        """
        Prepare the input XYZ file for paremetrization

        Parameters:

        Class variables:
        """

        self.xyz = XYZ(self.db_file, self.db_metal_file)
        self.xyz.readXYZ(self.input)
        self.xyz.calculateDistanceMatrix()
        self.xyz.generateAdjacencyMatrix()
        self.xyz.generateLinkList()
        self.xyz.connectedCompponents()
        self.xyz.assignChain()
        self.xyz.createPDB()
        self.xyz.writePDBFiles(self.outpath)
        self.solventPDB = self.xyz.filenames[0]
        self.solventmol2 = np.array(self.xyz.ligands)[0]


    def runAntechamber(self):
        """
        Read RESP charges generated by Multiwfn and saved in orca_opt.molden.chg

        Parameters:

        Class variables:
            - self.charges = contains the RESP charges generated by Multiwfn
        """

        self.amber = amberInterface(self.outpath)
        self.amber.antechamber(self.solventPDB.split('.')[0], charge = 0) # solvent should be neutral
        self.amber.runParmchk2(self.solventmol2)


    def getMultiwfnCharges(self):
        """
        Read RESP charges generated by Multiwfn and saved in orca_opt.molden.chg

        Parameters:

        Class variables:
            - self.charges = contains the RESP charges generated by Multiwfn
        """
        f = open(self.outpath + 'solvopt/orca_opt.molden.chg', 'r')
        self.charges = []
        iterator = 0
        for line in f:
            iterator += 1
            self.charges.append([line.split()[-1]]) # this maintains compatibility with the multiple fragment approach
        f.close()


    def changeCharges(self):
        fin = self.inputpath + '/A.mol2'  # read mol2 files corresponding to the solvent pdb
        fout = self.outpath + '/SLV.mol2'
        self.chargeChanger = ChargeChanger()
        self.chargeChanger.change(fin, fout, 'SLV', self.charges)
    def changeCharges_old(self):
        """
                Create final mol2 files for tleap, using the charges from the RESP calculation

                Parameters:
                    :param string path: folder where files are located/created
                Class variables:

                """
        iterator = 0  # keeps track of current atom id
        switch = 0

        fin = open(self.inputpath + '/A.mol2', 'r')  # read mol2 files corresponding to the solvent pdb
        fout = open(self.outpath + '/SLV.mol2', 'w')  # create new mol2 file
        for line in fin:
            if 'USER_CHARGES' in line:
                fout.write('RESP Charge\n')
            elif '@<TRIPOS>ATOM' in line:
                switch = 1
                fout.write(line)
            elif '@<TRIPOS>BOND' in line:
                fout.write(line)
                switch = 2
            elif '@<TRIPOS>SUBSTRUCTURE' in line:
                fout.write(line)
                switch = 3
            elif '@<TRIPOS>MOLECULE' in line:
                fout.write(line)
                switch = 4
            else:
                if switch == 1:

                    tmp_line = line.split()[:-1]  # remove old charge

                    tmp_line.append(float(self.charges[iterator]))  # add resp charge
                    tmp_line[-2] = 'SLV'
                    fout.write(
                        '{:>7} {:<7}    {:>7}    {:>7}    {:>7}   {:<2}    {}  {:>2} {:>9.6f}\n'.format(*tmp_line))
                    iterator += 1
                elif switch == 3:
                    fout.write('     1 {}         1 TEMP              0 ****  ****    0 ROOT\n'.format('SLV'))
                    switch = 0
                elif switch == 4:
                    fout.write('SLV\n')
                    switch = 0
                else:
                    if 'bcc' in line:
                        fout.write(line.replace('bcc', 'RESP Charge'))
                    else:
                        fout.write(line.replace('ar', '1'))
        fout.close()
        fin.close()

    def runORCA(self, epsilon:str, refrac:str, method: str = 'PBE0', basis: str = 'def2-SVP', DSP: str = 'D4', CPU: int = 12):
        try:
            os.mkdir(self.outpath + '/solvopt/')
        except:
            print('Previous calculations found, checking...\n')
        os.chdir(self.outpath + '/solvopt/')
        shutil.copyfile(self.input, self.outpath + '/solvopt/input.xyz')
        self.calc = Calculation(self.outpath)
        self.calc.checkpath()

        self.method = '! {} {} {}'.format(method,basis,DSP)
        f = open(self.outpath + '/solvopt/orca_opt.inp', 'w')
        f.write(self.orca_inp.format(self.method, '! OPT', CPU, epsilon, refrac))
        f.close()

        self.calc.calculate(calctype='opt', calcpath = '/solvopt')
        command = 'orca_2mkl orca_opt -molden'
        calc = subprocess.run([command], shell=True)

        os.chdir(self.outpath)

    def runMultiwfn(self):
        os.chdir(self.outpath + '/solvopt/')
        self.multiwfn = MultiWfnInterface(self.outpath+'/solvopt/', orcaname = 'orca_opt')
        self.multiwfn.run()

        os.chdir(self.outpath)
    def run(self, epsilon:str, refrac:str, method: str = 'PBE0', basis: str = 'def2-SVP', DSP: str = 'D4', CPU: int = 12):
        if os.path.exists(self.inputpath + '../done'):
            print('Previously performed parametrization is complete, using solvent\n')
            return
        self.parseXYZ()
        self.runORCA(epsilon, refrac, method, basis, DSP, CPU)
        self.runMultiwfn()
        self.runAntechamber()
        self.getMultiwfnCharges()
        self.changeCharges()
        dest = shutil.copyfile(self.inputpath + '/A.frcmod', self.inputpath + '/SLV.frcmod')
        print(Color.GREEN + 'Solvent parametrization complete!\n' + Color.END)
        f = open(self.inputpath + '../done', 'w')
        f.close()