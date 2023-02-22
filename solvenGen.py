import os
import shutil

import numpy as np
import subprocess

from .calculate import Calculation
from .multiWFN import MultiWfnInterface
from .inputparser import XYZ
from .amber import amberInterface


class solventParametrizer:
    def __init__(self, structurePath : str, outpath : str):
        self.db_file = os.path.split(__file__)[0] + '/db/atom-radius.txt'
        self.db_metal_file = os.path.split(__file__)[0] + '/db/metal-radius.txt'

        self.input = structurePath
        self.outpath = outpath
        self.solventPDB = None
        self.solventmol2 = None
        self.orca_inp = '''{}
        {}
        {}

        %PAL NPROCS {} END

        %scf
        maxiter 350
        end

        * xyzfile 0 1 input.xyz
        '''

    def parseXYZ(self):

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

        self.amber = amberInterface(self.outpath)
        self.amber.antechamber(self.solventPDB)
        self.amber.runParmchk2(self.solventmol2)


    def getMultiwfnCharges(self):
        pass

    def changeCharges(self):
        pass

    def runORCA(self, method: str = 'PBE0', basis: str = 'def2-SVP', DSP: str = 'D4', CPCM: str = 'Water', CPU: int = 12):
        os.mkdir(self.outpath + '/solvopt/')
        os.chdir(self.outpath + '/solvopt/')
        shutil.copyfile(self.input, self.outpath + '/solvopt/input.xyz')
        self.calc = Calculation(self.outpath)
        self.calc.checkpath()

        self.method = '! {} {} {}'.format(method,basis,DSP)
        f = open(self.outpath + '/solvopt/orca_opt.inp', 'w')
        f.write(self.orca_inp.format(self.method, '! OPT', '! CPCM({})'.format(CPCM), CPU))
        f.close()

        self.calc.calculate(calctype='opt', calcpath = '/solvopt')
        command = 'orca_2mkl orca_freq -molden'
        calc = subprocess.run([command], shell=True)

        os.chdir(self.outpath)

    def runMultiwfn(self):
        os.chdir(self.outpath + '/solvopt/')
        self.multiwfn = MultiWfnInterface(self.outpath+'/solvopt/', orcaname = 'orca_opt')
        self.multiwfn.run()

        os.chdir(self.outpath)
    def run(self):
        pass