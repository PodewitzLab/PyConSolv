import shutil
import os
from tkinter import *
import numpy as np

from .interfaces.cpptraj import CPPtraj
from .interfaces.mdengines import MDEngine
from .interfaces.parmed import Parmed
from .misc.coordinateCheck import XYZMapper
from .misc.counterion import Counterion
from .misc.counterionGen import counterionParametrizer
from .misc.frcmod import frcmodParser
from .misc.ions import ionlib
from .misc.mol2 import mol2Parser
from .misc.parameterChecker import ParameterChecker
from .utils.charge import ChargeChanger
from .misc.solvenGen import solventParametrizer
from .misc.solvent import Solvent
from .utils.colorgen import Color
from .interfaces.amber import amberInterface
from .interfaces.calculate import Calculation
from .misc.filestructure import Setup
from .misc.inputparser import XYZ
from .interfaces.multiWFN import MultiWfnInterface
from .misc.outgen import Faker
from .misc.ui import GUI
from .misc.restart import RestartFile
from .utils.copier import Copier


def printlog(message, file):
    print(message)
    f = open(file)
    f.write(message)
    f.close()


def error(step):
    """
    Print out an error message

    Parameters:
        - step: step at which the error happened

    Class variables:
    """

    print(Color.RED + 'Something went wrong, please check your input/output' + Color.END)
    print(Color.RED + '''
    ############################################
    ######             WARNING            ######
    ######      Calculation failed at:    ######
    ######{:^32}######
    ############################################
    '''.format(step) + Color.END)


class PyConSolv:
    """
    Run conformer generation in explicit solvent via MD simulations. This package relies on
    Ambertools(MCPB.py/antechamber/tleap), ORCA 5 and Multiwfn 3.8

    Parameters:
        - path: location of XYZ file

    Class variables:
        - self.engine - MD engine to be used for the equilibration and simulation
        - self.MDEngine - MD engine instance
        - self.parameterChecker - instance of frcmod parameter error checker
        - self.version - program version
        - self.hasMetal - True when a metal is part of the structure, False otherwise
        - self.restarter - object for reading restart files (see restart.py)
        - self.inputpath - full path to the folder where the input.xyz file is located
        - self.path - full path to the input.xyz file is located
        - self.status - variable used for error checking. 0 means error, 1 means all is well
        - self.restart - variable used to check the restart point
        - self.db_file - location of the atom-radius.txt file which contains ionic radius information
        - self.db_metal_file - location of the atom-radius.txt file which contains ionic radius information for metals
        - self.amber - Object containing an amberInterface (see amber.py)
        - self.MCPB - full path to the MCPB_setup folder
        - self.xyz - Object containing an XYZ (see inputparser.py)
        - self.solventsImplemented - list of supported solvents
        - self.refrac - refractive index of custom solvent
        - self.epsilon - permittivity of custom solvent
        - self.addSolvent - if true, the optimized and parametrized custom solvent will be stored
        - self.solventAbb - abbreviation for new solvent
        - self.solventParamPath - path to location of solvent XYZ file
        - self.solventPath - path to pre-parametrized solvents
        - self.counterIon - counterion to be used
        - self.counterIonsImplemented - list of supported counterions
        - self.bondsRT - list of bonds to be added for the restrained simulation
        - self.strengthRT - list of harmonic potential offsets
        - self.forceConstants - list of harmonic potential force constants
        - self.restraintWidth - list of restraint widths
        - self.map - map of atom id mapping to original file
    """

    def __init__(self, path):
        self.engine = None
        self.MDEngine = None
        self.parameterChecker = None
        path = os.path.abspath(path)
        self.counterIon = ''
        self.refrac = None
        self.epsilon = None
        self.solventParamPath = None
        self.version = '1.0.6.3'
        self.metals = ['LI', 'BE', 'NA', 'MG', 'AL', 'SI', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN', 'FE',
                       'CO', 'NI', 'CU', 'ZN',
                       'GA', 'GE', 'AS', 'SE', 'BR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG',
                       'CD', 'IN', 'SN', 'SB',
                       'TE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER',
                       'TM', 'YB', 'YB', 'LU', 'HF',
                       'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB', 'BI', 'PO', 'AT', 'FR', 'RA', 'AC',
                       'TH', 'PA', 'U', 'NP',
                       'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', 'MD', 'NO',
                       'LR']
        self.hasMetal = None
        self.restarter = None
        self.path = path
        self.inputFile = path.split('/')[-1]
        self.inputpath = '/'.join(path.split('/')[:-1])
        if self.inputFile != 'input.xyz':
            print('Copying {} to input.xyz\n'.format(self.inputFile))
            shutil.copyfile(self.path,self.inputpath + '/input.xyz')
            self.inputFile = 'input.xyz'
        if '.xyz' not in self.inputFile:
            error('Initialization... make sure the input is an xyz file')
        self.status = 0
        self.restart = 0
        self.db_file = os.path.split(__file__)[0] + '/db/atom-radius.txt'
        self.db_metal_file = os.path.split(__file__)[0] + '/db/metal-radius.txt'
        self.solventPath = os.path.split(__file__)[0] + '/solvents/'
        self.ionPath = os.path.split(__file__)[0] + '/counterions/'
        self.amber = None
        self.MCPB = self.inputpath + '/MCPB_setup'
        self.xyz = None
        ions = ionlib()
        amberions = list(ions.ionsinAmber.keys())
        self.counterIonsImplemented = amberions +['OTf-','BF4-', 'BARF-', 'PFC-', 'ScF6-', 'ClO4-', 'BPh4-', 'custom']
        self.solventsImplemented = ['Water', 'WaterOPC', 'WaterTIP4PEW','Acetonitrile', 'Acetone', 'Benzene', 'Cyclohexane', 'Chloroform', 'CCl4',
                                    'CH2Cl2', 'DMF', 'DMSO', 'Ethanol', 'Hexane', 'Methanol', 'Ammonia', 'Octanol',
                                    'THF', 'Toluene', 'custom']
        self.watermodels = ['Water', 'WaterOPC', 'WaterTIP4PEW']
        self.metalCheck()
        self.addSolvent = False
        self.solventAbb = ''
        self.boxsize = 20
        self.bondsRT = []
        self.strenghtRT = []
        self.forceConstants = []
        self.restraintWidth = []
        self.map = None


        self.startInfo()

    def startInfo(self):
        print(Color.BLUE + r'''

                  _____        _____             _____       _       
                 |  __ \      / ____|           / ____|     | |      
                 | |__) |   _| |     ___  _ __ | (___   ___ | |_   __
                 |  ___/ | | | |    / _ \| '_ \ \___ \ / _ \| \ \ / /
                 | |   | |_| | |___| (_) | | | |____) | (_) | |\ V / 
                 |_|    \__, |\_____\___/|_| |_|_____/ \___/|_| \_/  
                         __/ |                                       
                        |___/                                        
                                        Ver {}

        Welcome to PyConSolv, your friendly neighbourhood conformer generator

        Calculations will be set up in:

        {}

                '''.format(self.version, self.inputpath) + Color.END)
    def checkRestart(self):
        """
        Check for the existence of a restart file (pyconsolv.restart) and get last completed step

        Parameters:

        Class variables:
            - self.restarter - RestartFile object used for checking and writing restart file (see restart.py)
        """
        if os.path.exists(self.inputpath + '/pyconsolv.restart'):
            self.restarter = RestartFile(self.inputpath)
            self.restart = self.restarter.getstate()
            print(Color.PURPLE + 'Restart file found!\n' + Color.END)
            os.remove(self.inputpath + '/pyconsolv.restart')

    def prepareSolvent(self, method, basis, dsp, cpu):
        print('You have selected a custom solvent, please input the path to the solvent xyz file:\n')
        self.solventParamPath = os.path.abspath(input())

        print('Please enter the epsilon value for your custom solvent:\n')
        self.epsilon = input()

        print('Please enter the refractive index value for your custom solvent:\n')
        self.refrac = input()

        # print('Do you want to add this custom solvent to the list of available solvents?[y/n]\n')
        # add = str(input())
        # while add.lower() not in ['y','n']:
        #     print('Type "y" if you want to add this custom solvent to the list of available solvents or "n" if not:')
        #     add = str(input())
        # if add.lower() == 'y':
        #     self.addSolvent = True
        #     print('Give this solvent a unique abbreviations consisting of exactly 3 characters:\n')
        #     abb = str(input())
        #     while len(abb) != 3:
        #         print('It is mandatory to give exactly 3 characters! If you changed your mind please enter "quit"')
        #         abb=str(input())
        #         if abb.lower() == "quit":
        #             self.addSolvent = False
        #             print('Solvent files will not be added to pyConSolv!')
        #             break
        #     self.solventAbb = abb
        #     print('You chose {} as abbreviation!'.format(self.solventAbb))

        print(Color.GREEN + '''

        ###################
        Starting solvent parametrization in:
        {}
        ####################

        '''.format('/'.join(self.solventParamPath.split('/')[:-1]) + '/solv_param') + Color.END)
        self.solventParam = solventParametrizer(self.solventParamPath)
        self.solventCPCM = '''

        %CPCM       EPSILON      {}
                    REFRAC       {}
        END
        '''.format(self.epsilon, self.refrac)
        self.solventParam.run(method, basis, dsp, cpu, self.solventCPCM, 0)

    def prepareCounterion(self, method, basis, dsp, cpu, charge, cpcm):

        answered = False
        while not answered:
            answer = input('You system is charged, would you like to select an appropriate counterion? \n y/n: ')
            if answer == 'y':
                print('Please choose one of the following counterions:\n')
                print(', '.join(self.counterIonsImplemented) + '\n')
                while not answered:
                    ion = input()
                    if ion not in self.counterIonsImplemented:
                        print('Selected ion not supported, try again...\n')
                    else:
                        print('Please enter the number of counterions to be used:\n')
                        amount = int(input())
                        answered = True
                        f = open(self.inputpath + '/counterion', 'w')
                        f.write('{} {}'.format(ion, amount))
                        f.close()
            elif answer == 'n':
                ion = ''
                answered = True
            else:
                print('Wrong input, try again\n')


        if ion == 'custom':
            print('You have chosen a custom counterion for your system, which needs to be parametrized\n' +
                  'Please provide a path to the location of an XYZ file containing your ion\n')
            self.ionParamPath = input()
            chargeIon = input('Please enter the charge for your ion: \n')
            cIon = counterionParametrizer(self.ionParamPath)
            if cpcm == 'custom':
                cIon.run(method, basis, dsp, cpu, self.solventCPCM, int(chargeIon))
            else:
                cIon.run(method, basis, dsp, cpu, cpcm, int(chargeIon))
        self.counterIon = ion



    def setup(self, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP', dsp: str = 'D4',
              cpcm: str = 'Water', cpu: int = 12, multiplicity:int  = 1, opt: bool = True) -> int:
        """
        Run setup for creating the appropriate folders and parse XYZ file

        Parameters:
            :param int charge: charge of the complete system
            :param string method: ORCA 5 method line
            :param string basis: Basis set for ORCA calculations
            :param string dsp: Dispersion corrections
            :param string cpcm: CPCM solvation model solvent
            :param int cpu: number of CPU cores to be used
            :param int multiplicity: multiplicity for the system

        Class variables:
        """
        if cpcm == 'custom':
            self.prepareSolvent(method, basis, dsp, cpu)
        if cpcm not in self.solventsImplemented:
            print(Color.RED + 'Selected solvent is not yet implemented\n' + Color.END)
            return 0
        if cpcm in self.watermodels:
            cpcmname = 'Water'
        else:
            cpcmname = cpcm

        if charge != 0:
            if not os.path.exists(self.inputpath+'/counterion'):
                self.prepareCounterion(method, basis, dsp, cpu, charge, cpcmname)

        if self.restart == 0:
            self.xyz = XYZ(self.db_file, self.db_metal_file)
            self.xyz.prepareInput(self.inputpath + '/input.xyz')
            self.xyz = None
            setup = Setup(self.inputpath + '/' + self.inputFile, charge=charge, multi = multiplicity, opt = opt)
            setup.Method(method, basis, dsp, cpcmname, cpu, self.epsilon, self.refrac)
            self.status = setup.run()
            if self.status == 0:
                error('Setup')
                return 0

        if cpcmname == 'custom':
            shutil.copyfile('/'.join(self.solventParamPath.split('/')[:-1]) + '/SLV_param/SLV.frcmod',
                            self.MCPB + '/SLV.frcmod')
            shutil.copyfile('/'.join(self.solventParamPath.split('/')[:-1]) + '/SLV_param/SLV.mol2',
                            self.MCPB + '/SLV.mol2')
        else:
            solvent = Solvent()
            solvname = solvent.solventDict[cpcmname]
            if cpcmname != 'Water':
                shutil.copyfile(self.solventPath + '/{}.frcmod'.format(solvname),
                                self.MCPB + '/{}.frcmod'.format(solvname))
                shutil.copyfile(self.solventPath + '{}.mol2'.format(solvname),
                                self.MCPB + '/{}.mol2'.format(solvname))

        if self.counterIon == 'custom':
            shutil.copyfile('/'.join(self.ionParamPath.split('/')[:-1]) + '/CTI_param/CTI.frcmod',
                            self.MCPB + '/CTI.frcmod')
            shutil.copyfile('/'.join(self.ionParamPath.split('/')[:-1]) + '/CTI_param/CTI.mol2',
                            self.MCPB + '/CTI.mol2')
        elif self.counterIon == '':
            pass

        else:
            cIon = Counterion()
            ions = ionlib()
            if self.counterIon in list(ions.ionsinAmber.keys()):
                ionname = ions.ionsinAmber[self.counterIon]
            else:
                ionname = cIon.counterionDict[self.counterIon]
                shutil.copyfile(self.ionPath + '/{}.frcmod'.format(ionname),
                                self.MCPB + '/{}.frcmod'.format(ionname))
                shutil.copyfile(self.ionPath + '{}.mol2'.format(ionname),
                                self.MCPB + '/{}.mol2'.format(ionname))

        f = open(self.inputpath + '/simulation/solvent', 'w')
        f.write(cpcmname)
        f.close()

        print(Color.GREEN + 'Setup is complete, moving on to ORCA calculations...\n' + Color.END)

        self.map = self.mapper()

        self.restarter = RestartFile(self.inputpath)

        return 1

    def orca(self, opt: bool = True):
        """
        Run ORCA optimization and frequency calculations

        Parameters:
            :param opt bool: if set to False, a single point calculation will be performed
        Class variables:
        """
        try:
            self.restarter.write('setup')
        except AttributeError:
            return 0

        calculation = Calculation(self.inputpath + '/orca_calculations')
        if self.hasMetal:
            self.status = calculation.run(opt = opt)
        else:
            self.status = calculation.run(freq=False, opt = opt)
            shutil.copyfile(self.inputpath + '/orca_calculations/opt/orca_opt.xyz', self.inputpath + '/MCPB_setup/input.xyz')
        if self.status == 0:
            error('ORCA Calculations')
            return 0

        print(Color.GREEN + 'ORCA Calculations complete, moving on to MCPB setup...' + Color.END)
        self.restarter.write('orca')
        return 1

    def antechamber(self):
        """
        Run antechamber for each generated fragment

        Parameters:

        Class variables:
        """
        self.restarter.write('orca')
        shutil.copyfile(self.inputpath + '/orca_calculations/opt/orca_opt.xyz', self.inputpath + '/MCPB_setup/input.xyz')
        self.xyz = XYZ(self.db_file, self.db_metal_file)
        self.xyz.prepareMCPB(self.inputpath + '/MCPB_setup/input.xyz')
        pdbs = self.xyz.filenames

        print('The following pdb files were created: \n {}\n'.format(' '.join(pdbs)))
        print('Please enter fragment charges:\n')

        # get charges from user
        window = Tk()
        window.title('Fragment charge assignment')
        start = GUI(window, self.MCPB, pdbs)
        window.mainloop()
        window.destroy()

        # create mol2 files and run antechamber
        self.xyz.writeMol2Files()
        antechamberFiles = self.xyz.molNotCreated
        ligands = np.array(self.xyz.ligands)
        self.amber = amberInterface(self.MCPB)
        for filename in antechamberFiles:
            self.status = self.amber.antechamber(*filename)
            if self.status == 0:
                error('antechamber for {}'.format(filename))
                return 0

        print('Generating frcmod files for ligands:\n')
        for filename in ligands:
            self.status = self.amber.runParmchk2(filename)
            if self.status == 0:
                error('parmchk2 for {}'.format(filename))
                return 0
        self.hasMetal = self.xyz.hasMetal

        if self.hasMetal:  # if metal is detected, proceed with MCPB.py
            metals = self.xyz.metals
            self.xyz.writeMetalConnections(self.MCPB)  # write out metal connections file
            self.xyz.writeConnections(self.MCPB)  # write out connections file
            self.amber.inputFileGenerator(metals[0][1], ligands[:, 1])
            self.restarter.write('frcmod')
            return 1

        else:  # if no metal, proceed with tleap, presumes only 1 ligand
            # os.chdir(self.inputpath + '/equilibration/')
            print(antechamberFiles)
            shutil.copyfile(self.MCPB + '/' + str(antechamberFiles[0][0]) + '.mol2',
                            self.MCPB + '/LIG.mol2')
            shutil.copyfile(self.MCPB + '/' + antechamberFiles[0][0] + '.frcmod',
                            self.MCPB + '/LIG.frcmod')
            self.amber.tleapNoMetalSolv(self.MCPB)
            self.restarter.write('frcmod')
            return 1

    def multiwfn(self, cores: int) -> int:
        """
        Run Multiwfn charge calculations

        Parameters:
            :param int cores: number of cpu cores for Multiwfn

        Class variables:
        """
        self.restarter.write('frcmod')
        print(Color.GREEN + 'Fragments have been prepared, running MultiWfn task...\n\n' + Color.END)
        if self.xyz is None:
            self.xyz = XYZ(self.db_file, self.db_metal_file)
            self.xyz.hasMetal = self.hasMetal
            self.xyz.readFilenames(self.MCPB)
        if not self.hasMetal:
            print('filenames')
            print(self.xyz.filenames)
            residues = []
            for elem in self.xyz.filenames:
                residues.append(elem.replace('.pdb','').replace('\n',''))

            multiwfn = MultiWfnInterface(self.inputpath + '/orca_calculations/opt/', orcaname='orca_opt')
            self.status = multiwfn.run(cores)
            self.xyz.hasMetal = False
            self.xyz.readRESP(self.inputpath + '/orca_calculations/')
            chargeChanger = ChargeChanger()

            if len(residues) > 1:
                for residue in residues:
                    chargeChanger.change(self.MCPB + '/{}.mol2'.format(residue),
                                         self.MCPB + '/{}x.mol2'.format(residue), residue, self.xyz.charges)
                mol2parser = mol2Parser(self.MCPB,['{}x'.format(x) for x in residues])
                mol2parser.writeCombinedMol2()
                frcmodparser = frcmodParser(self.MCPB, residues)
                frcmodparser.writeCombinedFrcmod()
            else:
                 chargeChanger.change(self.MCPB + '/A.mol2',
                                         self.MCPB + '/LIG.mol2', 'A', self.xyz.charges)
        else:
            multiwfn = MultiWfnInterface(self.inputpath + '/orca_calculations/freq/')
            self.status = multiwfn.run(cores)

        if self.status == 0:
            error('MultiWfn Calculations')
            return 0
        self.restarter.write('multiwfn')
        return 1

    def MCPB_script(self):
        """
        Run MCPB.py for the system

        Parameters:

        Class variables:
        """

        if not self.hasMetal:
            self.restarter.write('mcpb')
            return 1

        self.restarter.write('multiwfn')

        print(Color.GREEN + 'Converting ORCA output to MCPB.py compatible input...\n' + Color.END)

        faker = Faker(self.inputpath + '/orca_calculations/freq/')
        faker.fakecrds()
        # faker.fakeesp()
        faker.fakeforce()

        shutil.copyfile(self.inputpath + '/orca_calculations/freq/fakechk.fchk',
                        self.inputpath + '/MCPB_setup/LIG_small_opt.fchk')
        shutil.copyfile(self.inputpath + '/orca_calculations/freq/fakelog.log',
                        self.inputpath + '/MCPB_setup/LIG_small_fc.log')


        print(Color.GREEN + 'Proceeding with MCPB steps...\n' + Color.END)

        if self.amber is None:
            self.amber = amberInterface(self.MCPB)

        self.status = self.amber.runMCPB('1')

        if self.status == 0:
            error('MCPB step 1')
            return 0

        if not self.amber.checkMCPBBonds(self.MCPB):
            self.status = self.amber.runMCPB('1')

            if self.status == 0:
                error('MCPB step 1')
                return 0

        self.status = self.amber.runMCPB('2')

        if self.status == 0:
            error('MCPB step 2')
            return 0

        if self.xyz is None:
            self.xyz = XYZ(self.db_file, self.db_metal_file)
            self.xyz.path = self.inputpath + '/MCPB_setup/'
            self.xyz.hasMetal = self.hasMetal

        if self.xyz.path is None:
            self.xyz.path = self.inputpath + '/MCPB_setup/'
            self.xyz.hasMetal = self.hasMetal
        self.xyz.createFinalMol2(self.inputpath)

        self.status = self.amber.runMCPB('4')

        if self.status == 0:
            error('MCPB step 4')
            return 0

        self.parameterChecker = ParameterChecker(self.MCPB)
        self.parameterChecker.run()

        self.restarter.write('mcpb')
        return 1

    def tleap(self, solvent: str, boxsize: int = 10) -> int:
        """
        Run tleap

        Parameters:
            :param string solvent: solvent for your box.

        Class variables:
        """
        self.restarter.write('mcpb')
        self.boxsize = boxsize

        os.chdir(self.MCPB)
        print('Solvent of choice is: {}'.format(solvent))
        if solvent not in self.solventsImplemented:
            print(Color.RED + 'Selected solvent is not yet implemented\n' + Color.END)
            return 0

        if self.amber is None:
            self.amber = amberInterface(self.MCPB)

        if self.hasMetal:
            solutename = 'mol'
        else:
            solutename = 'LIG'
        solv = Solvent()
        solv.applyItem(solvent, self.MCPB + '/LIG_tleap.in',
                          self.MCPB + '/LIG_tleap.in', self.MCPB, solutename)
        if os.path.exists(self.inputpath+'/counterion'):
            f = open(self.inputpath+'/counterion','r')
            for line in f:
                amount = line.split()[1]
                self.counterion = line.split()[0]
            f.close()
        else:
            self.counterion = ''
            amount = 0
        if self.counterIon != '':
            ion = Counterion()
            ionsAmber = ionlib()
            if self.counterIon in list(ionsAmber.ionsinAmber.keys()):
                ion.applyItem(self.counterIon, self.MCPB + '/LIG_tleap.in',
                              self.MCPB + '/LIG_tleap.in', self.MCPB, solutename, amount)
            else:
                ion.applyItem(self.counterIon, self.MCPB + '/LIG_tleap.in',
                                  self.MCPB + '/LIG_tleap.in', self.MCPB, solutename, amount)
        self.amber.tleapChecker(self.MCPB)
        self.amber.changeBoxSize(self.MCPB + '/LIG_tleap.in',self.boxsize- self.amber.defaultbox)
        self.status = self.amber.runTleap()

        if self.status == 0:
            error('Tleap')
            return 0

        self.restarter.write('tleap')

        print('Checking topology files for inconsistencies...')
        self.checkTop()
        print('Done!\n')

        print('Parametrization complete!\n')
        print('Solvent of choice is: {}\n'.format(solvent))
        return 1

    def checkTop(self):
        parmed = Parmed(self.MCPB)
        parmed.checkTop('LIG_solv')
        parmed.checkTop('LIG_dry')

    def equilibration(self, cpu: int = 12, engine = 'amber', cart: str = None, cartstr: int = 100):
        """
        Run system equilibration

        Parameters:
            :param int cpu = number of cpus to be used for the equilibration
            :param string engine = MD engine to be used for equilibration
            :param str cart: labels of cartesian restraints
            :param int cartstr: strength of cartesian restrains

        Class variables:
        """
        self.engine = engine
        self.restarter.write('tleap')

        print(Color.GREEN + 'Setting up system equilibration...\n' + Color.END)
        print('Writing input files in the equilibration folder...\n')


        if self.xyz is None:
            self.xyz = XYZ(self.db_file, self.db_metal_file)
            self.xyz.readFilenames(self.MCPB)
            if len(self.xyz.filenames) == 1:
                self.hasMetal = False
        self.checkTop()
        shutil.copyfile(self.MCPB + '/LIG_solv.inpcrd', self.inputpath + '/equilibration/00.rst7')
        shutil.copyfile(self.MCPB + '/LIG_solv.inpcrd', self.inputpath + '/equilibration/LIG_solv.inpcrd')
        shutil.copyfile(self.MCPB + '/LIG_solv.prmtop', self.inputpath + '/equilibration/LIG_solv.prmtop')

        restrain = None
        if self.bondsRT == []:
            if os.path.exists(self.inputpath + '/restraints'):
                with open(self.inputpath + '/restraints', 'r') as f:
                    for line in f:
                        self.bondsRT.append(line.split('-')[:-2])
                        self.strenghtRT.append(line.split('-')[-2])
                        self.forceConstants.append(line.split('-')[-1])

        if self.bondsRT != []:
            restrain = self.setupRestraint()

        if cart is not None:
            print('Restraining ids {} to cartesian coordinates with the strength of {}'.format(cart,cartstr))
            if cart == 'all':
                cart = '1-{}'.format(len(self.xyz.files))
        self.MDEngine = MDEngine(self.MCPB, engine = self.engine)
        print('Done!\n')
        print(Color.GREEN + 'Starting equilibration...' + Color.END)
        self.status = self.MDEngine.run(self.inputpath, cpus = cpu, restrain=restrain, cart = cart, cartstr = cartstr)


        if self.status == 0:
            error('Equilibration')
            print('If error stems from grid changes, you can modify the input file for the offending step and restart '
                  'the equilibration with the provided restart script in\n {}\n'.format(self.inputpath + '/equilibration'))
            return 0
        self.restarter.write('equilibration')
        return 1

    def setupRestraint(self):
        '''
        Create the restraint files needed for a simulation with amber #TODO modify to work with GROMACS
        :return: restraint string to be added to simulation
        '''

        restrain = '\n&wt TYPE=\'END\' /\nDISANG=disang.r\n'

        print('Generating restraints file\n')
        restraintTemplate = '''parm LIG_dry.prmtop
        reference LIG_dry.pdb
        rst {} reference offset {} rk2 {} rk3 {} out disang.{} width {}
        run
        quit'''
        restraintFile = []
        cpptraj = CPPtraj()
        counter = 0

        for restraint in self.bondsRT:
            atoms = ''

            for unit in restraint:
                atoms = atoms + ':1@{} '.format(unit)
            with open(self.MCPB + '/restraint.in', 'w') as f:
                f.write(restraintTemplate.format(atoms, self.strenghtRT[counter],self.forceConstants[counter],self.forceConstants[counter],counter, self.restraintWidth[counter]))
            cpptraj.run('restraint')
            with open(self.MCPB + '/disang.{}'.format(counter), 'r') as f:
                for line in f:
                    restraintFile.append(line)
            counter += 1

        with open(self.MCPB + '/disang.r', 'w') as f:
            for line in restraintFile:
                f.write(line)
        shutil.copyfile(self.MCPB + '/disang.r', self.inputpath + '/equilibration/disang.r')
        shutil.copyfile(self.MCPB + '/disang.r', self.inputpath + '/simulation/disang.r')
        return restrain

    def prepareSimulation(self, solvent: str, engine: str = 'amber', cart: str = None, cartstr: int = 100 ):
        """
        Function to copy scripts and inputs needed for running and analysing the simulation

        Parameters:
            :param string solvent = 3 letter keyword for the solvent used for the simulation
            :param string cart: id of the residues/atoms that should be restrained
            :param int cartstr: strength of the cartesian restraints in kcal/mol


        """
        self.restarter.write('equilibration')
        print('Preparing simulation...\n')

        sourceloc = os.path.split(__file__)[0]
        try:
            # PyConSolv / scripts_and_inputs
            shutil.copyfile(self.MCPB + '/LIG_dry.prmtop', self.inputpath + '/simulation/LIG_dry.prmtop')
            shutil.copyfile(self.MCPB + '/LIG_solv.prmtop', self.inputpath + '/simulation/LIG_solv.prmtop')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/align_dry.in', self.inputpath + '/simulation/align_dry.in')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/dry_sim.in', self.inputpath + '/simulation/dry_sim.in')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/align.in', self.inputpath + '/simulation/align.in')
            # shutil.copyfile(sourceloc + '/scripts_and_inputs/cluster_kmeans.in', self.inputpath + '/simulation/cluster_kmeans.in')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/dry.vmd', self.inputpath + '/simulation/dry.vmd')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/solv.vmd', self.inputpath + '/simulation/solv.vmd')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/dry_aligned.vmd', self.inputpath + '/simulation/dry_aligned.vmd')
            # shutil.copyfile(sourceloc + '/scripts_and_inputs/strip.sh', self.inputpath + '/simulation/strip.sh')
        except:
            print('Failed to copy files into simulation folder')
            return 1

        if engine == 'gromacs':
            if self.MDEngine == None:
                self.MDEngine = self.MDEngine = MDEngine(self.MCPB, engine = engine)
                self.MDEngine.MD.checkpath()

            executable = self.MDEngine.MD.executable

            shutil.copyfile(self.MCPB + '/LIG_dry.top', self.inputpath + '/simulation/LIG_dry.top')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/simulation_gro.mdp', self.inputpath + '/simulation/simulation_gro.mdp')
            shutil.copyfile(self.inputpath + '/equilibration/LIG_solv.top', self.inputpath + '/simulation/LIG_solv.top')
            shutil.copyfile(self.inputpath + '/equilibration/index.ndx', self.inputpath + '/simulation/index.ndx')
            shutil.copyfile(self.inputpath + '/equilibration/npt.gro', self.inputpath + '/simulation/npt.gro')
            shutil.copyfile(self.inputpath + '/equilibration/npt.cpt', self.inputpath + '/simulation/npt.cpt')

            f = open(self.inputpath + '/simulation/run_simulation_gro.sh', 'w')
            f.write('{} grompp -f simulation_gro.mdp -n index.ndx -c npt.gro -t npt.cpt -p LIG_solv.top -o sim-01.tpr\n'.format(executable))
            f.write('{} mdrun -deffnm sim-01 -nb gpu\n'.format(executable))
            f.close()

        else:
            shutil.copyfile(self.inputpath + '/equilibration/21.rst7', self.inputpath + '/simulation/eq.rst7')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/run_simulation.sh', self.inputpath + '/simulation/run-simulation.sh')
            if self.bondsRT != []:
                shutil.copyfile(sourceloc + '/scripts_and_inputs/simulation_restraint.in',
                                self.inputpath + '/simulation/simulation.in')
            elif cart is not None:
                if self.xyz is None:
                    self.xyz = XYZ(self.db_file, self.db_metal_file)
                    self.xyz.readFilenames(self.MCPB)
                if cart == 'all':
                    cart = '1-{}'.format(len(self.xyz.filenames))
                shutil.copyfile(sourceloc + '/scripts_and_inputs/simulation_cart.in',
                                self.inputpath + '/simulation/simulation.in')
                tmp = []
                with open(self.inputpath + '/simulation/simulation.in', 'r') as f:
                    ln = '   restraint_wt = {}, restraintmask = "!@H=&:{}",\n'.format(cartstr, cart)
                    for line in f:
                        if 'here' in line:
                            tmp.append(ln)
                        else:
                            tmp.append(line)
                with open(self.inputpath + '/simulation/simulation.in', 'w') as f:
                    for line in tmp:
                        f.write(line)
            else:
                shutil.copyfile(sourceloc + '/scripts_and_inputs/simulation.in', self.inputpath + '/simulation/simulation.in')

        solv = Solvent()
        if solvent in self.watermodels:
            solvent = 'Water'
        solvID = solv.solventDict[solvent]
        self.modifyDryScript(self.inputpath + '/simulation/dry_sim.in', solvID)
        self.modifyDryScript(self.inputpath + '/simulation/solv.vmd', solvID)
        self.restarter.write('DONE')
        print('Simulation setup complete, please execute the run_simulation.sh script in:\n {}\n to '
              'begin a 100ns cmd production run.\n'.format(self.inputpath + '/simulation'))
        print('A quick analysis of the simulation run can be performed using the \"pyconsolv sim-01 -a\" command, in your simulation folder\n\n')

        if self.addSolvent == True:
            Copier(self.inputpath+'/Solvent/solv_param/SLV.frcmod',
                   self.solventPath+'/{}.frcmod'.format(self.solventAbb)).copy()
            Copier(self.inputpath + '/Solvent/solv_param/SLV.mol2',
                   self.solventPath + '/{}.mol2'.format(self.solventAbb)).copy()

        print(Color.GREEN + 'My job here is done!' + Color.END)

        return 0



    def modifyDryScript(self, path, solvent):
        """
                Modify the dry script to account for the correct solvent.

                Parameters:
                    :param str path: location of dry_sim.in cpptraj script
                    :param str solvent: new solvent label
        """
        f = open(path, 'r')
        tmp = []
        for line in f:
            if 'SLV' in line:
                line = line.replace('SLV', solvent)
            tmp.append(line)
        f.close()

        f = open(path, 'w')
        for line in tmp:
            f.write(line)
        f.close()

    def metalCheck(self):
        """
        Function to check whether a metal is present in the input structure

        """
        f = open(self.path, 'r')
        for line in f:
            if not line.split():
                continue
            elif line.split()[0].upper() in self.metals:
                self.hasMetal = True
                break
        f.close()

    def checkRT(self):
        '''
        This function will require input from the user to define restraints to be used in the simulation.
        :return:
        '''
        Mapper = self.mapper()
        self.bondsRT = []
        self.strenghtRT = []
        self.forceConstants = []
        self.restraintWidth = []
        stop = False
        print('You have selected to perform a restrained simulation. In this approach, the position of the atoms you '
              'selected will be held in place via the utilization of virtual bonds, angles or dihedrals between them,'
              ' with a very high potential')

        while not stop:
            bond = input(Color.CYAN + 'Please input the virtual bonds in the format "atomid 1 - atomid 2" e.g. 1-2:\n' + Color.END)
            if (len(bond.split('-')) < 2 or len(bond.split('-')) > 4) and bond != 'n':
                print(Color.RED + 'Input is incorrect! Try again\n' + Color.END)
                continue

            if bond != 'n':
                try:
                    strength = float(input('Value to offset distance/angle/torsion in reference by (default 0):\n') or "0")
                except:
                    print(Color.RED + 'Input is incorrect! Try again\n' + Color.END)
                    continue
                try:
                    force = float(input('Please enter distance/angle/torsion restraint force constant (kcal/mol, default 30):\n') or "30")
                except:
                    print(Color.RED + 'Input is incorrect! Try again\n' + Color.END)
                    continue
                try:
                    width = float(input('Please enter distance/angle/torsion restraint width (Angstrom, default 0):\n') or "0")
                except:
                    print(Color.RED + 'Input is incorrect! Try again\n' + Color.END)
                    continue
                tmp = []
                for atom in bond.split('-'):
                    tmp.append(str(Mapper.mapReference(int(atom))))
                self.bondsRT.append(tmp)
                self.strenghtRT.append(strength)
                self.forceConstants.append(force)
                self.restraintWidth.append(width)
                print('Added restraint, enter "n" to stop\n')
            else:
                stop = True
        with open(self.inputpath + '/restraints', 'w') as f:
            for i in range(len(self.bondsRT)):
                f.write('{}-{}-{}\n'.format('-'.join(self.bondsRT[i]), self.strenghtRT[i],self.forceConstants[i]))

    def mapper(self) -> XYZMapper:
        '''
        Map optimized XYZ file back to the original XYZ file
        :return: XYZMapper object to help remap the MCPB.py optimized structure back to original input XYZ file atom ordering
        '''
        Mapper = XYZMapper(self.inputpath + '/input.xyz.original')
        Mapper.mapXYZ(self.inputpath + '/input.xyz')
        return Mapper




    def run(self, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP', dsp: str = 'D4', cpu: int = 12,
            solvent: str = 'Water', multiplicity: int = 1, engine: str = 'amber', opt: bool = True, box: int = 20, rst: bool = False,
            cart: str = None, cartstr: int = 100):
        """
        Run the conformer generation

        Parameters:
            :param str solvent: solvent to be used for MD simulation, check solvent list for available options
            :param int charge: charge of the complete system
            :param string method: ORCA 5 method line
            :param string basis: Basis set for ORCA calculations
            :param string dsp: Dispersion corrections
            :param int cpu: number of CPU cores to be used
            :param str solvent: solvent to be used
            :param int multiplicity: multiplicity of the system
            :param str engine: MD engine to be used for equilibration/simulation
            :param bool opt : if set to False, no geometry optimization will be performed
            :param int box : set box size for amber tleap
            :param bool rst : set if the simulation is of a transition state
            :param str cart : which coordinates should be restrained at a cartesian level
            :param int cartstr : strength of the cartesian coordinate restraints

        Class variables:
        """
        print(Color.GREEN + 'Entering initial setup...\n\n' + Color.END)
        print(cart,cartstr)


        self.checkRestart()
        self.setup(charge, method, basis, dsp, solvent, cpu, multiplicity, opt)
        if rst:
            self.checkRT()

        if self.restart < 2:
            if self.orca(opt = opt) == 0:
                return

        if self.restart < 3:
            if self.antechamber() == 0:
                return
            if not self.hasMetal:
                print('No metal has been found in your file, switching to simple parametrization')
                self.multiwfn(cpu)
                self.restart = 6  # skip mcpb and multiwfn if no metal

        if self.restart < 5:
            if self.multiwfn(cpu) == 0:
                return

        if self.restart < 6:
            if self.MCPB_script() == 0:
                return

        if self.restart < 7:
            if self.tleap(solvent, box) == 0:
                return

        if self.restart < 8:
            if self.equilibration(cpu, engine, cart=cart, cartstr = cartstr) == 0:
                return

        if self.restart < 9:
            if self.prepareSimulation(solvent, engine, cart = cart, cartstr=cartstr) == 0:
                return
        if self.restart == 9:
            print('This structure has already completed parametrization, please make sure you are using the correct input\n')