import shutil
import os
from tkinter import *
import numpy as np

from .interfaces.mdengines import MDEngine
from .misc.counterion import Counterion
from .misc.counterionGen import counterionParametrizer
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
    """

    def __init__(self, path):
        self.engine = None
        self.MDEngine = None
        path = os.path.abspath(path)
        self.counterIon = ''
        self.refrac = None
        self.epsilon = None
        self.solventParamPath = None
        self.version = '0.9.3.0'
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
        self.counterIonsImplemented = ['OTf-','BF4-', 'BARF-', 'PFC-', 'ScF6-', 'ClO4-', 'BPh4-', 'custom']
        self.solventsImplemented = ['Water', 'Acetonitrile', 'Acetone', 'Benzene', 'Cyclohexane', 'Chloroform', 'CCl4',
                                    'CH2Cl2', 'DMF', 'DMSO', 'Ethanol', 'Hexane', 'Methanol', 'Ammonia', 'Octanol',
                                    'THF', 'Toluene', 'custom']

        self.metalCheck()
        self.addSolvent = False
        self.solventAbb = ''


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
                        answered = True
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
              cpcm: str = 'Water', cpu: int = 12, multiplicity:int  = 1) -> int:
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

        if charge != 0:
            self.prepareCounterion(method, basis, dsp, cpu, charge, cpcm)

        if self.restart == 0:
            self.xyz = XYZ(self.db_file, self.db_metal_file)
            self.xyz.prepareInput(self.inputpath + '/input.xyz')
            self.xyz = None
            setup = Setup(self.inputpath + '/' + self.inputFile, charge=charge, multi = multiplicity)
            setup.Method(method, basis, dsp, cpcm, cpu, self.epsilon, self.refrac)
            self.status = setup.run()
            if self.status == 0:
                error('Setup')
                return 0

        if cpcm == 'custom':
            shutil.copyfile('/'.join(self.solventParamPath.split('/')[:-1]) + '/solv_param/SLV.frcmod',
                            self.MCPB + '/SLV.frcmod')
            shutil.copyfile('/'.join(self.solventParamPath.split('/')[:-1]) + '/solv_param/SLV.mol2',
                            self.MCPB + '/SLV.mol2')
        else:
            solvent = Solvent()
            solvname = solvent.solventDict[cpcm]
            if cpcm != 'Water':
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
            ionname = cIon.counterionDict[self.counterIon]
            shutil.copyfile(self.ionPath + '/{}.frcmod'.format(ionname),
                            self.MCPB + '/{}.frcmod'.format(ionname))
            shutil.copyfile(self.ionPath + '{}.mol2'.format(ionname),
                            self.MCPB + '/{}.mol2'.format(ionname))

        f = open(self.inputpath + '/simulation/solvent', 'w')
        f.write(cpcm)
        f.close()

        print(Color.GREEN + 'Setup is complete, moving on to ORCA calculations...\n' + Color.END)

        self.restarter = RestartFile(self.inputpath)

        return 1

    def orca(self):
        """
        Run ORCA optimization and frequency calculations

        Parameters:

        Class variables:
        """
        try:
            self.restarter.write('setup')
        except AttributeError:
            return 0

        calculation = Calculation(self.inputpath + '/orca_calculations')
        if self.hasMetal:
            self.status = calculation.run()
        else:
            self.status = calculation.run(freq=False)
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
            self.xyz.readFilenames(self.MCPB)  # todo this might not be needed
        if not self.hasMetal:
            multiwfn = MultiWfnInterface(self.inputpath + '/orca_calculations/opt/', orcaname='orca_opt')
            self.status = multiwfn.run(cores)
            self.xyz.hasMetal = False
            self.xyz.readRESP(self.inputpath + '/orca_calculations/')
            chargeChanger = ChargeChanger()
            chargeChanger.change(self.MCPB + '/A.mol2', self.MCPB + '/LIG.mol2', 'A', self.xyz.charges)
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

        self.restarter.write('mcpb')
        return 1

    def tleap(self, solvent: str) -> int:
        """
        Run tleap

        Parameters:
            :param string solvent: solvent for your box.

        Class variables:
        """
        self.restarter.write('mcpb')

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
        if self.counterIon != '':
            ion = Counterion()
            ion.applyItem(self.counterIon, self.MCPB + '/LIG_tleap.in',
                              self.MCPB + '/LIG_tleap.in', self.MCPB, solutename)
        self.amber.tleapChecker(self.MCPB)
        self.status = self.amber.runTleap()

        if self.status == 0:
            error('Tleap')
            return 0

        self.restarter.write('tleap')

        print('Parametrization complete!\n')
        print('Solvent of choice is: {}\n'.format(solvent))
        return 1

    def equilibration(self, cpu: int = 12, engine = 'amber'):
        """
        Run system equilibration

        Parameters:
            :param int cpu = number of cpus to be used for the equilibration
            :param string engine = MD engine to be used for equilibration

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

        #old amber only equilibration
        # shutil.copyfile(self.MCPB + '/LIG_solv.inpcrd', self.inputpath + '/equilibration/00.rst7') # different for different engines
        # shutil.copyfile(self.MCPB + '/LIG_solv.prmtop', self.inputpath + '/equilibration/LIG_solv.prmtop')
        #
        # if self.amber is None:
        #     self.amber = amberInterface(self.MCPB)
        # self.amber.equil(self.inputpath)
        #
        # print('Done!\n')
        # print(Color.GREEN + 'Starting equilibration...' + Color.END)
        #
        # self.status = self.amber.equilibrate(cpus = cpu)

        shutil.copyfile(self.MCPB + '/LIG_solv.inpcrd', self.inputpath + '/equilibration/00.rst7')
        shutil.copyfile(self.MCPB + '/LIG_solv.inpcrd', self.inputpath + '/equilibration/LIG_solv.inpcrd')
        shutil.copyfile(self.MCPB + '/LIG_solv.prmtop', self.inputpath + '/equilibration/LIG_solv.prmtop')
        self.MDEngine = MDEngine(self.MCPB, engine = self.engine)
        print('Done!\n')
        print(Color.GREEN + 'Starting equilibration...' + Color.END)
        self.status = self.MDEngine.run(self.inputpath, cpus = cpu)


        if self.status == 0:
            error('Equilibration')
            print('If error stems from grid changes, you can modify the input file for the offending step and restart '
                  'the equilibration with the provided restart script in\n {}\n'.format(self.inputpath + '/equilibration'))
            return 0
        self.restarter.write('equilibration')
        return 1

    def prepareSimulation(self, solvent: str):
        """
        Function to copy scripts and inputs needed for running and analysing the simulation

        Parameters:
            :param string solvent = 3 letter keyword for the solvent used for the simulation


        """
        self.restarter.write('equilibration')
        print('Preparing simulation...\n')
        sourceloc = os.path.split(__file__)[0]
        try:
            # PyConSolv / scripts_and_inputs
            shutil.copyfile(self.MCPB + '/LIG_dry.prmtop', self.inputpath + '/simulation/LIG_dry.prmtop')
            shutil.copyfile(self.MCPB + '/LIG_solv.prmtop', self.inputpath + '/simulation/LIG_solv.prmtop')
            shutil.copyfile(self.inputpath + '/equilibration/21.rst7', self.inputpath + '/simulation/eq.rst7')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/align_dry.in', self.inputpath + '/simulation/align_dry.in')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/run_simulation.sh', self.inputpath + '/simulation/run-simulation.sh')
            shutil.copyfile(sourceloc + '/scripts_and_inputs/simulation.in', self.inputpath + '/simulation/simulation.in')
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
        solv = Solvent()
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
            if line.split()[0].upper() in self.metals:
                self.hasMetal = True
                break
        f.close()

    def run(self, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP', dsp: str = 'D4', cpu: int = 12,
            solvent: str = 'Water', multiplicity: int = 1, engine: str = 'amber'):
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

        Class variables:
        """
        print(Color.GREEN + 'Entering initial setup...\n\n' + Color.END)

        self.checkRestart()
        self.setup(charge, method, basis, dsp, solvent, cpu, multiplicity)

        if self.restart < 2:
            if self.orca() == 0:
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
            if self.tleap(solvent) == 0:
                return

        if self.restart < 8:
            if self.equilibration(cpu, engine) == 0:
                return

        if self.restart < 9:
            if self.prepareSimulation(solvent) == 0:
                return
        if self.restart == 9:
            print('This structure has already completed parametrization, please make sure you are using the correct input\n')