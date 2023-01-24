import shutil
import os
from tkinter import *
import numpy as np

from .colorgen import Color
from .amber import amberInterface
from .calculate import Calculation
from .filestructure import Setup
from .inputparser import XYZ
from .multiWFN import MultiWfnInterface
from .outgen import Faker
from .ui import GUI
from .restart import RestartFile


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
    """

    def __init__(self, path):
        self.hasMetal = None
        self.restarter = None
        self.path = path
        self.inputpath = '/'.join(path.split('/')[:-1])
        self.status = 0
        self.restart = 0
        self.db_file = os.path.split(__file__)[0] + '/db/atom-radius.txt'
        self.db_metal_file = os.path.split(__file__)[0] + '/db/metal-radius.txt'
        self.amber = None
        self.MCPB = self.inputpath + '/MCPB_setup'
        self.xyz = None

        print(Color.BLUE + r'''

          _____        _____             _____       _       
         |  __ \      / ____|           / ____|     | |      
         | |__) |   _| |     ___  _ __ | (___   ___ | |_   __
         |  ___/ | | | |    / _ \| '_ \ \___ \ / _ \| \ \ / /
         | |   | |_| | |___| (_) | | | |____) | (_) | |\ V / 
         |_|    \__, |\_____\___/|_| |_|_____/ \___/|_| \_/  
                 __/ |                                       
                |___/                                        
                                Ver 0.1.1
                    
Welcome to PyConSolv, your friendly neighbourhood conformer generator

Calculations will be set up in:

{}

        '''.format(self.inputpath) + Color.END)

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

    def setup(self, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP', dsp: str = 'D4',
              cpcm: str = 'Water', cpu: int = 12) -> int:
        """
        Run setup for creating the appropriate folders

        Parameters:
            :param int charge: charge of the complete system
            :param string method: ORCA 5 method line
            :param string basis: Basis set for ORCA calculations
            :param string dsp: Dispersion corrections
            :param string cpcm: CPCM solvation model solvent
            :param int cpu: number of CPU cores to be used

        Class variables:
        """
        if self.restart == 0:
            setup = Setup(self.path, charge=charge)
            setup.Method(method, basis, dsp, cpcm, cpu)
            self.status = setup.run()
            if self.status == 0:
                error('Setup')
                return 0
            print(Color.GREEN + 'Setup is complete, moving on to ORCA calculations...\n' + Color.END)

        self.restarter = RestartFile(self.inputpath)
        return 1

    def orca(self):
        """
        Run ORCA optimization and frequency calculations

        Parameters:

        Class variables:
        """
        self.restarter.write('setup')
        calculation = Calculation(self.inputpath + '/orca_calculations')
        if self.hasMetal:
            self.status = calculation.run()
        else:
            self.status = calculation.run(freq=False)
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
        shutil.copyfile(self.inputpath + '/orca_calculations/freq/input.xyz', self.inputpath + '/MCPB_setup/input.xyz')

        self.xyz = XYZ(self.db_file, self.db_metal_file)
        self.xyz.prepareMCPB(self.inputpath + '/MCPB_setup/input.xyz')
        pdbs = self.xyz.filenames

        print('The following pdb files were created: \n {}\n'.format(' '.join(pdbs)))
        print('Please enter fragment charges:\n')

        # get charges from user
        window = Tk()
        start = GUI(window, self.MCPB, pdbs)
        window.mainloop()

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
            os.chdir(self.inputpath + '/equilibration/')
            shutil.copyfile(self.MCPB + '/' + str(antechamberFiles[0][0]) + '.mol2',
                            self.inputpath + '/equilibration/LIG.mol2')
            shutil.copyfile(self.MCPB + '/' + antechamberFiles[0][0] + '.frcmod',
                            self.inputpath + '/equilibration/LIG.frcmod')
            self.amber.tleapNoMetalSolv(self.inputpath + '/equilibration/')
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
            self.xyz.path = (self.inputpath + '/MCPB_setup/')
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
            :param string solvent: solvent for your box. Currently only water is available

        Class variables:
        """
        self.restarter.write('mcpb')

        print('Solvent of choice is: {}'.format(solvent))

        if self.amber is None:
            self.amber = amberInterface(self.MCPB)

        if not self.hasMetal:
            self.amber.tleapNoMetalSolv(self.MCPB)  # change to get values from ligands
        self.amber.tleapChecker(self.MCPB)
        self.status = self.amber.runTleap()

        if self.status == 0:
            error('Tleap')
            return 0

        self.restarter.write('tleap')

        print('Parametrization complete!\n')
        print('Solvent of choice is: {}'.format(solvent))
        return 1

    def equilibration(self):
        """
        Run system equilibration

        Parameters:

        Class variables:
        """
        self.restarter.write('tleap')

        print(Color.GREEN + 'Setting up system equilibration...\n' + Color.END)
        print('Writing input files in the equilibration folder...\n')
        shutil.copyfile(self.inputpath + '/MCPB_setup/LIG_solv.inpcrd', self.inputpath + '/equilibration/00.rst7')
        shutil.copyfile(self.inputpath + '/MCPB_setup/LIG_solv.prmtop',
                        self.inputpath + '/equilibration/LIG_solv.prmtop')

        if self.amber is None:
            self.amber = amberInterface(self.MCPB)
        self.amber.equil(self.inputpath)

        print('Done!\n')
        print(Color.GREEN + 'Starting equilibration...' + Color.END)

        self.status = self.amber.equilibrate()

        if self.status == 0:
            error('Equilibration')
            return 0
        self.restarter.write('equilibration')
        return 1

    def run(self, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP', dsp: str = 'D4',
            cpcm: str = 'Water', cpu: int = 12, solvent: str = 'Water'):
        """
        Run the conformer generation

        Parameters:
            :param str solvent: solvent to be used for MD simulation, check solvent list for available options
            :param int charge: charge of the complete system
            :param string method: ORCA 5 method line
            :param string basis: Basis set for ORCA calculations
            :param string dsp: Dispersion corrections
            :param string cpcm: CPCM solvation model solvent
            :param int cpu: number of CPU cores to be used

        Class variables:
        """
        print(Color.GREEN + 'Entering initial setup...\n\n' + Color.END)

        self.checkRestart()
        self.setup(charge, method, basis, dsp, cpcm, cpu)

        if self.restart < 2:
            if self.orca() == 0:
                return

        if self.restart < 3:
            if self.antechamber() == 0:
                return
            if not self.hasMetal:
                print('No metal has been found in your file, switching to simple parametrization')
                self.restart = 6# skip mcpb and multiwfn if no metal

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
            if self.equilibration() == 0:
                return
