import shutil
import os
from tkinter import *
import numpy as np

from .colorgen import color
from .amber import amberInterface
from .calculate import Calculation
from .filestructure import Setup
from .inputparser import XYZ
from .multiWFN import MultiWfnInterface
from .outgen import Faker
from .ui import GUI
from .restart import RestartFile

class PyConSolv():
    def __init__(self, path):
        self.inputpath = '/'.join(path.split('/')[:-1])
        self.path = path
        self.status = 0
        self.restart = 0
        self.db_file = os.path.split(__file__)[0] + '/db/atom-radius.txt'
        self.db_metal_file = os.path.split(__file__)[0] + '/db/metal-radius.txt'
        self.amber = None
        self.MCPB = self.inputpath+'/MCPB_setup'
        self.xyz = None

        print(color.BLUE + '''

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

        '''.format(self.inputpath) + color.END)
        
    
    def error(self, step):
        print(color.RED + 'Something went wrong, please check your input/output' + color.END)
        print(color.RED + '''
        ############################################
        ######             WARNING            ######
        ######      Calculation failed at:    ######
        ######{:^32}######
        ############################################
        '''.format(step) + color.END)
    
    def run(self, cores = 8, solvent = 'water', charge = 0):
        print(color.GREEN + 'Entering initial setup...\n\n' + color.END)

        #####################
        ########SETUP########
        #####################
        if os.path.exists(self.inputpath+'/pyconsolv.restart'):
            restarter = RestartFile(self.inputpath)
            self.restart = restarter.getstate()
            print(color.PURPLE + 'Restart file found!\n' + color.END)
            os.remove(self.inputpath+'/pyconsolv.restart')
        
        if self.restart == 0:
            setup = Setup(self.path, charge = charge)
            self.status = setup.run()
            if self.status == 0:
                self.error('Setup')
                return
            print(color.GREEN + 'Setup is complete, moving on to ORCA calculations...\n' + color.END)
        
        restarter = RestartFile(self.inputpath)

        #####################
        ########ORCA#########
        #####################
        if self.restart < 2:
            restarter.write('setup')
            calculation = Calculation(self.inputpath+'/orca_calculations')
            self.status = calculation.run()
            if self.status == 0:
                self.error('ORCA Calculations')
                return
        
            print(color.GREEN + 'ORCA Calculations complete, moving on to MCPB setup...' + color.END)
            restarter.write('orca')

        #####################
        ##GUI + Antechamber##
        #####################
        if self.restart < 3:
            restarter.write('orca')
            shutil.copyfile(self.inputpath+'/orca_calculations/freq/input.xyz',self.inputpath+'/MCPB_setup/input.xyz')
        
            self.xyz = XYZ(self.db_file, self.db_metal_file)
            self.xyz.prepareMCPB(self.inputpath+'/MCPB_setup/input.xyz')
            pdbs = self.xyz.filenames
        
            print('The following pdb files were created: \n {}\n'.format(' '.join(pdbs)))
            print('Please enter fragment charges:\n')
        
            # get charges from user
            window = Tk()
            start = GUI (window, self.MCPB, pdbs)
            window.mainloop()
        
            # create mol2 files and run antechamber
            self.xyz.writeMol2Files(self.path)
        
            antechamberFiles = self.xyz.molNotCreated
            metals = self.xyz.metals
            ligands = np.array(self.xyz.ligands)
            self.xyz.writeMetalConnections(self.MCPB)  # write out metal connections file

            self.amber = amberInterface(self.MCPB)
            for filename in antechamberFiles:
                self.status = self.amber.antechamber(*filename)
                if self.status == 0:
                    self.error('antechamber for {}'.format(filename))
                    return

            print('Generating frcmod files for ligands:\n')
            for filename in ligands:
                self.status = self.amber.runParmchk2(filename)
                if self.status == 0:
                    self.error('parmchk2 for {}'.format(filename))
                    return

            self.amber.inputFileGenerator(metals[0][1], ligands[:, 1])
            restarter.write('frcmod')

        #####################
        #######Multiwfn######
        #####################
        if self.restart < 5:
            restarter.write('frcmod')
            print(color.GREEN + 'Fragments have been prepared, running MultiWfn task...\n\n' + color.END)
        
            multiwfn = MultiWfnInterface(self.inputpath+'/orca_calculations/freq/')
            self.status = multiwfn.run(cores)
        
            if self.status == 0:
                self.error('MultiWfn Calculations')
                return
            restarter.write('multiwfn')


        #####################
        #######MCPB.py#######
        #####################
        if self.restart < 6:
            restarter.write('multiwfn')

            print(color.GREEN + 'Converting ORCA output to MCPB.py compatible input...\n' + color.END)
        
            faker = Faker(self.inputpath+'/orca_calculations/freq/')
            faker.fakecrds()
            # faker.fakeesp()
            faker.fakeforce()
        
            shutil.copyfile(self.inputpath+'/orca_calculations/freq/fakechk.fchk',self.inputpath+'/MCPB_setup/LIG_small_opt.fchk')
            shutil.copyfile(self.inputpath+'/orca_calculations/freq/fakelog.log',self.inputpath+'/MCPB_setup/LIG_small_fc.log')
        
            print(color.GREEN + 'Proceeding with MCPB steps...\n' + color.END)

            if self.amber == None:
                self.amber = amberInterface(self.MCPB)

            self.status = self.amber.runMCPB('1')
        
            if self.status == 0:
                self.error('MCPB step 1')
                return

            if self.amber.checkMCPBBonds('/home/rat/PyPer/mo-cat/MCPB_setup') == False:
                self.status = self.amber.runMCPB('1')

                if self.status == 0:
                    self.error('MCPB step 1')
                    return

            self.status = self.amber.runMCPB('2')
        
            if self.status == 0:
                self.error('MCPB step 2')
                return

            if self.xyz == None:
                self.xyz = XYZ(self.db_file, self.db_metal_file)
            self.xyz.createFinalMol2(self.inputpath)
        
        
            self.status = self.amber.runMCPB('4')
        
            if self.status == 0:
                self.error('MCPB step 4')
                return
        
            restarter.write('mcpb')

        #####################
        #######Tleap#########
        #####################
        if self.restart < 7:
            
            restarter.write('mcpb')

            print('Solvent of choice is: {}'.format(solvent))

            if self.amber == None:
                self.amber = amberInterface(self.MCPB)

            self.status = self.amber.runTleap()
        
            if self.status == 0:
                self.error('Tleap')
                return
        
            restarter.write('tleap')
        
            print('Parametrization complete!\n')
            print('Solvent of choice is: {}'.format(solvent))
        
        
        if self.restart < 8:
            restarter.write('tleap')

            print(color.GREEN + 'Setting up system equilibration...\n' + color.END)
            print('Writing input files in the equlibration folder...\n')
            shutil.copyfile(self.inputpath+'/MCPB_setup/LIG_solv.inpcrd',self.inputpath+'/equilibration/00.rst7')
            shutil.copyfile(self.inputpath+'/MCPB_setup/LIG_solv.prmtop',self.inputpath+'/equilibration/LIG_solv.prmtop')
            self.amber.equil(self.inputpath)
        
            print('Done!\n')
            print(color.GREEN + 'Starting equilibration...' + color.END)


            if self.amber == None:
                self.amber = amberInterface(self.MCPB)

            self.status = self.amber.equilibrate()

            if self.status == 0:
                self.error('Equilibration')
                return
            restarter.write('equilibration')
