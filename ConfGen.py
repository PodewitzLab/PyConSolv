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
        #self.db_file = r'/home/rat/PyPer/Nov_22/atom-radius.txt' # to be changed to the path of the python package
        
        print(color.BLUE + '''

          _____        _____             _____       _       
         |  __ \      / ____|           / ____|     | |      
         | |__) |   _| |     ___  _ __ | (___   ___ | |_   __
         |  ___/ | | | |    / _ \| '_ \ \___ \ / _ \| \ \ / /
         | |   | |_| | |___| (_) | | | |____) | (_) | |\ V / 
         |_|    \__, |\_____\___/|_| |_|_____/ \___/|_| \_/  
                 __/ |                                       
                |___/                                        
                                Ver 0.1
                    
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
    
    def run(self, cores = 8, solvent = 'water'):
        print(color.GREEN + 'Entering initial setup...\n\n' + color.END)
        
        if os.path.exists(self.inputpath+'/pyconsolv.restart'):
            restarter = RestartFile(self.inputpath)
            self.restart = restarter.getstate()
            print(color.PURPLE + 'Restart file found!\n' + color.END)
            os.remove(self.inputpath+'/pyconsolv.restart')
        
        if self.restart == 0:
            setup = Setup(self.path)
            self.status = setup.run()
            if self.status == 0:
                self.error('Setup')
                return
            print(color.GREEN + 'Setup is complete, moving on to ORCA calculations...\n' + color.END)
        
        restarter = RestartFile(self.inputpath)
        if self.restart < 2:
            restarter.write('setup')
            calculation = Calculation(self.inputpath+'/orca_calculations')
            self.status = calculation.run()
            if self.status == 0:
                self.error('ORCA Calculations')
                return
        
            print(color.GREEN + 'ORCA Calculations complete, moving on to MCPB setup...' + color.END)
            restarter.write('orca')
        
        if self.restart < 3:
            restarter.write('orca')
            shutil.copyfile(self.inputpath+'/orca_calculations/freq/input.xyz',self.inputpath+'/MCPB_setup/input.xyz')
        
            xyz = XYZ(self.db_file)
            xyz.prepareMCPB(self.inputpath+'/MCPB_setup/input.xyz')
            pdbs = xyz.filenames
        
            print('The following pdb files were created: \n {}\n'.format(' '.join(pdbs)))
        
            MCPB = self.inputpath+'/MCPB_setup'
        
            print('Please enter fragment charges:\n')
        
            # get charges from user
            window= Tk()
            start= GUI (window, MCPB, pdbs)
            window.mainloop()
        
            # create mol2 files and run antechamber
            xyz.writeMol2Files(self.path)
        
            antechamberFiles = xyz.molNotCreated
            metals = xyz.metals
            ligands = np.array(xyz.ligands)

            amber = amberInterface(MCPB)
            for filename in antechamberFiles:
                self.status = amber.antechamber(*filename)
                if self.status == 0:
                    self.error('antechamber for {}'.format(filename))
                    return
            restarter.write('antechamber')
        
        if self.restart < 4:
            restarter.write('antechamber')
            print('Generating frcmod files for ligands:\n')
            for filename in ligands:
                self.status = amber.runParmchk2(filename)
                if self.status == 0:
                    self.error('parmchk2 for {}'.format(filename))
                    return
        
            restarter.write('frcmod')
        
        if self.restart < 5:
            restarter.write('frcmod')
            amber.inputFileGenerator( metals[0][1], ligands[:,1])
            print(color.GREEN + 'Fragments have been prepared, running MultiWfn task...\n\n' + color.END)
        
            multiwfn = MultiWfnInterface(self.inputpath+'/orca_calculations/freq/')
            self.status = multiwfn.run(cores)
        
            if self.status == 0:
                self.error('MultiWfn Calculations')
                return
            restarter.write('multiwfn')        

        if self.restart < 6:
            restarter.write('multiwfn')

            print(color.GREEN + 'Converting ORCA output to MCPB.py compatible input...\n' + color.END)
        
            faker=Faker(self.inputpath+'/orca_calculations/freq/')
            faker.fakecrds()
            # faker.fakeesp()
            faker.fakeforce()
        
            shutil.copyfile(self.inputpath+'/orca_calculations/freq/fakechk.fchk',self.inputpath+'/MCPB_setup/LIG_small_opt.fchk')
            shutil.copyfile(self.inputpath+'/orca_calculations/freq/fakelog.log',self.inputpath+'/MCPB_setup/LIG_small_fc.log')
        
            print(color.GREEN + 'Proceeding with MCPB steps...\n' + color.END)
        
            self.status = amber.runMCPB('1')
        
            if self.status == 0:
                self.error('MCPB step 1')
                return
        
            self.status = amber.runMCPB('2')
        
            if self.status == 0:
                self.error('MCPB step 2')
                return
        
            xyz.createFinalMol2(self.inputpath)
        
        
            self.status = amber.runMCPB('4')
        
            if self.status == 0:
                self.error('MCPB step 4')
                return
        
            restarter.write('mcpb')
            
        if self.restart < 7:
            
            restarter.write('mcpb')

            print('Solvent of choice is: {}'.format(solvent))
            self.status = amber.runTleap()
        
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
            amber.equil(self.inputpath)
        
            print('Done!\n')
            print(color.GREEN + 'Starting equilibration...' + color.END)
        
            self.status = amber.equilibrate()     
            if self.status == 0:
                self.error('Equilibration')
                return
            restarter.write('equilibration')
