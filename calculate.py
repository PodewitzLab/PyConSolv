import os
import shutil
import subprocess

from .colorgen import color


class Calculation:
    def __init__(self, path):  #
        """
        Run ORCA calculations

        Parameters:
            - path = path that contains the input file. This will act as the root directory for the calculations.
                     Geometry optimizations will be performed in path/opt, frequency calculations in path/freq

        Class variables:
            - self.orcapath = path to call ORCA executable
            - self.path = root path for calculations
            - self.original_wd = current working directory when calculations are started
            - self.status = status of the calculation. 0 means an error occured and everything should be stopped
        """
        self.orcapath = None
        self.path = path
        self.original_wd = os.getcwd()
        self.status = 0

    def checkpath(self):
        """
        Check if ORCA is available in PATH

        Parameters:

        Class variables:
        """
        self.status = 0
        self.orcapath = ''
        systemPATH = os.environ.get('PATH', '').split(os.pathsep)
        for el in systemPATH:
            if 'orca' in el or 'ORCA' in el:
                print('Found ORCA in: ' + el)
                self.orcapath = el + '/orca '
                self.status = 1
                break
        if self.orcapath == '':
            print(color.RED + 'ORCA was not found on your system... is it in your PATH?' + color.END)
            self.status = 0
            return
        else:
            return

    def calculate(self, calctype):
        """
        Run orca calculations

        Parameters:
            - calctype = type of calculation, opt or freq

        Class variables:
        """
        if calctype == 'opt':
            output = 'orca_opt.out'
            inputfile = 'orca_opt.inp'
            os.chdir(self.path + '/opt')
            print('Running geometry optimization in ' + os.getcwd())
        elif calctype == 'freq':
            output = 'orca_freq.out'
            inputfile = 'orca_freq.inp'
            os.chdir(self.path + '/freq')
            print('Running frequency calculation in ' + os.getcwd())
            shutil.copyfile(self.path + '/opt/orca_opt.xyz', self.path + '/freq/input.xyz')
        else:
            print(color.RED + 'Unrecognized keyword for calculation!' + color.END)
            self.status = 0
            return

        command = self.orcapath + inputfile + ' > ' + output
        # command = [self.orcapath,inputfile]
        print(command)
        f = open('run_calc.sh', 'w')
        f.write('#!/bin/bash\n')
        f.write(command)
        f.close()
        script = './run_calc.sh'
        subprocess.run(['chmod u+x run_calc.sh'], shell=True)
        calc = subprocess.run([script], stdin=None)
        if calc.returncode == 0:
            print('''
            
Calculation completed successfully!
Moving on!

            ''')
            if calctype == 'freq':
                print('Generating molden input file from calculation..\n')
                command = 'orca_2mkl orca_freq -molden'
                calc = subprocess.run([command], shell=True)
                if calc.returncode == 0:
                    self.status = 1
                else:
                    print(color.RED + 'Could not create molden input file...\n' + color.END)
                    self.status = 0
                    os.chdir(self.original_wd)
                    return
            self.status = 1
            os.chdir(self.original_wd)
            return
        else:
            print(
                color.RED + 'Something went wrong with the ORCA calculation, please check output files in '
                + os.getcwd() + color.END)
            os.chdir(self.original_wd)
            self.status = 0
            return

    def run(self):
        """
        Run all ORCA calculations

        Parameters:

        Class variables:

        """
        self.checkpath()
        if self.status == 0:
            print(color.RED + 'Aborting calculation!' + color.END)
            return self.status
        self.calculate(calctype='opt')
        if self.status == 0:
            print(color.RED + 'Aborting calculation!' + color.END)
            return self.status
        self.calculate(calctype='freq')
        if self.status == 0:
            print(color.RED + 'Aborting calculation!' + color.END)
            return self.status
        print('ORCA calculations completed successfully!\n')
        os.chdir(self.original_wd)
        return self.status
