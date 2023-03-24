import os
import shutil
import subprocess

from ..utils.colorgen import Color


class MultiWfnInterface:
    def __init__(self, path: str, orcaname: str = 'orca_freq'):
        """
        Interface to MultiWfn

        Parameters:
            :param string path: path to orca calculation

        Class variables:
            - self.path = path to orca calculation
            - self.original_wd = working directory where original calculation was started from
            - self.status = checks if calculations failed or succeeded
        """
        self.path = path
        self.original_wd = os.getcwd()
        os.chdir(self.path)
        self.status = 0
        self.orcaname = orcaname

    def ECPcheck(self):
        """
        Checks ORCA output file for use of ECPs and creates a new molden.input file with the correct number of electrons
        This is needed to ensure RESP charge calculations are correct

        Parameters:

        Class variables:
        """
        ecp = []
        electrons = []
        skip = True

        f = open(self.path + '/' + self.orcaname + '.out', 'r')

        for line in f:
            if 'CARTESIAN COORDINATES (A.U.)' in line:
                skip = False
                continue
            if '----------------------------' in line:
                continue
            if 'ZA' in line:
                continue
            if 'core charge reduced due to ECP' in line:
                break
            if not skip:
                if not line.strip():
                    break # break loop if line is empty
                electrons.append(line.split()[2])
                if '*  ' in line:
                    ecp.append(line.split()[0:3])
                elif '--------------------------------' in line:
                    break
        f.close()
        if len(ecp) > 0:
            skip = True
            template = '{:<2}{:>4}{:>4}{:>22}{:>21}{:>21}\n'
            print('ECPs found, correcting molden input file\n')
            print('Atoms with ECP:')
            for at in ecp:
                print('{} with {} electrons\n'.format(at[1], int(float(at[2].replace('*', '')))))
            shutil.copyfile(self.path + '/' + self.orcaname + '.molden.input',
                            self.path + '/' + self.orcaname + '.molden.input_no_ECP')
            fin = open(self.path + '/' + self.orcaname + '.molden.input_no_ECP', 'r')
            fout = open(self.path + '/' + self.orcaname + '.molden.input', 'w')
            iterator = 0
            for line in fin:
                if '[Atoms] AU' in line:
                    skip = False
                    fout.write(line)
                    continue
                if not skip:
                    l = line.split()
                    l[2] = int(float(electrons[iterator].replace('*', '')))
                    fout.write(template.format(*l))
                    iterator += 1
                    if iterator >= len(electrons):
                        skip = True
                else:
                    fout.write(line)
            fin.close()
            fout.close()

    def generateMultiwfnInput(self):
        """
        Creates an input file for Multiwfn, for RESP charge calculation

        Parameters:
            - inputfile = full path to xyz inputfile

        Class variables:
        """
        template = '''{}
7
18
1




y
0
0
q
'''
        f = open(self.path + './multiwfn.input', 'w')
        f.write(template.format('./' + self.orcaname + '.molden.input'))
        f.close()

    def run(self, threads: int = 12) -> int:
        """
        Runs the MultiWfn RESP charge calculation

        Parameters:
            :param int threads = number of threads to use for calculation

        Class variables:
        """
        self.ECPcheck()
        self.generateMultiwfnInput()
        multiwfnCommand = 'Multiwfn -nt {} -silent < multiwfn.input > multiwfn.out'.format(threads)
        calc = subprocess.run([multiwfnCommand], shell=True)
        if calc.returncode == 0:
            print('RESP charge calculation completed successfully!\n')
            self.status = 1
        else:
            print(Color.RED + 'Something went wrong when running MultiWfn, please check multiwfn.out' + Color.END)
            self.status = 0
        return self.status
