import shutil
import os
from .colorgen import Color


class Setup:
    def __init__(self, path: str, charge: int = 0, multi: object = 1) -> NotImplemented:
        """
        Class for setting up all the necessary input files and folders for the generation of the parameters

        Parameters:
            - path = path to the input file in XMOL format
            - charge = charge of the structure, default is 0
            - multi = multiplicity, default is 1

        Class variables:
            - self.db_file = path to database file which contains atom radius information
            - self.path = folder containing the input file given as an argument
            - self.inputfile = name of the input file
            - self.charge = system charge
            - self.multi = system multiplicity
            - self.orca_inp = ORCA calculation input file template
            - self.check = if there are other files in the folder, warning will be displayed and operations stopped
        """
        self.db_file = os.path.split(__file__)[0] + '/db/atom-radius.txt'
        self.path = '/'.join(path.split('/')[:-1])
        self.inputFile = path.split('/')[-1]
        self.charge = charge
        self.multi = multi
        self.check = 1
        # ORCA input file template
        self.orca_inp = '''{}
{}
{}

%PAL NPROCS {} END

%scf
maxiter 350
end

* xyzfile {} {} input.xyz
'''

    def createFolders(self):
        """
        Create folder structure for calculations

        Parameters:

        Class variables:
            - self.calculation_folder = folder where all orca calculations will be performed
            - self.frequency_folder = folder where orca frequency calculations will be performed
            - self.optimization_folder = folder where orca geometry optimization calculations will be performed
            - self.MCPB_folder = folder where MCPB setup calculations will be performed
            - self.equilibration_folder = folder where system equilibration will be performed
            - self.simulation_folder = folder where a CMD will be performed
        """

        self.calculation_folder = self.path + '/' + 'orca_calculations'
        self.equilibration_folder = self.path + '/' + 'equilibration'
        self.simulation_folder = self.path + '/' + 'simulation'
        self.frequency_folder = self.calculation_folder + '/' + 'freq'
        self.optimization_folder = self.calculation_folder + '/' + 'opt'
        self.MCPB_folder = self.path + '/' + 'MCPB_setup'

        if len(os.listdir(self.path)) > 2:
            self.check = 0
            print(
                Color.RED + 'Your calculation directory includes unexpected folders/files. Please make sure it only contains your input xyz file!\n' + Color.END)
        else:
            print('Creating folders...')
            os.mkdir(self.calculation_folder)
            os.mkdir(self.equilibration_folder)
            os.mkdir(self.simulation_folder)
            os.mkdir(self.frequency_folder)
            os.mkdir(self.optimization_folder)
            os.mkdir(self.MCPB_folder)

    def createFiles(self):
        """
        Create input files for calculations

        Parameters:

        Class variables:
        """
        #####
        print('Creating necessary files...')
        f = open(self.optimization_folder + '/orca_opt.inp', 'w')
        f.write(self.orca_inp.format(self.method, '! OPT', self.solvent, self.cores, self.charge, self.multi))
        f.close()

        f = open(self.frequency_folder + '/orca_freq.inp', 'w')
        f.write(self.orca_inp.format(self.method, '! FREQ', self.solvent, self.cores, self.charge, self.multi))
        f.close()

        shutil.copy(self.path + '/' + self.inputFile, self.optimization_folder + '/' + self.inputFile)

    def Method(self, method: str = 'PBE0', basis: str = 'def2-SVP', DSP: str = 'D4', CPCM: str = 'Water',
               CPU: int = 12, epsilon: str = None, refrac: str = None):
        """
        Create folder structure for calculations

        Parameters:
            :param epsilon: permittivity for custom solvent
            :param refrac: refraction index for custom solvent
            :param int CPU: number of CPU cores, default is 12
            :param string CPCM: implicit solvent for orca calculations, default is Water
            :param string DSP: dispersion corrections, default is D4
            :param string basis:  ORCA basis set, default is def2-SVP
            :param string method: chosen QM method in ORCA input file format e.g. BP86 or MP2, default is PBE0

        Class variables:
            - self.method = ORCA method
            - self.solvent = solvent to be used for calculations, default is water
            - self.cores = number of CPU cores
        """
        self.method = '! {} {} {} '.format(method, basis, DSP)
        if CPCM == 'custom':
            self.solvent = '''
%CPCM       EPSILON      {}
            REFRAC       {}
END\n'''.format(epsilon, refrac)
        else:
            self.solvent = '! CPCM({})'.format(CPCM)
        self.cores = str(CPU)

    def run(self):
        """
        Run setup

        Parameters:

        Class variables:
        """
        status = 0
        self.createFolders()
        if self.check == 1:
            self.createFiles()

            print('Initial setup complete !')
            print('''

Set-up in folder {}
            
Your calculation will run with:
Method:     {}
Basis :     {}
Dispersion: {}
Solvent:    {}
CPU Cores:  {}
            '''.format(self.path, *self.method.replace('!', '').split(), self.solvent.replace('END\n', ''),
                       self.cores))
            status = 1
        return status
