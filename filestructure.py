import shutil
import os
from .colorgen import color


class Setup():
    def __init__(self, path, charge = 0, multi = 1):
        '''
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
        '''
        self.db_file = os.path.split(__file__)[0] + '/db/atom-radius.txt'
        #self.db_file =  r'/home/rat/PyPer/atom-radius.txt'#change to actual file
        self.path = '/'.join(path.split('/')[:-1])
        self.inputFile =  path.split('/')[-1]
        self.charge = charge
        self.multi = multi
        self.check = 1
        #### ORCA input file template
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
        '''
        Create folder structure for calculations
        
        Parameters:
            
        Class variables:
            - self.calculation_folder = folder where all orca calculations will be performed
            - self.frequency_folder = folder where orca frequency calculations will be performed
            - self.optimization_folder = folder where orca geometry optimization calculations will be performed
            - self.MCPB_folder = folder where MCPB setup calculations will be performed
            - self.equilibration_folder = folder where system equilibration will be performed
            - self.simulation_folder = folder where a CMD will be performed
        '''
        
        
        self.calculation_folder = self.path+'/'+ 'orca_calculations'
        self.equilibration_folder = self.path+'/'+ 'equilibration'
        self.simulation_folder = self.path+'/'+ 'simulation'
        self.frequency_folder = self.calculation_folder + '/' + 'freq'
        self.optimization_folder = self.calculation_folder + '/' + 'opt'
        self.MCPB_folder = self.path+'/'+ 'MCPB_setup'
        
        if len(os.listdir(self.path)) > 1:
            self.check = 0
            print(color.RED + 'Your calculation directory includes other folders/files. Please make sure it only contains your input xyz file!\n' + color.END)
        else:
            print('Creating folders...')
            os.mkdir(self.calculation_folder)
            os.mkdir(self.equilibration_folder)
            os.mkdir(self.simulation_folder)
            os.mkdir(self.frequency_folder)
            os.mkdir(self.optimization_folder)
            os.mkdir(self.MCPB_folder)
        
    def createFiles(self):
        '''
        Create input files for calculations
        
        Parameters:
            
        Class variables:
        '''
        #####
        print('Creating necessary files...')
        f = open(self.optimization_folder + '/orca_opt.inp','w')
        f.write(self.orca_inp.format(self.method, '! OPT', self.solvent, self.cores, self.charge, self.multi))
        f.close()
        
        f = open(self.frequency_folder + '/orca_freq.inp','w')
        f.write(self.orca_inp.format(self.method, '! FREQ', self.solvent, self.cores, self.charge, self.multi))
        f.close()
        
        shutil.copy(self.path + '/' + self.inputFile, self.optimization_folder + '/'  + self.inputFile)
        shutil.copy(self.path + '/' + self.inputFile, self.frequency_folder + '/'  + self.inputFile) # this should be changed to wait until optimization is run
        
    def Method(self, method='PBE0', basis = 'def2-SVP', DSP = 'D4', CPCM = 'Water', CPU = 8):
        '''
        Create folder structure for calculations
        
        Parameters:
            - method = chosen QM method in ORCA input file format e.g. BP86 or MP2, default is PBE0
            - basis = ORCA basis set, default is def2-SVP
            - DSP = dispersion corrections, default is D4
            - CPU = number of CPU cores, default is 8
        Class variables:
            - self.method = ORCA method
            - self.solvent = solvent to be used for calculations, default is water
            - self.cores = number of CPU cores
        '''
        self.method = '! {} {} {} '.format(method, basis, DSP)
        self.solvent = '! CPCM({})'.format(CPCM)
        self.cores = str(CPU)
    
    def run(self):
        '''
        Run setup
        
        Parameters:
            
        Class variables:
        '''
        status = 0
        self.Method()
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
            '''.format(self.path,*self.method.replace('!','').split(),self.solvent.replace('! CPCM',''), self.cores))
            status = 1
        return status
