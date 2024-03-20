import os
import shutil
import subprocess

from ..utils.colorgen import Color
from ..misc.converter import Converter


class gromacsInterface:
    def __init__(self, path: str):
        """
        Interfaces to the GROMACS program package to set up the simulation.

        Parameters:
            :param string path: full path to the folder where parametrization will take place

        Class variables:
            - self.path - location of input files
            - self.executable `- gromacs executable name
            - self.original_wd - initial working directory where class was created from
            - self.status - status of program (0 for error, 1 for success)
            - self.atoms - number of atoms in solute
            - self.system_size - number of atoms in system
        """

        self.executable = None
        self.path = path
        self.status = 1
        self.original_wd = os.getcwd()
        os.chdir(self.path)
        self.atoms = None
        self.system_size = None

        self.energy_minimize='''; em.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
'''

        self.nvt_equilibration = '''title                   = PyConSolv Equilibration
define                  = -DPOSRES  ; position restrain the solute
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale ; modified Berendsen thermostat
tc-grps                 = Solute SOL   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1 ; time constant, in ps
ref_t                   = 300     300 ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 300       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
'''

        self.npt_equilibration = '''title                   = PyConSolv NPT equilibration
define                  = -DPOSRES  ; position restrain the solute
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = yes       ; Restarting after NVT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Solute SOL   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = no        ; Velocity generation is off
'''

    def checkpath(self):
        """
        Check if GROMACS is available in PATH

        Parameters:

        Class variables:
        """
        path = ''
        self.status = 0
        path = shutil.which('gmx')
        if path is None:
            path = shutil.which('gmx_mpi')
        if path is None:
            print(Color.RED + 'GROMACS is not in your PATH!!!' + Color.END)
            self.status = 0
        else:
            self.executable = path.split('/')[-1]
            self.status = 1

    def equil(self, path: str, restrained_residue: str='1'):
        """
        Create input files for a stepwise equilibration

        Parameters:
            :param string path: full path to the location of the solvated structure
            :param string restrained_residue: residue to be restrained

        Class variables:
        """
        f = open('em.mdp', 'w')
        f.write(self.energy_minimize)
        f.close()
        f = open('nvt.mdp', 'w')
        f.write(self.nvt_equilibration)
        f.close()
        f = open('npt.mdp', 'w')
        f.write(self.npt_equilibration)
        f.close()

    def createRestraints(self, atoms):
        """
        Create restraint files for a stepwise equilibration

        Parameters:
            :param string atoms: ids of atoms to be restrained

        Class variables:
        """
        f = open('posre.itp', 'w')
        f.write('''; position restraints for system of GROningen MAchine for Chemical Simulation

[ position_restraints ]
;  i funct       fcx        fcy        fcz\n''')
        for i in range(atoms):
            f.write('{:>4}    1       1000       1000       1000\n'.format(i+1))
        f.close()

    def createIndex(self, atoms, system_size):

        f = open('index.ndx', 'w')

        template = '{:>4} '

        # create system
        f.write('[ System ]\n')
        line = ''
        for i in range(system_size):
            if (i + 1) % 15 == 0:
                f.write(line + template.format(i+1) + '\n')
                line = ''
            else:
                line = line + template.format(i + 1)
        if line != '':
            f.write(line + '\n')

        # create Other
        f.write('[ Other ]\n')
        line = ''
        for i in range(atoms):
            if (i + 1) % 15 == 0:
                f.write(line + template.format(i + 1) + '\n')
                line = ''
            else:
                line = line + template.format(i + 1)
        if line != '':
            f.write(line + '\n')

        #create solute
        f.write('[ Solute ]\n')
        line = ''
        for i in range(atoms):
            if (i+1) % 15 == 0:
                f.write(line + template.format(i+1) + '\n')
                line = ''
            else:
                line = line + template.format(i+1)
        if line != '':
            f.write(line + '\n')

        #create SOL
        f.write('[ SOL ]\n')
        line = ''
        for i in range(atoms, system_size):
            if (i + 1 - atoms) % 15 == 0:
                f.write(line + template.format(i + 1) + '\n')
                line = ''
            else:
                line = line + template.format(i + 1)
        if line != '':
            f.write(line + '\n')
        f.close()
    def runCommand(self, command, message):
        calc = subprocess.run(command, shell=True)
        if calc.returncode == 0:
            print('{}\n'.format(message))
            self.status = 1
        else:
            self.status = 0

    def equilibrate(self, cpus: int =8):
        os.chdir(self.path + '/../equilibration')
        command = '{} grompp -f em.mdp -c LIG_solv.gro -p LIG_solv.top -o em.tpr > em_prep.log'.format(self.executable)
        self.runCommand(command, 'Generated em.tpr')
        print('Minimizing\n')
        command = '{} mdrun -v -deffnm em > em_rum.log'.format(self.executable)
        self.runCommand(command, 'Energy minimization done!')

        command = '{} grompp -f nvt.mdp -c em.gro -r em.gro -p LIG_solv.top -o nvt.tpr -n index.ndx > nvt_prep.log'.format(self.executable)
        self.runCommand(command, 'Generated nvt.tpr')
        print('Equilibrating in NVT\n')
        command = '{} mdrun -v -deffnm nvt > nvt_run.log'.format(self.executable)
        self.runCommand(command, 'NVT equilibration done!')

        command = '{} grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p LIG_solv.top -o npt.tpr -n index.ndx> npt_prep.log'.format(self.executable)
        self.runCommand(command, 'Generated npt.tpr')
        print('Equilibrating in NPT\n')
        command = '{} mdrun -v -deffnm npt > npt_run.log'.format(self.executable)
        self.runCommand(command, 'Equilibration done!')

        return self.status

    def addModifyTop(self, file):
        tmp = []
        f = open(file, 'r')
        start = 0
        for line in f:

            if '[ moleculetype ]' in line:
                start += 1
            if start == 2:
                tmp.append('''; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif
''')
                start = 3
            tmp.append(line)
        f.close()

        f = open(file,'w')
        for el in tmp:
            f.write(el)
        f.close()
    def prepare(self, path = None, restrain: str = None, cart: str = None, cartstr: int = 100): #TODO implement restraints for gromacs
        os.chdir(self.path)
        self.checkpath()
        self.converter = Converter(self.path, intype='amber', outtype='gromacs')
        self.converter.getAtoms()
        self.converter.getSystemSize()
        self.atoms = self.converter.atoms
        self.system_size = self.converter.system_size
        self.converter.convert()
        self.addModifyTop('LIG_solv.top')
        shutil.copyfile(self.path+'/LIG_solv.top', path + '/equilibration/LIG_solv.top')
        shutil.copyfile(self.path + '/LIG_solv.gro', path + '/equilibration/LIG_solv.gro')
        os.chdir(path + '/equilibration')
        self.createIndex(self.atoms,self.system_size)
        self.createRestraints(self.atoms)
        self.equil(self.path)
        os.chdir(self.path)
