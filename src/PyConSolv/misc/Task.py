import os
import sys
from tkinter import Tk

from ..ConfGen import PyConSolv
from ..interfaces.calculate import Calculation
from .fragmenting import Fragmentor
from .ui import GUI
from ..utils.colorgen import Color


class Task:
    def __init__(self, inputfilepath):
        '''
        Create parametrization tasks
        '''
        self.inputfilepath = inputfilepath
        self.conf = PyConSolv(inputfilepath)
        self.conf.startInfo()

    def parametrize(self, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP', dsp: str = 'D4', cpu: int = 12,
            solvent: str = 'Water', multiplicity: int = 1, engine: str = 'amber', opt: bool = True, box: int = 20,
            rst: bool = False):
        '''
        Run parametrization task
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
        :param bool fragment : set if the structure should be fragmented
        :return:
        '''

        self.conf.run( charge, method, basis, dsp, cpu,
            solvent, multiplicity, engine, opt, box,
            rst)

    def fragment(self, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP', dsp: str = 'D4', cpu: int = 12,
            solvent: str = 'Water', multiplicity: int = 1, engine: str = 'amber', opt: bool = True, box: int = 20,
            rst: bool = False):
        '''
        Run parametrization task
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
        :param bool fragment : set if the structure should be fragmented
        :return:
        '''
        print(Color.CYAN + 'Substructure parametrization mode is on.\n' + Color.END)

        # Main structure
        # self.conf = PyConSolv(inputfilepath)
        self.conf.setup(charge, method, basis, dsp, solvent, cpu, multiplicity, opt=opt)
        self.conf.orca(opt=opt)
        self.setup(charge, method, basis, dsp, solvent, cpu, multiplicity, opt)

        if rst:
            self.conf.checkRT()

        if self.conf.restart < 2:
            if self.conf.orca(opt=opt) == 0:
                return

        if self.conf.restart < 3:
            if self.conf.antechamber() == 0:
                return
            if not self.conf.hasMetal:
                print('No metal has been found in your file, switching to simple parametrization')
                self.conf.multiwfn(cpu)
                self.conf.restart = 6  # skip mcpb and multiwfn if no metal

        if self.conf.restart < 5:
            if self.conf.multiwfn(cpu) == 0:
                return

        ##Substructure parametrization
        try:
            os.mkdir('/'.join(self.inputfilepath.split('/')[:-1] + ['substructure']))
            os.mkdir('/'.join(self.inputfilepath.split('/')[:-1] + ['substructure/prepare']))
        except:
            print('Could not create substructure folders')
            sys.exit()
        radius = float(
            input('Please input the radius around the metal to be included in the substructure (default 3):') or 3.0)
        frag = Fragmentor(self.inputfilepath, radius=radius)
        frag.run(filename='/substructure/substructure.xyz')
        inputfilepath = '/'.join(self.inputfilepath.split('/')[:-1] + ['/substructure/substructure.xyz'])

        hydrogens = []
        with open(inputfilepath) as f:
            next(f)
            next(f)
            for i, line in enumerate(f):
                if len(line.split())<4:
                    continue
                if line.split()[0] == 'H':
                    hydrogens.append(i)
        print(hydrogens)
        '''
        TODO:
        Optimize structure
        display in 3d to the user
        allow for editing
        '''

        # window = Tk()
        # window.title('Fragment display')
        # start = GUI(window, '/'.join(inputfilepath.split('/')[:-1]), pdbs)
        # window.mainloop()
        # window.destroy()

        customOrcaInput = '''{}
{}
{}

%PAL NPROCS {} END
%geom optimizehydrogens true
end

%scf
maxiter 350
end

* xyzfile {} {} input.xyz
'''
        self.conf = PyConSolv(inputfilepath)
        self.conf.setup( charge, method, basis, dsp, solvent, cpu, multiplicity, opt=True, customOrcaInput=customOrcaInput)
        self.conf.orca(opt = True)
        if self.conf.antechamber() == 0:
            return
        if not self.conf.hasMetal:
            print('No metal has been found in your file, switching to simple parametrization')
            # self.conf.multiwfn(cpu)

        if self.conf.multiwfn(cpu) == 0:
            return

        if self.conf.MCPB_script() == 0:
            return

