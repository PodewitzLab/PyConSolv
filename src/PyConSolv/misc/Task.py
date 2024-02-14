import os
import sys
from tkinter import Tk

from ..ConfGen import PyConSolv
from ..interfaces.calculate import Calculation
from .fragmenting import Fragmentor
from .ui import GUI
from ..utils.colorgen import Color


class Task:
    def __init__(self):
        '''
        Create parametrization tasks
        '''
        pass

    def parametrize(self, inputfilepath, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP', dsp: str = 'D4', cpu: int = 12,
            solvent: str = 'Water', multiplicity: int = 1, engine: str = 'amber', opt: bool = True, box: int = 20,
            rst: bool = False, fragment: bool = False):
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
        conf = PyConSolv(inputfilepath)
        conf.run( charge, method, basis, dsp, cpu,
            solvent, multiplicity, engine, opt, box,
            rst, fragment)

    def fragment(self, inputfilepath, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP', dsp: str = 'D4', cpu: int = 12,
            solvent: str = 'Water', multiplicity: int = 1, engine: str = 'amber', opt: bool = True, box: int = 20,
            rst: bool = False, fragment: bool = False):
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
        if fragment:
            print(Color.CYAN + 'Substructure parametrization mode is on.\n' + Color.END)

        try:
            os.mkdir('/'.join(inputfilepath.split('/')[:-1] + ['substructure']))
            os.mkdir('/'.join(inputfilepath.split('/')[:-1] + ['substructure/prepare']))
        except:
            print('Could not create substructure folders')
            sys.exit()
        radius = float(
            input('Please input the radius around the metal to be included in the substructure (default 3):') or 3.0)
        frag = Fragmentor(inputfilepath, radius=radius)
        frag.run(filename='/substructure/substructure.xyz')
        inputfilepath = '/'.join(inputfilepath.split('/')[:-1] + ['/substructure/substructure.xyz'])

        hydrogens = []
        with open(inputfilepath) as f:
            next(f)
            next(f)
            for i, line in enumerate(f):
                if line.split() == 'H':
                    hydrogens.append(i)
        print(hydrogens)
        optmize = Calculation('substructure/prepare')
        # inputfilepath = '/'.join(inputfilepath.split('/')[:-1] + ['/substructure/substructure.xyz'])
        '''
        TODO:
        Optimize structure
        display in 3d to the user
        allow for editing
        
        '''

        window = Tk()
        window.title('Fragment display')
        start = GUI(window, '/'.join(inputfilepath.split('/')[:-1]), pdbs)
        window.mainloop()
        window.destroy()