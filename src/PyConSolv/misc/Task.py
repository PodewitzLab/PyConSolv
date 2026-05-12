import os
import sys
import shutil
from pathlib import Path
from tkinter import Tk

from ..ConfGen import PyConSolv
from .fragmenting import Fragmentor
from .ui import FragmentReviewGUI
from ..utils.colorgen import Color


class Task:
    def __init__(self, inputfilepath: str):
        '''
        Create parametrization tasks
        :param str inputfilepath: path to the input XYZ file
        '''
        self.inputfilepath = inputfilepath
        self.conf = PyConSolv(inputfilepath)

    def parametrize(self, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP',
                    dsp: str = 'D4', cpu: int = 12, solvent: str = 'Water',
                    multiplicity: int = 1, engine: str = 'amber', opt: bool = True,
                    box: int = 20, rst: bool = False, memory: int = 3000,
                    cart: str = None, cartstr: int = 100):
        '''
        Run standard parametrization task (non-fragment mode)
        '''
        self.conf.run(charge, method, basis, dsp, cpu, memory,
                      solvent, multiplicity, engine, opt, box,
                      rst, cart, cartstr)

    def fragment(self, charge: int = 0, method: str = 'PBE0', basis: str = 'def2-SVP',
                 dsp: str = 'D4', cpu: int = 12, solvent: str = 'Water',
                 multiplicity: int = 1, engine: str = 'amber', opt: bool = True,
                 box: int = 20, rst: bool = False, memory: int = 3000,
                 radius: float = 4.0, cart: str = None, cartstr: int = 100,
                 forcefield: str = 'amber', charge_method: str = 'resp'):
        '''
        Run fragment-based parametrization:
        1. Extract fragment around metal
        2. Parametrize fragment (ORCA + MCPB) to get metal parameters
        3. Parametrize full ligands with GAFF
        4. Combine parameters
        5. Build system, equilibrate, prepare simulation
        '''
        # Print header
        self.conf.startInfo()
        print(Color.CYAN + 'Substructure parametrization mode is on.\n' + Color.END)

        base_path = Path(self.inputfilepath).parent

        # === STEP 1: Extract fragment ===
        print(Color.GREEN + '\n=== Step 1: Extracting fragment ===' + Color.END)

        substructure_dir = base_path / 'substructure'
        try:
            substructure_dir.mkdir(exist_ok=True)
        except Exception as e:
            print(f'Could not create substructure folder: {e}')
            sys.exit()

        print(f'Using radius of {radius} Angstroms for substructure extraction')

        frag = Fragmentor(self.inputfilepath, radius=radius)
        if not frag.run(filename='substructure.xyz'):
            print('Fragmentation failed')
            return

        # Move substructure file to correct location
        substructure_file = base_path / 'substructure.xyz'
        target_file = substructure_dir / 'substructure.xyz'
        if substructure_file.exists():
            substructure_file.rename(target_file)

        substructure_inputpath = str(target_file)

        # === STEP 2: User reviews fragment (GUI) ===
        print(Color.GREEN + '\n=== Step 2: Fragment review ===' + Color.END)
        print(Color.CYAN + 'Opening fragment review window...' + Color.END)

        window = Tk()
        window.title('Fragment Review - PyConSolv')
        review_gui = FragmentReviewGUI(window, substructure_inputpath, radius)
        window.mainloop()
        window.destroy()

        if not review_gui.is_confirmed():
            print(Color.RED + 'Fragment review cancelled by user. Aborting.' + Color.END)
            return

        print(Color.GREEN + 'Fragment confirmed. Proceeding with parametrization...' + Color.END)

        # === STEP 3: Parametrize fragment (to get metal parameters) ===
        print(Color.GREEN + '\n=== Step 3: Parametrizing fragment (metal parameters) ===' + Color.END)

        # Build solvent string for ORCA
        if solvent and solvent.lower() not in ['none', 'water']:
            solvent_line = f'! CPCM({solvent})'
        elif solvent and solvent.lower() == 'water':
            solvent_line = '! CPCM(Water)'
        else:
            solvent_line = ''

        # Custom ORCA input for fragment - optimize only hydrogens (capping atoms).
        # CHARMM path needs a Hessian for Seminario/FFTK, so append FREQ.
        opt_line = '! OPT FREQ' if forcefield == 'charmm' else '! OPT'
        customOrcaInput = '''! {method} {basis} {dsp}
{opt_line}
{solvent}

%PAL NPROCS {cpu} END
%maxcore {memory}
%geom optimizehydrogens true
end

%scf
maxiter 350
end

* xyzfile {charge} {multiplicity} input.xyz
'''.format(
            method=method,
            basis=basis,
            dsp=dsp,
            opt_line=opt_line,
            solvent=solvent_line,
            cpu=cpu,
            memory=memory,
            charge=charge,
            multiplicity=multiplicity
        )

        # Create PyConSolv for fragment and run through MCPB
        fragment_conf = PyConSolv(substructure_inputpath)
        fragment_conf.startInfo()
        fragment_conf.checkRestart()

        fragment_conf.setup(charge, method, basis, dsp, solvent, cpu, multiplicity,
                           memory=memory, opt=True, customOrcaInput=customOrcaInput)

        if fragment_conf.restart < 2:
            print(Color.GREEN + 'Running ORCA optimization for fragment...' + Color.END)
            if fragment_conf.orca(opt=True) == 0:
                print(Color.RED + 'Fragment ORCA calculation failed!' + Color.END)
                return

        # CHARMM fragment path: the fragment ORCA run already produced opt+freq
        # because customOrcaInput requested FREQ. Hand the Hessian straight
        # to ConfGen.runCharmmFromFragment() and stop — no MCPB, no GAFF.
        if forcefield == 'charmm':
            frag_hess = os.path.join(fragment_conf.inputpath,
                                     'orca_calculations/freq/orca.hess')
            frag_opt_xyz = os.path.join(fragment_conf.inputpath,
                                        'orca_calculations/opt/orca_opt.xyz')
            if not os.path.isfile(frag_hess):
                print(Color.RED + 'Fragment Hessian missing; expected {} '
                      '(ORCA FREQ must have run on the fragment).'
                      .format(frag_hess) + Color.END)
                return
            print(Color.GREEN + '\n=== CHARMM full-structure parametrization '
                  '(fragment Hessian) ===' + Color.END)
            self.conf = PyConSolv(self.inputfilepath)
            self.conf.runCharmmFromFragment(
                fragment_hess_file=frag_hess,
                fragment_xyz_file=frag_opt_xyz,
                charge=charge, method=method, basis=basis, dsp=dsp,
                cpu=cpu, memory=memory, solvent=solvent,
                multiplicity=multiplicity, opt=opt, box=box,
                charge_method=charge_method)
            print(Color.GREEN + '\n' + '=' * 50 + Color.END)
            print(Color.GREEN + 'Fragment-based CHARMM parametrization complete!'
                  + Color.END)
            print(Color.GREEN + '=' * 50 + Color.END)
            return

        if fragment_conf.restart < 3:
            print(Color.GREEN + 'Running antechamber for fragment...' + Color.END)
            if fragment_conf.antechamber() == 0:
                print(Color.RED + 'Fragment antechamber failed!' + Color.END)
                return

        if fragment_conf.restart < 5:
            print(Color.GREEN + 'Running MultiWfn for fragment (RESP charges)...' + Color.END)
            if fragment_conf.multiwfn(cpu) == 0:
                print(Color.RED + 'Fragment MultiWfn failed!' + Color.END)
                return

        if fragment_conf.restart < 6:
            print(Color.GREEN + 'Running MCPB.py for fragment (metal parameters)...' + Color.END)
            if fragment_conf.MCPB_script() == 0:
                print(Color.RED + 'Fragment MCPB failed!' + Color.END)
                return

        # Save the fragment's MCPB output (contains metal parameters)
        fragment_mcpb_path = substructure_dir / 'MCPB_setup'

        # === STEP 4: Parametrize full structure with GAFF ===
        print(Color.GREEN + '\n=== Step 4: Parametrizing full structure (GAFF) ===' + Color.END)

        # Create a new PyConSolv for the full structure
        self.conf = PyConSolv(self.inputfilepath)
        self.conf.startInfo()
        self.conf.checkRestart()

        self.conf.setup(charge, method, basis, dsp, solvent, cpu, multiplicity,
                       memory=memory, opt=opt)

        if rst:
            self.conf.checkRT()

        # Run ORCA on full structure (or skip if user wants)
        if self.conf.restart < 2:
            print(Color.GREEN + 'Running ORCA for full structure...' + Color.END)
            if self.conf.orca(opt=opt) == 0:
                print(Color.RED + 'Full structure ORCA calculation failed!' + Color.END)
                return

        # Run antechamber to get GAFF parameters for ligands
        if self.conf.restart < 3:
            print(Color.GREEN + 'Running antechamber for full structure (GAFF)...' + Color.END)
            if self.conf.antechamber() == 0:
                print(Color.RED + 'Full structure antechamber failed!' + Color.END)
                return

        # Run MultiWfn for charges
        if self.conf.restart < 5:
            print(Color.GREEN + 'Running MultiWfn for full structure...' + Color.END)
            if self.conf.multiwfn(cpu) == 0:
                print(Color.RED + 'Full structure MultiWfn failed!' + Color.END)
                return

        # === STEP 5: Combine parameters ===
        print(Color.GREEN + '\n=== Step 5: Combining parameters ===' + Color.END)

        # Copy metal-specific frcmod from fragment to full structure's MCPB folder
        self._combine_parameters(fragment_mcpb_path, Path(self.conf.MCPB))

        # Run MCPB for full structure (uses combined parameters)
        if self.conf.restart < 6:
            print(Color.GREEN + 'Running MCPB.py for full structure...' + Color.END)
            if self.conf.MCPB_script() == 0:
                print(Color.RED + 'Full structure MCPB failed!' + Color.END)
                return

        # === STEP 6: Build system with tleap ===
        if self.conf.restart < 7:
            print(Color.GREEN + '\n=== Step 6: Building system (tleap) ===' + Color.END)
            if self.conf.tleap(solvent, box) == 0:
                print(Color.RED + 'tleap failed!' + Color.END)
                return

        # === STEP 7: Equilibration ===
        if self.conf.restart < 8:
            print(Color.GREEN + '\n=== Step 7: Equilibration ===' + Color.END)
            if self.conf.equilibration(cpu, engine, cart=cart, cartstr=cartstr) == 0:
                print(Color.RED + 'Equilibration failed!' + Color.END)
                return

        # === STEP 8: Prepare simulation ===
        if self.conf.restart < 9:
            print(Color.GREEN + '\n=== Step 8: Preparing simulation ===' + Color.END)
            if self.conf.prepareSimulation(solvent, engine, cart=cart, cartstr=cartstr) == 0:
                print(Color.RED + 'Simulation preparation failed!' + Color.END)
                return

        print(Color.GREEN + '\n' + '='*50 + Color.END)
        print(Color.GREEN + 'Fragment-based parametrization complete!' + Color.END)
        print(Color.GREEN + '='*50 + Color.END)

    def _combine_parameters(self, fragment_mcpb_path: Path, full_mcpb_path: Path):
        '''
        Combine parameters from fragment MCPB with full structure.
        Copies metal-specific frcmod files from fragment to full structure.

        :param fragment_mcpb_path: Path to fragment's MCPB_setup folder
        :param full_mcpb_path: Path to full structure's MCPB_setup folder
        '''
        print(f'Combining parameters from fragment: {fragment_mcpb_path}')
        print(f'Into full structure: {full_mcpb_path}')

        # Look for metal frcmod files in fragment
        # MCPB.py typically creates files like: LIG_mcpbpy.frcmod
        fragment_frcmod_files = list(fragment_mcpb_path.glob('*mcpbpy*.frcmod'))

        if not fragment_frcmod_files:
            # Try other patterns
            fragment_frcmod_files = list(fragment_mcpb_path.glob('*.frcmod'))

        for frcmod_file in fragment_frcmod_files:
            dest = full_mcpb_path / f'fragment_{frcmod_file.name}'
            print(f'  Copying {frcmod_file.name} -> {dest.name}')
            try:
                shutil.copy(frcmod_file, dest)
            except Exception as e:
                print(f'  Warning: Could not copy {frcmod_file}: {e}')

        # Also copy any metal mol2 files that might have special parameters
        metal_mol2_files = list(fragment_mcpb_path.glob('*[A-Z][A-Z]*.mol2'))  # Metal names are typically 2 uppercase letters
        for mol2_file in metal_mol2_files:
            # Check if it's likely a metal file (FE, CU, ZN, etc.)
            name = mol2_file.stem
            if len(name) <= 3 and name.isupper():
                dest = full_mcpb_path / mol2_file.name
                if not dest.exists():
                    print(f'  Copying metal mol2: {mol2_file.name}')
                    try:
                        shutil.copy(mol2_file, dest)
                    except Exception as e:
                        print(f'  Warning: Could not copy {mol2_file}: {e}')

        print('Parameter combination complete.')
