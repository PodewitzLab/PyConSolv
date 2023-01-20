import os
import subprocess

from .colorgen import Color


class amberInterface:
    def __init__(self, path):
        """
        Runs main PDB generation functionality of the XYZ class.

        Parameters:
            - inputfile = full path to xyz inputfile

        Class variables:
        """

        self.path = path
        self.status = 0
        self.original_wd = os.getcwd()
        os.chdir(self.path)

    def checkpath(self):
        """
        Check if Amber is available in PATH

        Parameters:

        Class variables:
        """
        self.status = 0
        self.amberhome = os.getenv('AMBERHOME')
        if self.amberhome == '':
            print(Color.RED + 'Amber is not in your PATH!!!' + Color.END)
            self.status = 0
        else:
            self.status = 1

    def inputFileGenerator(self, metals, ligands):
        """
        Creates an input file for MCPB.py.

        Parameters:
            - metals = basename of the metal containing files - string
            - ligands = basename of the ligand containing files - list of strings

        Class variables:
            - self.inputfile = MCPB.py input file
        """
        # TODO: change input file to handle multiple iod IDs
        self.inputfile = '''original_pdb Full_PDB.pdb
group_name LIG
cut_off 2.8
ion_ids 1
software_version g16
ion_mol2files {}.mol2
naa_mol2files {}.mol2
frcmod_files {}.frcmod\n'''.format(metals, '.mol2 '.join(ligands), '.frcmod '.join(ligands))

        f = open(self.path + '/input.in', 'w')
        f.write(self.inputfile)
        f.close()

    def antechamber(self, name, charge):
        """
        Run Amber antechamber

        Parameters:
            - name = basename of the structure
            - charge = charge of the structure

        Class variables:
        """
        print('Running antechamber for {}.pdb'.format(name))
        antechamberCommand = 'antechamber -fi pdb -fo mol2 -i {}.pdb -o {}.mol2 -c bcc -pf y -nc {} > antechamber_{}.out'.format(
            name, name, int(float(charge)), name)
        calc = subprocess.run([antechamberCommand], shell=True)
        if calc.returncode == 0:
            print('Antechamber completed successfully for {}.pdb\n'.format(name))
            self.status = 1
        else:
            self.status = 0

        return self.status

    def runMCPB(self, step):
        """
        Run MCPB.py

        Parameters:
            - step = string, representing the step of MCPB.py to run

        Class variables:
        """
        os.chdir(self.path)
        MCPBCommand = 'MCPB.py -i input.in -s {} > MCPB_s{}.out'.format(step, step)
        calc = subprocess.run([MCPBCommand], shell=True)
        if calc.returncode == 0:
            print('MCPB.py step {} completed successfully\n'.format(step))
            self.status = 1
        else:
            self.status = 0

    def equilibrate(self, cpus=8):

        os.chdir(self.path + '/../equilibration')
        equilibrateCommand = [
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 1
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 2
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 3
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 4
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 5
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 6
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 7
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 8
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 9
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 10
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 11
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 12
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 13
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 14
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 15
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 16
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 17
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 18
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 19
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7',
            # 20
            'pmemd.cuda -O -i {:02d}.in -o {:02d}.out -p LIG_solv.prmtop -c {:02d}.rst7 -r {:02d}.rst7 -ref {:02d}.rst7'
            # 21
            ]

        for i in range(1, 22):
            if 'pmemd.cuda' in equilibrateCommand[i]:
                calc = subprocess.run(equilibrateCommand[i].format(i, i, i - 1, i, i - 1), shell=True)
            else:
                print(equilibrateCommand[i].format(cpus, i, i, i - 1, i, i - 1))
                calc = subprocess.run(equilibrateCommand[i].format(cpus, i, i, i - 1, i, i - 1), shell=True)
            if calc.returncode == 0:
                print('Equlibration step {} completed successfully\n'.format(i))
                self.status = 1
            else:
                self.status = 0
                print('Equlibration step {} failed!\n'.format(i))
                return self.status

    def runTleap(self):
        """
        Run tleap

        Parameters:

        Class variables:
        """

        tleapCommand = 'tleap -f LIG_tleap.in > tleap.log'
        calc = subprocess.run([tleapCommand], shell=True)
        if calc.returncode == 0:
            print('Tleap completed successfully\n')
            self.status = 1
        else:
            print('Warning, something went wrong while running Tleap!')
            self.status = 0
        return self.status

    def runParmchk2(self, name):
        """
        Run parmchk2

        Parameters:
            - name = basename of the structure

        Class variables:
        """
        parmchk2Command = 'parmchk2 -i {}.mol2 -o {}.frcmod -f mol2'.format(name[1], name[1])
        calc = subprocess.run([parmchk2Command], shell=True)
        if calc.returncode == 0:
            print('Generated {}.frcmod...\n'.format(name[1]))
            self.status = 1
        else:
            self.status = 0
            return self.status

    def equil(self, path, restrained_residue='1'):
        """
        Create input files for a stepwise equilibration

        Parameters:
            - path = full path to the location of the solvated structure
            - restrained_residue = residue to be restrained in steps 1-19

        Class variables:
        """
        minseq = [1000, 500, 200, 100, 50, 20, 10, 5, 4, 3, 2, 1, 0.5]

        # default amber settings
        minheader = "minimization\n&cntrl\n\timin=1,ncyc=500,maxcyc=1000,ntr=1,cut=8.0,\n\t"
        nvtheader = "equilibration\n&cntrl\n\tntb=1,ntc=2,ntf=2,ntt=3,gamma_ln=2.0,cut=8.0,\n\t"
        nptheader = "pressure equilibration\n&cntrl\n\tntb=2,ntp=1,pres0=1.0,tautp=2.0,\n\tntc=2,ntf=2,ntt=3,gamma_ln=2.0,\n\ttempi=300.0,temp0=300.0,\n\t"
        nvtheader2 = "equilibration\n&cntrl\n\tntb=1,ntc=2,ntf=2,ntt=3,ntr=1,gamma_ln=2.0,cut=8.0\n\t"

        restrained_residue = restrained_residue
        data = [[minheader, 'restraint_wt=1000.0,restraintmask="!@H="\n/\n'],  # 1
                [minheader, 'restraint_wt=1000.0,restraintmask="!@H=&:{}",\n/\n'.format(restrained_residue)],  # 2
                [nvtheader,
                 'ntr=1,restraint_wt=1000.0,restraintmask="!@H=&:{}",\n\tnstlim=100000,nmropt=1,dt=0.001,\n/\n'.format(
                     restrained_residue) + '&wt TYPE=\'TEMP0\', istep1=0, istep2=100000, value1=100.0,value2=300.0 /\n&wt TYPE=\'END\' /\n'],
                # 3
                [nptheader, 'ntr=1,restraint_wt=1000.0,restraintmask="!@H=&:{}",\n\tnstlim=30000,dt=0.001,\n/\n'.format(
                    restrained_residue)],  # 4s
                [nptheader, 'ntr=1,restraint_wt=1000.0,restraintmask="!@H=&:{}",\n\tnstlim=70000,dt=0.001,\n/\n'.format(
                    restrained_residue)],  # 4
                [nvtheader,
                 'ntr=1,restraint_wt=1000.0,restraintmask="!@H=&:{}",\n\tnstlim=100000,dt=0.001,nmropt=1,\n/\n'.format(
                     restrained_residue) + '&wt TYPE=\'TEMP0\', istep1=0, istep2=100000, value1=300.0,value2=100.0 /\n&wt TYPE=\'END\' /\n'],
                # 5
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[0], restrained_residue)],  # 6
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[1], restrained_residue)],  # 7
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[2], restrained_residue)],  # 8
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[3], restrained_residue)],  # 9
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[4], restrained_residue)],  # 10
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[5], restrained_residue)],  # 11
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[6], restrained_residue)],  # 12
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[7], restrained_residue)],  # 13
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[8], restrained_residue)],  # 14
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[9], restrained_residue)],  # 15
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[10], restrained_residue)],  # 16
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[11], restrained_residue)],  # 17
                [minheader,
                 'restraint_wt={:0.1f},restraintmask="!@H=&:{}",\n/\n'.format(minseq[12], restrained_residue)],  # 18
                [nvtheader.replace(',ntr=1', ''),
                 'nstlim=200000,dt=0.001,nmropt=1,\n/\n' + '&wt TYPE=\'TEMP0\', istep1=0, istep2=200000, value1=100.0,value2=300.0 /\n&wt TYPE=\'END\' /\n'],
                # 19
                [nptheader, 'nstlim=50000,dt=0.001,\n/\n']]  # 20

        # write input files
        for i in range(21):
            f = open(path + '/equilibration/{:>02}.in'.format(i + 1), 'w')
            f.write(data[i][0])
            f.write(data[i][1])
            f.close()

    def readMetalBonds(self, path):
        """
        Read atom-metal bond file

        Parameters:
            - path = full path to the folder containing the metalConnections file

        Returns:
            - metalbonds = [a] list of strings formatted as 'metalID @atomElementatomID atomID'
        Class variables:
        """

        metalbonds = []
        f = open(path + '/metalConnections', 'r')
        for line in f:
            metalbonds.append('{}-{} {}'.format(line.split()[0], line.split()[2], line.split()[1]))
        f.close()
        return metalbonds

    def readConnections(self, path):
        """
        Read atom-atom bond file

        Parameters:
            - path = full path to the folder containing the Connections file

        Returns:
            - connections = [a] list of atom ids, in the new pdb order
        Class variables:
        """

        connections = []
        f = open(path + '/Connections', 'r')
        for line in f:
            connections = connections + line.split()
        f.close()
        return connections

    def checkMCPBBonds(self, path):
        """
        Check whether all the bonds have been detected by MCPB.py. If not, the input.in file is modified to contain the additional bonds that are misssing.
        This overcomes the lack of support for metal-C bonds by default in MCPB.py

        Parameters:
            - path = full path to the folder containing the Connections file

        Returns:
            - Bool: False when connections are missing; True when all connections are present
        Class variables:
        """

        flag = 0
        residues = []
        f = open(path + '/MCPB_s1.out', 'r')
        for line in f:
            if 'Metal Site Information' in line:
                flag = 1
                continue
            if flag == 1:
                if 'The following residues are in the Metal Site' in line:
                    break
                if '@' in line:
                    residues.append(line.split()[0])
        f.close()

        metalbonds = self.readMetalBonds(path)
        connections = self.readConnections(path)

        inputAppend = []

        for bond in metalbonds:
            flag = 0
            for el in residues:
                if bond.split()[-1] in el:
                    flag = 1
            if flag == 0:
                a = int(bond.split()[0].split('-')[0]) + 1  # MCPB.py counts from 1
                tmp = bond.split()[0].split('-')[1]
                b = int(connections.index(tmp)) + 1  # MCPB.py counts from 1
                inputAppend.append('add_bonded_pairs {}-{}\n'.format(a, b))

        if inputAppend == []:
            return True
        else:
            print('Not all connections were autodetected by MCPB.py\n')
            f = open(path + '/input.in', 'a')
            for el in inputAppend:
                print('Adding : {}\n'.format(el.split()[-1]))
                f.write(el)
            f.close()
            return False

    def tleapChecker(self, path):
        """
        Remove duplicate lines that can sometimes appear within tleap

        Parameters:
            - path = full path to the folder containing the LIG_tleap.in file

        Class variables:
        """

        tmp = []
        f = open(path + '/LIG_tleap.in', 'r')
        for line in f:
            if line in tmp:
                continue
            else:
                tmp.append(line)
        f.close()

        f = open(path + '/LIG_tleap.in', 'w')
        for line in tmp:
            f.write(line)
        f.close()
        # todo add options for changing tleap input

    def tleapNoMetal(self, path, name):
        file = '''source oldff / leaprc.ff99SB
        source leaprc.gaff
        LIG = loadmol2 LIG.mol2
        loadamberparams LIG.frcmod
        check LIG
        saveoff LIG LIG.lib
        saveamberparm LIG LIG.prmtop LIP.rst7
        quit
        '''.replace('LIG', name)
        f = open(path + '/LIG_tleap.in', 'w')
        f.write(file)
        f.close()

    def tleapNoMetalSolv(self, path, name):
        file = '''source oldff / leaprc.ff99SB
        source leaprc.gaff
        source leaprc.water.tip3p
        LIG = loadmol2 LIG.mol2 
        loadamberparams LIG.frcmod
        savepdb LIG LIG_dry.pdb
        saveamberparm LIG LIG_dry.prmtop LIG_dry.inpcrd
        solvatebox LIG TIP3PBOX 20 iso
        savepdb LIG LIG_solv.pdb
        saveamberparm LIG LIG_solv.prmtop LIG_solv.inpcrd
        quit
        '''.replace('LIG', name)
        f = open(path + '/LIG_tleap.in', 'w')
        f.write(file)
        f.close()
