import os
import shutil
import subprocess
import sys

from .clustering import Cluster
from .solvent import Solvent
from ..interfaces.calculate import Calculation
from ..interfaces.cpptraj import CPPtraj


class Analysis:
    def __init__(self, path: str, alignMask: str = ''):
        """
        Analyse trajectory created from the PyConSolv simulation

        Parameters:
            - path = path that contains the simulation output file

        Class variables:
            - self.pyconsolv = check if the simulation was created using pyconsolv
            - self.orcapath = path to call ORCA executable
            - self.path = root path for calculations
            - self.homefolder = current working directory when calculations are started
            - self.simname = base name for the simulation
            - self.alignMask = atom mask to align the simulation to
            - self.orcafile = name of the orca input file for the single point calculations
            - self.reps = list of the cluster representatives
            - self.status = status of the calculation. 0 means an error occurred and everything should be stopped
            - self.aligndryFile = input file for cpptraj to align a dried solvent using the provided atommask
            - self.dryfile = input file for cpptraj to dry a simulation
            - self.alignfile = input file for cpptraj to aling a simulation to the provided atommask with solvent
            - self.solvent = solvent used for the simulation
        """
        self.pyconsolv = None
        self.solvent = None
        self.cpptraj = CPPtraj()
        self.rank = []
        self.path = os.path.abspath(path)
        self.homefolder = '/'.join(self.path.split('/')[:-1])
        self.simname = self.path.split('/')[-1].split('.')[0]
        self.alignMask = alignMask
        self.reps = []
        self.orcafile = 'orca_sp.inp'
        os.chdir(self.homefolder)

        self.alignfile ='''parm {}.prmtop
trajin {}
autoimage
align @{} first
trajout solv_aligned.nc
run
quit'''
        self.dryfile = '''parm {}.prmtop
trajin {}
strip :{}
autoimage
trajout dry.nc
run
quit'''
        self.aligndryFile = '''parm {}.prmtop
trajin {}
autoimage
align @{} first
trajout dry_aligned.nc
run
quit'''

    def alignSolv(self):
        """
        Align solvated trajectory

        Parameters:

        Class variables:
        """
        self.cpptraj.run('align')

    def align(self):
        """
        Align dry trajectory

        Parameters:

        Class variables:
        """
        self.cpptraj.run('align_dry')

    def dry(self):
        """
        Remove solvent from trajectory

        Parameters:

        Class variables:
        """
        self.cpptraj.run('dry_sim')

    def singlePoint(self, name: str):
        """
        Create folder and run single point calculation using orca. Folders are automatically created using the pdb file
        basename

        Parameters:
            :param string name: name of the pdb file for which a calculation shall be performed
        Class variables:
        """
        command = '{} {} > {}'.format(self.orcapath, self.orcafile, 'orca_sp.out')
        dirname = name.replace('.pdb', '')
        os.mkdir(dirname)
        shutil.copyfile(self.orcafile, dirname + '/' + self.orcafile)
        self.convertToXYZ(dirname)
        shutil.copyfile(dirname + '.xyz', dirname + '/input.xyz')
        os.chdir(dirname)
        calc = subprocess.run(command, shell=True)
        self.rank.append([dirname, float(self.getEnergy())])
        os.chdir(self.homefolder)

    def getReps(self):
        """
        Find all cluster representative pdb files generated by the clustering

        Parameters:

        Class variables:
        """
        files = os.listdir(self.homefolder)
        for file in files:
            if 'rep.c' in file and '.pdb' in file:
                self.reps.append(file)

    def convertToXYZ(self, basename):
        """
        Convert pdb files to XMOL format xyz files

        Parameters:
            :param string basename: basename for the pdb file to be converted
        Class variables:
        """
        data = []
        f = open(basename + '.pdb', 'r')
        for line in f:
            if 'ATOM' in line or 'ATM' in line:
                l = line.split()
                atomname = l[2]
                if len(atomname) > 1 and atomname[1].isdigit() == False:
                    ln = atomname[0] + atomname[1].lower()
                else:
                    ln = atomname[0]
                data.append(' '.join([ln, l[5], l[6], l[7]]))
        f.close()
        w = open(basename + '.xyz', 'w')
        w.write('{}\n'.format(str(len(data))))
        w.write('{} converted from pdb\n'.format(basename))
        for line in data:
            w.write(line + '\n')
        w.close()

    def getEnergy(self):
        """
        Get energy from orca single point calculation output

        Parameters:

        Class variables:
        """
        f = open('orca_sp.out', 'r')
        pattern = 'FINAL SINGLE POINT ENERGY'
        for line in f:
            if pattern in line:
                energy = line.split()[-1]
                break
        f.close()
        return energy

    def cluster(self, clustering):
        """
        Cluster the dry_aligned trajectory

        Parameters:
            :param string clustering: type of clustering to be performed, for options see ..misc.clustering.py
        Class variables:
        """
        cluster = Cluster(clustering)
        runtypes = cluster.runtypes
        if clustering in runtypes:
            print('clustering using {}\n'.format(clustering))
        else:
            print('clustering method not supported, please choose one of the following:\n')
            print(' '.join(runtypes) + '\n')
            return
        basename = 'cluster_{}'.format(clustering)
        inputstring = cluster.createInput()
        self.createInput(inputstring, basename)
        self.cpptraj.run(basename)
        print('Clustering done!\n')

    def createInput(self, inputstring, filename):
        """
        Create input file from string

        Parameters:
            :param string inputstring: input string to be written to file
            :param string filename: basename for the file to be written
        Class variables:
        """
        f = open('{}.in'.format(filename), 'w')
        f.write(inputstring)
        f.close()

    def writeFile(self,name, template, replacelist):
        f = open(name, 'w')
        for line in template.format(*replacelist):
            f.write(line)
        f.close()
    def useMask(self):
        self.writeFile('align.in', self.alignfile, ['LIG_solv', self.simname + '.nc', self.alignMask])
        self.writeFile('align_dry.in', self.aligndryFile, ['LIG_dry', 'dry.nc', self.alignMask])
        self.writeFile('dry_sim.in', self.dryfile, ['LIG_solv', self.simname + '.nc', self.solvent])



    def rankClusters(self):
        """
        Rank clusters according to energy

        Parameters:

        Class variables:
        """
        self.rank = sorted(self.rank, key=lambda x: x[1])
        f = open('cluster_ranking.dat', 'w')
        f.write('Cluster Energy\n')
        for el in self.rank:
            f.write('{} {} \n'.format(el[0],el[1]))
        f.close()

    def checkORCAFile(self) -> bool:
        """
        Check if orce file exists, if not, try to create one

        Parameters:

        Class variables:
        """
        found = False
        files = os.listdir(self.homefolder)
        for file in files:
            if file == 'orca_sp.inp':
                found = True

        if not found:
            try:
                tmp = []
                f = open(self.homefolder + '/../orca_calculations/opt/orca_opt.inp','r')
                for line in f:
                    if 'OPT' in line:
                        tmp.append(line.replace('OPT','SP'))
                        continue
                    tmp.append(line)
                f.close()
                f = open(self.homefolder + '/orca_sp.inp','w')
                for line in tmp:
                    f.write(line)
                f.close()
                found = True
            except:
                found = False
        return found

    def getSolvent(self):
        """
        Get the solvent used for the simulation

        Parameters:

        Class variables:
        """
        solv = Solvent()
        try:
            f = open('solvent','r')
            for line in f:
                s = line.replace('\n','')
            f.close()
            self.solvent = solv.solventDict[s]
        except:
            s = input('Unable to detect solvent used, please provide the solvent (e.g Water, CH2Cl2, custom):\n')
        self.solvent = solv.solventDict[s]

    def checkORCAPath(self):
        calc = Calculation(self.path)
        calc.checkpath()
        if calc.status == 0:
            print("ORCA was not found on your system.\n Aborting...\n")
            sys.exit()
        else:
            self.orcapath = calc.orcapath

    def checkPyConSolv(self) -> bool:
        if os.path.exists('solvent'):
            return True
        else:
            return False
    def run(self, clustering='kmeans', nosp = False):
        """
        Run clustering and ranking

        Parameters:
            :param string clustering: type of clustering to be performed, for options see ..misc.clustering.py

        Class variables:
        """
        self.pyconsolv = self.checkPyConSolv()
        if self.pyconsolv is False:
            print('This simulation was not created using PyConSolv, but can still be analyzed. Skipping to clustering...\n')
            self.useMask()
            self.align()
            self.cluster(clustering)
        else:
            self.getSolvent()
            self.useMask()
            self.alignSolv()
            self.dry()
            self.align()
            self.cluster(clustering)
        if not nosp:
            self.checkORCAPath()
            if not self.checkORCAFile():
                print('No calculation template file found, please make sure orca_sp.inp exists in this folder\n')
                return
            self.getReps()
            for rep in self.reps:
                print('Running single point for {}'.format(rep))
                self.singlePoint(rep)
            self.rankClusters()
            print(self.rank)
        else:
            print('nosp option detected, Clusters have not been ranked.\n')
