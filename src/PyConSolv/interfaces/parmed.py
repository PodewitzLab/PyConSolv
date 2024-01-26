import os
import parmed as pmd
import shutil

class Parmed():

    def __init__(self, path: str):
        '''
        Class to interface to parmed for various functions

        :param path: location of files which should be processed
        '''
        self.path = path
        self.original_dir = os.getcwd()

    def checkTop(self, name: str ='LIG_solv'):
        '''
        Function to check topology and input coordinates after generating the files with tleap. It is known that
        sometimes tleap does not generate fully correct geometries and topologies when multiple molecules are involved.
        :param string name: name of the files that should be checked, generally LIG_solv and LIG_dry
        :return:
        '''
        os.chdir(self.path)
        a = pmd.load_file('{}.prmtop'.format(name), xyz='{}.inpcrd'.format(name))
        a.rediscover_molecules(True)
        a.remake_parm()
        a.load_pointers()
        a.save('{}_checked.prmtop'.format(name))
        a.save('{}_checked.inpcrd'.format(name))
        a.save('{}_checked.pdb'.format(name))
        shutil.move('{}_checked.prmtop'.format(name),'{}.prmtop'.format(name))
        shutil.move('{}_checked.inpcrd'.format(name), '{}.inpcrd'.format(name))
        shutil.move('{}_checked.pdb'.format(name), '{}.pdb'.format(name))