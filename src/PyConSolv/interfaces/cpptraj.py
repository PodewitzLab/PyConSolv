from ..misc.pathchecker import PATH
import subprocess

class CPPtraj:
    def __init__(self):
        self.program = 'cpptraj'
        # self.status = PATH(self.program)

    def run(self, basename):
        command = 'cpptraj {}.in > {}.log'.format(basename,basename)
        calc = subprocess.run(command, shell=True)