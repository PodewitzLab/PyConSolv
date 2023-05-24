from ..utils.colorgen import Color
import os


class PATH:
    def __init__(self, program: str):
        self.status = 1
        self.program = program
        self.checkpath()
        return self.status

    def checkpath(self):
        """
        Check if program is available in PATH

        Parameters:

        Class variables:
        """
        self.status = 0
        self.programpath = ''
        systemPATH = os.environ.get('PATH', '').split(os.pathsep)
        for el in systemPATH:
            if self.program in el or self.program.upper() in el:
                print('Found {} in: '.format(self.program) + el)
                self.programpath = el + '/{} '.format(self.program)
                self.status = 1
                break
        if self.programpath == '':
            print(Color.RED + '{} was not found on your system... is it in your PATH?'.format(self.program) + Color.END)
            self.status = 0
            return
        else:
            return
