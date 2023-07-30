from .amber import amberInterface
from .gromacs import gromacsInterface


class MDEngine:
    def __init__(self, path: str, engine: str = 'amber'):
        self.status = 1
        match engine:
            case 'amber':
                self.MD = amberInterface(path)
            case 'gromacs':
                self.MD = gromacsInterface(path)

    def run(self, path, cpus):
        self.MD.prepare(path)
        self.status = self.MD.equilibrate(cpus)
        return self.status

