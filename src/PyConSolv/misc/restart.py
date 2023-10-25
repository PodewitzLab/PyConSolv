class RestartFile:
    def __init__(self, path: str):
        """
        Class for checking restart files and returning the last known state

        Parameters:
        :param string path: path to restart file (pyconsolv.restart)

        Class variables:
            - self.restart: path to restart file
            - self.state: calculation state, based on the restart file
        """
        self.restart = path + '/pyconsolv.restart'
        self.state = 0

    def getstate(self) -> int:
        f = open(self.restart, 'r')
        for line in f:
            match line.split()[0]:
                case 'setup':
                    self.state = 1
                case 'orca':
                    self.state = 2
                case 'antechamber':
                    self.state = 3
                case 'frcmod':
                    self.state = 4
                case 'multiwfn':
                    self.state = 5
                case 'mcpb':
                    self.state = 6
                case 'tleap':
                    self.state = 7
                case 'equilibration':
                    self.state = 8
                case 'DONE':
                    self.state = 9
                case other:
                    self.state = 0

        f.close()
        return self.state

    def write(self, state: str):
        f = open(self.restart, 'w')
        f.write(state)
        f.close()

    def parseInput(self):
        pass

    def writeInput(self):
        pass