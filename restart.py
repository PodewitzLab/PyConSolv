class RestartFile():
    def __init__(self,path):
        self.restart = path + '/pyconsolv.restart'
        self.state = 0
    def getstate(self):
        f = open(self.restart,'r')
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
                case other:
                    self.state = 0
            
        f.close()
        return self.state
    
    
    def write(self, state):
        f = open(self.restart,'w')
        f.write(state)
        f.close()
