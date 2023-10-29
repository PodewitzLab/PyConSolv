class FileParser:
    def __init__(self, path:str, filelist:list):
        """
        Run Multiwfn charge calculations

        Parameters:
            :param str path: path where mol2 files are located
            :param list mol2list: list od mol2file names, without the .mol2 extension

        Class variables:
        """
        self.filelist = filelist
        self.path = path + '/'
        self.files = []

    def readfiles(self, ftype:str):
        for name in self.filelist:
            tmp = []
            f = open(self.path + ftype.format(name), 'r')
            for line in f:
                tmp.append(line)
            f.close()
            self.files.append(tmp)

    def write(self, outfile: str, file: str):
        f = open(self.path + outfile, 'w')
        for line in file:
            f.write(line)
        f.close()