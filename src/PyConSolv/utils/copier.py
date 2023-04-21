import shutil
class Copier:
    def __init__(self, source, destiny):
        self.source = source
        self.destiny = destiny

    def copy(self):
        shutil.copyfile(self.source,self.destiny)