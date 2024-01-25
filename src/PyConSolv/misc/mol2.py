from .files import FileParser


class mol2Parser(FileParser):

    def joinMol2(self):
        template = '''@<TRIPOS>MOLECULE
A
    {}     {}     1     0     0
SMALL
bcc


@<TRIPOS>ATOM
{}
@<TRIPOS>BOND
{}
@<TRIPOS>SUBSTRUCTURE
     1 A           1 TEMP              0 ****  ****    0 ROOT
'''
        atoms = []
        connect = []
        atomsperfile = [0]
        for file in self.files:
            flag = 0
            for line in file:
                if '@<TRIPOS>ATOM' in line:
                    flag = 1
                    continue
                if '@<TRIPOS>BOND' in line:
                    flag = 2
                    continue
                if '@<TRIPOS>SUBSTRUCTURE' in line:
                    break
                if flag == 1 :
                    atoms.append(line.split())
                if flag == 2 :
                    connect.append(line.split())
            atomsperfile.append(len(atoms))

        atomstring = ''
        connectionstring = ''

        for i in range(len(atoms)):
            atoms[i][0] = str(i+1)
            atomstring = atomstring + ' '.join(atoms[i]) + '\n'

        counter = 0
        for i in range(len(connect)):
            if i > 0:
                if (atomsperfile[counter] + int(connect[i][0]))< int(connect[i-1][0]) :
                    counter += 1
            connect[i][0] = str(i+1)
            connect[i][1] = str(int(connect[i][1]) + atomsperfile[counter])
            connect[i][2] = str(int(connect[i][2]) + atomsperfile[counter])
            connectionstring = connectionstring + ' '.join(connect[i]) + '\n'

        return template.format(len(atoms), len(connect), atomstring, connectionstring)


    def writeCombinedMol2(self):
        self.readfiles('{}.mol2')
        mol2string = self.joinMol2()
        self.write('LIG.mol2', mol2string)



