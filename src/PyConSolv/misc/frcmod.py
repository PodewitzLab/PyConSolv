from .files import FileParser


class frcmodParser(FileParser):

    def joinFrcmod(self):
        template = '''Remark line goes here
MASS
{}

BOND
{}

ANGLE
{}

DIHE
{}

IMPROPER
{}

NONBON
{}
'''
        mass = []
        bond = []
        angle = []
        dihedral = []
        improper = []
        nonbonded = []
        for file in self.files:
            flag = 0
            for line in file:
                if 'MASS' in line:
                    flag = 1
                    continue
                if 'BOND' in line:
                    flag = 2
                    continue
                if 'ANGLE' in line:
                    flag = 3
                    continue
                if 'DIHE' in line:
                    flag = 4
                    continue
                if 'IMPROPER' in line:
                    flag = 5
                    continue
                if 'NONBON' in line:
                    flag = 6
                    continue
                if flag == 1 :
                    if len(line.split('  '))>1:
                        mass.append(line.split('  '))
                if flag == 2:
                    if len(line.split('  ')) > 1:
                        bond.append(line.split('  '))
                if flag == 3:
                    if len(line.split('  ')) > 1:
                        angle.append(line.split('  '))
                if flag == 4 :
                    if len(line.split('  '))>1:
                        dihedral.append(line.split('  '))
                if flag == 5 :
                    if len(line.split('  '))>1:
                        improper.append(line.split('  '))
                if flag == 6 :
                    if len(line.split('  '))>1:
                        nonbonded.append(line.split('  '))

        masstring = ''
        bondstring = ''
        anglestring = ''
        dihedralstring = ''
        improperstring = ''
        nonbondedstring = ''
        if len(mass) > 0:
            for i in range(len(mass)):
                if mass[i][0] not in masstring:
                    masstring = masstring + '  '.join(mass[i])
        else:
            masstring = '\n'
        if len(bond) > 0:
            for i in range(len(bond)):
                if bond[i][0] not in bondstring:
                    bondstring = bondstring + '  '.join(bond[i])
        else:
            bondstring = '\n'
        if len(angle) > 0:
            for i in range(len(angle)):
                if angle[i][0] not in anglestring:
                    anglestring = anglestring + '  '.join(angle[i])
        else:
            anglestring = '\n'
        if len(dihedral) > 0:
            for i in range(len(dihedral)):
                if dihedral[i][0] not in dihedralstring:
                    dihedralstring = dihedralstring + '  '.join(dihedral[i])
        else:
            dihedralstring = '\n'
        if len(improper) > 0:
            for i in range(len(improper)):
                if improper[i][0] not in improperstring:
                    improperstring = improperstring + '  '.join(improper[i])
        else:
            improperstring = '\n'
        if len(nonbonded) > 0:
            for i in range(len(nonbonded)):
                if nonbonded[i][0] not in nonbondedstring:
                    nonbondedstring = nonbondedstring + '  '.join(nonbonded[i])
        else:
            nonbondedstring = '\n'



        return template.format(masstring,bondstring,anglestring,dihedralstring,improperstring,nonbondedstring)


    def writeCombinedFrcmod(self):
        self.readfiles('{}.frcmod')
        frcmodString = self.joinFrcmod()
        self.write('LIG.frcmod', frcmodString)



