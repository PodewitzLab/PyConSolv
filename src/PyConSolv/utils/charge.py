class ChargeChanger:
    def __init__(self):
        self.iterator = 0

    def change(self, filein, fileout, resname, charges, fragmented = False, mapfile = None):
        print(filein)
        fin = open(filein, 'r')  # read mol2 files corresponding to the pdb
        fout = open(fileout, 'w')  # create new mol2 file
        switch = 0
        for line in fin:
            if 'USER_CHARGES' in line:
                fout.write('RESP Charge\n')
            elif '@<TRIPOS>ATOM' in line:
                switch = 1
                fout.write(line)
            elif '@<TRIPOS>BOND' in line:
                fout.write(line)
                switch = 2
            elif '@<TRIPOS>SUBSTRUCTURE' in line:
                fout.write(line)
                switch = 3
            elif '@<TRIPOS>MOLECULE' in line:
                fout.write(line)
                switch = 4
            else:
                if switch == 1:

                    tmp_line = line.split()[:-1]  # remove old charge
                    # if fragmented:
                    tmp_line.append(float(charges[self.iterator][-1]))  # add resp charge
                    # else:
                    # tmp_line.append(float(charges[self.iterator]))  # add resp charge
                    tmp_line[-2] = resname
                    if mapfile is not None:
                        tmp_line[5] = mapfile[self.iterator][-1]
                    fout.write(
                        '{:>7} {:<7}    {:>7}    {:>7}    {:>7}   {:<2}    {}  {:>2} {:>9.6f}\n'.format(*tmp_line))
                    self.iterator += 1
                elif switch == 3:
                    fout.write('     1 {}         1 TEMP              0 ****  ****    0 ROOT\n'.format(resname))
                    switch = 0
                elif switch == 4:
                    fout.write(resname + '\n')
                    switch = 0
                else:
                    if 'bcc' in line:
                        fout.write(line.replace('bcc', 'RESP Charge'))
                    else:
                        fout.write(line.replace('ar', '1'))
        fout.close()
        fin.close()