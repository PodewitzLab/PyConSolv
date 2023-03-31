import numpy as np


class Faker:
    def __init__(self, path: str, calc_label: str = 'orca_freq'):
        """
        Class containing methods for creating a dummy fchk and log file

        Parameters:
            :param string path: path of folder containing orca calculations
            :param string calc_label: label for orca calculations, default is orca_calc

        Class variables:

        """
        self.path = path
        self.calc_label = calc_label
        f = open(self.path + '/' + 'fakechk.fchk', 'w')
        f.close()
        f = open(self.path + '/' + 'fakelog.log', 'w')
        f.close()

    def writeout(self, list_f: list, title: str, factor: float = 1.):
        """
        Write list in gaussian fchk format (5 columns, 80 characters)

        Parameters:
            :param list list_f: 1D list of values to write out
            :param string title: header for information to be written out
            :param float factor: factor to multiply values by eg.:1.89 for Angstrom/Bohr

        Class variables:

        """
        f = open(self.path + '/' + 'fakechk.fchk', 'a')
        f.write(title)
        for i in range(len(list_f) // 5):
            for j in range(5):
                f.write(' {:14E}'.format(np.double(list_f[i * 5:(i + 1) * 5][j]) * factor))
            f.write('\n')
        if len(list_f) % 5 != 0:
            for i in range(len(list_f) % 5):
                f.write(' {:14E}'.format(np.double(list_f[-(i + 1)]) * factor))
        f.write('\n')
        f.close()

    def fakecrds(self, xyzfile: str = 'input.xyz'):
        """
        Write out optimized coordinates into fchk file

        Parameters:
            :param string xyzfile: name of the input xyz file, default is input.xyz

        Class variables:

        """
        # read most recent coordinates from optimized output structure
        f = open(self.path + '/' + xyzfile, 'r')
        count = 0
        file = []
        for line in f:
            count += 1
            if count <= 2:  # skip the first 2 lines in xyz file
                continue
            else:
                file.append(line.split()[1:])
        f.close()

        file = sum(file, [])  # flatten list to 1D
        self.writeout(file, 'Current cartesian coordinates    N=   ' + str(len(file)) + '\n', factor=1.88973)

    def fakeforce(self):
        """
        Write out forces into fchk file

        Parameters:

        Class variables:

        """
        f = open(self.path + '/' + self.calc_label + '.hess', 'r')
        index = 0
        size = -5
        hess = []
        for line in f:
            if '$hessian' in line:
                index = -1
                continue
            elif index == -1:
                size = int(line.split()[0])
            elif index == size + 1 and size > 0:
                index = 0
            elif index >= 1:
                hess.append(line.split()[1:])
            if "$" in line and size > 0:
                break
            if size < 0:
                index = 0
            else:
                index += 1
        f.close()
        del hess[-1]

        full_hess = []
        for i in range(size):
            tmp = []
            for j in range(int(len(hess) / size)):
                tmp = tmp + hess[i + size * j]
            full_hess.append(tmp)

        lin_hess = []
        for i in range(size):
            lin_hess = lin_hess + full_hess[i][:i + 1]
        self.writeout(lin_hess, 'Cartesian Force Constants                  R   N=         ' + str(
            int(((size + 1) * size) / 2)) + '\n', factor=1)
        # return lin_hess

    def fakeesp(self):
        ESP_data = []
        coords = []
        pointID = 0
        f = open(self.path + '/' + 'ESPfitpt.txt', 'r')
        for line in f:
            ESP_data.append(line.split())
        f.close()
        fit_points = float(ESP_data[0][0])
        del ESP_data[0]
        f = open(self.path + '/' + 'orca_freq.molden.chg', 'r')
        for line in f:
            coords.append(line.split())
        f.close()
        s1 = 'Atomic Center    1 is at  -0.855031 -0.000059  0.062990'
        s2 = 'ESP Fit Center   31 is at  -0.855031 -0.000059  2.035590'
        s3 = ''
        s4 = ''
        f = open(self.path + '/' + 'fakelog.log', 'w')
        f.write('''
 **********************************************************************

            Electrostatic Properties Using The SCF Density
            
 **********************************************************************\n\n'''
                )
        for i in range(len(coords)):
            pointID += 1
            f.write('       Atomic Center ' + '{:>4}'.format(pointID) + ' is at ' + '{:>10.6f}'.format(
                float(coords[i][1])) + '{:>10.6f}'.format(float(coords[i][2])) + '{:>10.6f}\n'.format(
                float(coords[i][3])))
        for i in range(len(ESP_data)):
            pointID += 1
            f.write('      ESP Fit Center ' + '{:>4}'.format(pointID) + ' is at ' + '{:>10.6f}'.format(
                float(ESP_data[i][0])) + '{:>10.6f}'.format(float(ESP_data[i][1])) + '{:>10.6f}\n'.format(
                float(ESP_data[i][2])))
        f.write('\n')
        f.write('''
 -----------------------------------------------------------------

              Electrostatic Properties (Atomic Units)

 -----------------------------------------------------------------
    Center     Electric         -------- Electric Field --------
               Potential          X             Y             Z
 -----------------------------------------------------------------
'''
                )
        pointID = 0
        for i in range(len(coords)):
            pointID += 1
            f.write(' {:>4}'.format(pointID) + ' Atom ' + '{:> .6f}\n'.format(float(coords[i][-1])))
        for i in range(len(ESP_data)):
            pointID += 1
            f.write(' {:>4}'.format(pointID) + ' Fit  ' + '{:> .6f}\n'.format(float(ESP_data[i][-1])))
        f.write('\n')
        f.close()

    def run(self):
        self.fakecrds()
        # self.fakeesp()
        self.fakeforce()
