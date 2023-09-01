class Cluster:
    def __init__(self, clustering, engine = 'amber'):
        """
        Class to generate cpptraj input files for clustering

        Parameters:
            :param str clustering: clustering method

        Class variables:
            - self.engine - simulation engine to be used
            - self.runtypes - clustering methods in cpptraj
            - self.clustering - chosen clustering method
            - self.kmeans - template input file for kmeans clustering
            - self.dpeaks - template input file for dpeaks clustering
            - self.dbscan - template input file for dbscan clustering
            - self.hierarchical - template input file for hierarchical clustering
        """

        self.engine = engine
        if engine == 'gromacs':
            self.suffix = 'xtc'
        else:
            self.suffix = 'nc'
        self.runtypes = ['kmeans', 'hierarchical', 'dbscan', 'dpeaks']
        self.clustering = clustering
        self.kmeans = r"""parm LIG_dry.prmtop
trajin dry_aligned.{}
cluster c1 \
 kmeans clusters {} randompoint maxit 500 \
 rms !@H= \
 sieve {} random \
 out cnumvtime.dat \
 summary summary.dat \
 info info.dat \
 cpopvtime cpopvtime.agr normframe \
 repout rep repfmt pdb \
 singlerepout singlerep.nc singlerepfmt netcdf \
 avgout avg avgfmt pdb
run
"""
        self.dbscan = r"""parm LIG_dry.prmtop
trajin dry_aligned.{}
cluster c1 \
 dbscan minpoints {} epsilon {} kdist {}\
 rms !@H= \
 sieve {} random \
 out cnumvtime.dat \
 summary summary.dat \
 info info.dat \
 cpopvtime cpopvtime.agr normframe \
 repout rep repfmt pdb \
 singlerepout singlerep.nc singlerepfmt netcdf \
 avgout avg avgfmt pdb
run
"""
        self.dpeaks = r"""parm LIG_dry.prmtop
trajin dry_aligned.{}
cluster c1 \
 dpeaks epsilon {} \
 rms !@H= \
 sieve {} random \
 out cnumvtime.dat \
 summary summary.dat \
 info info.dat \
 cpopvtime cpopvtime.agr normframe \
 repout rep repfmt pdb \
 singlerepout singlerep.nc singlerepfmt netcdf \
 avgout avg avgfmt pdb
run
"""
        self.hierarchical = r"""parm LIG_dry.prmtop
trajin dry_aligned.{}
cluster c1 \
 hieragglo clusters {} {} \
 rms !@H= \
 sieve {} random \
 out cnumvtime.dat \
 summary summary.dat \
 info info.dat \
 cpopvtime cpopvtime.agr normframe \
 repout rep repfmt pdb \
 singlerepout singlerep.nc singlerepfmt netcdf \
 avgout avg avgfmt pdb
run
"""

    def createInput(self):
        inputstr = None
        match self.clustering:
            case 'kmeans':
                clusters = input('Please enter the number of desired clusters: (default 10)\n') or '10'
                inputstr = self.kmeans.format(self.suffix, clusters, '{}')
            case 'dbscan':
                minpoints = input('Please enter the minimum number of points for a cluster: ()\n')
                epsilon = input('Please enter the epsilon value: ()\n')
                kdist = input('Please enter the kdist: ()\n')
                inputstr = self.dbscan.format(self.suffix, minpoints, epsilon, kdist, '{}')
            case 'dpeaks':
                epsilon = input('Please enter the epsilon value: ()\n')
                inputstr = self.dpeaks.format(self.suffix, epsilon, '{}')

            case 'hierarchical':
                clusters = input('Please enter the number of desired clusters: (default 10)\n') or '10'
                linkage = input('Please enter the linkage type (linkage, averagelinkage, complete): (Default averagelinkage)\n') or 'averagelinkage'
                match linkage:
                    case 'linkage':
                        inputstr = self.hierarchical.format(self.suffix, clusters, linkage, '{}')
                    case 'averagelinkage':
                        inputstr = self.hierarchical.format(self.suffix, clusters, linkage, '{}')
                    case 'complete':
                        inputstr = self.hierarchical.format(self.suffix, clusters, linkage, '{}')

        sieve = input('Please enter the value for the frame number sieve: (default 10)\n') or '10'
        return inputstr.format(sieve)
