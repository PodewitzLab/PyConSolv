class Cluster:
    def __init__(self, clustering):
        self.runtypes = ['kmeans', 'hierarchical', 'dbscan', 'dpeaks']
        self.clustering = clustering
        self.kmeans = r"""parm LIG_dry.prmtop
trajin dry_aligned.nc
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
trajin dry_aligned.nc
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
trajin dry_aligned.nc
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
trajin dry_aligned.nc
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
                inputstr = self.kmeans.format(clusters, '{}')
            case 'dbscan':
                minpoints = input('Please enter the minimum number of points for a cluster: ()\n')
                epsilon = input('Please enter the epsilon value: ()\n')
                kdist = input('Please enter the kdist: ()\n')
                inputstr = self.dbscan.format(minpoints, epsilon, kdist, '{}')
            case 'dpeaks':
                epsilon = input('Please enter the epsilon value: ()\n')
                inputstr = self.dbscan.format(epsilon, '{}')

            case 'hierarchical':
                clusters = input('Please enter the number of desired clusters: (default 10)\n') or '10'
                linkage = input('Please enter the linkage type: (Default average)\n') or 'average'
                match linkage:
                    case 'linkage':
                        inputstr = self.hierarchical.format(clusters, linkage)
                    case 'averagelinkage':
                        inputstr = self.hierarchical.format(clusters, linkage)
                    case 'complete':
                        inputstr = self.hierarchical.format(clusters, linkage)

        sieve = input('Please enter the value for the frame number sieve: (default 10)\n') or '10'
        return inputstr.format(sieve)
