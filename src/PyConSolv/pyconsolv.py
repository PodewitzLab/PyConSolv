#!/usr/bin/env python3
import sys
import os
import argparse

from PyConSolv.ConfGen import PyConSolv
from PyConSolv.misc.analysis import Analysis


def main():
    ver = '1.0.6.3'
    parser = argparse.ArgumentParser(prog='PyConSolv', description='Process commandline arguments for PyconSolv')
    parser.add_argument('input', help='input file in XYZ format')

    # parametrization
    parser.add_argument('-c', '--charge', nargs='?', default=0, type=int, help='charge of the system, default 0')
    parser.add_argument('-m', '--method', nargs='?', default='PBE0', type=str,
                        help='ORCA optimization/frequency calculations method of choice, default PBE0')
    parser.add_argument('-b', '--basis', nargs='?', default='def2-SVP', type=str,
                        help='basis set to be used for calculations, default def2-SVP')
    parser.add_argument('-d', '--dispersion', nargs='?', default='D4', type=str,
                        help='dispersion corrections, default = D4')
    parser.add_argument('-s', '--solvent', nargs='?', default='Water', type=str,
                        help='solvent to be used for MD simulations/ OM Calculations, default Water')
    parser.add_argument('-p', '--cpu', nargs='?', default=12, type=int,
                        help='number of cpu cores to be used for calculations, default 12')
    parser.add_argument('-mult', '--multiplicity', nargs='?', default=1, type=int,
                        help='multiplicity of the system, default 1')
    parser.add_argument('-noopt', '--noopt', action='store_false',
                        help='do not perform geometry optimization for parametrization')
    parser.add_argument('-box', '--box', nargs='?', default=10, type=int,
                        help='set the box size to use with ambertools, for solvating the system')
    parser.add_argument('-e', '--engine', nargs='?', default='amber', type=str,
                        help='MD engine to be used for equilibration and simulation')
    parser.add_argument('-rst', '--restraint', action='store_true',
                        help='set up system for a restrained simulation')
    parser.add_argument('-cart', '--cartesianrst', nargs='?', default=None, type=str,
                        help='set up system for a simulation with cartesian restraints')
    parser.add_argument('-cartstr', '--cartesianrststr', nargs='?', default=100, type=int,
                        help='strength of cartesian restraints in kcal/mol')

    # analysis
    parser.add_argument('-a', '--analyze', action='store_true', help='analyze a simulation')
    parser.add_argument('-nosp', '--nosp', action='store_true', help='do not run single point calculations')
    parser.add_argument('-mask', '--mask', nargs='?', default=0, type=str, help='atomid mask for clustering')
    parser.add_argument('-cluster', '--cluster', nargs='?', default=0, type=str, help='clustering method')
    parser.add_argument('-qmmm', '--qmmm', action='store_true',
                        help='use a qmmm approach to determine cluster energy ranking')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(ver))

    args = parser.parse_args()

    inputfilepath = os.path.abspath(args.input)

    if args.analyze:

        if not args.mask:
            print(
                'Warning, you have not provided an input mask for alignment, please provide a list of atom ids in the '
                'format: 1,2,3-10\n')
            mask = input()
        else:
            mask = args.mask
        if not args.cluster:
            print('Please provide a clustering method: [kmeans, dbscan, dpeaks, hierarchical]\n')
            cluster = input()
        else:
            cluster = args.cluster
        analysis = Analysis(path=inputfilepath, alignMask=mask)
        analysis.run(clustering=cluster, nosp=args.nosp, engine=args.engine, qmmm=args.qmmm)

    elif '.xyz' not in inputfilepath:
        print('Path does not contain a valid XYZ file\n')
        sys.exit()

    else:
        conf = PyConSolv(inputfilepath)
        conf.run(charge=args.charge, method=args.method, basis=args.basis, dsp=args.dispersion, cpu=args.cpu,
                 solvent=args.solvent, multiplicity=args.multiplicity, engine=args.engine, opt=args.noopt, box=args.box,
                 rst=args.restraint, cart = args.cartesianrst, cartstr = args.cartesianrststr)
    sys.exit()


if __name__ == '__main__':
    main()
