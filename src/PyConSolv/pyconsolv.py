#!/usr/bin/env python3
import sys
import os
import argparse

from PyConSolv.ConfGen import PyConSolv


def main():
    ver = '1.0.0.1'
    parser = argparse.ArgumentParser(prog = 'PyConSolv', description='Process commandline arguments for PyconSolv')
    parser.add_argument('input', help = 'input file in XYZ format')
    parser.add_argument('-c', '--charge',  nargs='?', default=0, type=int, help = 'charge of the system, default 0')
    parser.add_argument('-m', '--method', nargs='?', default='PBE0', type=str, help='ORCA optimization/frequency calculations method of choice, default PBE0')
    parser.add_argument('-b', '--basis', nargs='?', default='def2-SVP', type=str, help='basis set to be used for calculations, default def2-SVP')
    parser.add_argument('-d', '--dispersion', nargs='?', default='D4', type=str, help='dispersion corrections, default = D4')
    parser.add_argument('-s', '--solvent', nargs='?', default='Water', type=str, help='solvent to be used for MD simulations/ OM Calculations, default Water')
    parser.add_argument('-p', '--cpu', nargs='?', default=12, type=int, help='number of cpu cores to be used for calculations, default 12')
    parser.add_argument('-mult', '--multiplicity',  nargs='?', default=1, type=int, help = 'multiplicity of the system, default 1')
    parser.add_argument('-v', '--version', action = 'version', version = '%(prog)s {}'.format(ver))

    args = parser.parse_args()

    inputfilepath = os.path.abspath(args.input)

    if '.xyz' not in inputfilepath:
        print('Path does not contain a valid XYZ file\n')
        sys.exit()
    conf = PyConSolv(inputfilepath)
    conf.run(charge= args.charge , method = args.method, basis = args.basis , dsp = args.dispersion , cpu = args.cpu ,
            solvent = args.solvent, multiplicity = args.multiplicity )
    sys.exit()

if __name__ == '__main__':
    main()