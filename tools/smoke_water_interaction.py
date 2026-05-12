#!/usr/bin/env python3
"""End-to-end smoke test for the water-interaction charge-fitting pipeline.

Runs Phase A (polar-site detection + water placement) and Phase B (ORCA
single-point batch) on an XYZ structure, then prints a table of
interaction energies and dumps probe PDBs for visual review.

Examples
--------
    # Phase A only - no ORCA needed, just dump probe PDBs
    python tools/smoke_water_interaction.py methanol.xyz --no-orca

    # Full pipeline (requires ORCA in PATH, e.g. `module load orca/5.0.4`)
    python tools/smoke_water_interaction.py ligand.xyz -c 0 --cpu 1

    # Filter out polar atoms directly bonded to a metal (FFTK convention)
    python tools/smoke_water_interaction.py complex.xyz --skip-metal-bonded
"""
import argparse
import os
import sys
import time

# Run from repo root without installing the package.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src'))

import numpy as np

from PyConSolv.misc.inputparser import XYZ
from PyConSolv.misc.polar_sites import (
    detectPolarSites, placeWaters, writeAllProbes,
)
from PyConSolv.interfaces.water_interaction import WaterInteractionCharges
from PyConSolv.utils.colorgen import Color


DB_ROOT = os.path.join(HERE, '..', 'src', 'PyConSolv', 'db')


def load_structure(xyz_path: str):
    """Read XYZ + detect bonds via PyConSolv's existing connectivity logic.

    Returns: (coords, elements, bonds, metal_bonded_indices)
    """
    xyz = XYZ(os.path.join(DB_ROOT, 'atom-radius.txt'),
              os.path.join(DB_ROOT, 'metal-radius.txt'))
    xyz.readXYZ(xyz_path)
    xyz.calculateDistanceMatrix()
    xyz.generateAdjacencyMatrix()
    xyz.generateLinkList()

    coords = np.asarray(xyz.coords, dtype=float)
    elements = list(xyz.atoms)

    bonds = set()
    for i, nbrs in enumerate(xyz.linkList):
        for j in nbrs:
            bonds.add((min(i, j), max(i, j)))

    metal_bonded = set()
    for entry in xyz.metalBonds:
        parts = entry.split()
        a, b = int(parts[0]), int(parts[2])
        bonds.add((min(a, b), max(a, b)))
        metal_bonded.add(b)   # the ligand atom side
    return coords, elements, sorted(bonds), metal_bonded


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('xyz', help='input XYZ structure')
    parser.add_argument('-c', '--charge', type=int, default=0)
    parser.add_argument('-m', '--multiplicity', type=int, default=1)
    parser.add_argument('-o', '--out-dir', default='/tmp/pyconsolv_smoke')
    parser.add_argument('--level', default='HF 6-31G*',
                         help='ORCA SP level string (default: HF 6-31G*)')
    parser.add_argument('--cpu', type=int, default=1)
    parser.add_argument('--memory', type=int, default=1000)
    parser.add_argument('--no-orca', action='store_true',
                         help='only detect sites + dump probe PDBs')
    parser.add_argument('--skip-metal-bonded', action='store_true',
                         help='drop polar atoms directly bonded to a metal '
                              '(FFTK convention - they are not reliable probe sites)')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    coords, elements, bonds, metal_bonded = load_structure(args.xyz)
    print('Read {} atoms, {} bonds from {}'.format(
        len(elements), len(bonds), args.xyz))
    if metal_bonded:
        print('Polar atoms directly bonded to a metal: {}'.format(
            sorted(metal_bonded)))

    sites = detectPolarSites(coords, elements, bonds, atom_types=None)
    if args.skip_metal_bonded:
        before = len(sites)
        sites = [s for s in sites if s.heavy_idx not in metal_bonded]
        print('Filtered out {} metal-bonded site(s)'.format(before - len(sites)))

    n_donor = sum(1 for s in sites if s.role == 'donor')
    n_acc = sum(1 for s in sites if s.role == 'acceptor')
    print('\nPhase A: {} probe(s) -- {} donor, {} acceptor'.format(
        len(sites), n_donor, n_acc))
    for s in sites:
        print('  {}'.format(s.label))

    probes = placeWaters(sites, coords)
    probe_dir = os.path.join(args.out_dir, 'probes')
    writeAllProbes(probes, coords, elements, probe_dir)
    print('\nProbe PDBs: {}'.format(probe_dir))

    if args.no_orca or not probes:
        return 0

    print('\nPhase B: ORCA SP batch at {} ...'.format(args.level))
    wic = WaterInteractionCharges(args.out_dir,
                                    level=args.level,
                                    cpu=args.cpu, memory=args.memory)
    if not wic.checkpath():
        print(Color.RED + 'ORCA not available; aborting Phase B.' + Color.END)
        return 1

    t0 = time.time()
    results = wic.runBatch(probes, coords, elements,
                            charge=args.charge,
                            multiplicity=args.multiplicity)
    dt = time.time() - t0
    print('Done in {:.1f}s ({}/{} probes succeeded)'.format(
        dt, len(results), len(probes)))

    if not results:
        return 1
    print()
    print('  {:<40s} {:>12s} {:>12s} {:>8s}'.format(
        'label', 'E_int(kcal)', 'x1.16', 'r(A)'))
    print('  ' + '-' * 76)
    for r in results:
        print('  {:<40s} {:>12.3f} {:>12.3f} {:>8.3f}'.format(
            r.site_label, r.e_int_kcal, r.e_int_scaled, r.distance))
    return 0


if __name__ == '__main__':
    sys.exit(main())
