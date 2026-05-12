"""Polar-site detection and TIP3P probe placement for water-interaction
charge fitting.

Given a structure + connectivity + (optional) CGenFF atom types, return
per-polar-site descriptors and TIP3P water geometries placed in the
FFTK-canonical orientation for each site. Pure-functional; no QM, no
I/O beyond the optional PDB writer.

Geometry conventions follow FFTK:
  - Donor probes (X-H...O):   water O at 2.0 A from H along the X-H axis
  - Acceptor probes (X...H-O): water H at 1.85 A from the acceptor along
                                a lone-pair direction
  - TIP3P internal geometry:  O-H = 0.9572 A, H-O-H = 104.52 deg

Limitations:
  - Hybridization is read from CGenFF type prefixes (OG2*/NG2* = sp2;
    OG3*/NG3* = sp3) when atom_types are supplied. Without types the
    fallback is element + neighbour count, which mis-classifies aromatic
    or amide N-H (pyrrole, amide) as sp3 and may emit a spurious LP
    along an axis that lies in the aromatic plane. Pass CGenFF types
    whenever they are available.
  - Polar sites bonded directly to a metal are NOT filtered here; the
    caller should drop them from the list before running QM.
"""
import os
from dataclasses import dataclass
from typing import List, Optional

import numpy as np


TIP3P_OH = 0.9572
TIP3P_HOH_DEG = 104.52
DONOR_O_DIST = 2.0
ACCEPTOR_H_DIST = 1.85
TET_ANGLE_DEG = 109.471

POLAR_ELEMENTS = {'N', 'O', 'S', 'F'}

SP2_TYPE_PREFIXES = ('OG2', 'NG2')
SP3_TYPE_PREFIXES = ('OG3', 'NG3')


@dataclass
class PolarSite:
    heavy_idx: int
    role: str                       # 'donor' | 'acceptor'
    h_idx: Optional[int]
    direction: np.ndarray
    label: str


@dataclass
class Probe:
    site: PolarSite
    water_coords: np.ndarray        # (3, 3): [O, H1, H2]


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError('cannot normalize zero vector')
    return v / n


def _perpendicular(v: np.ndarray) -> np.ndarray:
    v = _unit(v)
    candidate = np.array([1.0, 0.0, 0.0])
    if abs(v[0]) > 0.9:
        candidate = np.array([0.0, 1.0, 0.0])
    return _unit(np.cross(v, candidate))


def _rotate(vec: np.ndarray, axis: np.ndarray, angle_rad: float) -> np.ndarray:
    axis = _unit(axis)
    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)
    return (vec * cos_a
            + np.cross(axis, vec) * sin_a
            + axis * np.dot(axis, vec) * (1 - cos_a))


def _buildNeighbours(n_atoms: int, bonds: list, elements: list) -> tuple:
    neighbours = [[] for _ in range(n_atoms)]
    for (i, j) in bonds:
        neighbours[i].append(j)
        neighbours[j].append(i)
    hydrogens = [[j for j in nbrs if elements[j].upper() == 'H']
                 for nbrs in neighbours]
    return neighbours, hydrogens


def _hybridization(atom_idx: int, atom_types: list, elements: list,
                    neighbours: list, hydrogens: list) -> str:
    if atom_idx < len(atom_types) and atom_types[atom_idx]:
        t = atom_types[atom_idx].upper()
        if any(t.startswith(p) for p in SP2_TYPE_PREFIXES):
            return 'sp2'
        if any(t.startswith(p) for p in SP3_TYPE_PREFIXES):
            return 'sp3'
    el = elements[atom_idx].upper()
    n_heavy = len([j for j in neighbours[atom_idx]
                    if elements[j].upper() != 'H'])
    n_h = len(hydrogens[atom_idx])
    if el == 'O':
        if n_heavy == 1 and n_h == 0:
            return 'sp2'
        return 'sp3'
    if el == 'N':
        if n_heavy + n_h == 3:
            return 'sp3'
        if n_heavy + n_h == 2:
            return 'sp2'
    return 'sp3'


def _spLonepairs(heavy_idx: int, coords: np.ndarray, elements: list,
                  neighbours: list, hybridization: str) -> list:
    """Return list of unit vectors along the lone pair directions."""
    centre = coords[heavy_idx]
    nbrs = neighbours[heavy_idx]
    if not nbrs:
        return []
    bond_dirs = [_unit(coords[j] - centre) for j in nbrs]

    if hybridization == 'sp2':
        if len(nbrs) == 1:
            c = nbrs[0]
            c_nbrs = [j for j in neighbours[c] if j != heavy_idx]
            if not c_nbrs:
                return []
            v1 = coords[c] - centre
            v2 = coords[c_nbrs[0]] - coords[c]
            normal = np.cross(v1, v2)
            if np.linalg.norm(normal) < 1e-8:
                return []
            normal = _unit(normal)
            front = -bond_dirs[0]
            return [_unit(_rotate(front, normal, np.radians(60))),
                    _unit(_rotate(front, normal, -np.radians(60)))]
        if len(nbrs) == 2:
            bisector = _unit(bond_dirs[0] + bond_dirs[1])
            return [-bisector]
        return []

    if hybridization == 'sp3':
        k = len(nbrs)
        if k == 3:
            return [_unit(-sum(bond_dirs, np.zeros(3)))]
        if k == 2:
            bisector = _unit(bond_dirs[0] + bond_dirs[1])
            perp = np.cross(bond_dirs[0], bond_dirs[1])
            if np.linalg.norm(perp) < 1e-8:
                return []
            perp = _unit(perp)
            half = np.radians(TET_ANGLE_DEG / 2)
            front = -bisector
            return [_unit(front * np.cos(half) + perp * np.sin(half)),
                    _unit(front * np.cos(half) - perp * np.sin(half))]
        if k == 1:
            front = -bond_dirs[0]
            perp = _perpendicular(bond_dirs[0])
            half = np.radians(TET_ANGLE_DEG / 2)
            return [_unit(front * np.cos(half) + perp * np.sin(half)),
                    _unit(front * np.cos(half) - perp * np.sin(half))]
    return []


def detectPolarSites(coords: np.ndarray, elements: list, bonds: list,
                      atom_types: list = None) -> list:
    """Identify donor (X-H) and acceptor (lone-pair) probe sites.

    coords     : (N, 3) Angstroms.
    elements   : list of N element symbols.
    bonds      : list of (i, j) 0-based pairs (one entry per bond).
    atom_types : optional list of N CGenFF atom-type strings; used for
                 hybridization. Without it, the fallback is element +
                 neighbour-count heuristics (see module docstring).

    Returns a list of PolarSite descriptors, one per probe direction.
    """
    if atom_types is None:
        atom_types = [''] * len(elements)
    neighbours, hydrogens = _buildNeighbours(len(elements), bonds, elements)

    sites = []
    for i, el in enumerate(elements):
        if el.upper() not in POLAR_ELEMENTS:
            continue
        for h in hydrogens[i]:
            sites.append(PolarSite(
                heavy_idx=i,
                role='donor',
                h_idx=h,
                direction=_unit(coords[h] - coords[i]),
                label='donor_{}{}_H{}'.format(el.upper(), i, h),
            ))
        hyb = _hybridization(i, atom_types, elements, neighbours, hydrogens)
        for k, lp in enumerate(_spLonepairs(i, coords, elements,
                                            neighbours, hyb)):
            sites.append(PolarSite(
                heavy_idx=i,
                role='acceptor',
                h_idx=None,
                direction=lp,
                label='acceptor_{}_{}{}_lp{}'.format(hyb, el.upper(),
                                                     i, k + 1),
            ))
    return sites


def placeWater(site: PolarSite, coords: np.ndarray) -> np.ndarray:
    """Return TIP3P water coords (3, 3) for a probe placed at ``site``."""
    half = np.radians(TIP3P_HOH_DEG / 2)
    if site.role == 'donor':
        H = coords[site.h_idx]
        axis = site.direction
        O = H + DONOR_O_DIST * axis
        front = axis
        perp = _perpendicular(front)
        H1 = O + TIP3P_OH * (front * np.cos(half) + perp * np.sin(half))
        H2 = O + TIP3P_OH * (front * np.cos(half) - perp * np.sin(half))
        return np.array([O, H1, H2])

    if site.role == 'acceptor':
        X = coords[site.heavy_idx]
        axis = site.direction
        H1 = X + ACCEPTOR_H_DIST * axis
        O = H1 + TIP3P_OH * axis
        h1_dir = -axis
        perp = _perpendicular(h1_dir)
        full = np.radians(TIP3P_HOH_DEG)
        h2_dir = h1_dir * np.cos(full) + perp * np.sin(full)
        H2 = O + TIP3P_OH * h2_dir
        return np.array([O, H1, H2])

    raise ValueError('unknown site role: {}'.format(site.role))


def placeWaters(sites: list, coords: np.ndarray) -> list:
    return [Probe(site=s, water_coords=placeWater(s, coords)) for s in sites]


def writeProbePDB(probe: Probe, ligand_coords: np.ndarray,
                   ligand_elements: list, out_path: str) -> str:
    """Write a single ligand+water PDB for visual inspection."""
    with open(out_path, 'w') as f:
        f.write('REMARK   {}\n'.format(probe.site.label))
        idx = 1
        for el, xyz in zip(ligand_elements, ligand_coords):
            name = '{}{}'.format(el, idx)[:4]
            f.write('ATOM  {:>5d} {:<4s} LIG  {:>4d}    '
                    '{:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00      LIG {:>2s}\n'
                    .format(idx, name, 1, xyz[0], xyz[1], xyz[2], el))
            idx += 1
        for name, el, xyz in zip(['OH2', 'H1', 'H2'], ['O', 'H', 'H'],
                                  probe.water_coords):
            f.write('ATOM  {:>5d} {:<4s} TIP3 {:>4d}    '
                    '{:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00      WAT {:>2s}\n'
                    .format(idx, name, 2, xyz[0], xyz[1], xyz[2], el))
            idx += 1
        f.write('END\n')
    return out_path


def writeAllProbes(probes: list, ligand_coords: np.ndarray,
                    ligand_elements: list, out_dir: str) -> list:
    """Dump per-probe PDBs into ``out_dir``. Returns list of file paths."""
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    paths = []
    for i, probe in enumerate(probes):
        path = os.path.join(out_dir,
                             'probe_{:03d}_{}.pdb'.format(i, probe.site.label))
        writeProbePDB(probe, ligand_coords, ligand_elements, path)
        paths.append(path)
    return paths