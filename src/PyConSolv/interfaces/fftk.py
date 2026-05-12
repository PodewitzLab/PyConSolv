"""FFTK-style bonded parameter derivation in native CHARMM36 format.

Phase A of the native-CHARMM pipeline: derives bond (Kb, b0) and angle
(Ktheta, theta0) parameters from an ORCA-optimized geometry and Cartesian
Hessian via the Seminario (1996) projection method. Emits CHARMM RTF
(MASS + minimal metal RESI) and PRM (BONDS + ANGLES) fragments that can
be merged with CGenFF ligand output.

Charges, dihedrals and LJ for metals are handled outside this module:
  - charges      : RESP via MultiWfn (Phase B, mirrors the AMBER path and
                   overwrites CGenFF atomic charges by atom name).
  - dihedrals    : Kchi=0 soft-coordination stubs for metal-crossing torsions
                   (parallel to MCPB.py zeroing in the AMBER path); CGenFF
                   values used elsewhere. No per-torsion QM PES fitting is
                   performed by design.
  - metal LJ     : looked up from db/charmm_metal_nonbonded.txt (broad
                   superset: alkali, alkaline earth, 3d/4d/5d transition
                   metals, main-group post-transition, lanthanides, actinides).
"""
import os

import numpy as np

from ..utils.colorgen import Color


HARTREE_TO_KCAL = 627.5094740631  # kcal/mol per Hartree
BOHR_TO_ANG = 0.52917721067       # A per Bohr
# Hessian conversion: Hartree/Bohr^2 -> kcal/mol/A^2
HA_BOHR2_TO_KCAL_ANG2 = HARTREE_TO_KCAL / (BOHR_TO_ANG ** 2)


# ---------------------------------------------------------------------------
# ORCA .hess parser
# ---------------------------------------------------------------------------

def parseOrcaHessian(hess_file: str) -> np.ndarray:
    """Read the $hessian block from an ORCA .hess file.

    ORCA stores the Hessian column-major in chunks of up to 5 columns.
    Returns a dense (3N, 3N) matrix in Hartree/Bohr^2.
    """
    if not os.path.isfile(hess_file):
        raise FileNotFoundError(hess_file)

    with open(hess_file, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines) and not lines[i].strip().startswith('$hessian'):
        i += 1
    if i >= len(lines):
        raise ValueError('No $hessian block found in {}'.format(hess_file))

    n = int(lines[i + 1].strip())
    H = np.zeros((n, n))
    i += 2
    col_start = 0
    while col_start < n:
        header = lines[i].split()
        ncols = len(header)
        i += 1
        for _ in range(n):
            parts = lines[i].split()
            row = int(parts[0])
            for k, val in enumerate(parts[1:]):
                H[row, col_start + k] = float(val)
            i += 1
        col_start += ncols
    return H


# ---------------------------------------------------------------------------
# Seminario core
# ---------------------------------------------------------------------------

def _crossBlock(H: np.ndarray, i: int, j: int) -> np.ndarray:
    """Return -H[3i:3i+3, 3j:3j+3] symmetrized."""
    block = -H[3 * i:3 * i + 3, 3 * j:3 * j + 3]
    return 0.5 * (block + block.T)


def _projectedSum(block: np.ndarray, unit_vec: np.ndarray) -> float:
    """Seminario projection: sum_i |lambda_i| * (v_i . u)^2."""
    eigvals, eigvecs = np.linalg.eigh(block)
    total = 0.0
    for k in range(3):
        proj = float(np.dot(eigvecs[:, k], unit_vec))
        total += abs(eigvals[k]) * proj * proj
    return total


def seminarioBond(H_bohr: np.ndarray, coords_ang: np.ndarray,
                  i: int, j: int) -> tuple:
    """Compute (Kb, b0) for bond i-j.

    Kb in kcal/mol/A^2 (CHARMM convention E = Kb (r - b0)^2, no 1/2).
    b0 in Angstroms.
    """
    rij_ang = coords_ang[j] - coords_ang[i]
    b0 = float(np.linalg.norm(rij_ang))
    u = rij_ang / b0
    # Hessian is in Hartree/Bohr^2; project and convert. The conversion from
    # Hartree/Bohr^2 to kcal/mol/A^2 already accounts for the length units.
    k_hartree_bohr2 = _projectedSum(_crossBlock(H_bohr, i, j), u)
    # CHARMM/AMBER form E = K (r-r0)^2 has d2E/dr2 = 2K, so divide by 2.
    Kb = 0.5 * k_hartree_bohr2 * HA_BOHR2_TO_KCAL_ANG2
    return Kb, b0


def seminarioAngle(H_bohr: np.ndarray, coords_ang: np.ndarray,
                   i: int, j: int, k: int) -> tuple:
    """Compute (Ktheta, theta0) for angle i-j-k (j is apex).

    Ktheta in kcal/mol/rad^2 (CHARMM E = Ktheta (theta - theta0)^2).
    theta0 in degrees.
    """
    rij = coords_ang[i] - coords_ang[j]
    rkj = coords_ang[k] - coords_ang[j]
    R_ij = float(np.linalg.norm(rij))
    R_kj = float(np.linalg.norm(rkj))
    u_ji = rij / R_ij
    u_jk = rkj / R_kj
    cos_theta = float(np.clip(np.dot(u_ji, u_jk), -1.0, 1.0))
    theta0 = float(np.degrees(np.arccos(cos_theta)))

    # Plane normal and in-plane perpendiculars
    n = np.cross(u_jk, u_ji)
    n_norm = float(np.linalg.norm(n))
    if n_norm < 1e-8:
        # Near-linear angle; Seminario degenerate. Return placeholder.
        return 0.0, theta0
    n /= n_norm
    u_PA = np.cross(n, u_ji)   # perpendicular to j-i, in-plane
    u_PC = np.cross(u_jk, n)   # perpendicular to j-k, in-plane

    kA_bohr2 = _projectedSum(_crossBlock(H_bohr, i, j), u_PA)
    kC_bohr2 = _projectedSum(_crossBlock(H_bohr, k, j), u_PC)

    # Convert k_A, k_C to kcal/mol/A^2 (length^-2 units)
    kA = kA_bohr2 * HA_BOHR2_TO_KCAL_ANG2
    kC = kC_bohr2 * HA_BOHR2_TO_KCAL_ANG2
    if kA <= 0.0 or kC <= 0.0:
        return 0.0, theta0
    inv = 1.0 / (R_ij * R_ij * kA) + 1.0 / (R_kj * R_kj * kC)
    Ktheta = 0.5 / inv  # factor 1/2 for CHARMM E = K(theta-theta0)^2
    return Ktheta, theta0


# ---------------------------------------------------------------------------
# Atom-type assignment
# ---------------------------------------------------------------------------

def metalAtomType(element: str) -> str:
    """Return a 4-char CHARMM atom type name for a metal element symbol.

    Convention: uppercase 'M' + uppercase element symbol, truncated to 4.
    e.g. Fe -> MFE, Cu -> MCU, Pt -> MPT, W -> MW.
    """
    sym = element.strip().upper()
    return ('M' + sym)[:4]


# ---------------------------------------------------------------------------
# Metal non-bonded lookup
# ---------------------------------------------------------------------------

def loadMetalNonbonded(db_path: str) -> dict:
    """Parse db/charmm_metal_nonbonded.txt into {element: (epsilon, Rmin/2)}.

    Values are in CHARMM conventions: epsilon in -kcal/mol, Rmin/2 in A.
    Lines starting with '#' are ignored.
    """
    table = {}
    if not os.path.isfile(db_path):
        return table
    with open(db_path, 'r') as f:
        for line in f:
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            element = parts[0].upper()
            eps = float(parts[1])
            rmin2 = float(parts[2])
            table[element] = (eps, rmin2)
    return table


# ---------------------------------------------------------------------------
# FFTKBondedInterface
# ---------------------------------------------------------------------------

class FFTKBondedInterface:
    """Seminario-based bonded parameter derivation for CHARMM36.

    Consumes an ORCA-optimized geometry + Hessian and produces native
    CHARMM RTF/PRM fragments suitable for merging with CGenFF output.

    Scope is option (c): any bond/angle that (a) touches a metal, OR
    (b) involves an atom whose CGenFF penalty exceeds `penalty_threshold`.
    Everything else is left to CGenFF.
    """

    def __init__(self, path: str, penalty_threshold: float = 50.0):
        self.path = path
        self.penalty_threshold = penalty_threshold
        self.status = 0
        if path and not os.path.isdir(path):
            os.makedirs(path, exist_ok=True)

    # ----- scope selection -----------------------------------------------

    def _shouldRefitBond(self, a: int, b: int, metal_indices: set,
                          penalties: dict) -> bool:
        if a in metal_indices or b in metal_indices:
            return True
        return (penalties.get(a, 0.0) > self.penalty_threshold or
                penalties.get(b, 0.0) > self.penalty_threshold)

    def _shouldRefitAngle(self, a: int, b: int, c: int, metal_indices: set,
                           penalties: dict) -> bool:
        if metal_indices.intersection({a, b, c}):
            return True
        return any(penalties.get(x, 0.0) > self.penalty_threshold
                   for x in (a, b, c))

    # ----- main entry point ----------------------------------------------

    def deriveParameters(self, xyz_coords_ang: np.ndarray,
                         elements: list,
                         hessian: np.ndarray,
                         bonds: list,
                         angles: list,
                         metal_indices: list,
                         atom_types: list,
                         penalties: dict = None,
                         out_basename: str = 'fftk') -> tuple:
        """Derive bonded params and write RTF + PRM fragments.

        Parameters
        ----------
        xyz_coords_ang : (N, 3) array of Cartesian coords in Angstroms.
        elements       : list of N element symbols (e.g. 'Fe', 'N', 'C', ...).
        hessian        : (3N, 3N) Cartesian Hessian in Hartree/Bohr^2.
        bonds          : list of (i, j) 0-based atom index pairs.
        angles         : list of (i, j, k) with j the apex.
        metal_indices  : indices of metal atoms.
        atom_types     : list of N strings; CGenFF type for ligand atoms,
                         anything for metals (will be overridden with MFE etc).
        penalties      : dict {atom_idx: cgenff_penalty} for scope selection.
        out_basename   : filename stem for the output .rtf / .prm files.

        Returns (rtf_path, prm_path).
        """
        penalties = penalties or {}
        metal_set = set(metal_indices)

        # Convert coords to Bohr for compatibility with Hessian units.
        coords_bohr = xyz_coords_ang / BOHR_TO_ANG
        # We implement Seminario with coords in Angstroms *and* convert the
        # Hessian conversion factor at the end. Pass coords_ang through.

        # Override atom types for metals with Mxx convention.
        types = list(atom_types)
        for m in metal_indices:
            types[m] = metalAtomType(elements[m])

        # ----- derive bond parameters -----
        bond_params = []  # list of (type_i, type_j, Kb, b0)
        for (i, j) in bonds:
            if not self._shouldRefitBond(i, j, metal_set, penalties):
                continue
            Kb, b0 = seminarioBond(hessian, xyz_coords_ang, i, j)
            bond_params.append((types[i], types[j], Kb, b0))

        # ----- derive angle parameters -----
        angle_params = []  # (type_i, type_j, type_k, Ktheta, theta0)
        for (i, j, k) in angles:
            if not self._shouldRefitAngle(i, j, k, metal_set, penalties):
                continue
            Kth, th0 = seminarioAngle(hessian, xyz_coords_ang, i, j, k)
            angle_params.append((types[i], types[j], types[k], Kth, th0))

        # ----- write RTF (metal MASS + minimal single-atom RESI blocks) -----
        rtf_path = os.path.join(self.path, out_basename + '.rtf')
        nb_db = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                             'db', 'charmm_metal_nonbonded.txt')
        nb_table = loadMetalNonbonded(nb_db)
        self._writeRTF(rtf_path, metal_indices, elements, types, nb_table)

        # ----- write PRM (BONDS + ANGLES sections, plus NONBONDED for metals) -----
        prm_path = os.path.join(self.path, out_basename + '.prm')
        self._writePRM(prm_path, bond_params, angle_params,
                       metal_indices, elements, types, nb_table)

        self.status = 1
        return rtf_path, prm_path

    # ----- writers --------------------------------------------------------

    def _writeRTF(self, rtf_path: str, metal_indices: list, elements: list,
                  types: list, nb_table: dict):
        lines = [
            '* Metal topology fragment produced by PyConSolv (Seminario/FFTK)\n',
            '*\n',
            '   36   1\n',
            '\n',
        ]
        # MASS entries for each unique metal type
        seen = set()
        for m in metal_indices:
            t = types[m]
            if t in seen:
                continue
            seen.add(t)
            sym = elements[m]
            mass = _elementMass(sym)
            lines.append('MASS  -1  {:<6s} {:8.4f} ! {} center\n'.format(
                t, mass, sym))
        lines.append('\n')
        # One single-atom RESI per metal type, charge placeholder 0.0
        # (real charge comes from the full-system RESP/CGenFF assignment;
        # this RESI exists only so ParmEd can find the atom-type/mass pair).
        for t in sorted(seen):
            resname = t[:4]
            lines.append('RESI {:<4s}       0.00\n'.format(resname))
            lines.append('GROUP\n')
            lines.append('ATOM {:<4s} {:<6s}  0.00\n'.format(resname, t))
            lines.append('\n')
        lines.append('END\n')
        with open(rtf_path, 'w') as f:
            f.writelines(lines)

    def _writePRM(self, prm_path: str,
                  bond_params: list, angle_params: list,
                  metal_indices: list, elements: list, types: list,
                  nb_table: dict):
        lines = [
            '* Metal + high-penalty bonded parameters from Seminario/FFTK\n',
            '*\n',
            '\n',
        ]
        if bond_params:
            lines.append('BONDS\n')
            lines.append('!atom_i atom_j   Kb       b0\n')
            for (ti, tj, Kb, b0) in bond_params:
                lines.append('{:<6s} {:<6s} {:10.3f} {:8.4f}\n'.format(
                    ti, tj, Kb, b0))
            lines.append('\n')

        if angle_params:
            lines.append('ANGLES\n')
            lines.append('!atom_i atom_j atom_k   Ktheta   theta0\n')
            for (ti, tj, tk, Kth, th0) in angle_params:
                lines.append('{:<6s} {:<6s} {:<6s} {:10.3f} {:10.4f}\n'.format(
                    ti, tj, tk, Kth, th0))
            lines.append('\n')

        # Soft dihedral stubs for any torsion crossing a metal:
        # CHARMM lookups happen by type, so emit a wildcard around metal types.
        # We don't know the full dihedral set at this stage; charmm_builder
        # will fall back to CGenFF types for non-metal torsions.
        seen_metal_types = sorted({types[m] for m in metal_indices})
        if seen_metal_types:
            lines.append('DIHEDRALS\n')
            lines.append('!wildcard soft torsions around metal centers\n')
            for mt in seen_metal_types:
                lines.append('X      {:<6s} X      X      0.0    1    0.0\n'.format(mt))
            lines.append('\n')

        # NONBONDED entries for metals, from the lookup table.
        if seen_metal_types:
            lines.append('NONBONDED\n')
            lines.append('!type            epsilon     Rmin/2\n')
            for m in metal_indices:
                t = types[m]
                sym = elements[m].upper()
                if sym not in nb_table:
                    print(Color.YELLOW + 'No CHARMM NONBONDED entry for {}; '
                                          'using generic 0.0 / 2.0. Expand '
                                          'db/charmm_metal_nonbonded.txt.'.format(sym)
                          + Color.END)
                    eps, rmin2 = 0.0, 2.0
                else:
                    eps, rmin2 = nb_table[sym]
                lines.append('{:<6s}   0.0   {:10.6f}   {:10.6f}\n'.format(
                    t, eps, rmin2))
            lines.append('\n')

        lines.append('END\n')
        with open(prm_path, 'w') as f:
            f.writelines(lines)


# ---------------------------------------------------------------------------
# Minimal element -> mass table for MASS lines
# ---------------------------------------------------------------------------

_ATOMIC_MASSES = {
    'H': 1.00800, 'C': 12.0107, 'N': 14.00700, 'O': 15.99900,
    'F': 18.99800, 'P': 30.97400, 'S': 32.06000, 'CL': 35.45300,
    'BR': 79.90400, 'I': 126.90400,
    'NA': 22.98977, 'MG': 24.30500, 'K': 39.09830, 'CA': 40.07800,
    'MN': 54.93805, 'FE': 55.84500, 'CO': 58.93320, 'NI': 58.69340,
    'CU': 63.54600, 'ZN': 65.38000,
    'RU': 101.07000, 'RH': 102.90550, 'PD': 106.42000, 'AG': 107.86820,
    'IR': 192.21700, 'PT': 195.08400, 'AU': 196.96657, 'HG': 200.59200,
    'W': 183.84000, 'MO': 95.95000, 'RE': 186.20700,
}


def _elementMass(element: str) -> float:
    return _ATOMIC_MASSES.get(element.strip().upper(), 0.0)
