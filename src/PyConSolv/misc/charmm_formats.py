"""Minimal parsers and writers for CHARMM RTF and PRM files.

Only the sections that PyConSolv actually needs to merge are modelled
(MASS, RESI, PRES for RTF; BONDS, ANGLES, DIHEDRALS, IMPROPERS, NONBONDED
for PRM). Everything else is preserved as raw text blocks.
"""
import os


RTF_SECTIONS = ('MASS', 'RESI', 'PRES')
PRM_SECTIONS = ('BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS', 'IMPROPER',
                'NONBONDED', 'NBFIX', 'CMAP', 'HBOND')


def _strip_comment(line: str) -> str:
    idx = line.find('!')
    return line if idx < 0 else line[:idx]


def parseRTF(path: str) -> dict:
    """Return {'mass': [lines], 'resi': {name: [lines]}, 'pres': {name: [lines]}, 'header': [lines]}."""
    out = {'mass': [], 'resi': {}, 'pres': {}, 'header': []}
    current = None
    current_name = None
    with open(path, 'r') as f:
        for line in f:
            token = _strip_comment(line).strip().split()
            if not token:
                if current in ('resi', 'pres'):
                    out[current][current_name].append(line)
                continue
            head = token[0].upper()
            if head == 'MASS':
                out['mass'].append(line)
                current = None
            elif head == 'RESI':
                current = 'resi'
                current_name = token[1]
                out['resi'][current_name] = [line]
            elif head == 'PRES':
                current = 'pres'
                current_name = token[1]
                out['pres'][current_name] = [line]
            elif head == 'END':
                current = None
            elif current in ('resi', 'pres'):
                out[current][current_name].append(line)
            else:
                out['header'].append(line)
    return out


def writeRTF(data: dict, path: str):
    with open(path, 'w') as f:
        for line in data.get('header', []):
            f.write(line)
        for line in data.get('mass', []):
            f.write(line)
        f.write('\n')
        for name, lines in data.get('resi', {}).items():
            for ln in lines:
                f.write(ln)
            f.write('\n')
        for name, lines in data.get('pres', {}).items():
            for ln in lines:
                f.write(ln)
            f.write('\n')
        f.write('END\n')


def parsePRM(path: str) -> dict:
    """Return {section_name: [lines]} for BONDS/ANGLES/.../NONBONDED, plus 'header'."""
    sections = {s: [] for s in PRM_SECTIONS}
    sections['header'] = []
    current = 'header'
    with open(path, 'r') as f:
        for line in f:
            token = _strip_comment(line).strip().split()
            if not token:
                sections[current].append(line)
                continue
            head = token[0].upper()
            if head in PRM_SECTIONS:
                current = head
                continue
            if head == 'END':
                current = 'header'
                continue
            sections[current].append(line)
    return sections


def writePRM(data: dict, path: str):
    order = ['BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPERS', 'IMPROPER',
             'NONBONDED', 'NBFIX', 'CMAP', 'HBOND']
    with open(path, 'w') as f:
        for line in data.get('header', []):
            f.write(line)
        for section in order:
            lines = data.get(section, [])
            if not lines:
                continue
            f.write('\n{}\n'.format(section))
            for ln in lines:
                f.write(ln)
        f.write('\nEND\n')


def mergeRTF(files: list, out_path: str) -> str:
    merged = {'mass': [], 'resi': {}, 'pres': {}, 'header': []}
    seen_mass = set()
    for path in files:
        data = parseRTF(path)
        if not merged['header']:
            merged['header'] = data['header']
        for line in data['mass']:
            key = tuple(line.split()[:3])  # MASS idx type mass
            if key in seen_mass:
                continue
            seen_mass.add(key)
            merged['mass'].append(line)
        merged['resi'].update(data['resi'])
        merged['pres'].update(data['pres'])
    writeRTF(merged, out_path)
    return out_path


def mergePRM(files: list, out_path: str) -> str:
    merged = {s: [] for s in PRM_SECTIONS}
    merged['header'] = []
    seen = {s: set() for s in PRM_SECTIONS}
    for path in files:
        data = parsePRM(path)
        if not merged['header']:
            merged['header'] = data['header']
        for section in PRM_SECTIONS:
            for line in data.get(section, []):
                tokens = _strip_comment(line).split()
                if not tokens:
                    merged[section].append(line)
                    continue
                key = tuple(tokens[:4 if section != 'NONBONDED' else 1])
                if key in seen[section]:
                    continue
                seen[section].add(key)
                merged[section].append(line)
    writePRM(merged, out_path)
    return out_path


def rewriteRTFCharges(rtf_path: str, charges: list,
                       atom_order: list = None,
                       target_total: float = None,
                       out_path: str = None) -> str:
    """Overwrite per-atom charges in a CGenFF-style RTF with RESP values.

    Walks every RESI/PRES block and replaces the 4th whitespace-separated
    field of each `ATOM name type charge ...` line. Non-ATOM lines are
    preserved verbatim so comments, IMPR/BOND records and section order
    are untouched.

    Parameters
    ----------
    rtf_path     : input RTF file (typically the CGenFF `splitSTR` output).
    charges      : list of new atom charges in the same order as the ATOM
                   records appear in the file. If `atom_order` is supplied,
                   `charges` is treated as a {name: q} mapping instead (via
                   atom_order[i] -> charges[i]).
    atom_order   : optional list of atom names aligning the XYZ order with
                   the RTF ATOM-record order; if provided, enables per-name
                   lookup and tolerates mismatched ordering.
    target_total : if set, round each charge to 4 decimals and absorb the
                   rounding drift onto the atom with the largest |q| so the
                   final sum equals target_total exactly (to 4 dp).
    out_path     : where to write the result. Defaults to overwriting
                   `rtf_path` in place.
    """
    out_path = out_path or rtf_path
    if atom_order is not None:
        name_to_charge = {name: float(q) for name, q in zip(atom_order, charges)}
    else:
        name_to_charge = None
        charges = [float(q) for q in charges]

    with open(rtf_path, 'r') as f:
        lines = f.readlines()

    atom_hits = []  # list of (line_idx, resolved_charge)
    in_resi = False
    charge_iter = 0
    for idx, line in enumerate(lines):
        stripped = line.lstrip()
        upper = stripped.upper()
        if upper.startswith('RESI') or upper.startswith('PRES'):
            in_resi = True
            continue
        if upper.startswith('END'):
            in_resi = False
            continue
        if not in_resi:
            continue
        if not upper.startswith('ATOM'):
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        name = parts[1]
        if name_to_charge is not None:
            if name not in name_to_charge:
                continue
            q = name_to_charge[name]
        else:
            if charge_iter >= len(charges):
                continue
            q = charges[charge_iter]
            charge_iter += 1
        atom_hits.append((idx, q))

    if target_total is not None and atom_hits:
        rounded = [(i, round(q, 4)) for i, q in atom_hits]
        drift = round(target_total - sum(q for _, q in rounded), 4)
        if abs(drift) >= 1e-5:
            k = max(range(len(rounded)), key=lambda j: abs(rounded[j][1]))
            rounded[k] = (rounded[k][0], round(rounded[k][1] + drift, 4))
        atom_hits = rounded

    for line_idx, q in atom_hits:
        raw = lines[line_idx]
        leading = raw[:len(raw) - len(raw.lstrip())]
        parts = raw.strip().split()
        parts[3] = '{:.4f}'.format(q)
        lines[line_idx] = leading + ' '.join(parts) + '\n'

    with open(out_path, 'w') as f:
        f.writelines(lines)
    return out_path


def readRTFAtoms(rtf_path: str) -> list:
    """Return [(atom_name, atom_type, charge), ...] in RTF order.

    Walks every RESI/PRES block and collects each `ATOM name type charge`
    line. Comments and BOND/IMPR/etc. records are ignored.
    """
    atoms = []
    in_resi = False
    with open(rtf_path, 'r') as f:
        for line in f:
            stripped = line.lstrip()
            upper = stripped.upper()
            if upper.startswith('RESI') or upper.startswith('PRES'):
                in_resi = True
                continue
            if upper.startswith('END'):
                in_resi = False
                continue
            if not in_resi or not upper.startswith('ATOM'):
                continue
            parts = stripped.split()
            if len(parts) >= 4:
                try:
                    atoms.append((parts[1], parts[2], float(parts[3])))
                except ValueError:
                    continue
    return atoms


def readNonbondedLJ(prm_path: str) -> dict:
    """Return {atom_type: (epsilon, Rmin/2)} parsed from NONBONDED section.

    CGenFF / CHARMM NONBONDED line format:
        type  ignored  eps  Rmin/2  [1-4 cols]
    1-4 entries (when present) are dropped; only the 1-2 LJ pair is kept.
    """
    data = parsePRM(prm_path)
    out = {}
    for line in data.get('NONBONDED', []):
        tokens = _strip_comment(line).split()
        if len(tokens) < 4:
            continue
        try:
            eps = float(tokens[2])
            rmin2 = float(tokens[3])
        except ValueError:
            continue
        out[tokens[0]] = (eps, rmin2)
    return out


def xyzToCharmmPDB(xyz_path: str, pdb_path: str, resname: str = 'LIG',
                   segid: str = 'LIG') -> str:
    """Convert a PyConSolv XYZ file to a CHARMM-compatible PDB.

    CHARMM truncates atom names to 4 characters and expects a segid in
    columns 73-76. Each atom gets a unique name (element + index).
    """
    from .inputparser import XYZ
    db_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'db')
    xyz = XYZ(os.path.join(db_dir, 'atom-radius.txt'),
              os.path.join(db_dir, 'metal-radius.txt'))
    xyz.readXYZ(xyz_path)

    with open(pdb_path, 'w') as f:
        f.write('REMARK   generated by PyConSolv\n')
        for i, (elem, x, y, z) in enumerate(zip(xyz.atoms,
                                                 xyz.coords[:, 0],
                                                 xyz.coords[:, 1],
                                                 xyz.coords[:, 2]), start=1):
            name = '{}{}'.format(elem, i)[:4]
            f.write('ATOM  {:>5d} {:<4s} {:<3s} {:>4d}    '
                    '{:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00      {:<4s}{:>2s}\n'.format(
                        i, name, resname[:3], 1, x, y, z, segid[:4], elem))
        f.write('END\n')
    return pdb_path
