import os
import shutil
import subprocess

from ..utils.colorgen import Color


class CGenFFInterface:
    """Wrapper around the local CGenFF program for small-molecule CHARMM36
    parametrization. Produces CHARMM stream (.str) files that contain both
    a topology (RTF) and parameter (PRM) section.
    """

    def __init__(self, cgenff_path: str = None):
        self.cgenff_path = cgenff_path or shutil.which('cgenff') or 'cgenff'
        self.status = 0
        self.last_str = None

    def checkpath(self) -> bool:
        if shutil.which(self.cgenff_path) is None and not os.path.isfile(self.cgenff_path):
            print(Color.RED + 'CGenFF executable not found (set cgenff_path or add to PATH)' + Color.END)
            self.status = 0
            return False
        self.status = 1
        return True

    def parametrize(self, mol2_file: str, output_dir: str = None) -> str:
        """Run CGenFF on a mol2 file. Returns path to the resulting .str file."""
        output_dir = output_dir or os.path.dirname(os.path.abspath(mol2_file))
        base = os.path.splitext(os.path.basename(mol2_file))[0]
        str_file = os.path.join(output_dir, base + '.str')

        cmd = '{} {} -o {}'.format(self.cgenff_path, mol2_file, str_file)
        calc = subprocess.run([cmd], shell=True, cwd=output_dir)
        if calc.returncode != 0:
            print(Color.RED + 'CGenFF failed for {}'.format(mol2_file) + Color.END)
            self.status = 0
            return ''

        self.status = 1
        self.last_str = str_file
        print('CGenFF parametrization complete: {}'.format(str_file))
        return str_file

    def parsePenalties(self, str_file: str) -> dict:
        """Extract per-atom and parameter penalty scores from a .str file.

        CGenFF flags assigned atom types/params with a 'penalty' comment that
        indicates how reliable the guess is. Returns {identifier: score}.
        """
        penalties = {}
        if not os.path.isfile(str_file):
            return penalties
        with open(str_file, 'r') as f:
            for line in f:
                if 'penalty' not in line.lower():
                    continue
                parts = line.split()
                try:
                    idx = [p.lower() for p in parts].index('penalty=')
                except ValueError:
                    try:
                        idx = [p.lower() for p in parts].index('penalty')
                    except ValueError:
                        continue
                try:
                    score = float(parts[idx + 1].rstrip(','))
                except (IndexError, ValueError):
                    continue
                key = parts[1] if len(parts) > 1 else line.strip()
                penalties[key] = score
        return penalties

    def parseAtomTypes(self, str_file: str) -> dict:
        """Return {atom_name: atom_type} extracted from the RESI block.

        CGenFF writes one ATOM line per atom with fields:
            ATOM  <name>  <type>  <charge>
        """
        types = {}
        in_rtf = False
        with open(str_file, 'r') as f:
            for line in f:
                lower = line.strip().lower()
                if lower.startswith('read rtf'):
                    in_rtf = True
                    continue
                if lower == 'end':
                    in_rtf = False
                    continue
                if not in_rtf:
                    continue
                parts = line.split()
                if len(parts) >= 4 and parts[0].upper() == 'ATOM':
                    types[parts[1]] = parts[2]
        return types

    def splitSTR(self, str_file: str) -> tuple:
        """Split a CGenFF .str file into (rtf_text, prm_text).

        The stream file is organised as `read rtf ... end` and
        `read para ... end` blocks. We copy each block verbatim.
        """
        rtf_lines = []
        prm_lines = []
        section = None
        with open(str_file, 'r') as f:
            for line in f:
                stripped = line.strip().lower()
                if stripped.startswith('read rtf'):
                    section = 'rtf'
                    continue
                if stripped.startswith('read para'):
                    section = 'prm'
                    continue
                if stripped == 'end':
                    section = None
                    continue
                if section == 'rtf':
                    rtf_lines.append(line)
                elif section == 'prm':
                    prm_lines.append(line)
        return ''.join(rtf_lines), ''.join(prm_lines)
