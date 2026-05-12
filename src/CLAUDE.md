# CLAUDE.md - PyConSolv Project Documentation

## Project Overview

PyConSolv is a Python package for automated parametrization and solvation of metal complexes for molecular dynamics simulations. It provides tools for:

- ORCA-based geometry optimization and electronic structure calculations
- MCPB.py-based metal force field parametrization
- GAFF parametrization for organic ligands via antechamber
- Automated solvation box setup with AMBER/tleap
- Equilibration and production MD preparation

## Project Structure

```
src/PyConSolv/
├── pyconsolv.py          # CLI entry point
├── ConfGen.py            # Main PyConSolv class (full workflow)
├── interfaces/
│   ├── amber.py          # AMBER/MCPB.py interface
│   └── orca.py           # ORCA interface
├── misc/
│   ├── inputparser.py    # XYZ file parsing, bond detection
│   ├── fragmenting.py    # Substructure extraction (Fragmentor class)
│   ├── Task.py           # Fragment workflow orchestration
│   └── ui.py             # GUI components (charge assignment, fragment review)
└── db/
    ├── atom-radius.txt   # Covalent radii for organic atoms
    └── metal-radius.txt  # Ionic radii for metal atoms
```

## Key Classes

### PyConSolv (ConfGen.py)
Main orchestrator for standard (non-fragment) parametrization workflow.
- `run()`: Full workflow - ORCA optimization, MCPB parametrization, solvation, equilibration
- `setup()`: Configure ORCA input with method, basis, dispersion, solvent
- `orca()`: Run ORCA geometry optimization
- `antechamber()`: Generate AMBER parameters via MCPB.py

### Fragmentor (misc/fragmenting.py)
Extracts substructures around metal centers for cheaper QM calculations.
- `findMetalIndices()`: Locate all metal atoms in structure
- `checkRadius()`: Find atoms within radius of ANY metal
- `checkBreakPoint()`: Ensure clean breaks, preserve rings
- `cutStructure()`: Generate XYZ with capping hydrogens at proper bond distances
- `prepareXYZ()`: Initialize metal-ligand bond connectivity

### Task (misc/Task.py)
Orchestrates the fragment-based parametrization workflow:
1. Extract fragment around metal(s)
2. GUI review of extracted substructure
3. Parametrize fragment (ORCA + MCPB) for metal parameters
4. Parametrize full ligands (GAFF) for organic parameters
5. Combine parameters
6. Run tleap for solvation
7. Equilibration
8. Prepare production files

### XYZ (misc/inputparser.py)
Parses XYZ files and detects bonding using covalent/ionic radii.
- `readXYZ()`: Parse XMOL format
- `calculateDistanceMatrix()`: Compute interatomic distances
- `generateAdjacencyMatrix()`: Determine bonding from radii
- `generateLinkList()`: Create connectivity lists
- `metalBonds`: List of metal-ligand bonds in format "metal_idx @ElementLigand_idx ligand_idx"

## CLI Usage

### Standard Parametrization
```bash
pyconsolv input.xyz -c 0 -m PBE0 -b def2-SVP -d D4 -s Water
```

### Fragment Mode (for large complexes)
```bash
pyconsolv input.xyz -c 0 -f -r 4.0 -mem 4000
```

### Key Arguments
- `-c, --charge`: System charge
- `-m, --method`: DFT method (default: PBE0)
- `-b, --basis`: Basis set (default: def2-SVP)
- `-d, --dispersion`: Dispersion correction (D3, D4, or N for none)
- `-s, --solvent`: Implicit solvent for ORCA
- `-f, --fragment`: Enable fragment mode
- `-r, --radius`: Extraction radius around metal (default: 4.0 Å)
- `-mem, --memory`: Memory per core in MB

## Bond Detection

Bond detection uses database files with atomic radii:
- `db/atom-radius.txt`: Covalent radii for C, N, O, S, P, H, etc.
- `db/metal-radius.txt`: Ionic radii for transition metals

A bond exists if: `distance < (radius_A + radius_B) * 1.15`

Metal-ligand bonds are stored separately in `xyz.metalBonds` for special handling in MCPB.py input generation.

## Multi-Metal Support

The code supports structures with multiple metal centers:

### Detection
- `Fragmentor.findMetalIndices()` returns a list of all metal atom indices
- `checkRadius()` includes atoms within radius of ANY metal in the structure

### Fragment Extraction
When extracting substructures:
- All metals and their coordination spheres are included
- Rings are preserved if any ring atom is within radius of any metal
- Capping hydrogens are placed at 1.09 Å from parent atoms

### MCPB.py Input
- `amber.inputFileGenerator()` accepts a list of metals
- Generates proper `ion_ids` and `ion_mol2files` entries for all metals

## Metal-Metal Bond Handling

For structures with direct metal-metal bonds (e.g., bimetallic complexes):

### Bond Detection (inputparser.py)
Metal-metal bonds are detected when both atoms are metals:
```python
if self.isMetal(a[0]) and self.isMetal(a[1]):
    # Use sum of both metal ionic radii
    dist = self.metalRadius.get(a[0].upper()) + self.metalRadius.get(a[1].upper())
```

### Duplicate Prevention
A `recorded_metal_bonds` set prevents the same M-M bond from being recorded twice:
```python
bond_key = tuple(sorted([i, j]))
if bond_key not in recorded_metal_bonds:
    self.metalBonds.append(bond_info)
    recorded_metal_bonds.add(bond_key)
```

### Parametrization Workflow
1. Both metals are included in fragment extraction
2. MCPB.py handles the M-M bond as a special case
3. Force field parameters generated for:
   - Individual metal coordination environments
   - M-M stretching potential
   - M-M-L angle bending terms
4. RESP charges calculated for entire fragment including both metals

### MCPB.py Input Generation
For multi-metal systems, the input file contains:
- Multiple entries in `ion_ids` (e.g., "1 2" for two metals)
- Multiple mol2 files in `ion_mol2files` (e.g., "Fe.mol2 Cu.mol2")

## Fragment Workflow Details

### Step 1: Fragment Extraction
- Fragmentor identifies all metals
- Atoms within specified radius of any metal are kept
- Rings are preserved to avoid cutting aromatic systems
- Capping hydrogens replace cut bonds at proper distances

### Step 2: GUI Review
- FragmentReviewGUI displays 2D structure
- Uses XYZ class for bond detection (not RDKit's buggy rdDetermineBonds)
- User can confirm or cancel

### Step 3: Fragment Parametrization
- ORCA optimization of fragment
- MCPB.py generates metal parameters
- RESP charges from fragment QM

### Step 4: Full Ligand Parametrization
- GAFF parameters via antechamber
- Handles organic portions of ligands

### Step 5: Parameter Combination
- Merge metal parameters from MCPB.py
- Merge organic parameters from GAFF
- Handle parameter conflicts

### Steps 6-8: MD Setup
- tleap solvation with specified box size
- Equilibration protocol
- Production file preparation

## Common Issues

### Dispersion Bug (Fixed)
The `-d N` flag (no dispersion) previously caused undefined variable. Fixed by adding explicit else clause.

### Capping Hydrogen Positions
Capping H atoms are repositioned to 1.09 Å from parent atom to ensure proper bond detection and QM convergence.

### RDKit Bond Detection
RDKit's automatic bond detection is unreliable for metal complexes. Always use the XYZ class with database radii files.

## Development Notes

- Branch `47-substructure-parametrization`: Fragment mode development
- Branch `1.0.7`: Stable release with memory parameter and fixes
- Use pathlib for all path operations
- Metal detection uses `xyz.isMetal()` method

---

# CHARMM36 Force Field Integration

## Status: Phase 1 + native FFTK bonded derivation complete

**Phase 1 is implemented and tested.** VMD was removed from the plan entirely; solvation is now done with Packmol + ParmEd. PSF generation uses `parmed.charmm.CharmmPsfFile` instead of `psfgen`. **easyPARM has also been removed**: metal-center and high-penalty bonded parameters are derived natively from the ORCA Hessian via Seminario (1996) projection in `interfaces/fftk.py`, emitting pure CHARMM36 RTF/PRM — no conversion layer, no external parametrization binary.

The FFTK-style pipeline is implemented in phases (see "FFTK-native parametrization" section further down):
- **Phase A — bonded (Seminario):** ✅ complete. Bond Kb/b0 and angle Kθ/θ0 derived from QM Hessian, scope = metal-touching + CGenFF-penalty > 50. Ground-truth diatomic recovery test in `tests/test_fftk.py`.
- **Phase B — charges (RESP via MultiWfn):** ✅ complete. Rather than reimplementing FFTK's water-interaction protocol, we reuse the existing MultiWfn RESP infrastructure from the AMBER path. After CGenFF produces the RTF, `ConfGen._charmmApplyRESP()` runs MultiWfn on the ORCA wavefunction (`freq` when a metal is present, `opt` otherwise — mirrors the AMBER path), reads the `.molden.chg` file, and `charmm_formats.rewriteRTFCharges()` overwrites the ATOM records by atom name. Rounding drift is absorbed onto the largest-|q| atom so the total matches the requested molecular charge exactly to 4 dp. Caveat: mild calibration mismatch with CHARMM's TIP3P (CGenFF was tuned against water-interaction data, RESP targets gas-phase ESP); acceptable tradeoff for metal complexes where CGenFF penalties around the metal are high anyway.
- **Phase C — dihedrals (QM PES scan fit):** ✅ skipped by design, matching the AMBER pipeline. CGenFF provides ligand dihedrals (parallel to GAFF's role in the AMBER path), and metal-crossing dihedrals are already emitted as soft wildcard torsions `X M<Sym> X X 0.0 1 0.0` by `fftk.py` (parallel to MCPB.py zeroing metal dihedrals). No custom per-torsion PES-scan fitting is performed.
- **Phase D — metal LJ table expansion:** ⏳ starter table shipped in `db/charmm_metal_nonbonded.txt` (18 entries).

Remaining work is in Phase 5 (integration runs against real CGenFF/Packmol binaries and validation vs NAMD/OpenMM) and Phase 6 (user-facing docs).

## Overview

Implement native CHARMM36 force field support alongside existing AMBER workflow. All tools run locally (no server dependencies).

## Target Workflow

```
Input XYZ
    ↓
ORCA optimization (existing)
    ↓
antechamber (AM1-BCC mol2)
    ↓
CGenFF (local) → Ligand RTF/PRM (via STR split)
    ↓
easyPARM (ORCA Hessian) → Metal RTF/PRM
    ↓
charmm_formats.merge{RTF,PRM} → combined params
    ↓
ParmEd CharmmPsfFile → PSF (metal-ligand bonds added manually)
    ↓
Packmol → Solvated box; ParmEd attachSolvatedBox merges coords+topology
    ↓
Output: PSF, PRM, PDB (ready for NAMD / OpenMM)
```

## Dependencies to Install

### Required External Tools
- [x] **CGenFF program** (local installation, requires license from SilcsBio)
- [x] **Packmol** binary (replaces VMD for solvation)
- [x] **ParmEd** (already used; now also handles PSF generation, replacing psfgen)
- [x] ~~VMD~~ – removed from dependencies
- [x] ~~easyPARM~~ – removed; metal bonded params derived natively via Seminario in `interfaces/fftk.py`

### Python Dependencies
- NumPy (already required) — used for the Hessian linear algebra
- Packmol is a standalone binary, not a Python package. MDAnalysis not required – ParmEd covers coordinate handling.

---

## Implementation Tasks

### Phase 1: Core Infrastructure — ✅ COMPLETE

#### Task 1.1: Create CHARMM interface module — ✅
- [x] Created `src/PyConSolv/interfaces/charmm.py`
- [ ] Implement `CharmmInterface` class (parallel to `amberInterface`)
- [ ] Methods needed:
  ```python
  class CharmmInterface:
      def __init__(self, path: str)
      def checkDependencies(self) -> bool  # Check VMD, CGenFF, easyPARM
      def generateRTF(self, mol2_file: str) -> str  # Topology file
      def generatePRM(self, mol2_file: str) -> str  # Parameter file
      def runPsfgen(self, rtf_files: list, pdb_file: str) -> str  # Generate PSF
      def solvate(self, psf: str, pdb: str, box_size: float) -> tuple
      def addIons(self, psf: str, pdb: str, concentration: float) -> tuple
  ```

#### Task 1.2: Create CGenFF interface — ✅
- [x] Created `src/PyConSolv/interfaces/cgenff.py`
- [ ] Implement local CGenFF execution
- [ ] Parse STR output files (topology + parameters)
- [ ] Extract penalty scores for parameter quality assessment
- [ ] Methods:
  ```python
  class CGenFFInterface:
      def __init__(self, cgenff_path: str = None)
      def parametrize(self, mol2_file: str, output_dir: str) -> str  # Returns STR file
      def parsePenalties(self, str_file: str) -> dict  # {atom: penalty_score}
      def splitSTR(self, str_file: str) -> tuple  # Returns (rtf_content, prm_content)
  ```

#### Task 1.3: Native Seminario/FFTK metal parametrization — ✅ (replaced easyPARM)
- [x] Created `src/PyConSolv/interfaces/fftk.py`
- [x] `parseOrcaHessian()` reads sparse column-major `$hessian` into a dense NumPy matrix
- [x] `seminarioBond(H, coords, i, j)` returns CHARMM-convention (Kb, b0) — validated against a synthetic diatomic in `tests/test_fftk.py`
- [x] `seminarioAngle(H, coords, i, j, k)` returns (Kθ, θ0); linear-angle degeneracy handled
- [x] `FFTKBondedInterface.deriveParameters()` writes a CHARMM RTF (MASS + minimal metal RESI) + PRM (BONDS + ANGLES + soft-metal DIHEDRALS + metal NONBONDED) fragment
- [x] Scope = option (c): metal-touching + CGenFF-penalty > configurable threshold
- [x] Metal atom types follow `M<Sym>` convention (MFE, MCU, MZN, MPT, MW, …)
- [x] Metal LJ from `db/charmm_metal_nonbonded.txt` (Won/Roux/Babu starter values)
- [ ] Wrapper for easyPARM Python package
- [ ] Convert ORCA Hessian to easyPARM format
- [ ] Handle multi-metal systems
- [ ] Methods:
  ```python
  class EasyPARMInterface:
      def __init__(self, path: str)
      def convertOrcaHessian(self, hess_file: str) -> str  # Convert to Cartesian Hessian
      def parametrizeMetal(self, xyz: str, hessian: str, charges: list) -> tuple
      def generateCharmmFiles(self, output_dir: str) -> tuple  # (rtf, prm)
  ```

#### Task 1.4: PSF generation — ✅ (replaced psfgen with ParmEd)
- [x] Created `src/PyConSolv/interfaces/charmm_builder.py`
- Uses `parmed.charmm.CharmmParameterSet` + `CharmmPsfFile`
- Metal-ligand bonds appended manually as `pmd.Bond` entries (no TCL patches needed)
- `attachSolvatedBox()` re-loads the Packmol output and merges coordinates + topology into the PSF

#### Task 1.5: Solvation — ✅ (replaced VMD with Packmol)
- [x] Created `src/PyConSolv/interfaces/packmol.py`
- `SOLVENT_PROPERTIES` dict with density + MW for water/TIP3P, acetonitrile, methanol, DMSO, DMF, chloroform, …
- Count math: `density_g_per_ml * volume_A3 * 1e-24 / mw * AVOGADRO`
- Box dimensions parsed from PDB CRYST1 (columns 30:38, 38:46, 46:54)
- Template PDBs: `db/water_tip3p.pdb`, `db/ion_sod.pdb`, `db/ion_cla.pdb`
- User-supplied custom solvent PDBs accepted as long as density + MW are supplied

---

### Phase 2: File Format Handling — ✅ COMPLETE

#### Task 2.1: CHARMM file parsers/writers — ✅
- [x] `src/PyConSolv/misc/charmm_formats.py`
- [x] RTF parser/writer (MASS, RESI, PRES sections)
- [x] PRM parser/writer (BONDS, ANGLES, DIHEDRALS, IMPROPERS, NONBONDED, NBFIX, CMAP, HBOND)
- [x] PSF handled via ParmEd (`CharmmPsfFile`) – no custom parser needed
- [x] STR split handled in `cgenff.py` (`read rtf … end` / `read para … end` blocks)

#### Task 2.2: Coordinate format conversions — ✅
- [x] `charmm_formats.xyzToCharmmPDB()` – truncates names to 4 chars, writes segid in columns 73-76, unique atom names from element + index
- [x] Residue naming parameterized via `resname` / `segid` arguments

#### Task 2.3: Atom type mapping — ✅
- [x] `src/PyConSolv/db/charmm_atomtypes.txt`

---

### Phase 3: Workflow Integration — ✅ (standard mode); ⏳ fragment mode pending

Note: the CLI flag landed as `-ff/--forcefield`, NOT `-e/--engine`. The existing `-e/--engine` flag is reserved for MD engine selection (amber/gromacs) and is independent of the force field.

#### Task 3.1: Modify ConfGen.py — ✅
- [x] Added `runCharmm()` method parallel to existing AMBER workflow
- [x] ORCA / charge / antechamber steps reused from the AMBER path

#### Task 3.2: Modify pyconsolv.py CLI — ✅
- [x] Added `-ff, --forcefield` argument: choices=['amber', 'charmm'], default='amber'
- [x] Dispatches `conf.runCharmm()` vs `conf.run()` based on `args.forcefield`

#### Task 3.3: Create CHARMM workflow orchestrator — ✅
- [x] `ConfGen.runCharmm()` orchestrates ORCA → CGenFF → easyPARM → merge → PSF → Packmol
- [ ] Workflow steps:
  1. ORCA optimization (existing)
  2. RESP charge calculation (existing via Multiwfn)
  3. CGenFF for organic ligands
  4. easyPARM for metal centers
  5. Merge RTF/PRM files
  6. psfgen to build PSF
  7. VMD solvate
  8. VMD ionize

#### Task 3.4: Fragment mode CHARMM support — ✅
- [x] `Task.fragment(..., forcefield='charmm')` plumbed through from the CLI (`-f -ff charmm`).
- [x] Fragment ORCA call appends `FREQ` when `forcefield=='charmm'` so a Hessian is available for FFTK.
- [x] New `ConfGen.runCharmmFromFragment()` runs ORCA opt on the full structure, CGenFF + RESP on the full ligand, and FFTK/Seminario on the **fragment's** Hessian (atoms/bonds indexed to the fragment). Type-based bonded terms (M<Sym> + CGenFF types) merge cleanly into the full-structure PRM.
- [x] Metal-bonded parameters derived once from the cheaper fragment QM; full-structure freq is avoided.

---

### Phase 4: Parameter Merging — ✅ baseline done

#### Task 4.1: RTF merger — ✅
- [x] `charmm_formats.mergeRTF` deduplicates MASS by (idx, type, mass) and preserves RESI/PRES

#### Task 4.2: PRM merger — ✅
- [x] `charmm_formats.mergePRM` merges BONDS/ANGLES/DIHEDRALS/IMPROPERS/NONBONDED/NBFIX/CMAP/HBOND
- [x] Duplicate key prevention (first-wins)

#### Task 4.3: Cross-term handling — ✅ minimal
- [x] Metal-ligand bonds added directly into PSF via `charmm_builder.addMetalBonds()`
- [ ] Pending: richer conflict reporting when CGenFF and easyPARM both emit the same atom type

---

### Phase 5: Testing & Validation — ✅ unit tests; ⏳ integration/validation pending

#### Task 5.1: Unit tests — ✅
- [x] `tests/test_cgenff.py`, `tests/test_easyparm.py`, `tests/test_packmol.py`, `tests/test_charmm_interface.py`, `tests/test_charmm_formats.py`
- [x] 82 unit tests across 14 modules; subprocess calls mocked via `unittest.mock`
- [x] `tests/helpers.py` provides `TempDir` (restores cwd, since `easyparm.__init__` does `os.chdir`), XYZ fixtures, and `radii_files()`
- [x] `tests/__init__.py` injects `src/` into `sys.path`; legacy tests updated to `import tests` shim

#### Task 5.2: Integration tests — ⏳ pending
- [ ] End-to-end run on a simple organic molecule (CGenFF only, real binaries)
- [ ] Real metal complex (CGenFF + easyPARM)
- [ ] Multi-metal complex
- [ ] Fragment mode + CHARMM

#### Task 5.3: Validation — ⏳ pending
- [ ] Compare energies: AMBER vs CHARMM for the same structure
- [ ] Verify solvation box geometry (density, ion concentrations)
- [ ] Round-trip test through NAMD / OpenMM

---

### Phase 6: Documentation

#### Task 6.1: Update CLAUDE.md
- [ ] Document CHARMM workflow
- [ ] Document new CLI options
- [ ] Document file formats

#### Task 6.2: User documentation
- [ ] Installation guide for CHARMM dependencies
- [ ] Tutorial: CHARMM parametrization
- [ ] Troubleshooting guide

---

## File Structure After Implementation

```
src/PyConSolv/
├── pyconsolv.py              # CLI (add -e/--engine)
├── ConfGen.py                # Add runCharmm() workflow
├── interfaces/
│   ├── amber.py              # Existing
│   ├── orca.py               # Existing
│   ├── charmm.py             # CHARMM orchestrator (CGenFF + FFTK Seminario + Packmol)
│   ├── cgenff.py             # CGenFF wrapper – STR split into RTF/PRM, penalty + atom-type parsing
│   ├── fftk.py               # Native Seminario bonded derivation, pure CHARMM36 output
│   ├── packmol.py            # Packmol solvent box builder (replaces VMD solvate)
│   └── charmm_builder.py     # ParmEd-based PSF assembly + attachSolvatedBox
├── misc/
│   ├── inputparser.py        # Existing
│   ├── fragmenting.py        # Existing
│   ├── Task.py               # TODO: forcefield-aware fragment workflow
│   ├── ui.py                 # Existing
│   └── charmm_formats.py     # RTF/PRM parse/write/merge + xyzToCharmmPDB
└── db/
    ├── atom-radius.txt       # Existing
    ├── metal-radius.txt      # Existing
    ├── charmm_atomtypes.txt  # Element → CGenFF / metal fallback type mapping
    ├── charmm_metal_nonbonded.txt  # Element → (epsilon, Rmin/2) for FFTK metal LJ
    ├── water_tip3p.pdb       # Packmol water template
    ├── ion_sod.pdb           # Packmol Na+ template
    └── ion_cla.pdb           # Packmol Cl- template
```

Notes:
- `psfgen.py` and `vmd_solvate.py` were planned but never created — PSF generation is handled by `charmm_builder.py` (ParmEd), solvation by `packmol.py`.
- `easyparm.py` was implemented and subsequently removed — replaced by native Seminario/FFTK derivation in `fftk.py`, keeping all output in pure CHARMM36 format with no conversion layer.

---

## CLI (as implemented)

```bash
# AMBER workflow (default)
pyconsolv input.xyz -c 0 -m PBE0 -b def2-SVP

# CHARMM36 workflow (CGenFF + easyPARM + Packmol)
pyconsolv input.xyz -c 0 -m PBE0 -b def2-SVP -ff charmm

# Fragment mode with CHARMM (pending – Task.py not yet forcefield-aware)
pyconsolv input.xyz -c 0 -f -r 4.0 -ff charmm
```

`-e/--engine` remains reserved for MD-engine selection (amber/gromacs) and is independent of `-ff/--forcefield`.

---

## Key Technical Decisions

### CGenFF Local Execution
CGenFF requires a license from SilcsBio. The program reads MOL2 files and outputs STR (stream) files containing both topology and parameters.

Command: `cgenff mol2_file.mol2 -o output.str`

### easyPARM Integration
easyPARM is Python-based and can be imported directly:
```python
from easyparm import EasyParm
ep = EasyParm(xyz_file, hessian_file, charges)
ep.run()
ep.write_charmm('output')  # Writes RTF and PRM
```

### PSF generation (ParmEd, replaces psfgen)
```python
import parmed as pmd
from parmed.charmm import CharmmParameterSet, CharmmPsfFile

params = CharmmParameterSet(merged_rtf, merged_prm)
psf = CharmmPsfFile.from_structure(structure)
psf.load_parameters(params)
# Metal-ligand bonds added manually:
psf.bonds.append(pmd.Bond(psf.atoms[i], psf.atoms[j]))
psf.write_psf('output.psf')
```

### Solvation (Packmol, replaces VMD solvate)
```
tolerance 2.0
filetype pdb
output solvated.pdb

structure solute.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
end structure

structure water_tip3p.pdb
  number <computed from density * volume / mw * N_A>
  inside box 0. 0. 0. Lx Ly Lz
end structure
```
Ion counts derived from target concentration; `attachSolvatedBox()` re-loads the Packmol output and merges coordinates + topology into the PSF.

---

## Remaining Priority Order

1. **FFTK Phase D — metal LJ table expansion** — ✅ broad superset shipped (alkali, alkaline earth, full 3d/4d/5d TM rows, main-group post-transition, lanthanides, Ac/Th/U).
4. **Phase 3.4**: Fragment mode CHARMM support in `Task.py`.
5. **Phase 5.2–5.3**: Integration runs against real CGenFF/Packmol binaries and validation vs NAMD/OpenMM.
6. **Phase 6**: User-facing tutorial and troubleshooting docs.

---

## Notes

- easyPARM v4.0 supports multiple metals in same structure
- CGenFF penalty scores > 50 indicate parameters need QM refinement
- CHARMM uses all-atom representation (explicit hydrogens everywhere)
- Water models: TIP3P (CHARMM-modified), TIP4P, SPC/E
- Ion parameters from CHARMM36 ion topology files
