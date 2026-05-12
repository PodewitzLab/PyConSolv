# PyConSolv - Claude Context

## Project Overview

PyConSolv is a Python tool for automated parametrization of metal-containing complexes for molecular dynamics simulations. It integrates quantum mechanical calculations (ORCA) with molecular mechanics tools (AmberTools, MCPB.py) to generate force field parameters for metal complexes.

## Key Directories

```
src/PyConSolv/
├── pyconsolv.py          # CLI entry point
├── ConfGen.py            # Main PyConSolv class with full workflow
├── db/                   # Database files for atomic radii
│   ├── atom-radius.txt   # Ionic/covalent radii for elements
│   └── metal-radius.txt  # Metal ion radii
├── interfaces/
│   ├── amber.py          # AmberTools interface (antechamber, MCPB.py, tleap)
│   ├── orca.py           # ORCA quantum chemistry interface
│   ├── parmed.py         # ParmEd interface for topology manipulation
│   ├── charmm.py         # CHARMM orchestrator (CGenFF + native FFTK + Packmol)
│   ├── cgenff.py         # CGenFF wrapper – organic ligand RTF/PRM via STR
│   ├── fftk.py           # Native Seminario/FFTK bonded derivation (replaces easyPARM)
│   ├── packmol.py        # Packmol solvent box builder (replaces VMD solvate)
│   └── charmm_builder.py # ParmEd-based PSF assembly + solvated box merge
├── misc/
│   ├── inputparser.py    # XYZ class for structure parsing and bond detection
│   ├── fragmenting.py    # Fragmentor class for substructure extraction
│   ├── Task.py           # Task orchestration for fragment workflow
│   ├── ui.py             # GUI classes for visualization
│   ├── charmm_formats.py # RTF/PRM parsers/writers + mergers, xyzToCharmmPDB
│   ├── parameterChecker.py
│   └── tleapAdderInterface.py
└── utils/
    └── colorgen.py       # Terminal color output
```

## Core Classes

### PyConSolv (ConfGen.py)
Main workflow class. Key methods:
- `run()` - Full AMBER parametrization workflow
- `runCharmm()` - Full CHARMM36 workflow (ORCA → CGenFF → native Seminario/FFTK → merge → PSF → Packmol solvate)
- `runCharmmFromFragment()` - Fragment-mode CHARMM. Takes a pre-computed fragment ORCA opt+freq and runs CGenFF/RESP on the full structure + FFTK on the fragment Hessian; metal-bonded terms transfer by atom type (M<Sym>)
- `setup()` - Initialize folders, parse XYZ, create ORCA input
- `orca()` - Run ORCA optimization/frequency calculations
- `antechamber()` - Generate mol2/frcmod with GAFF
- `multiwfn()` - Calculate RESP charges
- `MCPB_script()` - Run MCPB.py for metal parameters
- `tleap()` - Build solvated system
- `equilibration()` - Run MD equilibration
- `prepareSimulation()` - Prepare production MD files

### XYZ (inputparser.py)
Structure parsing and bond detection:
- `readXYZ()` - Parse XMOL format files
- `calculateDistanceMatrix()` - Compute atom pair distances
- `generateAdjacencyMatrix()` - Detect bonds using ionic radii from database
- `generateLinkList()` - Build connectivity list (excludes metal bonds)
- `metalBonds` - List of metal-ligand bonds (format: "metal_idx @ElementLigand_idx ligand_idx")
- Bond detection uses: `distance <= (radius1 + radius2) * 0.6`

### Fragmentor (fragmenting.py)
Extract substructure around metal(s) for cheaper QM calculations:
- `findMetalIndices()` - Find ALL metal atoms in structure (returns list)
- `checkRadius()` - Find atoms within radius of ANY metal
- `checkBreakPoint()` - Handle ring systems, identify cut points
- `cutStructure()` - Generate XYZ with capping hydrogens at proper bond distances
- `prepareXYZ()` - Populate linkList with metal-ligand bonds for all metals
- `hydrogenate_list` - Indices of atoms converted to capping H
- `metal_indices` - List of all metal atom indices

### Task (Task.py)
Orchestrates fragment-based parametrization:
1. Extract fragment around metal
2. Show GUI for user review
3. Parametrize fragment (ORCA → MCPB) for metal parameters
4. Parametrize full structure with GAFF
5. Combine parameters
6. Build system, equilibrate, prepare simulation

## CLI Usage

```bash
# Standard parametrization
pyconsolv input.xyz -c 0 -m PBE0 -b def2-SVP

# Fragment-based parametrization
pyconsolv input.xyz -c 0 -f -r 4.0

# CHARMM36 parametrization (CGenFF + easyPARM + Packmol, no VMD)
pyconsolv input.xyz -c 0 -ff charmm

# Key flags:
# -f, --fragment     Enable fragment mode
# -r, --radius       Radius around metal for fragment (default 4.0 Å)
# -c, --charge       System charge
# -m, --method       ORCA method (default PBE0)
# -b, --basis        Basis set (default def2-SVP)
# -d, --dispersion   Dispersion correction (default D4, N for none)
# -mem, --memory     Memory per core for ORCA
# -p, --cpu          Number of CPU cores
# -ff, --forcefield  amber (default) | charmm – selects parametrization pipeline
# -e, --engine       MD engine (amber | gromacs) – independent of -ff
```

## CHARMM36 Workflow (`runCharmm`)

Parallel pipeline that produces PSF/PRM/PDB for NAMD or OpenMM without requiring VMD:
1. ORCA optimization + frequency (shared with AMBER path)
2. `antechamber` generates an AM1-BCC MOL2 used as CGenFF input
3. `CGenFFInterface.parametrize()` runs the local CGenFF binary, splits the resulting STR into RTF + PRM, and exposes penalty scores via `parsePenalties()`
4. `MultiWfnInterface.run()` computes RESP charges from the ORCA wavefunction (same protocol as the AMBER path — `freq` wavefunction when a metal is present, `opt` otherwise). `ConfGen._charmmApplyRESP()` overwrites CGenFF charges in the ligand RTF by atom name via `charmm_formats.rewriteRTFCharges()`, absorbing rounding drift onto the largest-|q| atom so the total charge lands exactly on target.
5. `FFTKBondedInterface.deriveParameters()` parses the ORCA `$hessian` block (sparse column-major → dense) and applies the Seminario (1996) projection method to derive bond Kb/b0 and angle Kθ/θ0 directly in CHARMM36 units. Scope: any bond/angle that (a) touches a metal OR (b) involves an atom whose CGenFF penalty exceeds 50. Metal atom types follow the `M<sym>` convention (MFE, MCU, MZN, …); metal LJ parameters come from `db/charmm_metal_nonbonded.txt` (Won/Roux/Babu values). No external parametrization binary required.
6. `charmm_formats.mergeRTF` / `mergePRM` deduplicate MASS, BONDS, ANGLES, DIHEDRALS, IMPROPERS, NONBONDED entries
7. `charmm_builder` loads the merged parameters through `parmed.charmm.CharmmParameterSet`, builds a `CharmmPsfFile`, and manually appends metal-ligand bonds as `pmd.Bond` entries
8. `PackmolInterface` computes solvent counts from density/MW, emits a `.inp` script, runs the `packmol` binary, and `attachSolvatedBox()` re-loads and merges the box into the PSF

### Solvation (Packmol path, replaces VMD)
- `interfaces/packmol.py` owns `SOLVENT_PROPERTIES` (water/TIP3P, acetonitrile, methanol, DMSO, DMF, chloroform, …) with density + MW
- Count math: `density_g_per_ml * volume_A3 * 1e-24 / mw * AVOGADRO`
- Box dimensions parsed from PDB CRYST1 (columns 30:38, 38:46, 46:54)
- Template PDBs shipped in `db/`: `water_tip3p.pdb`, `ion_sod.pdb`, `ion_cla.pdb`
- User-supplied custom solvent PDBs are accepted as long as density + MW are known

## Important Implementation Details

### Bond Detection
- Uses database files with ionic radii (`db/atom-radius.txt`)
- Column mapping: `df.values[:, 1]` = element symbol, `df.values[:, 5]` = radius
- Metal bonds stored separately in `xyz.metalBonds`
- `linkList` excludes metal-containing bonds
- Metal-ligand bonds: `distance <= (metalRadius + ligandRadius) * 1.0`
- Organic bonds: `distance <= (radius1 + radius2) * 0.6`

### Multi-Metal & Metal-Metal Bond Support
Supports bimetallic/polymetallic complexes including direct M-M bonds:

**Bond Detection (`generateAdjacencyMatrix`)**:
- Detects M-M bonds using sum of metal radii
- Detects M-L bonds regardless of which atom is listed first
- Uses `recorded_metal_bonds` set to avoid duplicate entries
- Format in `metalBonds`: `"metal_idx @Element+ligand_idx ligand_idx"`

**Fragment Extraction**:
- `findMetalIndices()` returns list of ALL metal indices
- `checkRadius()` includes atoms within radius of ANY metal
- Both metals and their coordination spheres included in fragment

**MCPB.py Integration**:
- Multiple metals listed as separate ions: `ion_ids 1 2 3`
- Separate mol2 files for each metal: `ion_mol2files FE.mol2 CU.mol2`
- `checkMCPBBonds()` adds missing bonds via `add_bonded_pairs`

**Parametrization Workflow for M-M Bonds**:
1. Fragment extraction includes both metals
2. ORCA optimizes entire fragment including M-M bond
3. Frequency calculation computes M-M stretch force constants
4. MCPB.py derives M-M bond parameters from QM Hessian
5. Force field includes M-M bond/angle/dihedral parameters

**Limitations**:
- MCPB.py expects metals in separate residues
- Very short M-M bonds (metal clusters) may need manual adjustment
- Charge distribution between metals requires careful consideration

### Capping Hydrogens (fragmenting.py)
When cutting structure, atoms outside the radius that bond to kept atoms become capping H:
- Positioned at 1.09 Å from parent atom (typical X-H bond distance)
- Direction preserved from original heavy atom position

### GUI Visualization (ui.py)
- `FragmentReviewGUI` - Shows extracted fragment for user confirmation
- Uses XYZ class for bond detection, then adds bonds to RDKit molecule
- Metal bonds parsed from `xyz.metalBonds` (format: 3 space-separated parts)

## Restart Mechanism
Uses restart levels to allow resuming interrupted jobs:
- Level 0: Start
- Level 2: After ORCA
- Level 3: After antechamber
- Level 5: After MultiWfn
- Level 6: After MCPB
- Level 7: After tleap
- Level 8: After equilibration
- Level 9: Complete

## External Dependencies
- ORCA 5.x (quantum chemistry)
- AmberTools (antechamber, MCPB.py, tleap, pmemd)
- MultiWfn (wavefunction analysis, RESP charges)
- RDKit (structure visualization)
- NumPy, Pandas, Matplotlib, ParmEd
- **CHARMM mode only:** CGenFF (SilcsBio binary), Packmol (binary). VMD and easyPARM are NOT required — metal bonded parameters are derived natively via Seminario/FFTK.

## Test Infrastructure
- `tests/__init__.py` injects `src/` into `sys.path`; all test modules do `import tests` then import `PyConSolv...` directly
- `tests/helpers.py` provides `TempDir`, `SAMPLE_XYZ`/`SAMPLE_METAL_XYZ` fixtures, and `radii_files()`
- Unit tests cover colorgen, copier, ions, restart, charge, inputparser, fragmenting, amber interface (mocked subprocess), MD engines, CHARMM format mergers, CGenFF, Packmol, the FFTK Seminario back-end (including ground-truth diatomic force-constant recovery), and the CharmmInterface orchestrator
- Subprocess calls in interface tests are mocked via `unittest.mock` – no external binaries needed

## Recent Changes (Branch 1.0.7 + substructure integration)
- Added `-f` flag for fragment mode
- Added `-r` flag for extraction radius
- Fixed dispersion bug (`dsp` undefined when not 'N')
- Added `customOrcaInput` parameter to `setup()`
- Moved `startInfo()` from `__init__` to `run()`
- Added `FragmentReviewGUI` for fragment visualization
- Capping H repositioned to proper bond distance
- Fixed metal bond parsing in UI (was checking `>= 4` parts, should be `>= 3`)
- **Multi-metal support:**
  - `Fragmentor.findMetalIndices()` finds ALL metals
  - `checkRadius()` includes atoms within radius of ANY metal
  - `prepareXYZ()` handles bonds for all metals
  - `inputFileGenerator()` accepts list of metals, generates proper `ion_ids`
  - `ConfGen.antechamber()` passes all metal names to MCPB input
- **Metal-metal bond handling:**
  - `generateAdjacencyMatrix()` detects M-M bonds using metal radii
  - Added detection for M-L bonds when metal is second atom
  - Duplicate bond prevention using `recorded_metal_bonds` set
  - M-M bond parameters derived from QM Hessian via MCPB.py
- **CHARMM36 support (Phase 1 complete):**
  - New `-ff/--forcefield {amber,charmm}` flag (distinct from `-e/--engine`)
  - `ConfGen.runCharmm()` orchestrates the pipeline
  - Interfaces: `cgenff.py`, `fftk.py`, `charmm.py`, `charmm_builder.py`, `packmol.py`
  - VMD dependency removed – solvation via Packmol + ParmEd
  - easyPARM dependency removed – metal/high-penalty bonded terms derived natively via Seminario/FFTK in pure CHARMM36 format
  - `misc/charmm_formats.py` for RTF/PRM parse/write/merge and `xyzToCharmmPDB`
  - `db/charmm_atomtypes.txt`, `db/charmm_metal_nonbonded.txt`, `db/water_tip3p.pdb`, `db/ion_sod.pdb`, `db/ion_cla.pdb`
- **CHARMM36 native FFTK phases:**
  - **Phase A (bonded Seminario):** ✅ `fftk.py` parses ORCA `$hessian`, derives Kb/b0 and Kθ/θ0 directly in CHARMM36 units. Scope: metal-touching + CGenFF-penalty > 50.
  - **Phase B (RESP charges via MultiWfn):** ✅ `ConfGen._charmmApplyRESP()` overwrites CGenFF ATOM charges by name; rounding drift absorbed onto largest-|q| atom to hit target total exactly.
  - **Phase C (dihedral fitting):** skipped by design — CGenFF provides ligand dihedrals (parallel to GAFF's role in AMBER) and metal-crossing torsions are already emitted as `X M<Sym> X X 0.0 1 0.0` wildcards (parallel to MCPB.py zeroing).
  - **Phase D (metal LJ table):** ✅ `db/charmm_metal_nonbonded.txt` broad superset — alkali, alkaline earth, full 3d/4d/5d TM rows, main-group post-transition, lanthanides, Ac/Th/U (Roux/Babu-Lim/Won/Li-Merz).
- **Fragment-mode CHARMM (Task 3.4):** ✅ `Task.fragment(..., forcefield='charmm')` + `ConfGen.runCharmmFromFragment()`. Fragment ORCA input appends `FREQ`; the fragment Hessian feeds FFTK while the full structure only needs ORCA opt (CGenFF + RESP). Metal-bonded params transfer via type-based merge (M<Sym> atom types).
- **Test suite:** 51/51 CHARMM-path tests pass (`test_fftk`, `test_charmm_interface`, `test_charmm_formats`, `test_cgenff`, `test_packmol`). 8 pre-existing legacy failures in `test_inputparser`/`test_amber`/`test_filestructure` are due to a stale `tests/Testfiles/input.xyz` and predate this branch's CHARMM work.

## CHARMM roadmap — what's left
- **Phase 5.2 — Integration tests against real binaries.** End-to-end runs with actual CGenFF + Packmol on (a) simple organic, (b) mono-metal complex, (c) bimetallic, (d) fragment mode.
- **Phase 5.3 — Validation.** Round-trip PSF/PRM/PDB through NAMD and OpenMM; compare single-point energies vs the AMBER path on the same geometry; verify solvation density / ion counts.
- **Phase 6 — User docs.** Installation guide (CGenFF license, Packmol), `-ff charmm` tutorial, troubleshooting notes (high-penalty warnings, missing metal LJ entries, fragment extraction radius tuning).
