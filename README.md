![License](https://img.shields.io/badge/license-GPL3-blue)
[![Latest Version](https://img.shields.io/badge/release-v.1.0.6.3-red)](https://pypi.org/project/PyConSolv/1.0.6.2/)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.3c00798-blue)](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00798)

# Latest:
Please update version 1.0.2 to 1.0.3.1+ to fix a bug with generating parameters. Versions before 1.0.2 should not be affected.

Changelog from v 1.0.1:

- More counterions are parametrized and supported by default
- Box size can now be specified using the `-box` parameter
- Support for multiple non-covalently bound solutes
- Support for QM/MM calculations as criteria for the energy ranking of generated conformers
- Support for restrained simulations involving transition states (Currently only for AmberMD)
- Support for cartesian coordinate restraints for residues, suitable for GIST analysis using `-cart` for the residues and `-cartstr` for the restraint strength

# PyConSolv

A python based interface for generation of conformers of transition metal complexes in explicit solvent.
The interface bridges to the well known MCPB.py package available within ambertools. The input required 
consists of only a simple xyz file and all required steps for parametrization are performed automatically,
with minimal user intervention.

Publication:

[PyConSolv: A Python Package for Conformer Generation of (Metal-Containing) Systems in Explicit Solvent
R. A. Talmazan and M. Podewitz
Journal of Chemical Information and Modeling 2023, 63, 17, 5400–5407
DOI: 10.1021/acs.jcim.3c00798](https://pubs.acs.org/doi/full/10.1021/acs.jcim.3c00798)



## Features
Utilizes freely available software, with high performance

18 predefined solvents and 6 counterions (along with 40+ single atom ions), with the ability to use any solvent or counterion

Automated molecule splitting for transition metal parametrization

Utilizes ORCA 5.0 for quantum mechanical optimizations/frequency calculations

Utilizes MultiWfn for the generation of the RESP charges

Automated equilibration of simulation box

Automated clustering


## Requirements

Python >= 3.10

AmberTools >= 20

ORCA >= 5.0

MultiWfn >= 3.8

## Installation

The creation of a new virtual environment is highly recommended:

using conda:
```
conda create -c conda-forge --name PyConSolv python=3.10 rdkit numpy pandas parmed
conda activate PyConSolv
pip install PyConSolv
```

using pip:
```
python3 -m venv env
source env/bin/activate
pip install numpy pandas rdkit parmed PyConSolv
```

## Usage

### Console:
```
pyconsolv [-h] [-c [CHARGE]] [-m [METHOD]] [-b [BASIS]] [-d [DISPERSION]] [-s [SOLVENT]] [-p [CPU]] [-mult [MULTIPLICITY]] [-noopt] [-a [ANALYZE]] [-mask [MASK]] [-cluster [CLUSTER]] [-nosp] [-v] input
```

**positional arguments:**  
input file in XYZ format

**options that affect simulation setup:**  
 
  -c [CHARGE], --charge [CHARGE] charge of the system, default 0  
  -m [METHOD], --method [METHOD] ORCA optimization/frequency calculations method of choice, default PBE0  
  -b [BASIS], --basis [BASIS] basis set to be used for calculations, default def2-SVP  
  -d [DISPERSION], --dispersion [DISPERSION] dispersion corrections, default = D4  
  -s [SOLVENT], --solvent [SOLVENT] solvent to be used for MD simulations/ OM Calculations, default Water  
  -p [CPU], --cpu [CPU] number of cpu cores to be used for calculations, default 12  
  -mult [MULTIPLICITY], --multiplicity [MULTIPLICITY] multiplicity of the system, default 1   
  -noopt perform a single point calculation instead of a geometry optimization  
  -box specify box size for your system  
  -e, --engine         choice of simulation engine  
  -rst, --restraint perform a restrained simulation, useful for transition states  
  
**options that affect analysis**:   
  -a , --analyze analyze a simulation  
  -mask [MASK], --mask [MASK] atomid mask for clustering  
  -cluster [CLUSTER], --cluster [CLUSTER] clustering method  
  -nosp skip single point calculations for clusters  
  -qmmm, --qmmm use a qmmm approach to determine cluster energy ranking


**general options:**
  -h, --help            show this help message and exit  
  -v, --version         show program's version number and exit  


see user manual for more details


### Jupyter Notebook

```
from PyConSolv import ConfGen

conf = ConfGen(path/to/input.xyz)

conf.run([options])
```



