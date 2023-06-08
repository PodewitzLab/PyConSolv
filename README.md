
# PyConSolv

A python based interface for generation of conformers of transition metal complexes in explicit solvent.
The interface bridges to the well known MCPB.py package available within ambertools. The input required 
consists of only a simple xyz file and all required steps for parametrization are performed automatically,
with minimal user intervention.



## Features
Utilizes freely available software, with high performance

18 predefined solvents and 6 counterions, with the ability to use any solvent or counterion

Automated molecule splitting for transition metal parametrization

Utilizes ORCA 5.0 for quantum mechanical optimizations/frequency calculations

Utilizes MultiWfn for the generation of the RESP charges

Automated equilibration of simulation box

Automated clustering


## Requirements

Python >=3.10

AmberTools 20+

ORCA 5.0+

MultiWfn 3.8+

## Installation

The creation of a new virtual environment is highly recommended:

using conda:
```
conda create -c conda-forge --name PyConSolv python=3.10 rdkit numpy pandas
conda activate PyConSolv
pip install PyConSolv
```

using pip:
```
python3 -m venv env
source env/bin/activate
pip install numpy pandas rdkit PyConSolv
```

## Usage

### Console:
```
pyconsolv [-h] [-c [CHARGE]] [-m [METHOD]] [-b [BASIS]] [-d [DISPERSION]] [-s [SOLVENT]] [-p [CPU]] [-mult [MULTIPLICITY]] [-a [ANALYZE]] [-mask [MASK]] [-cluster [CLUSTER]] [-nosp] [-v] input
```

positional arguments:  
input input file in XYZ format

options:  
  -h, --help            show this help message and exit  
  -c [CHARGE], --charge [CHARGE] charge of the system, default 0  
  -m [METHOD], --method [METHOD] ORCA optimization/frequency calculations method of choice, default PBE0  
  -b [BASIS], --basis [BASIS] basis set to be used for calculations, default def2-SVP  
  -d [DISPERSION], --dispersion [DISPERSION] dispersion corrections, default = D4  
  -s [SOLVENT], --solvent [SOLVENT] solvent to be used for MD simulations/ OM Calculations, default Water  
  -p [CPU], --cpu [CPU] number of cpu cores to be used for calculations, default 12  
  -mult [MULTIPLICITY], --multiplicity [MULTIPLICITY] multiplicity of the system, default 1  
  -nosp skip single point calculations for clusters
  -a , --analyze analyze a simulation  
  -mask [MASK], --mask [MASK] atomid mask for clustering  
  -cluster [CLUSTER], --cluster [CLUSTER] clustering method  
  -v, --version         show program's version number and exit  

see user manual for more details


### Jupyter Notebook

```
from PyConSolv import ConfGen

conf = ConfGen(path/to/input.xyz)

conf.run([options])
```



