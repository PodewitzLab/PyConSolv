
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

## Usage

### Console:
pyconsolv [options]

see user manual for console options


### Jupyter Notebook

python

from PyConSolv import ConfGen

conf = ConfGen(path/to/input.xyz)

conf.run()



