# urethane_top_md

# Dependencies #
numpy
parmed

# Usage: #
```bash
export PYTHONPATH=path/to/urethane_top_md # Files are imported, this solves any import issue
python3 path/to/urethane_top_md/polymer_manager.py path/to/urethane_polymer.pdb [ optional: path/to/gromacs.ndx ]
```

PDB has to have unique atomnames. ParmEd limitation.
No atom name renaming code snippet is written.

It has functionality to update GROMACS ndx file with a new atom selection. Now selection "Backobone atom + C=O and HN" is hard coded. 
Before using the optional functionality, a flag has to be changed.
