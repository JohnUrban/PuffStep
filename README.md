# PuffStep
PuffStep: HMM-based approach(es) to analyzing genomics data.

## Overview
PuffStep is a Python 3 toolkit for analyzing genomics data using Hidden Markov Models (HMMs). 

## Available Sub-commands
- **normalize** - Normalize bedGraph files using various protocols in preparation for puffcn HMM.
- **puffcn** - Copy number analysis using HMM to identify regions of different copy numbers.
- **summits** - Find summit bins within a set of regions based on the HMM state path (or other BED regions) and the normalized data (or other bedGraph).

## Requirements
- Python 3.6 or higher
- NumPy
- SciPy
- See `requirements.txt` for full dependencies

## Installation
```bash
pip install -r requirements.txt
python setup.py install
```

## Usage
```bash
# View available commands
python pufferfish/pufferfish_main.py --help

# Run puffcn analysis
python pufferfish/pufferfish_main.py puffcn -l late_stage.bedGraph -1 --outpfx output

# Normalize data
python pufferfish/pufferfish_main.py normalize -l data.bedGraph -1 -o normalized.bedGraph
```


## Ancestry
PuffStep descends from "pufferfish": https://github.com/JohnUrban/pufferfish

PuffStep is the modernized version, updated to:
- Use Python 3.6+
- Remove R/rpy2 dependencies (pure Python implementation)
- Streamline functionality to core analysis tools

PuffStep obtained a new name since "pufferfish" is the name of an index used elsewhere in bioinformatics (e.g. see the pufferfish index for Salmon quantification).

## Changes from pufferfish
- Migrated to Python 3
- Removed R/rpy2 dependencies - all HMM calculations now use pure Python (NumPy/SciPy)
- Removed deprecated sub-commands (mapreads, getcov, findpuffs, dump)
- Removed geckocorrection and switchblade modules
- Simplified and cleaned up codebase
