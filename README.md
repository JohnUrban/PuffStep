# PuffStep
PuffStep: HMM-based approach(es) to analyzing genomics data.

## Overview
PuffStep is a Python 3 toolkit for analyzing genomics data using Hidden Markov Models (HMMs). It is the modernized version of pufferfish, updated to:
- Use Python 3.6+
- Remove R/rpy2 dependencies (pure Python implementation)
- Streamline functionality to core analysis tools

## Available Sub-commands
- **normalize** - Normalize bedGraph files using various protocols
- **puffcn** - Copy number analysis using HMM to identify regions of different copy numbers
- **summits** - Find summits in normalized data
- **generate** - Generate synthetic HMM data for testing
- **filter** - Filter and process bedGraph files

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

## Changes from pufferfish
- Migrated to Python 3
- Removed R/rpy2 dependencies - all HMM calculations now use pure Python (NumPy/SciPy)
- Removed deprecated sub-commands (mapreads, getcov, findpuffs, dump)
- Removed geckocorrection and switchblade modules
- Simplified and cleaned up codebase
