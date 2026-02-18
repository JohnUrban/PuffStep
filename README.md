# PuffStep
PuffStep: HMM-based approach(es) to analyzing genomics data.

## Overview
PuffStep is a Python 3 toolkit for analyzing genomics data using Hidden Markov Models (HMMs). 

Original application:

	- Copy Number segmentation 

	- e.g. to identify DNA "Puff" sequences "for FISH" experiments in Sciara coprophila. Hence its original name (pufferfish).

Other applications:

	- Chromatin modification ChIP-seq segmentation (e.g. H3K27me3)

	- Transcript factor binding ChIP-seq/Cut-and-Run (e.g. EcR)

	- Any genomic enrichment data.



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
python puffStep.py --help

# Normalize data - a variety of options offered. Examples:
python puffStep.py normalize -l late.bedGraph -e early.bedGraph --protocol1 > RCN.protocol1.bedGraph

python puffStep.py normalize -l late.bedGraph -e early.bedGraph --protocol12 --log2 > log2RCN.protocol12.bedGraph

# Run puffcn analysis
python puffStep.py puffcn -i RCN.bedGraph > states.bedGraph

# Find summits within peaks / amplicons / enriched regions.
python puffStep.py summits --states states.bedGraph -i RCN.bedGraph
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




PuffStep - PufferFish Version Info:
============
- Updated to PuffStep version 1.1.20260218 (02/18/2026)

- Various updates since 2020.

- pufferfish - 0.1.20200925 (09/25/2020)
	- Many new features and updates including new normalization protocols, discrete HMMs, Viterbi Training, new utilities, segmented motif stuff, etc.

- pufferfish - 0.0.0b (06/01/2020)

	- Now its own Git Repo.

- pufferfish - 0.0.0 (02/15/2016)

	- early development was in subdir of sciaraTools: https://github.com/JohnUrban/sciara-project-tools (see deprecated subdirectory)





Please cite PuffStep (pufferfish) as:
-------------------------------------------

- John Urban's 2016 PhD Thesis from Brown University: 

	- "The genome and DNA puff sequences of the fungus fly, Sciara coprophila, and genome-wide methods for studying DNA replication."

	- https://repository.library.brown.edu/studio/item/bdr:733543/

	- See Chapter 4


