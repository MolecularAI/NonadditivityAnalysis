 **Please note: this repository is no longer being maintained.**

# NonadditivityAnalysis

Notebook for standardization of activity data, nonadditivity analysis and its evaluation.

A Python package and corresponding Jupyter notebook for:

1. Cleaning and standardizing ChEMBL activity data
2. Running nonadditivity analysis based on NAA code published by C. Kramer [1]
3. Evaluating the nonadditivity results

[1] Kramer C (2019) Nonadditivity Analysis. J Chem Inf Model 59:4034–4042. 
    https://doi.org/10.1021/acs.jcim.9b00631

## Installation

First, install a copy of the RDKit cheminformatics toolkit, available
from http://rdkit.org/. The easiest way is to install via PyPI with
`pip install rdkit-pypi`.

Install directly from source with:

```bash
$ pip install git+https://github.com/MolecularAI/NonadditivityAnalysis.git
```

Install the code in development mode with:

```bash
$ git clone git+https://github.com/MolecularAI/NonadditivityAnalysis.git
$ cd NonadditivityAnalysis
$ pip install -e .
```

## Usage

The Jupyter notebook can be run for any ChEMBL assay via the `get_processed_assay_df()` function.
While `ChEMBL1614027` is used as an example, any can be used by changing the `assay_chembl_id` variable.
