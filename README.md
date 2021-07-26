# NonadditivityAnalysis

Notebook for standardization of activity data, nonadditivity analysis and its evaluation.

A jupyter notebook for:
1. Cleaning and standardizing ChEMBL activity data
2. Running nonadditivity analysis based on NAA code published by C. Kramer [1]
3. Evaluating the nonadditivity results

[1] Kramer C (2019) Nonadditivity Analysis. J Chem Inf Model 59:4034–4042. 
    https://doi.org/10.1021/acs.jcim.9b00631


--------------------

## Requirements

Installation requirements are the same as for the published NAA code:

- A copy of the RDKit cheminformatics toolkit, available
from http://rdkit.org/. The easiest way is to install via PyPI with
  `pip install rdkit-pypi`.
- A running version of NAA, nonadditivity analysis code, available from
https://github.com/KramerChristian/NonadditivityAnalysis
    ```shell
    $ pip install git+https://github.com/KramerChristian/NonadditivityAnalysis.git
    ```

Apart from this, standard scientific python libraries like scipy and 
numpy are required as well as seaborn and matplotlib for plot generation in the analysis part.
```shell
conda install -c anaconda scipy
conda install -c anaconda numpy
conda install -c conda-forge matplotlib
conda install seaborn
```

-------------------

## Usage


The jupyter notebook can be run directly with gzipped activity data downloaded from ChEMBL, 
as an example the activity data for ChEMBL1614027 (ChEMBL Version 27) is included in this package.
'my_path' and 'my_name' has to be adjusted at the beginning of the jupyter notebook to 
reflect the user's specific path and output name.





