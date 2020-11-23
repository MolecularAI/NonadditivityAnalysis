# NonadditivityAnalysis
Notebook for standardization of actvity data, nonadditivity analysis and its evaluation.

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
from http://rdkit.org/ 

- A running version of mmpdb, a matched molecular pairs
database generation and analysis toolkit, available from
http://github.com/rdkit/mmpdb

- A running version of NAA, nonadditivity analysis code, available from
https://github.com/KramerChristian/NonadditivityAnalysis


In order to run the NAA code directly from jupyter notebook, 
you need to set an environmetal variable.
Therefore, on command line, generate directories and files:

```shell
cd $CONDA_PREFIX
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh
```

then edit ./etc/conda/activate.d/env_vars.sh as follows:
```shell
#!/bin/sh
export NAA='/path/to/naa/code/'
```

and edit ./etc/conda/deactivate.d/env_vars.sh as follows:
```shell
#!/bin/sh
unset NAA
```

Apart from this, standard scientific python libraries like scipy and 
numpy are required as well as seaborn for plot generation in the analysis part.
 

-------------------

## Usage


The jupyter notebook can be run directly with gzipped activity data downloaded from ChEMBL.
'my_path' and 'my_name' has to be adjusted at the beginning of the jupyter notebook to 
reflect the user's specific path and output name.





