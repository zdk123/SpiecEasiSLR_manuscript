# Notebook for _Disentangling microbial associations from hidden environmental and technical factors via latent graphical models_

This is a notebook containing data and scripts for reproducing the analysis.

The code and data, including source OTU biom files, mapping files, intermediate RData outputs,
as well as provenance relationships
are hosted as a [Synapse project](https://www.synapse.org/#!Synapse:syn20843558/)
and [google drive](https://drive.google.com/open?id=1ulTlSO4_FhLfv4nvmburD0_W6foWOhzN).
The code without the data is hosted at [github](https://github.com/zdk123/SpiecEasiSLR_manuscript).
These setup instructions are for Synapse, which requires an account and,
for command line access, and an API key.


## Installation ##

For running the network inference on your own data, base R and SpiecEasi is all you need.
This involves a few auxiliary R packages (devtools, phyloseq, Matrix, Rcpp and RcppArmadillo)
See the (SpiecEasi github page)[https://github.com/zdk123/SpiecEasi] for details.

For reproducing the manuscript or for accessing the input data or output results via Synapse,
you can set up the dependencies in one-shot in two different ways.

### docker ###
All stably-versioned dependencies are available from a Synapse-hosted docker image.

Mount the current directory as the working directory to use in interactive mode:
```sh
IM=docker.synapse.org/syn20843558/spieceasi-slr:d6cf2fb
docker pull ${IM}
docker tag ${IM} spieceasi-slr:latest
docker run -w /root -v ${PWD}:/root -ti spieceasi-slr:latest /bin/bash
```

Then this repo can be cloned from within the container:
```sh
git clone https://github.com/zdk123/SpiecEasiSLR_manuscript.git
cd SpiecEasiSLR_manuscript
```
### conda ###

Clone this repo, set up an conda environment and install SpiecEasi at the specified commit:
```sh
conda env create -n seslr -f docker/environment_minimal.yml
Rscript -e "devtools::install_github('zdk123/SpiecEasi', ref='d6cf2fb', upgrade='never')"
```

SpiecEasi is also (on bioconda)[https://bioconda.github.io/recipes/r-spieceasi/README.html],
but for version stability (i.e. to reproduce the manuscript's results,
stick to the devtools installation at the git commit).


## Project structure ##
This structure for this project is the same for the Synapse and Google Drive.

Main data and code files:
```
.
├── CompNet/ (simulated network eval)
│   ├── main.R (runner script)
│   └── simulator/ (helper simulator scripts)
├── QMP/ (Simulate data based on the quantitative microbiome project)
│   ├── working_dir/ (qiime working directory)
│   └── simulator/
│       ├── main.R (main script for running simulator protocols)
│       ├── model_functions.R (generate models from QMP counts)
│       ├── method_functions.R (network fitting functions)
│       └── eval_functions.R (evaluate functions)
├── public/ (public data and analyses)
│   ├── qiita/
│   │   └── XYZ/ (dataset name)
│   │       ├── README.txt
│   │       ├── BIOM/ (from qiita)
│   │       ├── mapping_files/ (from qiita)
│   │       ├── combined_mapping_file.txt
│   │       ├── process_biom.R (script 1)
│   │       ├── filter_phy.R (script 2)
│   │       ├── fit_nets.R (script 3)
│   │       └── robPCA.R (script 4)
│   ├── scripts/
│   │   ├── process_fun.R (helper scripts for 1)
│   │   ├── filter_fun.R (helper scripts for 2)
│   │   ├── analyze_net_funs.R (helper script for 3)
│   │   └── analyze_rob_pca_funs.R (helper scripts for 4)
│   ├── data/
│   │   ├── metadata.txt (public data info)
│   │   └── feature_cats.txt
│   ├── process_all.R (loop through public datasets. Calls script 1)
│   ├── filter_all.R (loop through public datasets. Calls script 2)
│   ├── fit_nets_all.R (loop through public datasets. Calls script 3)
│   ├── analyze_rob_pca_all.R (loop through public datasets. Calls script 4)
│   ├── analyze_compare_datasets.R
│   └── analyze_consensus_net.R
├── Demo/ (Latent variable graphical model)
│   └── make_demo.R (R script)
├── figures/
├── Synapse/
└── docker/
    ├── Dockerfile
    ├── environment.yml (conda for docker)
    └── environment_minimal.yml (conda for users)
```

Intermediate output files for public datasets:
```
.
└── public/ (public data and analyses)
    ├── qiita/
    │   └── XYZ/ (dataset name)
    │       ├── XYZphy.RData (R data output from 1)
    │       ├── XYZphyfilt.RData (R data output from 2)
    │       ├── XYZNetFits.RData (R data output from 3)
    │       └── rob_pca.RData (R data output from 4)
    ├── all_nets_list.RData (intermediate R data output)
    └── signed_nets_list.RData (intermediate R data output)
```

## Sync from Synapse ##

[SynapseSync](https://github.com/zdk123/SynapseSync) is a python lightweight library
for syncing Google Drive-backed files from Synapse to a local directory
(developed to support this repo, but should work for general use).
This is preinstalled in the docker/conda environment.

## Usage ##

To recreate figures & intermediate datasets

### Synthetic compositional datasets ###
R package dependencies:
* simulator
* SpiecEasi
* latex2exp
* ggplot2
* dplyr
* Matrix

```sh
cd CompNet
Rscript main.R
```

### Public datasets ###
R package dependencies:
* SpiecEasi
* pulsar
* dplyr
* ggplot2

Each of the the qiita datasets subdirectory contains processing scripts (1-4)
designed to process/filter/network the data. These can be run individually
from within each subdirectory.

The top level scripts will loop through the dataset subdirectories, call the
inner scripts and collect the results to a common data structure:
```sh
cd public
Rscript process_all.R
Rscript filter_all.R
```

Network inference step should be skipped if limited computational resources:
```sh
Rscript fit_nets_all.R
```

Analysis scripts:
```sh
cd public
Rscript analyze_compare_methods.R
Rscript analyze_rob_pca_all.R
Rscript analyze_compare_datasets.R
Rscript analyze_consensus_net.R
```

## Demo ##
Generate the composite figures for the toy example by running.
This requires the xkcd package for funky fonts/axes.

```sh
Rscript make_demo.R
```
