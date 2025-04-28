# AMSD (Aggregate Mutation Spetrum Distance) method and code to accompany manuscript identifying carcinogen and ancestry effects on cancer mutation spectra

*Figure 1 from manuscript here to explain figure*

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Reproducing manuscript analysis](#reproducing-manuscript-analysis)
- [Citation](#citation)

## Installation

### Prerequisites for AMSD
- R			# minimum version
- devtools	# for installation
- tidyverse	# minimum version

### Install AMSD R package

```
library(devtools)
# install_github("sfhart33/AMSD_cancer_mutation_spectra") # Not yet operation as an R package
```

## Usage

## Reproducing manuscript analysis

```
source("amsd_simulations.R") # Generates simulations from Supp Fig 1

source("amsd_mouse_carcinogens.R") # Runs AMSD on mouse carcinogen groupings from [Riva et al.](https://www.nature.com/articles/s41588-020-0692-4)
source("amsd_mouse_plotting.R") # Plots the results

source("amsd_asbestos.R") # Runs AMSD on asbestos data set from [Mangiante et al.](https://www.nature.com/articles/s41588-023-01321-1)

source("amsd_tcga_ancestry.R") # Runs AMSD on TCGA cancer spectra (mutation calls from [Ellrott et al.](https://www.sciencedirect.com/science/article/pii/S2405471218300966?via%3Dihub)), grouped by cancer type and ancestry (ancestry calls from [Carrot-Zhang et al.](https://www.sciencedirect.com/science/article/pii/S1535610820302117?via%3Dihub))
source("amsd_tcga_plotting.R") # Plots the results
```



### Prerequisites for reproducing manuscript analysis
- 

## Citation

If you use any of these methods or data for your own research, please use the following citation:

*Not yet posted or published... coming soon*
