# AMSD (Aggregate Mutation Spectrum Distance) method to detect significant differences in mutation spectra

<br/>

![AMSD method](outputs/Figure1_2025-04-28_AMSD.png)

- Input mutation spectra for each sample in a cohort divided into two groups. As an example here we have a group of tumors unexposed (blue) or exposed (gold) to carcinogen.
- AMSD aggregates mutation spectra for each group and calculates the cosine distance between the aggregate spectra (green).
- AMSD also randomly reshuffles group labels to calculate the cosine distance between randomly sampled tumors (grey), repeating 1000+ times to create a null distribution expectation.
- AMSD then calculates a p-value from the fraction of random samplings that are greater than or equal to the observed distance between the two groups, as visualized in the histogram.
- To interpret what mutational mechanisms may be behind a significant difference, we recommend applying mutational signature fitting to the aggregate spectra and/or individual samples and compared between the two groups.

<br/>

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Reproducing manuscript analysis](#reproducing-manuscript-analysis)
- [Citation](#citation)

<br/>



## Installation

**TO DO: Not yet operation as an R package**

```
# install.packages("devtools")
# library(devtools)

install_github("sfhart33/AMSD_cancer_mutation_spectra")
```

This package was built and tested with R version 4.3.1. Internally, this package uses tidyverse functions. Users do not need to load tidyverse separately; required functions are imported as needed.

<br/>


## Usage

**TO DO: Walk through example analysis from paper**


## Reproducing manuscript analysis

### Dependencies
*versions indicate those used for original analysis*
- tidyverse     v2.0.0
- ggpubr        v0.6.0   (for multi-panel plotting)
- RColorBrewer  v1.1-3   (color palettes)
- svglite       v2.1.3   (output scalable vector graphics)
- sigfit        v2.2     (signature fitting tool)

### Run scripts to generate figures
```
source("amsd_simulations.R")         # Generates simulations from Supp Fig 1
source("amsd_mouse_carcinogens.R")   # Runs AMSD on mouse tumors grouped by carcinogen exposure
source("amsd_mouse_plotting.R")      # Plots the mouse carcinogen results
source("amsd_asbestos.R")            # Runs AMSD on mesothelioma grouped by professional asbestos exposure
source("amsd_tcga_ancestry.R")       # Runs AMSD on TCGA cancer types grouped by genetic ancestry 
source("amsd_tcga_plotting.R")       # Plots the TCGA ancestry results
```

<br/>

## Citation

If you use this R package or results for your own research, please use the following citation:

- *Not yet posted or published... coming soon*

All input data is from publically available datasets. If you used any of the results, we also recommend citing the cooresponding data source:

- Mouse carcinogen data is from [Riva et al.](https://www.nature.com/articles/s41588-020-0692-4)
- Asbestos data set from [Mangiante et al.](https://www.nature.com/articles/s41588-023-01321-1)
- TCGA ancestry calls from [Carrot-Zhang et al.](https://www.sciencedirect.com/science/article/pii/S1535610820302117?via%3Dihub)
- TCGA mutation calls from [Ellrott et al.](https://www.sciencedirect.com/science/article/pii/S2405471218300966?via%3Dihub)

Finally, if you build off the AMSD method for new applications, we recommend citing the origional method implementation, developed by [Sasani et al.](https://elifesciences.org/articles/89096)

