# AMSD (Aggregate Mutation Spectrum Distance) method and code to accompany manuscript

# AMSD overview

![alt text](./outputs/Figure1_2025-04-28_AMSD)

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Reproducing manuscript analysis](#reproducing-manuscript-analysis)
- [Citation](#citation)




## Installation

**Not yet operation as an R package**

```
# install.packages("devtools")
# library(devtools)

install_github("sfhart33/AMSD_cancer_mutation_spectra")
```

This package was built and tested with R version 4.3.1. Internally, this package uses tidyverse functions. Users do not need to load tidyverse separately; required functions are imported as needed.




## Usage




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

## Citation

If you use this R package or results for your own research, please use the following citation:

- *Not yet posted or published... coming soon*

All input data is from publically available datasets. If you used any of the results, we also recommend citing the cooresponding data source:

- Mouse carcinogen data is from [Riva et al.](https://www.nature.com/articles/s41588-020-0692-4)
- Asbestos data set from [Mangiante et al.](https://www.nature.com/articles/s41588-023-01321-1)
- TCGA ancestry calls from [Carrot-Zhang et al.](https://www.sciencedirect.com/science/article/pii/S1535610820302117?via%3Dihub)
- TCGA mutation calls from [Ellrott et al.](https://www.sciencedirect.com/science/article/pii/S2405471218300966?via%3Dihub)

Finally, if you build off the AMSD method for new applications, we recommend citing the origional method implementation, developed by [Sasani et al.](https://elifesciences.org/articles/89096)

