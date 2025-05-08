# Code to accompany manuscript: *A signature-agnostic test for differences between tumor mutation spectra reveals carcinogen and ancestry effects*

<br/>

## Table of Contents
- [Installation](#installation)
- [Reproducing manuscript analysis](#reproducing-manuscript-analysis)
- [Usage](#usage)
- [Citation](#citation)

<br/>

## Installation

Install package "[mutspecdist](https://github.com/sfhart33/mutspecdist)" to run AMSD

```
# install.packages("devtools")
# library(devtools)

devtools::install_github("sfhart33/mutspecdist")
```

<br/>

## Reproducing manuscript analysis

### Dependencies
*versions indicate those used for original analysis*
- tidyverse     v2.0.0
- ggpubr        v0.6.0   (for multi-panel plotting)
- ggrepel       v0.9.6   (for figure labels)
- RColorBrewer  v1.1-3   (color palettes)
- svglite       v2.1.3   (output scalable vector graphics)
- sigfit        v2.2     (signature fitting tool)

### Run scripts to generate figures
```
source("amsd_simulations.R")         # Generates simulations from Supp Fig 1
source("amsd_mouse_carcinogens.R")   # Runs AMSD on mouse tumors grouped by carcinogen exposure
source("amsd_mouse_plotting.R")      # Plots the mouse carcinogen results
source("amsd_asbestos.R")            # Runs AMSD on mesothelioma grouped by professional asbestos exposure and plots results
source("amsd_tcga_ancestry.R")       # Runs AMSD on TCGA cancer types grouped by genetic ancestry 
source("amsd_tcga_plotting.R")       # Plots the TCGA ancestry results
```

<br/>

## Usage

See example analysis in the [mutspecdist README](https://github.com/sfhart33/mutspecdist) for how to run AMSD on your own data

![AMSD method](outputs/Figure1_2025-04-29.png)

- Input mutation spectra for each sample in a cohort divided into two groups. As an example here we have a group of tumors unexposed (blue) or exposed (gold) to carcinogen.
- AMSD aggregates mutation spectra for each group and calculates the cosine distance between the aggregate spectra (green).
- AMSD also randomly reshuffles group labels to calculate the cosine distance between randomly sampled tumors (grey), repeating 1000+ times to create a null distribution expectation.
- AMSD then calculates a p-value from the fraction of random samplings that are greater than or equal to the observed distance between the two groups, as visualized in the histogram.
- To interpret what mutational mechanisms may be behind a significant difference, we recommend applying mutational signature fitting to the aggregate spectra and/or individual samples and comparing between the two groups.

<br/>

## Citation

If you use this R package or results for your own research, please use the following citation:

- Hart S.F.M., Alcala N., Feder A.F., Harris K. (under review) *A signature-agnostic test for differences between tumor mutation spectra reveals carcinogen and ancestry effects*.

All input data is from publically available datasets. If you used any of the results, we also recommend citing the cooresponding data source:

- Mouse carcinogen data is from [Riva et al.](https://www.nature.com/articles/s41588-020-0692-4)
- Asbestos data set from [Mangiante et al.](https://www.nature.com/articles/s41588-023-01321-1)
- TCGA ancestry calls from [Carrot-Zhang et al.](https://www.sciencedirect.com/science/article/pii/S1535610820302117?via%3Dihub)
- TCGA mutation calls from [Ellrott et al.](https://www.sciencedirect.com/science/article/pii/S2405471218300966?via%3Dihub)

Finally, if you build off the AMSD method for new applications, we recommend citing the origional method implementation, developed by [Sasani et al.](https://elifesciences.org/articles/89096)

