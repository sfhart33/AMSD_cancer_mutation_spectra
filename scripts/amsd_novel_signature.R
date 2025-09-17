
library(tidyverse)
library(sigfit)
data("cosmic_signatures_v3.2")

exp_signatures <- read.delim("../inputs/human_sbs96_filtered_v1_0.txt", 
           sep = "\t", 
           header = TRUE, 
           check.names = FALSE) %>%
    separate(MutationType, into = c("a","b"), sep = "\\[", remove = FALSE) %>%
    separate(b, into = c("b","c"), sep = "\\]") %>%
    arrange(b,a,c) %>%
    select(-a,-b,-c) %>%
    column_to_rownames(var = "MutationType") %>%
    t() %>%
    as.data.frame()

plot_spectrum(exp_signatures, pdf_path = "exposure_signatures_filtered.pdf")
exp_signatures

# match_signatures(exp_signatures[1:10,], cosmic_signatures_v3.2)
# exp_signatures*1000

mcmc_samples_fit <- fit_signatures(exp_signatures*1000, 
                                     cosmic_signatures_v3.2,
                                     iter = 2000, 
                                     warmup = 1000, 
                                     chains = 1, 
                                     seed = 1756)
plot_exposures(mcmc_samples = mcmc_samples_fit, pdf_path = "exposure_signatures_filtered_fraction.pdf")
plot_reconstruction(mcmc_samples = mcmc_samples_fit, pdf_path = "exposure_signatures_filtered_reconstruction.pdf")

################
exp_signatures2 <- read.delim("../inputs/human_sbs96_unfiltered_v1_0_TkNkoXo.txt", 
                             sep = "\t", 
                             header = TRUE, 
                             check.names = FALSE) %>%
  separate(MutationType, into = c("a","b"), sep = "\\[", remove = FALSE) %>%
  separate(b, into = c("b","c"), sep = "\\]") %>%
  arrange(b,a,c) %>%
  select(-a,-b,-c) %>%
  column_to_rownames(var = "MutationType") %>%
  t() %>%
  as.data.frame()
exp_signatures2

plot_spectrum(exp_signatures2, pdf_path = "exposure_signatures_unfiltered.pdf")


mcmc_samples_fit2 <- fit_signatures(exp_signatures2*1000, 
                                   cosmic_signatures_v3.2,
                                   iter = 2000, 
                                   warmup = 1000, 
                                   chains = 1, 
                                   seed = 1756)
plot_exposures(mcmc_samples = mcmc_samples_fit2, pdf_path = "exposure_signatures_unfiltered_fraction.pdf")
plot_reconstruction(mcmc_samples = mcmc_samples_fit2, pdf_path = "exposure_signatures_unfiltered_reconstruction.pdf")
