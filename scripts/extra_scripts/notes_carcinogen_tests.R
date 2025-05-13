library(tidyverse)

# load data
  setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github/AMSD_cancer_mutation_spectra/scripts")
  
  mouse_amsd_output <- readRDS("../outputs/mouse_amsd_output.rds")
  mouse_amsd_output$padj_BH <- p.adjust(mouse_amsd_output$pvalues, method="BH")
  mouse_amsd_output$padj_BONF <- p.adjust(mouse_amsd_output$pvalues, method="bonferroni")
  # mexposuresig <- readRDS("../inputs/mouse_exposuresig.rds")
  # mexposure <- readRDS("../inputs/mouse_exposure.rds")
  # mouse_carcinogen_counts <- readRDS("../inputs/mouse_carcinogen_spectra.rds") # counts
  # mouse_carcinogen_spectra <- mouse_carcinogen_counts/rowSums(mouse_carcinogen_counts) # spectra
  
# merge with metadata
  list <- read_delim("../inputs/mouse_carcinogens.txt") %>%
    arrange(Chemicals)
  livers <- filter(mouse_amsd_output, tissue == "LIVER") %>%
    arrange(exposure) %>%
    cbind(list)
  
  lung_ex <- c("ANTIMONY TRIOXIDE","COBALT METAL", "ISOBUTYL NITRITE", "NICKEL OXIDE",
    "NICKEL SUBSULFIDE", "NICKEL SULFATE HEXAHYDRATE", "SODIUM TUNGSTATE DIHYDRATE",
    "VANADIUM PENTOXIDE", "VINYLIDENE CHLORIDE")
  lungs <- filter(mouse_amsd_output, tissue == "LUNG") %>%
    arrange(exposure) %>%
    cbind(filter(list, Chemicals %in% lung_ex))
  lungs %>%
      arrange(pvalues)

  
  
  both <- rbind(livers, lungs) %>%
    arrange(pvalues)
  both
  sig_clear <- filter(both, pvalues < 0.05,
                      `NTP Bioassay Result**` == "Clear evidence") %>%
    nrow()
  sig_nclear <- filter(both, pvalues < 0.05,
                       `NTP Bioassay Result**` != "Clear evidence") %>%
    nrow()
  nsig_clear <- filter(both, pvalues > 0.05,
                       `NTP Bioassay Result**` == "Clear evidence") %>%
    nrow()
  nsig_nclear <- filter(both, pvalues > 0.05,
                        `NTP Bioassay Result**` != "Clear evidence") %>%
    nrow()  
  
  
  sig_clearB <- filter(both, pvalues < 0.05,
         `Ames Test**` == "Positive") %>%
    nrow()
  sig_nclearB <- filter(both, pvalues < 0.05,
         `Ames Test**` != "Positive") %>%
    nrow()
  nsig_clearB <- filter(both, pvalues > 0.05,
         `Ames Test**` == "Positive") %>%
    nrow()
  nsig_nclearB <- filter(both, pvalues > 0.05,
         `Ames Test**` != "Positive") %>%
    nrow()  
  
  fisher.test(matrix(c(sig_clear, sig_nclear, nsig_clear, nsig_nclear), nrow = 2)) # NTP
  fisher.test(matrix(c(sig_clearB, sig_nclearB, nsig_clearB, nsig_nclearB), nrow = 2)) # Ames
  
############
  
  sig_clear <- filter(both, padj_BH < 0.05,
                      `NTP Bioassay Result**` == "Clear evidence") %>%
    nrow()
  sig_nclear <- filter(both, padj_BH < 0.05,
                       `NTP Bioassay Result**` != "Clear evidence") %>%
    nrow()
  nsig_clear <- filter(both, padj_BH > 0.05,
                       `NTP Bioassay Result**` == "Clear evidence") %>%
    nrow()
  nsig_nclear <- filter(both, padj_BH > 0.05,
                        `NTP Bioassay Result**` != "Clear evidence") %>%
    nrow()    
  sig_clearB <- filter(both, padj_BH < 0.05,
                                   `Ames Test**` == "Positive") %>%
    nrow()
  sig_nclearB <- filter(both, padj_BH < 0.05,
                        `Ames Test**` != "Positive") %>%
    nrow()
  nsig_clearB <- filter(both, padj_BH > 0.05,
                        `Ames Test**` == "Positive") %>%
    nrow()
  nsig_nclearB <- filter(both, padj_BH > 0.05,
                         `Ames Test**` != "Positive") %>%
    nrow()  
  
  fisher.test(matrix(c(sig_clear, sig_nclear, nsig_clear, nsig_nclear), nrow = 2)) # NTP
  fisher.test(matrix(c(sig_clearB, sig_nclearB, nsig_clearB, nsig_nclearB), nrow = 2)) # Ames
  
  
  sig_clear <- filter(both, padj_BONF < 0.05,
                      `NTP Bioassay Result**` == "Clear evidence") %>%
    nrow()
  sig_nclear <- filter(both, padj_BONF < 0.05,
                       `NTP Bioassay Result**` != "Clear evidence") %>%
    nrow()
  nsig_clear <- filter(both, padj_BONF > 0.05,
                       `NTP Bioassay Result**` == "Clear evidence") %>%
    nrow()
  nsig_nclear <- filter(both, padj_BONF > 0.05,
                        `NTP Bioassay Result**` != "Clear evidence") %>%
    nrow()  
  fisher.test(matrix(c(sig_clear, sig_nclear, nsig_clear, nsig_nclear), nrow = 2)) # NTP  