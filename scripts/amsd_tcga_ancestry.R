library(tidyverse)
source("amsd_functions.R")

# Ancestry calls from: https://gdc.cancer.gov/about-data/publications/CCG-AIM-2020

# ANCESTRY
  anc_calls <- read.table("../inputs/tcga_ancestry_calls.txt", 
                          header=TRUE,
                          comment.char="",
                          sep = "\t") %>%
    select(IID = patient, tumor_type, consensus_ancestry) %>%
    mutate(anc2 = consensus_ancestry, anc3 = consensus_ancestry,) %>%
    mutate(across(anc2, str_replace, 'afr_admix', 'admix')) %>%
    mutate(across(anc2, str_replace, 'eas_admix', 'admix')) %>%
    mutate(across(anc2, str_replace, 'sas_admix', 'admix')) %>%
    mutate(across(anc2, str_replace, 'eur_admix', 'admix')) %>%
    mutate(across(anc3, str_replace, 'afr_admix', 'afr')) %>%
    mutate(across(anc3, str_replace, 'eas_admix', 'eas')) %>%
    mutate(across(anc3, str_replace, 'sas_admix', 'sas')) %>%
    mutate(across(anc3, str_replace, 'eur_admix', 'eur'))

# SPECTRA
  tcga_3mer <- read.table("../inputs/tcga_mutation_spectra.txt",
                          sep="\t",
                          header = TRUE) %>%
    separate(ID, into = c("a","b","c",NA,NA,NA,NA), sep = "-") %>%
    mutate(IID = paste(a,b,c,sep = "-")) %>%
    select(-a,-b,-c) %>%
    select(IID, everything())
  tcga_3mer_spectra <- tcga_3mer
  tcga_3mer_spectra[,2:97] <- tcga_3mer[,2:97]/rowSums(tcga_3mer[,2:97])
  tcga_3mer_spectra$mut_counts <- rowSums(select(tcga_3mer, -IID))
  tcga_3mer$mut_counts <- rowSums(select(tcga_3mer, -IID))
# merge
  anc_spectra <- inner_join(anc_calls, tcga_3mer_spectra)
  anc_counts <- inner_join(anc_calls, tcga_3mer)

  
######## TO DO #############v
  
  
# anc_calls
# tcga_3mer
# anc_spectra
# anc_counts
# 
# nrow(tcga_3mer)
# tcga_3mer$IID %>% unique() %>% length()
# 
#   column_to_rownames(var = "IID") %>%
#   select(-mut_counts)
