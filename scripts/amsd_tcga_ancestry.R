library(tidyverse)
source("amsd_functions.R")

# Ancestry calls from: https://gdc.cancer.gov/about-data/publications/CCG-AIM-2020

# ANCESTRY (anc2 excludes admixed, anc3 includes with closest ancestry group)
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
  saveRDS(anc_spectra, "../outputs/ancestry_spectra.rds")
  
# Run AMSD
 min_mutations <- 10
 min_sample <- 5
  
  anc_tumor_counts <- anc_counts %>% 
    filter(mut_counts >= min_mutations) %>%
    count(tumor_type, anc3) %>%
    filter(anc3 %in% c("afr", "eas","eur")) %>%
    filter(n >= min_sample) %>%
    select(-n)
  
  anc_tumor_comparisons <- select(anc_tumor_counts, tumor_type, ancestry1=anc3) %>%
    full_join(select(anc_tumor_counts, tumor_type, ancestry2=anc3),
              relationship = "many-to-many") %>%
    mutate(comparison = paste0(ancestry1, "_", ancestry2)) %>%
    filter(comparison %in% c("afr_eas","afr_eur","eas_eur"))
  anc_counts

  # Blank variables for outputs
  pvalues <- c()
  cosines <- c()
  reps <- 100000
  perms <- data.frame(rep = 1:reps)

  # run for each carcinogen
  for(count in 1:nrow(anc_tumor_comparisons)){

    # set variables
    tissue_type <- anc_tumor_comparisons[count,"tumor_type"]
    ancestry1 <- anc_tumor_comparisons[count,"ancestry1"]
    ancestry2 <- anc_tumor_comparisons[count,"ancestry2"]

    # Sample groupings
    samples_anc1 <- anc_spectra %>%
      filter(tumor_type == tissue_type,
             anc3 == ancestry1,
             mut_counts >= min_mutations) %>%
      select(-IID, -tumor_type, -consensus_ancestry, -anc2, -anc3, -mut_counts)
    samples_anc2 <- anc_spectra %>%
      filter(tumor_type == tissue_type,
             anc3 == ancestry2,
             mut_counts >= min_mutations) %>%
      select(-IID, -tumor_type, -consensus_ancestry, -anc2, -anc3, -mut_counts)
    
    # Run AMSD
    amsd_output <- amsd(samples_anc1,
                        samples_anc2,
                        mean_or_sum = "mean",
                        n_sim = reps,
                        seed = 123)
    
    # Save output
    pvalues <- c(pvalues, amsd_output$p)
    cosines <- c(cosines, amsd_output$cosine)
    perms1 <- data.frame(perms = amsd_output$sims)
    colnames(perms1) <- paste0(tissue_type, ".", ancestry1,"_", ancestry2)
    perms <- cbind(perms, perms1)
    print(paste(count, "of", nrow(anc_tumor_comparisons)))
    
  }

  
  # save outputs
  ancestry_amsd_output <- anc_tumor_comparisons
  ancestry_amsd_output$pvalues <- pvalues
  ancestry_amsd_output$cosines <- cosines

# number of samples per run
  n_anc_samples <- c()
  for(count in 1:nrow(anc_tumor_comparisons)){
    tissue_type <- anc_tumor_comparisons[count,"tumor_type"]
    ancestry1 <- anc_tumor_comparisons[count,"ancestry1"]
    ancestry2 <- anc_tumor_comparisons[count,"ancestry2"]
    n_ancestry1 <- filter(anc_counts,
                          tumor_type == tissue_type,
                          anc3 == ancestry1) %>%
      nrow()
    n_ancestry2 <- filter(anc_counts,
                          tumor_type == tissue_type,
                          anc3 == ancestry2) %>%
      nrow()
    n_anc_samples <- c(n_anc_samples, min(n_ancestry1,n_ancestry2))
  }
  ancestry_amsd_output$min_anc_n <- n_anc_samples
# save
  saveRDS(ancestry_amsd_output, "../outputs/ancestry_amsd_output.rds")
  saveRDS(perms, "../outputs/ancestry_amsd_perms.rds")
