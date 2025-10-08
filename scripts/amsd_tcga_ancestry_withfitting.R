library(tidyverse)
# source("amsd_functions.R")
library(mutspecdist)
library(reticulate)
# use_python("C:/Users/sfhar/AppData/Local/Programs/Python/Python313/python.exe", required = TRUE)
library(SigProfilerAssignmentR)
library(sigfit)
data("cosmic_signatures_v3.2")

# Ancestry calls from: https://gdc.cancer.gov/about-data/publications/CCG-AIM-2020

# # function to convert your colnames into SigProfiler SBS96 format
# convert_to_SBS96 <- function(df) {
#   mat <- as.matrix(df)
#   rownames(mat) <- NULL
#   
#   # transpose so samples = columns
#   mat_t <- t(mat)
#   
#   # Build fake SBS96 names (SigProfiler expects them in the COSMIC-defined order)
#   # Here we just number them SBS1 ... SBS96 unless you want exact context labels
#   colnames(mat_t) <- paste0("Sample", seq_len(ncol(mat_t)))
#   rownames(mat_t) <- paste0("Channel", seq_len(nrow(mat_t)))
#   
#   return(mat_t)
# }
make_sigprofiler_input <- function(mat) {
  # helper: convert your column names -> SBS96 format
  convert_to_sbs96 <- function(colnames_vec) {
    sapply(colnames_vec, function(name) {
      parts <- strsplit(name, "\\.")[[1]]
      sub <- parts[1]     
      ctx <- parts[2]     
      
      sub_parts <- strsplit(sub, "_")[[1]]
      ref <- sub_parts[1]  
      alt <- sub_parts[2]  
      
      left  <- substr(ctx, 1, 1)
      mid   <- substr(ctx, 2, 2)
      right <- substr(ctx, 3, 3)
      
      if (mid != ref) {
        warning(paste("Ref base mismatch in", name, "— using context mid:", mid))
      }
      
      paste0(left, "[", ref, ">", alt, "]", right)
    })
  }
  
  # step 1: convert names
  new_names <- convert_to_sbs96(colnames(mat))
  
  # step 2: replace column names
  converted <- mat
  colnames(converted) <- new_names
  
  # step 3: transpose so mutation types are rows
  converted_t <- t(converted)
  
  # step 4: add rownames as first column
  out_df <- data.frame(MutationType = rownames(converted_t), converted_t, 
                       check.names = FALSE, stringsAsFactors = FALSE) %>%
    arrange(MutationType) %>%
    select(-MutationType)
  out_df
}




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
###################### NEW    THIS DOESNT WORK YET
    # --- Run SigProfiler fitting for ancestry group 1 (sample-level) ---
    sp_anc1 <- anc_counts %>%
      filter(tumor_type == tissue_type,
             anc3 == ancestry1,
             mut_counts >= min_mutations) %>%
      select(-IID, -tumor_type, -consensus_ancestry, -anc2, -anc3, -mut_counts) %>%
      make_sigprofiler_input()  # rows = samples, cols = SBS96
    outdir_anc1 <- paste0("../outputs/sigprofiler_", tissue_type, "_", ancestry1)
    
    "//gs-ddn2/gs-vol1/home/sfhart/github/AMSD_cancer_mutation_spectra/outputs/figure_3_revised.png"(
      samples = sp_anc1,
      output = outdir_anc1,
      input_type = "matrix",
      context_type = "96",
      collapse_to_SBS96 = TRUE,
      exome = TRUE
    )
    
    exp_anc1 <- read.delim(
      file.path(outdir_anc1, "Assignment_Solution", "Activities", "Assignment_Solution_Activities.txt"),
      sep = "\t", check.names = FALSE
    )
    
    # keep sample IDs from the colnames of sp_anc1
    exp_anc1 <- exp_anc1 %>%
      mutate(
        sample_id = colnames(sp_anc1),
        ancestry = ancestry1,
        tissue_type = tissue_type
      )
    
    # --- Run SigProfiler fitting for ancestry group 2 (sample-level) ---
    sp_anc2 <- anc_counts %>%
      filter(tumor_type == tissue_type,
             anc3 == ancestry2,
             mut_counts >= min_mutations) %>%
      select(-IID, -tumor_type, -consensus_ancestry, -anc2, -anc3, -mut_counts) %>%
      make_sigprofiler_input()
    outdir_anc2 <- paste0("../outputs/sigprofiler_", tissue_type, "_", ancestry2)
    
    cosmic_fit(
      samples = sp_anc2,
      output = outdir_anc2,
      input_type = "matrix",
      context_type = "96",
      collapse_to_SBS96 = TRUE,
      exome = TRUE
    )
    
    exp_anc2 <- read.delim(
      file.path(outdir_anc2, "Assignment_Solution", "Activities", "Assignment_Solution_Activities.txt"),
      sep = "\t", check.names = FALSE
    )
    
    exp_anc2 <- exp_anc2 %>%
      mutate(
        sample_id = colnames(sp_anc2),
        ancestry = ancestry2,
        tissue_type = tissue_type
      )
    
    # --- Combine exposures ---
    all_exposures <- bind_rows(exp_anc1, exp_anc2)
    
    # Save combined exposures for this comparison
    saveRDS(
      all_exposures,
      paste0("../outputs/sigprofiler_exposures_", tissue_type,
             "_", ancestry1, "_vs_", ancestry2, ".rds")
    )
    
    # (Optional) clean up output dirs
    unlink(outdir_anc1, recursive = TRUE)
    unlink(outdir_anc2, recursive = TRUE)
    ###################### NEW 
    
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
