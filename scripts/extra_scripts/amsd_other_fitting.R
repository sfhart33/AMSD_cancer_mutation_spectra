library(tidyverse)
library(sigfit)
library(ggrepel)
library(parallel)
data("cosmic_signatures_v3.2")
library(mutspecdist)
library(reticulate)
# For home computer
# use_python("C:/Users/sfhar/AppData/Local/Programs/Python/Python313/python.exe", required = TRUE)
library(SigProfilerAssignmentR)



# testing with real data
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

UCEC_eur_names <- anc_calls %>%
  filter(tumor_type == "UCEC", anc3 == "eur") %>%
  pull(IID)
UCEC_eur <- tcga_3mer %>%
  filter(IID %in% UCEC_eur_names) %>%
  select(-IID, -mut_counts)

UCEC_eas_names <- anc_calls %>%
  filter(tumor_type == "UCEC", anc3 == "eas") %>%
  pull(IID)
UCEC_eas <- tcga_3mer %>%
  filter(IID %in% UCEC_eas_names) %>%
  select(-IID, -mut_counts)

UCEC_afr_names <- anc_calls %>%
  filter(tumor_type == "UCEC", anc3 == "afr") %>%
  pull(IID)
UCEC_afr <- tcga_3mer %>%
  filter(IID %in% UCEC_afr_names) %>%
  select(-IID, -mut_counts)




#--- helper: convert sigfit-style spectra to SigProfiler SBS96 input
convert_sigfit_to_SigProfiler <- function(sigfit_matrix) {
  mat_t <- t(sigfit_matrix)
  convert_label <- function(x) {
    ref   <- substr(x, 2, 2)
    alt   <- substr(x, 6, 6)
    left  <- substr(x, 1, 1)
    right <- substr(x, 3, 3)
    paste0(left, "[", ref, ">", alt, "]", right)
  }
  rownames(mat_t) <- sapply(rownames(mat_t), convert_label)
  
  sbs96_channels <- c(
    "A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T",
    "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T",
    "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T",
    "C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T",
    "C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T",
    "C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T",
    "G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T",
    "G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T",
    "G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T",
    "T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T",
    "T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
    "T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"
  )
  
  mat_ordered <- mat_t[sbs96_channels, , drop = FALSE]
  df <- as.data.frame(mat_ordered)
  df <- cbind(MutationType = rownames(df), df)
  rownames(df) <- NULL
  return(df)
}

UCEC_afr
test <- UCEC_afr
colnames(test) <- sbs96_channels 
test

cosmic_fit(
  samples = test,
  output = "test",
  input_type = "matrix",
  context_type = "96",
  collapse_to_SBS96 = TRUE#,
  # cosmic_version = 3.3,
  # genome_build = "GRCh38"
)



#--- main comparison function
compare_spectra_sigprofiler_top <- function(n_samples = 5,
                                            n_mutations = 50,
                                            sig_probs = c(SBS1 = 0.3, SBS5 = 0.6, SBS18 = 0.1),
                                            additional_sig = "SBS2",
                                            frac_extra = 0.2,
                                            n_sim = 1000,
                                            seed = 123#,
                                            #cosmic_version = 3.3,
                                            #genome_build = "GRCh38"
) {
  set.seed(seed)
  
  # --- Simulate spectra ---
  no_exp <- simulate_spectra(n_samples, n_mutations, sig_probs, cosmic_signatures_v3.2)
  no_exp <- as.data.frame(do.call(rbind, no_exp))
  
  n_extra <- round(n_mutations * frac_extra)   # single number
  
  if (n_extra > 0) {
    with_exp <- simulate_spectra(
      n_samples, n_mutations, sig_probs, cosmic_signatures_v3.2,
      additional_sig,
      n_extra = rep(n_extra, n_samples)
    )
  } else {
    # When n_extra = 0, do NOT pass additional_sig or n_extra
    with_exp <- simulate_spectra(
      n_samples, n_mutations, sig_probs, cosmic_signatures_v3.2
    )
  }
  with_exp <- as.data.frame(do.call(rbind, with_exp))
  
  # --- Convert to SigProfiler input ---
  sp_no  <- convert_sigfit_to_SigProfiler(no_exp)
  sp_yes <- convert_sigfit_to_SigProfiler(with_exp)
  
  # --- Define output directories ---
  out_no  <- paste0("Assignment_no_", seed)
  out_yes <- paste0("Assignment_yes_", seed)
  
  # --- Fit with SigProfiler ---
  cosmic_fit(
    samples = sp_no,
    output = out_no,
    input_type = "matrix",
    context_type = "96",
    collapse_to_SBS96 = TRUE#,
    #cosmic_version = cosmic_version,
    #genome_build = genome_build
  )
  
  cosmic_fit(
    samples = sp_yes,
    output = out_yes,
    input_type = "matrix",
    context_type = "96",
    collapse_to_SBS96 = TRUE#,
    #cosmic_version = cosmic_version,
    #genome_build = genome_build
  )
  
  # --- Read exposures ---
  exp_no <- read.delim(file.path(out_no, "Assignment_Solution",
                                 "Activities", "Assignment_Solution_Activities.txt"),
                       header = TRUE, sep = "\t", check.names = FALSE)
  
  exp_yes <- read.delim(file.path(out_yes, "Assignment_Solution",
                                  "Activities", "Assignment_Solution_Activities.txt"),
                        header = TRUE, sep = "\t", check.names = FALSE)
  
  exp_no$group <- "A"
  exp_yes$group <- "B"
  df <- bind_rows(exp_no, exp_yes)
  
  # --- Cleanup output dirs ---
  unlink(out_no, recursive = TRUE)
  unlink(out_yes, recursive = TRUE)
  
  # --- Statistical tests ---
  results <- df %>%
    pivot_longer(
      cols = -c(Samples, group),
      names_to = "signature",
      values_to = "exposure"
    ) %>%
    group_by(signature) %>%
    summarise(
      mean_A = mean(exposure[group == "A"]),
      mean_B = mean(exposure[group == "B"]),
      p_ttest = t.test(exposure ~ group)$p.value,
      p_wilcox = wilcox.test(exposure ~ group)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      # Multiple testing corrections
      p_ttest_Bonf = p.adjust(p_ttest, "bonferroni"),
      p_ttest_BH   = p.adjust(p_ttest, "BH"),
      p_wilcox_Bonf = p.adjust(p_wilcox, "bonferroni"),
      p_wilcox_BH   = p.adjust(p_wilcox, "BH"),
      n_samples = n_samples,
      n_mutations = n_mutations,
      additional_sig = additional_sig,
      frac_extra = frac_extra,
      seed = seed
    ) %>%
    arrange(p_ttest)
  
  # --- AMSD significance test ---
  amsd_p <- amsd(no_exp, with_exp, n_sim = n_sim, seed = seed, mean_or_sum = "sum")$p
  
  results <- results %>%
    mutate(amsd_p = amsd_p)
  
  # --- Return only top row (most significant by t-test) ---
  top_result <- results %>% slice(1)
  
  return(top_result)
}

