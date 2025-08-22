library(tidyverse)
library(sigfit)
library(ggrepel)
library(parallel)
data("cosmic_signatures_v3.2")
library(mutspecdist)
library(reticulate)
use_python("C:/Users/sfhar/AppData/Local/Programs/Python/Python313/python.exe", required = TRUE)
library(SigProfilerAssignmentR)

n_samples = 125
n_mutations = 50
sig_probs <- c(SBS1 = 0.3, SBS5 = 0.6, SBS18 = 0.1)
signatures <- cosmic_signatures_v3.2
additional_sig <- "SBS40"
n_extra = rep(n_mutations*0.15,n_samples)
n_sim = 1000
seed = 123
#set.seed(123)

# ss <- c("SBS2","SBS40") # spiky and flat signatures to test
# ns <- c(5,25,125,625) # samples per exposure to test
# ms <- c(50,2500) # mutations per sample to test (approx WES and WGS)
# fs <- c(0.2,0.1,0.05,0.02)# extra mutations - fraction of total mutations
# xs <- 1:100 # replicates to test each parameter set

# Run on a base set, then with exposures
no_exposure_test <- simulate_spectra(n_samples = n_samples,
                                     n_mutations = n_mutations,
                                     sig_probs = sig_probs,
                                     signatures = signatures)
no_exposure_test <- as.data.frame(do.call(rbind, no_exposure_test))
no_exposure_test2 <- no_exposure_test/rowSums(no_exposure_test) # convert to mutaion fractions
with_exposure_test <- simulate_spectra(n_samples = n_samples,
                                       n_mutations = n_mutations,
                                       sig_probs = sig_probs,
                                       signatures = signatures,
                                       additional_sig = additional_sig,
                                       n_extra = n_extra)
with_exposure_test <- as.data.frame(do.call(rbind, with_exposure_test))
with_exposure_test2 <- with_exposure_test/rowSums(with_exposure_test) # convert to mutaion fractions

amsd(no_exposure_test, with_exposure_test, n_sim = n_sim, seed = seed, 'mean_or_sum' = 'sum')$p
# amsd(no_exposure_test2, with_exposure_test2, n_sim = n_sim, seed = seed)


# fitting
# no_exposure_test
# with_exposure_test
exposures <- fit_signatures(counts = no_exposure_test, 
                              signatures = cosmic_signatures_v3.2,
                              iter = 2000, 
                              warmup = 1000, 
                              chains = 1, 
                              seed = 1756) %>%
  retrieve_pars(par = "exposures") #, hpd_prob = 0.90)


exposures2 <- fit_signatures(counts = with_exposure_test, 
                              signatures = cosmic_signatures_v3.2,
                              iter = 2000, 
                              warmup = 1000, 
                              chains = 1, 
                              seed = 1756) %>%
  retrieve_pars(par = "exposures") #, hpd_prob = 0.90)


# Convert matrices to dataframes with group labels
df1 <- as.data.frame(exposures$mean) %>% mutate(group = "A")
df2 <- as.data.frame(exposures2$mean) %>% mutate(group = "B")

df <- bind_rows(df1, df2)

# Run t-tests across all signatures
results <- df %>%
  pivot_longer(-group, names_to = "signature", values_to = "exposure") %>%
  group_by(signature) %>%
  summarise(
    p_value = t.test(exposure ~ group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj_Bonf = p.adjust(p_value, method = "bonferroni"),
         p_adj_BH = p.adjust(p_value, method = "BH"),
         p_adj_BY = p.adjust(p_value, method = "BY")) %>%
  arrange(p_value)

# Most significant signature
# best <- results %>% slice(1)
# 
# best


amsd(no_exposure_test, with_exposure_test, n_sim = n_sim, seed = seed, 'mean_or_sum' = 'sum')$p
head(results)




######
no_exposure_testB <- simulate_spectra(n_samples = n_samples,
                                     n_mutations = n_mutations,
                                     sig_probs = sig_probs,
                                     signatures = signatures)
no_exposure_testB <- as.data.frame(do.call(rbind, no_exposure_testB))
no_exposure_testB
no_exposure_testC <- simulate_spectra(n_samples = n_samples,
                                     n_mutations = n_mutations,
                                     sig_probs = sig_probs,
                                     signatures = signatures)
no_exposure_testC <- as.data.frame(do.call(rbind, no_exposure_testC))
no_exposure_testC
amsd(no_exposure_testB, no_exposure_testC, n_sim = n_sim, seed = seed, 'mean_or_sum' = 'sum')$p

exposures <- fit_signatures(counts = no_exposure_testB, 
                            signatures = cosmic_signatures_v3.2,
                            iter = 2000, 
                            warmup = 1000, 
                            chains = 1, 
                            seed = 1756) %>%
  retrieve_pars(par = "exposures") #, hpd_prob = 0.90)


exposures2 <- fit_signatures(counts = no_exposure_testC, 
                             signatures = cosmic_signatures_v3.2,
                             iter = 2000, 
                             warmup = 1000, 
                             chains = 1, 
                             seed = 1756) %>%
  retrieve_pars(par = "exposures") #, hpd_prob = 0.90)


# Convert matrices to dataframes with group labels
df1 <- as.data.frame(exposures$mean) %>% mutate(group = "A")
df2 <- as.data.frame(exposures2$mean) %>% mutate(group = "B")

df <- bind_rows(df1, df2)

# Run t-tests across all signatures
results <- df %>%
  pivot_longer(-group, names_to = "signature", values_to = "exposure") %>%
  group_by(signature) %>%
  summarise(
    p_value = t.test(exposure ~ group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj_Bonf = p.adjust(p_value, method = "bonferroni"),
         p_adj_BH = p.adjust(p_value, method = "BH"),
         p_adj_BY = p.adjust(p_value, method = "BY")) %>%
  arrange(p_value)


amsd(no_exposure_testB, no_exposure_testC, n_sim = n_sim, seed = seed, 'mean_or_sum' = 'sum')$p
head(results)


compare_spectra_ttest <- function(countsA, countsB, signatures = cosmic_signatures_v3.2, 
                                  iter = 2000, warmup = 1000, chains = 1, seed = 1234) {
  
  # Fit signatures
  exposuresA <- fit_signatures(
    counts = countsA,
    signatures = signatures,
    iter = iter, warmup = warmup, chains = chains, seed = seed
  ) %>% retrieve_pars(par = "exposures")
  
  exposuresB <- fit_signatures(
    counts = countsB,
    signatures = signatures,
    iter = iter, warmup = warmup, chains = chains, seed = seed
  ) %>% retrieve_pars(par = "exposures")
  
  # Convert to long dataframe
  dfA <- as.data.frame(exposuresA$mean) %>% mutate(group = "A")
  dfB <- as.data.frame(exposuresB$mean) %>% mutate(group = "B")
  
  df <- bind_rows(dfA, dfB)
  
  # Run t-tests across all signatures
  results <- df %>%
    pivot_longer(-group, names_to = "signature", values_to = "exposure") %>%
    group_by(signature) %>%
    summarise(
      p_value = t.test(exposure ~ group)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      p_adj_Bonf = p.adjust(p_value, method = "bonferroni"),
      p_adj_BH   = p.adjust(p_value, method = "BH"),
      p_adj_BY   = p.adjust(p_value, method = "BY")
    ) %>%
    arrange(p_value)
  
  # Return full results + best hit
  results
}

# compare_spectra_ttest_exome <- function(countsA, countsB, signatures = cosmic_signatures_v3.2, 
#                                   iter = 2000, warmup = 1000, chains = 1, seed = 1234) {
#   
#   # Fit signatures
#   exposuresA <- fit_signatures(
#     counts = countsA,
#     signatures = signatures,
#     opportunities = "human-exome",
#     iter = iter, warmup = warmup, chains = chains, seed = seed
#   ) %>% retrieve_pars(par = "exposures")
#   
#   exposuresB <- fit_signatures(
#     counts = countsB,
#     signatures = signatures,
#     opportunities = "human-exome",
#     iter = iter, warmup = warmup, chains = chains, seed = seed
#   ) %>% retrieve_pars(par = "exposures")
#   
#   # Convert to long dataframe
#   dfA <- as.data.frame(exposuresA$mean) %>% mutate(group = "A")
#   dfB <- as.data.frame(exposuresB$mean) %>% mutate(group = "B")
#   
#   df <- bind_rows(dfA, dfB)
#   
#   # Run t-tests across all signatures
#   results <- df %>%
#     pivot_longer(-group, names_to = "signature", values_to = "exposure") %>%
#     group_by(signature) %>%
#     summarise(
#       p_value = t.test(exposure ~ group)$p.value,
#       .groups = "drop"
#     ) %>%
#     mutate(
#       p_adj_Bonf = p.adjust(p_value, method = "bonferroni"),
#       p_adj_BH   = p.adjust(p_value, method = "BH"),
#       p_adj_BY   = p.adjust(p_value, method = "BY")
#     ) %>%
#     arrange(p_value)
#   
#   # Return full results + best hit
#   results
# }

compare_spectra_ttest_exome <- function(countsA, countsB, signatures = cosmic_signatures_v3.2, 
                                        iter = 2000, warmup = 1000, chains = 1, seed = 1234) {
  
  # Fit signatures
  exposuresA <- suppressMessages(
    fit_signatures(
      counts = countsA,
      signatures = signatures,
      opportunities = "human-exome",
      iter = iter, warmup = warmup, chains = chains, seed = seed, verbose = FALSE
    ) %>% retrieve_pars(par = "exposures")
  )
  
  exposuresB <- suppressMessages(
    fit_signatures(
      counts = countsB,
      signatures = signatures,
      opportunities = "human-exome",
      iter = iter, warmup = warmup, chains = chains, seed = seed, verbose = FALSE
    ) %>% retrieve_pars(par = "exposures")
  )
  
  # Convert to long dataframe
  dfA <- as.data.frame(exposuresA$mean) %>% mutate(group = "A", sample = row_number())
  dfB <- as.data.frame(exposuresB$mean) %>% mutate(group = "B", sample = row_number())
  df  <- bind_rows(dfA, dfB)
  
  df_long <- df %>%
    pivot_longer(cols = -c(group, sample), names_to = "signature", values_to = "exposure")
  
  # Run t-tests across all signatures, add group means
  results <- df_long %>%
    group_by(signature) %>%
    summarise(
      mean_A = mean(exposure[group == "A"]),
      mean_B = mean(exposure[group == "B"]),
      p_value = t.test(exposure ~ group)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      p_adj_Bonf = p.adjust(p_value, method = "bonferroni"),
      p_adj_BH   = p.adjust(p_value, method = "BH"),
      p_adj_BY   = p.adjust(p_value, method = "BY")
    ) %>%
    arrange(p_value)
  
  # Return summary + per-sample exposures
  return(list(
    results = results,
    fitted  = df_long
  ))
}



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

compare_spectra_ttest_exome(no_exposure_testB, no_exposure_testC)


eas_afr <- compare_spectra_ttest_exome(UCEC_eas, UCEC_afr)

sigs <- eas_afr$results %>%
  filter(p_value < 0.05,
         mean_A + mean_B > 0.02) %>%
  pull(signature)
eas_afr$fitted

eas_afr$results %>%
  ggplot(aes(x = (mean_A - mean_B), y = -log10(p_value), label = signature)) +
    geom_point()+
    geom_text_repel()+
  geom_hline(yintercept = -log10(0.05))+
  geom_hline(yintercept = -log10(0.05/78))+
  xlab("Difference in mean exposure: EAS vs AFR UCEC")+
  theme_classic()

eas_afr$fitted %>%
  mutate(group = case_when(group == "A" ~ "EAS",
                           group == "B"  ~ "AFR")) %>%
  filter(signature %in% sigs) %>%
  ggplot(aes(signature, exposure, fill = group)) +
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

t.test(pull(filter(eas_afr$fitted, group == "A", signature == "SBS2"), exposure),
       pull(filter(eas_afr$fitted, group == "B", signature == "SBS2"), exposure))
t.test(pull(filter(eas_afr$fitted, group == "A", signature == "SBS10a"), exposure),
       pull(filter(eas_afr$fitted, group == "B", signature == "SBS10a"), exposure))
wilcox.test(pull(filter(eas_afr$fitted, group == "A", signature == "SBS2"), exposure),
            pull(filter(eas_afr$fitted, group == "B", signature == "SBS2"), exposure))
wilcox.test(pull(filter(eas_afr$fitted, group == "A", signature == "SBS10a"), exposure),
            pull(filter(eas_afr$fitted, group == "B", signature == "SBS10a"), exposure))

eas_eur <- compare_spectra_ttest_exome(UCEC_eur, UCEC_eas)
eas_eur$results
eas_eur$fitted


sigs <- eas_eur$results %>%
  filter(p_value < 0.05,
         mean_A + mean_B > 0.02) %>%
  pull(signature)

eas_eur$results %>%
  ggplot(aes(x = (mean_A - mean_B), y = -log10(p_value), label = signature)) +
  geom_point()+
  geom_text_repel()+
  geom_hline(yintercept = -log10(0.05))+
  geom_hline(yintercept = -log10(0.05/78))+
  xlab("Difference in mean exposure: EUR vs EAS UCEC")+
  theme_classic()

eas_eur$fitted %>%
  mutate(group = case_when(group == "A" ~ "EUR",
                           group == "B"  ~ "EAS")) %>%
  filter(signature %in% c("SBS2","SBS13","SBS10a","SBS10b")) %>%
  ggplot(aes(signature, exposure, fill = group)) +
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


####### simulation as a function

compare_spectra <- function(n_samples = 5,
                            n_mutations = 2500,
                            sig_probs = c(SBS1 = 0.3, SBS5 = 0.6, SBS18 = 0.1),
                            signatures = cosmic_signatures_v3.2,
                            additional_sig = "SBS2",
                            frac_extra = 0.02,
                            n_sim = 1000,
                            seed = 123,
                            iter = 2000,
                            warmup = 1000,
                            chains = 1) {
  
  # Number of extra mutations per sample
  n_extra <- rep(n_mutations * frac_extra, n_samples)
  
  # Simulate without exposure
  no_exposure <- simulate_spectra(
    n_samples = n_samples,
    n_mutations = n_mutations,
    sig_probs = sig_probs,
    signatures = signatures
  )
  no_exposure <- as.data.frame(do.call(rbind, no_exposure))
  
  # Simulate with exposure
  with_exposure <- simulate_spectra(
    n_samples = n_samples,
    n_mutations = n_mutations,
    sig_probs = sig_probs,
    signatures = signatures,
    additional_sig = additional_sig,
    n_extra = n_extra
  )
  with_exposure <- as.data.frame(do.call(rbind, with_exposure))
  
  # --- AMSD test ---
  amsd_p <- amsd(no_exposure, with_exposure,
                 n_sim = n_sim, seed = seed, mean_or_sum = "sum")$p
  
  # --- Fit signatures ---
  exposures1 <- fit_signatures(counts = no_exposure,
                               signatures = signatures,
                               iter = iter, warmup = warmup,
                               chains = chains, seed = seed) %>%
    retrieve_pars(par = "exposures")
  exposures2 <- fit_signatures(counts = with_exposure,
                               signatures = signatures,
                               iter = iter, warmup = warmup,
                               chains = chains, seed = seed) %>%
    retrieve_pars(par = "exposures")
  
  df1 <- as.data.frame(exposures1$mean) %>% mutate(group = "A")
  df2 <- as.data.frame(exposures2$mean) %>% mutate(group = "B")
  df <- bind_rows(df1, df2)
  
  # --- Run t-tests ---
  results <- df %>%
    pivot_longer(-group, names_to = "signature", values_to = "exposure") %>%
    group_by(signature) %>%
    summarise(
      p_value = t.test(exposure ~ group)$p.value,
      .groups = "drop"
    ) %>%
    mutate(p_adj_Bonf = p.adjust(p_value, method = "bonferroni"),
           p_adj_BH   = p.adjust(p_value, method = "BH"),
           p_adj_BY   = p.adjust(p_value, method = "BY")) %>%
    arrange(p_value)
  
  best <- results %>% slice(1)
  
  # --- Return summary row ---
  tibble(
    n_samples = n_samples,
    n_mutations = n_mutations,
    additional_sig = additional_sig,
    frac_extra = frac_extra,
    amsd_p = amsd_p,
    best_signature = best$signature,
    best_p = best$p_value,
    best_p_adj_Bonf = best$p_adj_Bonf,
    best_p_adj_BH = best$p_adj_BH,
    best_p_adj_BY = best$p_adj_BY
  )
}

compare_spectra()

# --- Example usage ---
out1 <- compare_spectra()
out2 <- compare_spectra(n_samples = 5, frac_extra = 0.05, additional_sig = "SBS40", n_mutations = 50)

results_all <- bind_rows(out1, out2)
results_all
bind_rows(results_all, out2)
bind_rows(results, out2)
c(rep(50,10))


results <-data.frame()
count <- 0
for(sam in rep(5,10)){
  for(ex in c(0,0.1)){
    out <- compare_spectra(n_samples = sam, frac_extra = ex, additional_sig = "SBS2", n_mutations = 50)
    results <- bind_rows(results,out)
    count <- count + 1
    print(count)
  }
}
results1 <- results
results
for(sam in rep(25,10)){
  for(ex in c(0,0.2)){
    out <- compare_spectra(n_samples = sam, frac_extra = ex, additional_sig = "SBS40", n_mutations = 50)
    results <- bind_rows(results,out)
    count <- count + 1
    print(count)
  }
}
results2 <- results
results
for(sam in rep(5,10)){
  for(ex in c(0,0.1)){
    out <- compare_spectra(n_samples = sam, frac_extra = ex, additional_sig = "SBS40", n_mutations = 2500)
    results <- bind_rows(results,out)
    count <- count + 1
    print(count)
  }
}
results[1:20,] %>%
  ggplot(aes(-log10(amsd_p), -log10(best_p_adj_Bonf), color = frac_extra))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05))+
  geom_vline(xintercept = -log10(0.05))+
  ggtitle("n=5,mut=50,sig=2")

results[21:40,] %>%
  ggplot(aes(-log10(amsd_p), -log10(best_p_adj_Bonf), color = frac_extra))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05))+
  geom_vline(xintercept = -log10(0.05))+
  ggtitle("n=25,mut=50,sig=40")

results[41:60,] %>%
  ggplot(aes(-log10(amsd_p), -log10(best_p_adj_Bonf), color = frac_extra))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05))+
  geom_vline(xintercept = -log10(0.05))+
  ggtitle("n=5,mut=2500,sig=40")


convert_sigfit_to_SigProfiler <- function(sigfit_matrix) {
  # Load library (for downstream use)
  library(SigProfilerAssignmentR)
  
  # 1. Transpose so rows = mutation types, cols = samples
  mat_t <- t(sigfit_matrix)
  
  # 2. Convert mutation labels from "ACA>AAA" → "A[C>A]A"
  convert_label <- function(x) {
    ref <- substr(x, 2, 2)       # central base
    alt <- substr(x, 6, 6)       # mutated base
    left <- substr(x, 1, 1)      # left context
    right <- substr(x, 3, 3)     # right context
    paste0(left, "[", ref, ">", alt, "]", right)
  }
  rownames(mat_t) <- sapply(rownames(mat_t), convert_label)
  
  # 3. Ready-made standard SBS96 order
  sbs96_channels <- c(
    "A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T",
    "A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T",
    "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T",
    "A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T",
    "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T",
    "A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T",
    "C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T",
    "C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T",
    "C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T",
    "C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T",
    "C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T",
    "C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T",
    "G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T",
    "G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T",
    "G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T",
    "G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T",
    "G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T",
    "G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T",
    "T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T",
    "T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T",
    "T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
    "T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
    "T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T",
    "T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T"
  )
  
  # 4. Reorder rows to match SigProfiler SBS96 order
  mat_ordered <- mat_t[sbs96_channels, , drop = FALSE]
  
  # 5. Return as dataframe with MutationType column
  df <- as.data.frame(mat_ordered)
  df <- cbind(MutationType = rownames(df), df)
  rownames(df) <- NULL
  return(df)
}


sp_input <- convert_sigfit_to_SigProfiler(with_exposure_test)
Assignment <- cosmic_fit(
  samples = sp_input,       # can also use "SigProfilerAssignment_input.txt"
  output = "AssignmentResults",
  input_type = "matrix",    # since we passed a data frame / matrix
  context_type = "96",      # SBS96
  collapse_to_SBS96 = TRUE, # ensures 96-channel fitting
  cosmic_version = 3.3,     # latest COSMIC v3
  genome_build = "GRCh38"   # or "GRCh37" depending on your data
)