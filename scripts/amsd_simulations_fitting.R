n_samples = 5
n_mutations = 2500
sig_probs <- c(SBS1 = 0.3, SBS5 = 0.6, SBS18 = 0.1)
signatures <- cosmic_signatures_v3.2
additional_sig <- "SBS2"
n_extra = rep(n_mutations*0.02,n_samples)
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

compare_spectra_ttest_exome <- function(countsA, countsB, signatures = cosmic_signatures_v3.2, 
                                  iter = 2000, warmup = 1000, chains = 1, seed = 1234) {
  
  # Fit signatures
  exposuresA <- fit_signatures(
    counts = countsA,
    signatures = signatures,
    opportunities = "human-exome",
    iter = iter, warmup = warmup, chains = chains, seed = seed
  ) %>% retrieve_pars(par = "exposures")
  
  exposuresB <- fit_signatures(
    counts = countsB,
    signatures = signatures,
    opportunities = "human-exome",
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



# testing with real data

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

compare_spectra_ttest(no_exposure_testB, no_exposure_testC)


eas_afr <- compare_spectra_ttest(UCEC_eas, UCEC_afr)
eas_afr
eas_eur <- compare_spectra_ttest(UCEC_eur, UCEC_eas)
eas_eur

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
