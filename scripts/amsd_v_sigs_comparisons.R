library(ggrepel)

############### MICE
# after amsd_mouse_plotting

compute_sig_tests_loop <- function(df,
                                   tissues = NULL,    # optional character vector, e.g. c("LIVER","LUNG")
                                   min_n = 2,         # minimum sample count in each group to run tests
                                   top_by = c("wilcox", "ttest")) {
  
  top_by <- match.arg(top_by)
  
  df <- dplyr::as_tibble(df)
  if (!is.null(tissues)) df <- df %>% filter(tissue %in% tissues)
  
  tissues_vec <- unique(df$tissue)
  if (length(tissues_vec) == 0) stop("No tissues found in the data (after optional filtering).")
  
  results_list <- vector("list", length = 0)
  rowid <- 1L
  
  for (t in tissues_vec) {
    df_t <- df %>% filter(tissue == t)
    
    # need a SPONTANEOUS group to compare against
    if (!any(df_t$exposure == "SPONTANEOUS")) {
      warning("Tissue '", t, "' has no SPONTANEOUS group — skipping.")
      next
    }
    
    exposures_to_test <- sort(unique(df_t$exposure))
    exposures_to_test <- exposures_to_test[exposures_to_test != "SPONTANEOUS"]
    if (length(exposures_to_test) == 0) {
      message("Tissue '", t, "' has no non-SPONTANEOUS exposures — skipping.")
      next
    }
    
    sigs <- unique(df_t$name)
    
    for (exp in exposures_to_test) {
      for (sig in sigs) {
        spont_vals <- df_t %>% filter(exposure == "SPONTANEOUS", name == sig) %>% pull(value)
        exp_vals   <- df_t %>% filter(exposure == exp,           name == sig) %>% pull(value)
        
        n_spont <- length(spont_vals)
        n_exp   <- length(exp_vals)
        mean_spont <- if (n_spont > 0) mean(spont_vals, na.rm = TRUE) else NA_real_
        mean_exp   <- if (n_exp > 0)   mean(exp_vals, na.rm = TRUE)   else NA_real_
        
        p_ttest <- NA_real_
        p_wilcox <- NA_real_
        
        if (n_spont >= min_n && n_exp >= min_n) {
          p_ttest  <- tryCatch(t.test(exp_vals, spont_vals)$p.value, error = function(e) NA_real_)
          p_wilcox <- tryCatch(wilcox.test(exp_vals, spont_vals, exact = FALSE)$p.value, error = function(e) NA_real_)
        }
        
        results_list[[rowid]] <- tibble(
          tissue     = t,
          exposure   = exp,
          name       = sig,
          n_spont    = n_spont,
          n_exp      = n_exp,
          mean_spont = mean_spont,
          mean_exp   = mean_exp,
          mean_diff  = mean_exp - mean_spont,
          p_ttest    = p_ttest,
          p_wilcox   = p_wilcox
        )
        rowid <- rowid + 1L
      } # sig
    } # exposure
  } # tissue
  
  results_df <- bind_rows(results_list)
  
  # multiple-testing corrections per tissue × exposure
  results_df <- results_df %>%
    group_by(tissue, exposure) %>%
    mutate(
      p_ttest_Bonf  = ifelse(all(is.na(p_ttest)), NA_real_, p.adjust(p_ttest, method = "bonferroni")),
      p_ttest_BH    = ifelse(all(is.na(p_ttest)), NA_real_, p.adjust(p_ttest, method = "BH")),
      p_wilcox_Bonf = ifelse(all(is.na(p_wilcox)), NA_real_, p.adjust(p_wilcox, method = "bonferroni")),
      p_wilcox_BH   = ifelse(all(is.na(p_wilcox)), NA_real_, p.adjust(p_wilcox, method = "BH"))
    ) %>%
    ungroup()
  
  # Top hit per tissue / exposure (choose which test to rank by)
  if (top_by == "wilcox") {
    top_hits <- results_df %>%
      group_by(tissue, exposure) %>%
      slice_min(order_by = p_wilcox, n = 1, with_ties = FALSE) %>%
      ungroup()
  } else {
    top_hits <- results_df %>%
      group_by(tissue, exposure) %>%
      slice_min(order_by = p_ttest, n = 1, with_ties = FALSE) %>%
      ungroup()
  }
  
  return(list(all = results_df, top = top_hits))
}


mexposuresig3 <- mexposuresig2 %>% filter(tissue %in% c("LIVER","LUNG"))

res <- compute_sig_tests_loop(mexposuresig3, tissues = c("LIVER","LUNG"), min_n = 2, top_by = "wilcox")
res$top %>% print(n=29)
# All per-signature tests:
test_count <- filter(res$all, !is.na(p_ttest)) %>%
  nrow()

# Top (most significant) signature per tissue×exposure:
merged_amsdsig <- mouse_amsd_output %>%
  mutate(exposure = str_replace(exposure, "1,2,3_TRICHLOROPROPANE", "TCP")) %>%
  arrange(mouse_amsd_output, tissue, exposure) %>%
  full_join(select(arrange(res$top, tissue, exposure), -p_ttest_Bonf, -p_ttest_BH, -p_wilcox_Bonf, -p_wilcox_BH))
merged_amsdsig %>%
  ggplot(aes(x = -log10(pvalues), y = -log10(p_ttest), color = tissue, label = exposure))+
    geom_point()+
  geom_label_repel() +
  geom_hline(yintercept = -log10(0.05))+
  geom_hline(yintercept = -log10(0.05/test_count))+
  geom_vline(xintercept = -log10(0.05))+
  geom_vline(xintercept = -log10(0.05/29))+
    xlab("AMSD p-value")+
    ylab("Top signature pvaue (ttest)")
merged_amsdsig %>%
  ggplot(aes(x = -log10(pvalues), y = -log10(p_wilcox), color = tissue, label = exposure))+
  geom_point()+
  geom_label_repel() +
  geom_hline(yintercept = -log10(0.05))+
  geom_hline(yintercept = -log10(0.05/test_count))+
  geom_vline(xintercept = -log10(0.05))+
  geom_vline(xintercept = -log10(0.05/29))+
  xlab("AMSD p-value")+
  ylab("Top signature pvaue (wilcox rank sum)")

mexposuresig3 %>%
  filter(exposure %in% c("VANADIUM_PENTOXIDE", "SPONTANEOUS"), tissue == "LUNG",name == "mSBS40") %>%
  ggplot(aes(x=exposure, y = value, color = name))+
  geom_boxplot()+
  geom_jitter()

#############################
#TCGA

# # 1. List files
# files <- list.files("../outputs", pattern = "^sigprofiler_exposures_.*\\.rds$", full.names = TRUE)
# 
# analyze_file <- function(file) {
#   dat <- readRDS(file)
#   
#   # 2. Convert counts to fractions per row
#   sig_cols <- grep("^SBS", names(dat), value = TRUE)
#   dat <- dat %>%
#     mutate(row_sum = rowSums(across(all_of(sig_cols)))) %>%
#     mutate(across(all_of(sig_cols), ~ .x / row_sum)) %>%
#     select(-row_sum)
#   
#   # Get which ancestries are present
#   ancestries <- unique(dat$ancestry)
#   if (length(ancestries) != 2) return(NULL)  # skip if not exactly 2
#   
#   # 3. Run wilcox test for each SBS
#   results <- map_dfr(sig_cols, function(sig) {
#     test <- wilcox.test(dat[[sig]] ~ dat$ancestry)
#     tibble(
#       signature = sig,
#       ancestry1 = ancestries[1],
#       ancestry2 = ancestries[2],
#       p_value = test$p.value,
#       median1 = median(dat[[sig]][dat$ancestry == ancestries[1]]),
#       median2 = median(dat[[sig]][dat$ancestry == ancestries[2]]),
#       file = basename(file),
#       tumor_type = unique(dat$tissue_type)
#     )
#   })
#   
#   results
# }
#########

# 1. List files
files <- list.files("../outputs", pattern = "^sigprofiler_exposures_.*\\.rds$", full.names = TRUE)

analyze_file <- function(file) {
  dat <- readRDS(file)
  
  # 2. Convert counts to fractions per row
  sig_cols <- grep("^SBS", names(dat), value = TRUE)
  dat <- dat %>%
    mutate(row_sum = rowSums(across(all_of(sig_cols)))) %>%
    mutate(across(all_of(sig_cols), ~ .x / row_sum)) %>%
    select(-row_sum)
  
  # add metadata columns
  dat <- dat %>%
    mutate(file = basename(file),
           tumor_type = unique(dat$tissue_type))
  
  # Get which ancestries are present
  ancestries <- unique(dat$ancestry)
  if (length(ancestries) != 2) {
    return(list(results = NULL, raw = dat))
  }
  
  # 3. Run wilcox test for each SBS
  results <- map_dfr(sig_cols, function(sig) {
    test <- wilcox.test(dat[[sig]] ~ dat$ancestry)
    #test <- t.test(dat[[sig]] ~ dat$ancestry)
    tibble(
      signature = sig,
      ancestry1 = ancestries[1],
      ancestry2 = ancestries[2],
      p_value = test$p.value,
      mean1 = mean(dat[[sig]][dat$ancestry == ancestries[1]]),
      mean2 = mean(dat[[sig]][dat$ancestry == ancestries[2]]),
      median1 = median(dat[[sig]][dat$ancestry == ancestries[1]]),
      median2 = median(dat[[sig]][dat$ancestry == ancestries[2]]),
      file = basename(file),
      tumor_type = unique(dat$tissue_type)
    )
  })
  
  list(results = results, raw = dat)
}

# 4. Apply to all files
all_outputs <- map(files, analyze_file)

# Separate into results and raw data
all_results <- map_dfr(all_outputs, "results")   # comparison results
all_raw     <- map_dfr(all_outputs, "raw")       # per-sample fractions

###########

tcga_test_count <- filter(all_results, !is.na(p_value)) %>% nrow()

# 5. just top hits
top_hits <- all_results %>%
  filter(!is.na(p_value), is.finite(p_value)) %>%   # remove NA/NaN p-values
  group_by(file) %>%
  slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

full_join(ancestry_amsd_output, top_hits) %>%
  #filter(min_anc_n > 10) %>%
  #ggplot(aes(x = -log10(pvalues), y = -log10(p_value), color = comparison, label = tumor_type, size = min_anc_n))+
  ggplot(aes(x = -log10(pvalues), y = -log10(p_value), color = comparison, label = tumor_type))+
  geom_point()+
  geom_label_repel() +
  geom_hline(yintercept = -log10(0.05))+
  geom_hline(yintercept = -log10(0.05/tcga_test_count))+
  geom_vline(xintercept = -log10(0.05))+
  geom_vline(xintercept = -log10(0.05/67))+
  xlab("AMSD p-value")+
  ylab("Top signature pvalue (wilcox)")


all_raw  %>%
  filter(tumor_type == "KIRP") %>% 
  select(-file) %>%
  unique() %>%
  ggplot(aes(x=ancestry,y=SBS85, color = tumor_type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.6)
all_raw  %>%
  select(-file) %>%
  unique() %>%
  filter(tumor_type == "SARC") %>%
  ggplot(aes(x=ancestry,y=SBS54, color = tumor_type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.6)
all_raw  %>%
  select(-file) %>%
  unique() %>%
  filter(tumor_type == "OV") %>%
  ggplot(aes(x=ancestry,y=SBS38, color = tumor_type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.6)

all_raw %>%
  filter(tumor_type == "OV", ancestry == "eas") %>%
  select(-file) %>%
  distinct() %>%
  # keep only SBS columns that are not all zeros
  select(where(~ !all(. == 0))) %>%
  pivot_longer(
    cols = starts_with("SBS"),
    names_to = "signature",
    values_to = "exposure"
  ) %>%
  filter(signature %in% c("SBS38","SBS53","SBS49","SBS7c", "SBS89","SBS3","SBS39")) %>%
  ggplot(aes(x = signature, y = exposure, color = sample_id)) +
    geom_jitter(height = 0,width = 0.2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

# quick check
all_raw_long %>%
  filter(tumor_type == "OV", ancestry == "eas") %>%
  head()


plot_sig_volcano <- function(anc1, anc2, tumor){
  p <- filter(ancestry_amsd_output,
                   tumor_type == tumor,
                   ancestry1 == anc1,
                   ancestry2 == anc2) %>%
  pull(pvalues)
  subset_tcga <- all_results %>%
    mutate(ancestries = paste0(ancestry1, "_v_", ancestry2)) %>%
    filter(tumor_type == tumor, 
           ancestries == paste0(anc1, "_v_", anc2), 
           !is.na(p_value))
  subset_tcga %>%
    ggplot(aes(x = (mean1 - mean2), y = -log10(p_value), label = signature))+
    geom_point()+
    geom_hline(yintercept = -log10(0.05))+
    geom_hline(yintercept = -log10(0.05/nrow(subset_tcga)))+
    geom_hline(yintercept = -log10(0.05/tcga_test_count))+
    geom_label_repel()+
    labs(title = paste0(tumor,": ",anc1, " vs ", anc2, ", AMSD pvalue = ", p),
         x = "Difference in means")
}
plot_sig_volcano("afr", "eas","LUAD")
plot_sig_volcano("eas", "eur","UCEC")
plot_sig_volcano("eas", "eur","KIRP")
plot_sig_volcano("eas", "eur","SARC")

pdf("../outputs/tcga_sigantures_volcano_plots.pdf")
for(i in 1:nrow(ancestry_amsd_output)){
  sample <- full_join(ancestry_amsd_output, top_hits) %>%
    arrange(pvalues)
  plot1 <- plot_sig_volcano(sample[i,"ancestry1"], sample[i,"ancestry2"],sample[i,"tumor_type"])
  print(plot1)
}
dev.off()

###########################
# load permutations
ancestry_amsd_perms <- readRDS(file = "../outputs/ancestry_amsd_perms.rds")
ancestry_amsd_perms
ggplot(ancestry_amsd_perms, aes(SARC.eas_eur))+
  geom_histogram()+
  geom_vline(xintercept = ancestry_amsd_output[53,"cosines"])+
  labs(title = "SARC: eas vs eur", x = "cosine distance",y = "permutation count")
ggplot(ancestry_amsd_perms, aes(KIRP.eas_eur))+ 
  geom_histogram()+
  geom_vline(xintercept = ancestry_amsd_output[26,"cosines"])+
  labs(title = "KIRP: eas vs eur", x = "cosine distance",y = "permutation count")
ggplot(ancestry_amsd_perms, aes(OV.eas_eur))+ 
  geom_histogram()+
  geom_vline(xintercept = ancestry_amsd_output[42,"cosines"])+
  labs(title = "OV: eas vs eur", x = "cosine distance",y = "permutation count")
ggplot(ancestry_amsd_perms, aes(LUAD.afr_eas))+ 
  geom_histogram()+
  geom_vline(xintercept = ancestry_amsd_output[34,"cosines"])+
  labs(title = "LUAD: afr vs eas", x = "cosine distance",y = "permutation count")
ggplot(ancestry_amsd_perms, aes(UCEC.eas_eur))+ 
  geom_histogram()+
  geom_vline(xintercept = ancestry_amsd_output[66,"cosines"])+
  labs(title = "UCEC: eas vs eur", x = "cosine distance",y = "permutation count")
