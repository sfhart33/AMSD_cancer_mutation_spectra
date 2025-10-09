library(ggrepel)
library(ggpubr)
# library(ggtext)
# library(patchwork)
#source("amsd_weighted_v_unweighted.R")
# setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github/AMSD_cancer_mutation_spectra/scripts")
source("amsd_vs_sigvar.R")

############### MICE
# load data
  mouse_carcinogen_counts <- readRDS("../inputs/mouse_carcinogen_spectra.rds") # counts
  mouse_carcinogen_spectra <- mouse_carcinogen_counts/rowSums(mouse_carcinogen_counts) # spectra
  # list of samples and attributes
  sample_table <- mouse_carcinogen_spectra %>%
    rownames_to_column(var = "label") %>%
    select(label) %>%
    separate(label, c("tissue","rest"), sep = "_", extra = "merge", remove = FALSE) %>%
    separate(rest, 
             into = c("exposure", "rep"), 
             sep = "_(?=\\d+$)", # Regex to split based on underscores near the end
             extra = "merge",            # Handle extra columns by merging into one
             fill = "right")             # Handle missing parts gracefully
  mexposuresig <- readRDS("../inputs/mouse_exposuresig.rds")
  mexposuresigX <- mexposuresig/rowSums(mexposuresig) 
  mexposuresig2 <- mexposuresigX %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    pivot_longer(cols = -sample) %>%
    mutate(label = str_replace_all(sample, " ", "_")) %>%
    left_join(sample_table) %>%
    mutate(exposure = (replace(exposure, exposure == "1,2,3_TRICHLOROPROPANE", "TCP")))
  mexposuresig3 <- mexposuresig2 %>% filter(tissue %in% c("LIVER","LUNG"))

# get p-values for all tests
#####################################################################
  compute_sig_tests_loop <- function(df,
                                     tissues = NULL,    # optional vector, e.g. c("LIVER","LUNG")
                                     min_n = 2,         # min samples per group
                                     top_by = c("wilcox", "ttest")) {
    
    top_by <- match.arg(top_by)
    df <- dplyr::as_tibble(df)
    if (!is.null(tissues)) df <- df %>% filter(tissue %in% tissues)
    
    tissues_vec <- unique(df$tissue)
    if (length(tissues_vec) == 0) stop("No tissues found in the data (after filtering).")
    
    results_list <- list()
    rowid <- 1L
    
    for (t in tissues_vec) {
      df_t <- df %>% filter(tissue == t)
      
      if (!any(df_t$exposure == "SPONTANEOUS")) {
        warning("Tissue '", t, "' has no SPONTANEOUS group — skipping.")
        next
      }
      
      exposures_to_test <- setdiff(unique(df_t$exposure), "SPONTANEOUS")
      if (length(exposures_to_test) == 0) next
      
      sigs <- unique(df_t$name)
      
      for (exp in exposures_to_test) {
        for (sig in sigs) {
          spont_vals <- df_t %>% filter(exposure == "SPONTANEOUS", name == sig) %>% pull(value)
          exp_vals   <- df_t %>% filter(exposure == exp,           name == sig) %>% pull(value)
          
          n_spont <- length(spont_vals)
          n_exp   <- length(exp_vals)
          mean_spont <- if (n_spont > 0) mean(spont_vals, na.rm = TRUE) else NA_real_
          mean_exp   <- if (n_exp > 0)   mean(exp_vals, na.rm = TRUE)   else NA_real_
          
          p_ttest  <- NA_real_
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
    
    # Bonferroni correction *within each tissue × exposure*
    results_df <- results_df %>%
      group_by(tissue, exposure) %>%
      mutate(
        n_tests = sum(!is.na(p_ttest)),  # number of valid tests (signatures)
        # p_ttest_Bonf  = ifelse(all(is.na(p_ttest)), NA_real_, pmin(p_ttest * n_tests, 1)),
        # p_wilcox_Bonf = ifelse(all(is.na(p_wilcox)), NA_real_, pmin(p_wilcox * n_tests, 1))
      ) %>%
      ungroup()
    
    # Optional: top hit per tissue/exposure
    top_hits <- results_df %>%
      group_by(tissue, exposure) %>%
      slice_min(order_by = if (top_by == "wilcox") p_wilcox else p_ttest, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    # Return both full results and top hits
    return(list(all = results_df, top = top_hits))
  }
  
#################################################################  
  
  res <- compute_sig_tests_loop(mexposuresig3, tissues = c("LIVER","LUNG"), min_n = 2, top_by = "wilcox")
  res$top %>% print(n=29)
  # All per-signature tests:
  test_count <- filter(res$all, !is.na(p_ttest)) %>%
    nrow()
  
  # Top (most significant) signature per tissue×exposure:
  merged_amsdsig <- mouse_amsd_output %>%
    mutate(exposure = str_replace(exposure, "1,2,3_TRICHLOROPROPANE", "TCP")) %>%
    arrange(mouse_amsd_output, tissue, exposure) %>%
    full_join(arrange(res$top, tissue, exposure)) %>%
    mutate(p_wilcox_Bonf = pmin(p_wilcox * n_tests, 1))

# plot 
  # merged_amsdsig %>%
  #   ggplot(aes(x = -log10(pvalues), y = -log10(p_ttest), color = tissue, label = exposure))+
  #     geom_point()+
  #   geom_label_repel() +
  #   geom_hline(yintercept = -log10(0.05))+
  #   geom_hline(yintercept = -log10(0.05/test_count))+
  #   geom_vline(xintercept = -log10(0.05))+
  #   geom_vline(xintercept = -log10(0.05/29))+
  #     xlab("AMSD p-value")+
  #     ylab("Top signature p-value (ttest)")
  merged_amsdsig %>%
    ggplot(aes(x = -log10(pvalues), y = -log10(p_wilcox), color = tissue, label = exposure))+
    geom_point()+
    geom_smooth(method = "lm",
                inherit.aes = FALSE,
                aes(-log10(pvalues),
                    -log10(p_wilcox)),
                color = "black")+
    theme_classic()+
    geom_label_repel() +
    geom_hline(yintercept = -log10(0.05))+
    geom_hline(yintercept = -log10(0.05/test_count))+
    geom_vline(xintercept = -log10(0.05))+
    geom_vline(xintercept = -log10(0.05/29))+
    xlab("AMSD p-value")+
    ylab("Top signature p-value \n(wilcoxon rank sum)")

  mouse_amsdsig_plot <- merged_amsdsig %>%
    ggplot(aes(x = -log10(pvalues), y = -log10(p_wilcox_Bonf), color = tissue, label = exposure))+
    geom_point()+
    geom_smooth(method = "lm",
                inherit.aes = FALSE,
                aes(-log10(pvalues),
                    -log10(p_wilcox_Bonf)),
                color = "black")+
    theme_classic()+
    geom_label_repel() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = -log10(0.05))+
    geom_hline(yintercept = -log10(0.05/29))+
    geom_vline(xintercept = -log10(0.05))+
    geom_vline(xintercept = -log10(0.05/29))+
    xlab("AMSD p-value")+
    ylab("Top signature adj p-value \n(wilcoxon rank sum)")
  mouse_amsdsig_plot
  #outlier
  mexposuresig3 %>%
    filter(exposure %in% c("VANADIUM_PENTOXIDE", "SPONTANEOUS"), tissue == "LUNG",name == "mSBS40") %>%
    ggplot(aes(x=exposure, y = value, color = name))+
    geom_boxplot()+
    geom_jitter()

# spectra comparison
  ox <- mouse_carcinogen_spectra[pull(filter(sample_table, tissue == "LIVER", exposure == "OXAZEPAM"), label),]
  ox
  tcp <- mouse_carcinogen_spectra[pull(filter(sample_table, tissue == "LIVER", exposure == "1,2,3_TRICHLOROPROPANE"), label),]
  tcp
  spon <- mouse_carcinogen_spectra[pull(filter(sample_table, tissue == "LIVER", exposure == "SPONTANEOUS"), label),]
  spon
  
  ox_long <- ox %>%
    as.data.frame() %>%
    mutate(sample = paste0("ox_", row_number()), group = "ox") %>%
    pivot_longer(-c(sample, group), names_to = "trinuc", values_to = "freq")
  
  spon_long <- spon %>%
    as.data.frame() %>%
    mutate(sample = paste0("spon_", row_number()), group = "spon") %>%
    pivot_longer(-c(sample, group), names_to = "trinuc", values_to = "freq")
  
  tcp_long <- tcp %>%
    as.data.frame() %>%
    mutate(sample = paste0("tcp_", row_number()), group = "tcp") %>%
    pivot_longer(-c(sample, group), names_to = "trinuc", values_to = "freq")
  
  dat_long <- bind_rows(ox_long, spon_long, tcp_long)
  
  # Run Wilcoxon per trinucleotide
  results_ox <- dat_long %>%
    group_by(trinuc) %>%
    summarise(
      pval = wilcox.test(freq[group == "ox"], freq[group == "spon"])$p.value,
      mean_ox = mean(freq[group == "ox"], na.rm = TRUE),
      mean_spon = mean(freq[group == "spon"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      effect = mean_ox - mean_spon,            # raw difference
      log2fc = log2((mean_ox + 1e-6) / (mean_spon + 1e-6)),  # fold change
      padj = p.adjust(pval, method = "BH")     # FDR correction
    )
  results_ox
  results_tcp <- dat_long %>%
    group_by(trinuc) %>%
    summarise(
      pval = wilcox.test(freq[group == "tcp"], freq[group == "spon"])$p.value,
      mean_tcp = mean(freq[group == "tcp"], na.rm = TRUE),
      mean_spon = mean(freq[group == "spon"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      effect = mean_tcp - mean_spon,            # raw difference
      log2fc = log2((mean_tcp + 1e-6) / (mean_spon + 1e-6)),  # fold change
      padj = p.adjust(pval, method = "BH")     # FDR correction
    )
  results_tcp
  # Volcano plot (effect size vs -log10 p-value)
  ggplot(results_ox, aes(x = effect, y = -log10(pval), label = trinuc)) +
    geom_point(aes(color = padj < 0.1)) +
    geom_label_repel()+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() +
    labs(
      x = "Mean difference (oxazempam - spontaneous)",
      y = "-log10(p-value)",
      color = "FDR < 0.1"
    )
  ggplot(results_tcp, aes(x = effect, y = -log10(pval), label = trinuc)) +
    geom_point(aes(color = padj < 0.1)) +
    geom_label_repel()+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() +
    labs(
      x = "Mean difference (TCP - spontaneous)",
      y = "-log10(p-value)",
      color = "FDR < 0.1"
    )
#############################
#TCGA

# 1. List files
files <- list.files("../outputs", pattern = "^sigprofiler_exposures_.*\\.rds$", full.names = TRUE)

analyze_file <- function(file) {
  dat <- readRDS(file)
  
  # 2. Signature columns
  sig_cols <- grep("^SBS", names(dat), value = TRUE)
  
  # --- keep counts ---
  dat_counts <- dat %>%
    select(all_of(c(sig_cols, "ancestry", "tissue_type"))) %>%
    mutate(file = basename(file),
           tumor_type = unique(dat$tissue_type))
  
  # --- convert to fractions ---
  dat_frac <- dat %>%
    mutate(row_sum = rowSums(across(all_of(sig_cols)))) %>%
    mutate(across(all_of(sig_cols), ~ .x / row_sum)) %>%
    select(-row_sum) %>%
    mutate(file = basename(file),
           tumor_type = unique(dat$tissue_type))
  
  # Get which ancestries are present
  ancestries <- unique(dat$ancestry)
  if (length(ancestries) != 2) {
    return(list(results = NULL, raw_frac = dat_frac, raw_counts = dat_counts))
  }
  
  # 3a. Wilcox tests on FRACTIONS
  results_frac <- map_dfr(sig_cols, function(sig) {
    test <- wilcox.test(dat_frac[[sig]] ~ dat_frac$ancestry)
    tibble(
      signature = sig,
      ancestry1 = ancestries[1],
      ancestry2 = ancestries[2],
      p_value   = test$p.value,
      mean1     = mean(dat_frac[[sig]][dat_frac$ancestry == ancestries[1]]),
      mean2     = mean(dat_frac[[sig]][dat_frac$ancestry == ancestries[2]]),
      sd1       = sd(dat_frac[[sig]][dat_frac$ancestry == ancestries[1]]),
      sd2       = sd(dat_frac[[sig]][dat_frac$ancestry == ancestries[2]]),
      n1        = length(dat_frac[[sig]][dat_frac$ancestry == ancestries[1]]),
      n2        = length(dat_frac[[sig]][dat_frac$ancestry == ancestries[2]]),
      file      = basename(file),
      tumor_type = unique(dat$tissue_type),
      measure   = "fraction"
    )
  })
  # Add Bonferroni correction within this comparison (per file)
  results_frac <- results_frac %>%
    mutate(
      p_value_Bonf = p.adjust(p_value, method = "bonferroni"),
      n_tests = n()  # keep number of signatures for reference
    )
  
  # 3b. Wilcox tests on COUNTS
  results_count <- map_dfr(sig_cols, function(sig) {
    test <- wilcox.test(dat_counts[[sig]] ~ dat_counts$ancestry)
    tibble(
      signature = sig,
      ancestry1 = ancestries[1],
      ancestry2 = ancestries[2],
      p_value   = test$p.value,
      mean1     = mean(dat_counts[[sig]][dat_counts$ancestry == ancestries[1]]),
      mean2     = mean(dat_counts[[sig]][dat_counts$ancestry == ancestries[2]]),
      sd1     = sd(dat_counts[[sig]][dat_counts$ancestry == ancestries[1]]),
      sd2     = sd(dat_counts[[sig]][dat_counts$ancestry == ancestries[2]]),
      n1        = length(dat_counts[[sig]][dat_frac$ancestry == ancestries[1]]),
      n2        = length(dat_counts[[sig]][dat_frac$ancestry == ancestries[2]]),
      file      = basename(file),
      tumor_type = unique(dat$tissue_type),
      measure   = "count"
    )
  })
  # Add Bonferroni correction within this comparison (per file)
  results_count <- results_count %>%
    mutate(
      p_value_Bonf = p.adjust(p_value, method = "bonferroni"),
      n_tests = n()
    )
  
  # merge into one table
  results <- bind_rows(results_frac, results_count)
  
  list(results = results, raw_frac = dat_frac, raw_counts = dat_counts)
}

# 4. Apply to all files
all_outputs <- map(files, analyze_file)

# Big combined results table
all_results <- map_dfr(all_outputs, "results")

# Raw per-sample data
all_raw_frac  <- map_dfr(all_outputs, "raw_frac")
all_raw_count <- map_dfr(all_outputs, "raw_counts")

tcga_test_count <- filter(all_results, measure == "fraction", !is.na(p_value)) %>% nrow()

# 5. just top hits
top_hits <- all_results %>%
  filter(!is.na(p_value), is.finite(p_value), measure == "fraction") %>%   # remove NA/NaN p-values
  group_by(file) %>%
  slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

###########
## MAIN PLOT
  ancestry_amsd_output <- readRDS("../outputs/ancestry_amsd_output.rds")
  full_join(ancestry_amsd_output, top_hits) %>%
    #filter(min_anc_n > 10) %>%
    #ggplot(aes(x = -log10(pvalues), y = -log10(p_value), color = comparison, label = tumor_type, size = min_anc_n))+
    ggplot(aes(x = -log10(pvalues), y = -log10(p_value), color = comparison, label = tumor_type))+
    geom_point()+
    geom_smooth(method = "lm",
                inherit.aes = FALSE,
                aes(-log10(pvalues),
                    -log10(p_value)),
                color = "black")+
    theme_classic()+
    geom_label_repel() +
    geom_hline(yintercept = -log10(0.05))+
    geom_hline(yintercept = -log10(0.05/tcga_test_count))+
    geom_vline(xintercept = -log10(0.05))+
    geom_vline(xintercept = -log10(0.05/67))+
    xlab("AMSD p-value")+
    ylab("Top signature p-value \n(wilcoxen rank sum)")
  
  tcga_amsdsig_plot <- full_join(ancestry_amsd_output, top_hits) %>%
    #filter(min_anc_n > 10) %>%
    #ggplot(aes(x = -log10(pvalues), y = -log10(p_value), color = comparison, label = tumor_type, size = min_anc_n))+
    ggplot(aes(x = -log10(pvalues), y = -log10(p_value_Bonf), color = comparison, label = tumor_type))+
    geom_point()+
    geom_smooth(method = "lm",
                inherit.aes = FALSE,
                aes(-log10(pvalues),
                    -log10(p_value_Bonf)),
                color = "black")+
    theme_classic()+
    geom_label_repel() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    geom_hline(yintercept = -log10(0.05))+
    geom_hline(yintercept = -log10(0.05/67))+
    geom_vline(xintercept = -log10(0.05))+
    geom_vline(xintercept = -log10(0.05/67))+
    xlab("AMSD p-value")+
    ylab("Top signature p-value \n(wilcoxen rank sum)")
  
  ggarrange(mouse_amsdsig_plot, tcga_amsdsig_plot, nrow=2, labels = c("A", "B"))
  
  allplots <- ggarrange(mouse_amsdsig_plot, 
                        tcga_amsdsig_plot, 
                        mouse_plot, 
                        anc_plot, 
                        nrow=2, 
                        ncol=2,  
                        labels = c("A", "B", "C", "D"))
  allplots
  
  ggsave("../outputs/amsd_v_sigs_weighted_v_unweighted.png",
         plot = allplots,
         width = 14,
         height = 14,
         units = "in"
  )
###########
  
  # Prepare the data
  esca_data <- all_raw_frac %>%
    filter(tumor_type == "ESCA") %>% 
    select(-file) %>%
    distinct() %>%
    select(ancestry, SBS46, SBS88) %>%
    pivot_longer(cols = c("SBS46", "SBS88"), names_to = "exp")
  
  # Compute N per ancestry × exposure
  n_counts <- esca_data %>%
    group_by(ancestry, exp) %>%
    summarise(N = n(), .groups = "drop")
  
  # Plot with N labels
  ggplot(esca_data, aes(x = ancestry, y = value, color = exp)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.2, alpha = 0.6) +
    geom_text(
      data = n_counts,
      aes(x = ancestry, y = 0, label = paste0("N=", N), color = exp),
      position = position_dodge(width = 0.75),
      vjust = 1.5,
      size = 3,
      show.legend = FALSE
    ) +
    theme_minimal()
  
  
  all_raw_frac  %>%
    filter(tumor_type == "ESCA") %>% 
    select(-file) %>%
    unique() %>%
    select(ancestry, SBS46, SBS88) %>%
    pivot_longer(cols = c("SBS46", "SBS88"), names_to = "exp") %>%
    ggplot(aes(x=ancestry,y=value, color = exp))+
    geom_boxplot(outliers = FALSE)+
    geom_jitter(height = 0, width = 0.2, alpha = 0.6)
  
  all_raw_frac  %>%
    filter(tumor_type == "HNSC") %>% 
    select(-file) %>%
    unique() %>%
    ggplot(aes(x=ancestry,y=SBS56, color = tumor_type))+
    geom_boxplot(outliers = FALSE)+
    geom_jitter(height = 0, width = 0.2, alpha = 0.6)
  all_raw_frac  %>%
    filter(tumor_type == "HNSC") %>% 
    select(-file) %>%
    unique() %>%
    ggplot(aes(x=ancestry,y=SBS29, color = tumor_type))+
    geom_boxplot(outliers = FALSE)+
    geom_jitter(height = 0, width = 0.2, alpha = 0.6)
  
all_raw_frac  %>%
  filter(tumor_type == "KIRP") %>% 
  select(-file) %>%
  unique() %>%
  ggplot(aes(x=ancestry,y=SBS85, color = tumor_type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.6)
all_raw_frac  %>%
  select(-file) %>%
  unique() %>%
  filter(tumor_type == "SARC") %>%
  ggplot(aes(x=ancestry,y=SBS54, color = tumor_type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.6)
all_raw_frac  %>%
  select(-file) %>%
  unique() %>%
  filter(tumor_type == "OV") %>%
  ggplot(aes(x=ancestry,y=SBS38, color = tumor_type))+
  geom_boxplot(outliers = FALSE)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.6)

all_raw_frac %>%
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


# plot_sig_volcano <- function(anc1, anc2, tumor){
#   p <- filter(ancestry_amsd_output,
#                    tumor_type == tumor,
#                    ancestry1 == anc1,
#                    ancestry2 == anc2) %>%
#   pull(pvalues)
#   subset_tcga <- all_results %>%
#     mutate(ancestries = paste0(ancestry1, "_v_", ancestry2)) %>%
#     filter(tumor_type == tumor, 
#            ancestries == paste0(anc1, "_v_", anc2), 
#            !is.na(p_value))
#   subset_tcga %>%
#     ggplot(aes(x = (mean1 - mean2), y = -log10(p_value), label = signature))+
#     geom_point()+
#     geom_hline(yintercept = -log10(0.05))+
#     geom_hline(yintercept = -log10(0.05/nrow(subset_tcga)))+
#     geom_hline(yintercept = -log10(0.05/tcga_test_count))+
#     geom_label_repel()+
#     labs(title = paste0(tumor,": ",anc1, " vs ", anc2, ", AMSD pvalue = ", p),
#          x = "Difference in means")
# }
plot_sig_volcano <- function(anc1, anc2, tumor, measure = c("fraction", "count")) {
  measure <- match.arg(measure)
  
  # AMSD p-value for context
  p <- ancestry_amsd_output %>%
    filter(
      tumor_type == tumor,
      ancestry1 == anc1,
      ancestry2 == anc2
    ) %>%
    pull(pvalues)
  
  # Subset results by ancestry + tumor + measure
  subset_tcga <- all_results %>%
    mutate(ancestries = paste0(ancestry1, "_v_", ancestry2)) %>%
    filter(
      tumor_type == tumor,
      ancestries == paste0(anc1, "_v_", anc2),
      measure == !!measure,         # key filter: only use chosen measure
      !is.na(p_value)
    )
  
  # Number of tests for Bonferroni
  n_tests <- nrow(subset_tcga)
  
  # Plot
  ggplot(subset_tcga, aes(x = mean1 - mean2, y = -log10(p_value), label = signature)) +
    geom_point(color = ifelse(measure == "fraction", "darkorange", "steelblue")) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05 / n_tests), color = "blue", linetype = "dashed") +
    geom_label_repel(max.overlaps = 20) +
    labs(
      title = paste0(
        tumor, ": ", anc1, " vs ", anc2,
        ", AMSD p = ", signif(p, 3),
        " [", measure, "]"
      ),
      x = "Difference in means",
      y = "-log10(p-value)"
    )
}

plot_sig_volcano("afr", "eas","LUAD", measure = "count")
plot_sig_volcano("afr", "eas","LUAD", measure = "fraction")
plot_sig_volcano("eas", "eur","LUAD", measure = "count")
plot_sig_volcano("eas", "eur","LUAD", measure = "fraction")


plot_sig_volcano("eas", "eur","BLCA", measure = "count")
plot_sig_volcano("eas", "eur","BLCA", measure = "fraction")

plot_sig_volcano("eas", "eur","UCEC")
plot_sig_volcano("eas", "eur","KIRP")
plot_sig_volcano("eas", "eur","SARC")
plot_sig_volcano("afr", "eur","ESCA")
plot_sig_volcano("eas", "eur","HNSC")
plot_sig_volcano("afr", "eur","HNSC")

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


################# Looking comparison plots again
plot_sig_means <- function(anc1, anc2, tumor, measure = c("fraction", "count")) {
  measure <- match.arg(measure)
  
  # Subset results
  subset_tcga <- all_results %>%
    mutate(ancestries = paste0(ancestry1, "_v_", ancestry2)) %>%
    filter(
      tumor_type == tumor,
      ancestries == paste0(anc1, "_v_", anc2),
      measure == !!measure,
      !is.na(p_value)
    )
  
  ggplot(subset_tcga, aes(x = mean1, y = mean2, label = signature)) +
    geom_point(color = ifelse(measure == "fraction", "darkorange", "steelblue")) +
    # Error bars for ancestry1 (x-axis)
    geom_errorbarh(aes(xmin = mean1 - sd1/sqrt(n1), xmax = mean1 + sd1/sqrt(n1)),
                   height = 0, alpha = 0.5) +
    # Error bars for ancestry2 (y-axis)
    geom_errorbar(aes(ymin = mean2 - sd2/sqrt(n2), ymax = mean2 + sd2/sqrt(n2)),
                  width = 0, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_label_repel(max.overlaps = 20) +
    # ggrepel::geom_label_repel(
    #   data = subset_tcga %>%
    #     filter(p_value < 0.05,
    #            abs(mean1 - mean2) > 0.01,
    #            mean1/mean2 > 1.1 | mean2/mean1 > 1.1),
    #   aes(label = signature),
    #   max.overlaps = 20
    # ) +
    labs(
      title = paste0(
        tumor, ": ", anc1, " vs ", anc2, " (", measure, ")"
      ),
      x = paste0("Mean in ", anc1),
      y = paste0("Mean in ", anc2)
    ) +
    theme_minimal()
}
plot_sig_means("afr", "eas", "LUAD", measure = "fraction")
plot_sig_means("afr", "eas", "LUAD", measure = "count")
plot_sig_means("eas", "eur", "LUAD", measure = "fraction")
plot_sig_means("eas", "eur", "LUAD", measure = "count")
plot_sig_means("eas", "eur", "ESCA", measure = "fraction")
plot_sig_means("eas", "eur", "ESCA", measure = "count")


plot_sig_means("eas", "eur", "UCEC", measure = "fraction")
plot_sig_means("eas", "eur", "UCEC", measure = "count")
plot_sig_means("eas", "eur", "KIRP", measure = "fraction")
plot_sig_means("eas", "eur", "KIRP", measure = "count")

pdf("../outputs/tcga_sigantures_comparison_plots.pdf")
for(i in 1:nrow(ancestry_amsd_output)){
  sample <- full_join(ancestry_amsd_output, top_hits) %>%
    arrange(pvalues)
  plot1 <- plot_sig_means(sample[i,"ancestry1"], sample[i,"ancestry2"],sample[i,"tumor_type"], measure = "fraction")
  print(plot1)
  plot2 <- plot_sig_means(sample[i,"ancestry1"], sample[i,"ancestry2"],sample[i,"tumor_type"], measure = "count")
  print(plot2)
}
dev.off()


# Function to make multi-panel plots for top hits
plot_top_hits <- function(df, p_cutoff = 0.05, measure_type = "fraction") {
  top_hits <- df %>%
    filter(pvalues < p_cutoff) %>%
    arrange(pvalues)
  
  # Generate plots
  plots <- map(seq_len(nrow(top_hits)), function(i) {
    row <- top_hits[i, ]
    plot_sig_means(row$ancestry1, row$ancestry2, row$tumor_type, measure_type)
  })
  
  # Combine into a grid
  wrap_plots(plots, ncol = 5)
}
frac_plot <- plot_top_hits(ancestry_amsd_output, p_cutoff = 0.01, measure_type = "fraction")
count_plot <- plot_top_hits(ancestry_amsd_output, p_cutoff = 0.01, measure_type = "count")

frac_count_plots <- ggarrange(frac_plot, 
                              count_plot, 
                              nrow=2, 
                              ncol=1,  
                              labels = c("A", "B"))
frac_count_plots

ggsave("../outputs/signature_comparisons_sigAMSD.png",
       plot = frac_count_plots,
       width = 14,
       height = 14,
       units = "in"
)

#############################################################
##############################
# remake figure 3 with new sigprofiler fitting:
# Read the first 67 lines
lines <- readLines("amsd_tcga_plotting.R", n = 67)
# Evaluate those lines
eval(parse(text = lines))
ancestry_volcano


# new plotting function

plot_sig_means_clean <- function(anc1, anc2, tumor, tumor_name, signatures) {

  # Subset results
  subset_tcga <- all_results %>%
    mutate(ancestries = paste0(ancestry1, "_v_", ancestry2)) %>%
    filter(
      tumor_type == tumor,
      ancestries == paste0(anc1, "_v_", anc2),
      measure == "fraction",
      !is.na(p_value)
    )
  
  ggplot(subset_tcga, aes(x = mean1, y = mean2, label = signature)) +
    geom_point() +
    # Error bars for ancestry1 (x-axis)
    geom_errorbarh(aes(xmin = mean1 - sd1/sqrt(n1), xmax = mean1 + sd1/sqrt(n1)),
                   height = 0, alpha = 0.5) +
    # Error bars for ancestry2 (y-axis)
    geom_errorbar(aes(ymin = mean2 - sd2/sqrt(n2), ymax = mean2 + sd2/sqrt(n2)),
                  width = 0, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    # geom_label_repel(max.overlaps = 20) +
    ggrepel::geom_label_repel(
      data = subset_tcga %>%
        filter(signature %in% signatures),
      aes(label = signature),
      avoid_points = TRUE,   # <—
      box.padding = 0.5,
      point.padding = 1,
      force = 1,
      force_pull = 0.1,
      max.overlaps = 50,
      min.segment.length = 0,
      segment.color = "gray60"
    ) +
    labs(
      title = tumor_name,
      x = toupper(anc1),
      y = toupper(anc2)
    ) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))
}

esca <- plot_sig_means_clean("eas", "eur", "ESCA", "Esophageal", c("SBS1", "SBS5", "SBS16","SBS17a", "SBS17b"))+
  coord_flip()

lihc <- plot_sig_means_clean("eas", "eur", "LIHC", "Liver", c("SBS5", "SBS16","SBS22", "SBS24"))+
  coord_flip()
# plot_sig_means_clean("eas", "eur", "LIHC", "Liver", c("SBS1","SBS5", "SBS16","SBS22", "SBS24"))+
#   xlim(0,0.2)+ ylim(0,0.2)+ coord_flip()+ 
#   annotate("segment", x = 0.175, y = 0.18, xend = 0.195, yend = 0.2,
#            arrow = arrow(length = unit(0.02, "npc"))) + 
#   annotate("label", x = 0.175, y = 0.18, label = "SBS5", hjust = 0.5, fill = "white") 
  
lua1 <- plot_sig_means_clean("afr", "eas", "LUAD", "Lung Adeno.", c("SBS1", "SBS5", "SBS4", "SBS40"))+
  coord_flip()
lua2 <- plot_sig_means_clean("afr", "eur", "LUAD", "Lung Adeno.", c("SBS5", "SBS4", "SBS40"))+
  coord_flip()

uce1 <- plot_sig_means_clean("eas", "eur", "UCEC", "Uterine", c("SBS1", "SBS5", "SBS10a", "SBS10b"))+
  coord_flip()
uce2 <- plot_sig_means_clean("afr", "eas", "UCEC", "Uterine", c("SBS1", "SBS5", "SBS10a", "SBS10b"))
coa1 <- plot_sig_means_clean("eas", "eur", "COAD", "Colorectal", c("SBS1", "SBS5", "SBS10a", "SBS10b"))+
  coord_flip()
coa1
coa2 <- plot_sig_means_clean("afr", "eas", "COAD", "Colorectal", c("SBS1", "SBS5", "SBS10a", "SBS10b"))
coa2
sig_plots <- ggarrange(esca, lihc, lua1, lua2, uce1, uce2, coa1, coa2, nrow = 2, ncol = 4)
  

fig3 <- ggarrange(ancestry_volcano,
                  sig_plots,
                  # ggarrange(NULL,
                  #           sig_plots,
                  #           NULL,
                  #           ncol = 3,
                  #           widths = c(0.025,0.95,0.025)),
                  nrow = 2,
                  heights = c(0.6,0.4),
                  labels = c("A", "B"))

ggsave("../outputs/Figure3_revision.png",
       plot = fig3,
       width = 10,
       height = 10,
       units = "in")
ggsave("../outputs/Figure3_revision.svg",
       plot = fig3,
       width = 10,
       height = 10,
       units = "in")
