library(tidyverse)
library(mutspecdist)
library(SigProfilerAssignmentR)
library(furrr)
library(future)

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

# #--- main comparison function
# compare_spectra_sigprofiler <- function(n_samples = 5,
#                                         n_mutations = 50,
#                                         sig_probs = c(SBS1 = 0.3, SBS5 = 0.6, SBS18 = 0.1),
#                                         additional_sig = "SBS2",
#                                         frac_extra = 0.05,
#                                         n_sim = 1000,
#                                         seed = 123,
#                                         cosmic_version = 3.3,
#                                         genome_build = "GRCh38") {
#   set.seed(seed)
#   
#   # --- Simulate spectra ---
#   no_exp <- simulate_spectra(n_samples, n_mutations, sig_probs, cosmic_signatures_v3.2)
#   no_exp <- as.data.frame(do.call(rbind, no_exp))
#   
#   n_extra <- rep(n_mutations * frac_extra, n_samples)
#   with_exp <- simulate_spectra(n_samples, n_mutations, sig_probs, cosmic_signatures_v3.2,
#                                additional_sig = additional_sig, n_extra = n_extra)
#   with_exp <- as.data.frame(do.call(rbind, with_exp))
#   
#   # --- Convert to SigProfiler input ---
#   sp_no  <- convert_sigfit_to_SigProfiler(no_exp)
#   sp_yes <- convert_sigfit_to_SigProfiler(with_exp)
#   
#   # --- Define output directories ---
#   out_no  <- paste0("Assignment_no_", seed)
#   out_yes <- paste0("Assignment_yes_", seed)
#   
#   # --- Fit with SigProfiler ---
#   cosmic_fit(
#     samples = sp_no,
#     output = out_no,
#     input_type = "matrix",
#     context_type = "96",
#     collapse_to_SBS96 = TRUE,
#     cosmic_version = cosmic_version,
#     genome_build = genome_build
#   )
#   
#   cosmic_fit(
#     samples = sp_yes,
#     output = out_yes,
#     input_type = "matrix",
#     context_type = "96",
#     collapse_to_SBS96 = TRUE,
#     cosmic_version = cosmic_version,
#     genome_build = genome_build
#   )
#   
#   # --- Read exposures ---
#   exp_no <- read.delim(file.path(out_no, "Assignment_Solution",
#                                  "Activities", "Assignment_Solution_Activities.txt"),
#                        header = TRUE, sep = "\t", check.names = FALSE)
#   
#   exp_yes <- read.delim(file.path(out_yes, "Assignment_Solution",
#                                   "Activities", "Assignment_Solution_Activities.txt"),
#                         header = TRUE, sep = "\t", check.names = FALSE)
#   
#   exp_no$group <- "A"
#   exp_yes$group <- "B"
#   df <- bind_rows(exp_no, exp_yes)
#   
#   # --- Cleanup output dirs ---
#   unlink(out_no, recursive = TRUE)
#   unlink(out_yes, recursive = TRUE)
#   
#   # --- Statistical tests ---
#   results <- df %>%
#     pivot_longer(
#       cols = -c(Samples, group),   # <- explicitly exclude Samples and group
#       names_to = "signature",
#       values_to = "exposure"
#     ) %>%
#     group_by(signature) %>%
#     summarise(
#       mean_A = mean(exposure[group == "A"]),
#       mean_B = mean(exposure[group == "B"]),
#       p_ttest = t.test(exposure ~ group)$p.value,
#       p_wilcox = wilcox.test(exposure ~ group)$p.value,
#       .groups = "drop"
#     ) %>%
#     mutate(
#       # Multiple testing corrections
#       p_ttest_Bonf = p.adjust(p_ttest, "bonferroni"),
#       p_ttest_BH   = p.adjust(p_ttest, "BH"),
#       p_ttest_BY   = p.adjust(p_ttest, "BY"),
#       p_wilcox_Bonf = p.adjust(p_wilcox, "bonferroni"),
#       p_wilcox_BH   = p.adjust(p_wilcox, "BH"),
#       n_samples = n_samples,
#       n_mutations = n_mutations,
#       additional_sig = additional_sig,
#       frac_extra = frac_extra,
#       seed = seed
#     ) %>%
#     arrange(p_ttest)
#   
#   list(
#     results = results,
#     exposures = df
#   )
# }
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
      additional_sig = additional_sig,
      n_extra = n_extra
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

# compare_spectra_sigprofiler()

run_parameter_grid <- function(
    param_grid,
    n_reps = 5,
    output_file = "top_results.tsv"
) {
  all_results <- list()
  
  run_id <- 1
  
  for (i in seq_len(nrow(param_grid))) {
    # Extract parameter set
    params <- param_grid[i, ]
    
    for (rep in seq_len(n_reps)) {
      cat("Running param set", i, "rep", rep, "...\n")
      
      # Call your modified comparison function
      res <- compare_spectra_sigprofiler_top(
        n_samples     = params$n_samples,
        n_mutations   = params$n_mutations,
        sig_probs     = eval(parse(text = params$sig_probs)),  # if stored as string
        additional_sig= params$additional_sig,
        frac_extra    = params$frac_extra,
        n_sim         = params$n_sim,
        seed          = as.integer(params$seed) + rep # vary seed per replicate
      )
      
      # Tag with metadata
      res$top_result$replicate <- rep
      res$top_result$param_set <- i
      res$top_result$run_id    <- run_id
      
      all_results[[run_id]] <- res$top_result
      run_id <- run_id + 1
    }
  }
  
  # Combine all into a single dataframe
  df_all <- bind_rows(all_results)
  
  # Save to file
  write.table(df_all, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(df_all)
}


run_parameter_grid_parallel <- function(param_grid, n_reps = 5, output_file = "top_results.tsv") {
  plan(multisession, workers = max(1, parallel::detectCores() - 1))
  
  tasks <- expand.grid(
    i = seq_len(nrow(param_grid)),
    rep = seq_len(n_reps),
    KEEP.OUT.ATTRS = FALSE
  )
  
  results <- future_map(
    seq_len(nrow(tasks)),
    function(task_id) {
      tryCatch({
        task <- tasks[task_id, ]
        i <- task$i
        rep <- task$rep
        params <- param_grid[i, ]
        
        message("Running task ", task_id, 
                " (param set ", i, ", rep ", rep, ")")
        
        res <- compare_spectra_sigprofiler_top(
          n_samples     = params$n_samples,
          n_mutations   = params$n_mutations,
          sig_probs     = eval(parse(text = params$sig_probs)),
          additional_sig= params$additional_sig,
          frac_extra    = params$frac_extra,
          n_sim         = params$n_sim,
          seed          = as.integer(params$seed) + rep
        )
        
        res$replicate <- rep
        res$param_set <- i
        res$run_id    <- task_id
        res
      },
      error = function(e) {
        message("❌ Error in task ", task_id, ": ", e$message)
        return(NULL)
      })
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )
  
  df_all <- bind_rows(Filter(Negate(is.null), results))
  write.table(df_all, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  return(df_all)
}


# Define a parameter grid (as a data frame)
param_grid <- expand.grid(
  n_samples = c(10,10),
  n_mutations = c(50,500),
  sig_probs = c("c(SBS1 = 0.3, SBS5 = 0.6, SBS18 = 0.1)"),
  additional_sig = c("SBS2", "SBS40"),
  frac_extra = c(0.1,0.2),
  n_sim = 1000,
  seed = 123
)

# Run
df_all <- run_parameter_grid_parallel(param_grid, n_reps = 3, output_file = "top_results.tsv")
df_all
