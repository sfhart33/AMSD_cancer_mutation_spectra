# load libraries
  library(tidyverse)

# cosine distance
cosine_dist <- function(A, B){
  cosine_similarity <- sum(A * B) / (sqrt(sum(A^2)) * sqrt(sum(B^2)))
  cosine_distance <- 1 - cosine_similarity
  return(cosine_distance)
}
  
# AMSD function
amsd <- function(set1,
                 set2,
                 mean_or_sum = "mean",  # or "sum"
                 n_sim = 1000,
                 seed = NULL) {   
  
  # Validate inputs
  if (!is.matrix(set1) && !is.data.frame(set1)) stop("set1 must be a matrix or data.frame")
  if (!is.matrix(set2) && !is.data.frame(set2)) stop("set2 must be a matrix or data.frame")
  
  if (!mean_or_sum %in% c("mean", "sum")) {
    stop("Argument 'mean_or_sum' must be either 'mean' or 'sum'")
  }
  
  # Define aggregation function
  aggragate_spectra <- if (mean_or_sum == "mean") {
    function(data) colMeans(data, na.rm = TRUE)
  } else {
    function(data) colSums(data, na.rm = TRUE)
  }
  
  
  # Compute observed distance
  spectra1 <- aggragate_spectra(set1)
  spectra2 <- aggragate_spectra(set2)
  observed_distance <- cosine_dist(spectra1, spectra2)[[1]]
  
  # Prepare permutation dataset
  combined_set <- rbind(set1, set2)

  # Run permutations
  set.seed(seed)
  dist_rands <- numeric(n_sim)
  group_size <- nrow(set1)
  
  for (k in seq_len(n_sim)) {
    indices <- sample(seq_len(nrow(combined_set)), group_size)
    spectra_group1 <- aggragate_spectra(combined_set[indices, , drop = FALSE])
    spectra_group2 <- aggragate_spectra(combined_set[-indices, , drop = FALSE])
    dist_rands[k] <- cosine_dist(spectra_group1, spectra_group2)[[1]]
  }
  
  # Calculate p-value
  p_value <- max(c(mean(dist_rands >= observed_distance), 1 / n_sim))
  
  # Return output
  return(list(cosine = observed_distance, p = p_value, sims = dist_rands))
}
    
# plot amsd output function
plot_amsd_histogram <- function(output){
  plot1 <- data.frame(sims = output$sims) %>%
    ggplot(aes(sims)) +
    geom_histogram()+
    geom_vline(xintercept = output$cosine, linetype = "dashed")+
    xlab("Cosine distance")+
    ggtitle(paste0("p = ", output$p)) +
    theme_classic()
  return(plot1)
}
    
    
    
# simulate spectra function
simulate_spectra <- function(n_samples,
                             n_mutations,
                             sig_probs,
                             signatures,
                             additional_sig = NULL,
                             n_extra = 0,
                             seed = NULL) {
  set.seed(seed)
  
# loop through samples and save 3mer counts
  spectras <- list()
  for (i in 1:n_samples) {
    
  # Loop through each signature and sample 3mer counts
    spectra <- c()
    for(s in 1:length(sig_probs)){
      
    # Which signature
      sig <- names(sig_probs)[s]
      sig_spectra <- signatures[sig,]
      
    # Sample mutations from given signatures
      mutations <- sample(names(sig_spectra),
                          size = (n_mutations) * sig_probs[sig],
                          replace = TRUE,
                          prob = sig_spectra)
      
      spectra[[s]] <- table(factor(mutations, levels = names(sig_spectra)))
    }
    
  # Add additional mutations if applicable
    if (!is.null(additional_sig) && n_extra[i] > 0) {
      extra_mutations <- sample(names(signatures[additional_sig,]),
                                size = n_extra[i],
                                replace = TRUE,
                                prob = signatures[additional_sig,])
      spectra[[length(sig_probs)+1]] <- table(factor(extra_mutations, levels = names(sig_spectra)))
    } 
    
    spectras[[i]] <- colSums(do.call(rbind, spectra))
  }
  return(spectras)
}

# simulate spectra and amsd together
simulate_spectra_amsd <- function(n_samples, n_mutations, sig_probs, signatures, additional_sig, n_extra, n_sim, seed) {
  
  # Run on a base set, then with exposures
  no_exposure_test <- simulate_spectra(n_samples = n_samples,
                                       n_mutations = n_mutations,
                                       sig_probs = sig_probs,
                                       signatures = signatures)
  no_exposure_test <- as.data.frame(do.call(rbind, no_exposure_test))
  with_exposure_test <- simulate_spectra(n_samples = n_samples,
                                         n_mutations = n_mutations,
                                         sig_probs = sig_probs,
                                         signatures = signatures,
                                         additional_sig = additional_sig,
                                         n_extra = n_extra)
  with_exposure_test <- as.data.frame(do.call(rbind, with_exposure_test))
  amsd(no_exposure_test, with_exposure_test, n_sim = n_sim, seed = seed) %>%
    return()
}
    
# # test out
#     library(sigfit) # For plotting spectra and cosmic signatures
#     data(cosmic_signatures_v3)
# 
#     # set parameters
#     sig_probs <- c(SBS1 = 0.2, SBS5 = 0.7, SBS18 = 0.1)
#     n_samples <- 10
#     n_mutations <- 1000
#     additional_sig <- "SBS2"
#     signatures = cosmic_signatures_v3
#     n_extra = rep(100,n_samples)
#     seed = NULL
#     n_sim = 2000
# 
#     test_p <- simulate_spectra_amsd(
#       n_samples = n_samples,
#       n_mutations = n_mutations,
#       sig_probs = sig_probs,
#       additional_sig = additional_sig,
#       signatures = cosmic_signatures_v3,
#       n_extra = n_extra,
#       n_sim = n_sim,
#       seed = seed)
#     test_p$p
# 
#     plot_amsd_histogram(test_p)
    