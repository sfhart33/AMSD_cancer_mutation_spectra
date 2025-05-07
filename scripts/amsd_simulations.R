library(tidyverse)
library(sigfit)
  data("cosmic_signatures_v3.2")
setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github/AMSD_cancer_mutation_spectra/scripts")
source("amsd_functions.R")

  
  
# test various situations
  set.seed(123)
  results <- list()
  rep <- 0
  seed = NULL
  n_sim = 1000
  sig_probs <- c(SBS1 = 0.3, SBS5 = 0.6, SBS18 = 0.1)
  signatures <- cosmic_signatures_v3.2
  ss <- c("SBS2","SBS40") # spiky and flat signatures to test
  ns <- c(5,25,125,625) # samples per exposure to test
  ms <- c(50,2500) # mutations per sample to test (approx WES and WGS)
  fs <- c(0.2,0.1,0.05,0.02)# extra mutations - fraction of total mutations
  xs <- 1:100 # replicates to test each parameter set
  total_tests <- length(ss)*length(ns)*length(ms)*length(fs)*length(xs)
  print(paste(rep,"of",total_tests, "tests"))
  for(s in ss){ 
    for(n in ns){ 
      for(m in ms){ 
        for(f in fs){ 
          for(x in xs){
            rep <- rep+1
            # set parameters
            n_samples <- n
            n_mutations <- m
            additional_sig <- s
            n_extra = rep(n_mutations*f,n_samples)
            # Run simluation
            p <- simulate_spectra_amsd(
              n_samples = n_samples,
              n_mutations = n_mutations,
              sig_probs = sig_probs,
              signatures = signatures,
              additional_sig = additional_sig,
              n_extra = n_extra,
              n_sim = n_sim,
              seed = seed)$p
            # Print to output list and update progress
            print(paste(rep,"of",total_tests, "tests"))
            #print(c(n,m,s,p,f))
            results[[rep]] <- c(n,m,s,p,f)
          }
        }
      } 
    } 
  }
  output <- do.call(rbind, results) %>%
    as.data.frame()
  colnames(output) <- c("n_samples","n_mutations","exposure","pvalue", "extra_muts")
  output <- output %>%
    mutate(group = paste0(extra_muts,"% ",exposure))
  
  saveRDS(output, "../outputs/amsd_simulation_output.rds")
  output <- readRDS("../outputs/amsd_simulation_output.rds")
  
  output
  output2 <- output %>%
    group_by(n_samples, n_mutations, exposure, extra_muts) %>%
    summarize(mean_p = mean(as.numeric(pvalue)),
              sd_p = sd(as.numeric(pvalue)),
              success_frac = sum(pvalue <= 0.05)/n())
  output2

  simulation_plot <- output2 %>%
    ggplot(aes(x = factor(n_samples, levels = c("5","25","125","625")),
               y = success_frac,
               color = factor(extra_muts, levels = c("0.02","0.05","0.1","0.2")),
               group = factor(extra_muts, levels = c("0.02","0.05","0.1","0.2"))))+
    geom_point()+
    geom_line() +
    facet_grid(n_mutations ~ exposure) +
    guides(color = guide_legend(title = "Extra mutations per \nexposure sample (%)"))+
    xlab("Sample count (same # exposed and non-exposed)")+
    ylab("Difference detected \n(p<0.05, fraction of 100 simulations)")+
    # ggtitle("AMSD strength of detection in \nWES (50 mu/sample) and WGS (2500 mu/sample)")+
    theme_classic()
  simulation_plot
  ggsave("../outputs/amsd_simulations.png",
         plot = simulation_plot,
         width = 7,
         height = 5,
         units = "in"
         )
  