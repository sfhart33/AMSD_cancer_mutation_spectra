library(tidyverse)
source("amsd_functions.R")

# load data from mouse paper (processed into spectra)
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

  
# Run for all options
  carcinogen_table <- sample_table %>%
    select(tissue, exposure) %>%
    filter(tissue %in% c("LIVER", "LUNG"), 
           exposure != "SPONTANEOUS") %>%
    count(tissue, exposure) %>%
    arrange(tissue, desc(n))

# Blank variables for outputs
  pvalues <- c()
  cosines <- c()
  reps <- 100000
  perms <- data.frame(rep = 1:reps)
  
# run for each carcinogen
  for(count in 1:nrow(carcinogen_table)){
    
    # set variables
    tissue_type <- carcinogen_table[count,"tissue"]
    exp_type <- carcinogen_table[count,"exposure"]
    
    # Sample groupings
    samples_exp <- filter(sample_table, tissue == tissue_type, exposure == exp_type) %>% pull(label)
    samples_nonexp <- filter(sample_table, tissue == tissue_type, exposure == "SPONTANEOUS") %>% pull(label)
    
    # Run AMSD
    amsd_output <- amsd(mouse_carcinogen_spectra[samples_exp,],
         mouse_carcinogen_spectra[samples_nonexp,],
         mean_or_sum = "mean",
         n_sim = reps,
         seed = 1234)
    
    # Save output
    pvalues <- c(pvalues, amsd_output$p)
    cosines <- c(cosines, amsd_output$cosine)
    perms1 <- data.frame(perms = amsd_output$sims)
    colnames(perms1) <- paste0(tissue_type, ".", exp_type)
    perms <- cbind(perms, perms1)
    print(paste(count, "of", nrow(carcinogen_table)))
    
  }
  
# save outputs
  mouse_amsd_output <- carcinogen_table
  mouse_amsd_output$pvalues <- pvalues
  mouse_amsd_output$cosines <- cosines
  
  saveRDS(mouse_amsd_output, "../outputs/mouse_amsd_output.rds")
  saveRDS(perms, "../outputs/mouse_amsd_perms.rds")
  mouse_amsd_output <- readRDS("../outputs/mouse_amsd_output.rds")
  
  pdf("../outputs/mouse_amsd_output.pdf")
  # volcano plot summary of everything together
  mutate(mouse_amsd_output, log10pval = -log10(pvalues)) %>%
    ggplot()+
      geom_point(aes(x=cosine_dif, y = log10pval, color = tissue, size = n))+
      geom_hline(yintercept = -log10(0.05))+
      geom_hline(yintercept = -log10(0.05/nrow(mouse_amsd_output)))+
      geom_hline(yintercept = -log10(1/reps))+
      geom_text(aes(x=0.225, y = (-log10(0.05/nrow(mouse_amsd_output))+0.1)), label = "FDR=0.05")+
      geom_text(aes(x=0.225, y = (-log10(0.05)+0.1)), label = "p=0.05")+
      geom_text(aes(x=0.225, y = (-log10(1/reps)+0.1)), label = "theoretical max")+
      theme_classic()+
      xlim(0,0.25)+
      scale_size_continuous(range = c(3, 5))+
      xlab("Cosine distance")+
      ylab("-log10(p-value)")
  dev.off()