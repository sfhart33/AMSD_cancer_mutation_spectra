library(tidyverse)
library(ggpubr)
setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github/AMSD_cancer_mutation_spectra/scripts")
source("amsd_functions.R")

# load data
  samples <- read.table("../inputs/asbestos_sample_data.tsv", sep = "\t", header = TRUE)
  asbestos_sbs_spectra <- readRDS("../inputs/asbestos_sbs_spectra.rds") %>%
    rownames_to_column(var = "sample") %>%
    separate(sample, c("sample","seq"), sep = "\\.") %>% # drop extra label
    # filter(seq == "TN") %>%
    select(-seq) %>%
    column_to_rownames(var = "sample")
  CNV_matrix <- read.table("../inputs/asbestos_CNVs.CNV48.matrix.tsv", sep = "\t", header = TRUE) %>%
    column_to_rownames(var = "MutationType")
  SV_matrix <- read.table("../inputs/asbestos_SVs.SV32.matrix.tsv", sep = "\t", header = TRUE) %>%
    column_to_rownames(var = "MutationType")

samples
asbestos_sbs_spectra 
CNV_matrix
SV_matrix

# which are exposed
  samples %>%
    count(Professional.Asbestos)
  samples_exp <- filter(samples, Professional.Asbestos == "Exposed") %>%
    separate(Sample, c("Sample", NA, NA), sep = "_T") %>%
    pull(Sample)
  samples_nonexp <- filter(samples, Professional.Asbestos == "Non exposed") %>%
    separate(Sample, c("Sample", NA, NA), sep = "_T") %>%
    pull(Sample)
  
# run AMSD on SBS
  reps = 10000
  spectra <- asbestos_sbs_spectra/rowSums(asbestos_sbs_spectra)
  exp_spectra <- spectra[samples_exp[samples_exp %in% rownames(spectra)],]
  nonexp_spectra <- spectra[samples_nonexp[samples_nonexp %in% rownames(spectra)],]
  amsd_output_sbs_mean <- amsd(exp_spectra,
                      nonexp_spectra,
                      mean_or_sum = "mean",
                      n_sim = reps,
                      seed = 123)
  plot_amsd_histogram(amsd_output_sbs_mean)

  spectra <- asbestos_sbs_spectra
  exp_spectra <- spectra[samples_exp[samples_exp %in% rownames(spectra)],]
  nonexp_spectra <- spectra[samples_nonexp[samples_nonexp %in% rownames(spectra)],]
  amsd_output_sbs_sum <- amsd(exp_spectra,
                      nonexp_spectra,
                      mean_or_sum = "sum",
                      n_sim = reps,
                      seed = 123)
  plot_amsd_histogram(amsd_output_sbs_sum)
  
# AMSD on CNV (weighted equally)
  spectra <- CNV_matrix %>%
    apply(2,function(x){x/sum(x)}) %>% # for freq not counts
    t()
  samples_exp <- filter(samples, Professional.Asbestos == "Exposed") %>%
    separate(Sample, into = c(NA, NA, "replicate"), sep = "_", remove = FALSE) %>%
    filter(Sample %in% rownames(spectra)) %>%
    filter(replicate != "T2") %>% # don't replicate samples with two 
    pull(Sample)
  samples_nonexp <- filter(samples, Professional.Asbestos == "Non exposed") %>%
    separate(Sample, into = c(NA, NA, "replicate"), sep = "_", remove = FALSE) %>%
    filter(Sample %in% rownames(spectra)) %>%
    filter(replicate != "T2") %>% # don't replicate samples with two 
    pull(Sample)
  amsd_output_cnv_mean <- amsd(spectra[samples_exp,],
                      spectra[samples_nonexp,],
                      mean_or_sum = "mean",
                      n_sim = reps,
                      seed = 123)
  plot_amsd_histogram(amsd_output_cnv_mean)
# AMSD on CNV (NOT weighted equally)
  spectra <- CNV_matrix %>%
    t()
  amsd_output_cnv_sum <- amsd(spectra[samples_exp,],
                      spectra[samples_nonexp,],
                      mean_or_sum = "sum",
                      n_sim = reps,
                      seed = 123)
  plot_amsd_histogram(amsd_output_cnv_sum)
  
  
  total_cnv <- colSums(CNV_matrix) %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample")
  
  exclude_high_samples <- total_cnv %>%
    filter(. > 500) %>%
    pull(sample)
  
  total_cnv %>%
    arrange(desc(.))
  
  counts_exp <- total_cnv %>%
    filter(sample %in% samples_exp) %>%
    mutate(exp = "abs-exposed")
  
  counts_nonexp <- total_cnv %>%
    filter(sample %in% samples_nonexp) %>%
    mutate(exp = "non-exposed")
  
  high_outlier_plot <- rbind(counts_exp, counts_nonexp) %>%
    ggplot(aes(x = ., fill = exp))+
    geom_histogram()+
    xlab("CNV count in each tumor")+
    theme_classic()+
    theme(legend.position = c(0.75, 0.5))
  high_outlier_plot
  
# AMSD on CNV (NOT weighted equally, drop high-CNV outliers)
  spectra <- CNV_matrix %>%
    t()
  samples_exp <- filter(samples, Professional.Asbestos == "Exposed") %>%
    separate(Sample, into = c(NA, NA, "replicate"), sep = "_", remove = FALSE) %>%
    filter(Sample %in% rownames(spectra)) %>%
    filter(replicate != "T2") %>% # don't replicate samples with two 
    filter(!(Sample %in% exclude_high_samples)) %>%
    pull(Sample)
  samples_nonexp <- filter(samples, Professional.Asbestos == "Non exposed") %>%
    separate(Sample, into = c(NA, NA, "replicate"), sep = "_", remove = FALSE) %>%
    filter(Sample %in% rownames(spectra)) %>%
    filter(replicate != "T2") %>% # don't replicate samples with two 
    filter(!(Sample %in% exclude_high_samples)) %>%
    pull(Sample)
  amsd_output_cnv_sum2 <- amsd(spectra[samples_exp,],
                               spectra[samples_nonexp,],
                               mean_or_sum = "sum",
                               n_sim = reps,
                               seed = 123)
  plot_amsd_histogram(amsd_output_cnv_sum2)  
 # AMSD on SV (weighted equally)
  spectra <- SV_matrix %>%
    apply(2,function(x){x/sum(x)}) %>% # for freq not counts
    t()
  samples_exp <- filter(samples, Professional.Asbestos == "Exposed") %>%
    separate(Sample, into = c(NA, NA, "replicate"), sep = "_", remove = FALSE) %>%
    filter(Sample %in% rownames(spectra)) %>%
    filter(replicate != "T2") %>% # don't replicate samples with two 
    pull(Sample)
  samples_nonexp <- filter(samples, Professional.Asbestos == "Non exposed") %>%
    separate(Sample, into = c(NA, NA, "replicate"), sep = "_", remove = FALSE) %>%
    filter(Sample %in% rownames(spectra)) %>%
    filter(replicate != "T2") %>% # don't replicate samples with two 
    pull(Sample)
  amsd_output_sv_mean <- amsd(spectra[samples_exp,],
                               spectra[samples_nonexp,],
                               mean_or_sum = "mean",
                               n_sim = reps,
                               seed = 123)
  plot_amsd_histogram(amsd_output_sv_mean)
# AMSD on SV (NOT weighted equally)
  spectra <- SV_matrix %>%
    t()
  amsd_output_sv_sum <- amsd(spectra[samples_exp,],
                              spectra[samples_nonexp,],
                              mean_or_sum = "sum",
                              n_sim = reps,
                              seed = 123)
  plot_amsd_histogram(amsd_output_sv_sum)
  
# CNV spectra
  spectra1 <- rowSums(CNV_matrix[,samples_exp], na.rm=TRUE)
  spectra2 <- rowSums(CNV_matrix[,samples_nonexp], na.rm=TRUE)
  
  exp_cnv_plot <- as.data.frame(spectra1) %>%
    rownames_to_column(var = "MuType") %>%
    ggplot(aes(MuType, spectra1))+
    geom_col()+
    ylab("count") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle("Prof exp to asbestos (sum, n = 75 tumors)")
  
  non_exp_cnv_plot <- as.data.frame(spectra2) %>%
    rownames_to_column(var = "MuType") %>%
    ggplot(aes(MuType, spectra2))+
    geom_col()+
    ylab("count") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle("NOT prof exp to asbestos (sum, n = 28 tumors)")
  exp_cnv_plot
  non_exp_cnv_plot

# plots
  p1 <- plot_amsd_histogram(amsd_output_sbs_mean) +
    ggtitle("SBS spectra, means") +
    geom_label(x = amsd_output_sbs_mean$cosine, y = 1000, label = paste0("p=",amsd_output_sbs_mean$p), angle=90, hjust = 0)
  p2 <- plot_amsd_histogram(amsd_output_sbs_sum) +
    ggtitle("SBS spectra, sums") +
    geom_label(x = amsd_output_sbs_sum$cosine, y = 1000, label = paste0("p=",amsd_output_sbs_sum$p), angle=90, hjust = 0)
  p3 <- plot_amsd_histogram(amsd_output_cnv_mean) +
    ggtitle("CNV spectra, means") +
    geom_label(x = amsd_output_cnv_mean$cosine, y = 1000, label = paste0("p=",amsd_output_cnv_mean$p), angle=90, hjust = 0)
  p4 <- plot_amsd_histogram(amsd_output_cnv_sum) +
    ggtitle("CNV spectra, sums") +
    geom_label(x = amsd_output_cnv_sum$cosine, y = 1000, label = paste0("p=",amsd_output_cnv_sum$p), angle=90)
  p5 <- plot_amsd_histogram(amsd_output_cnv_sum2) +
    ggtitle("CNV spectra, sums,\nno outliers") +
    geom_label(x = amsd_output_cnv_sum2$cosine, y = 1000, label = paste0("p=",amsd_output_cnv_sum2$p), angle=90, hjust = 0)

  asbestos_fig1 <- ggarrange(p1,p3,
                             high_outlier_plot, 
                             p2,p4,p5,
                             nrow = 2,
                             ncol = 3,
                             labels = c("A", "C", "E", "B","D", "F")
  )
  asbestos_fig2 <- ggarrange(exp_cnv_plot,
                             non_exp_cnv_plot,
                             nrow = 2,
                             ncol = 1,
                             labels = c("G", "H"))
  asbestos_fig <- ggarrange(asbestos_fig1,
                            asbestos_fig2,
                            nrow = 2,
                            ncol = 1)
  #asbestos_fig
  ggsave("../outputs/asbestos_supp_fig.png",
         plot = asbestos_fig,
         width = 7, height = 10, units = "in")
  