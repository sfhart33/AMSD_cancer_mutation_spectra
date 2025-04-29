library(tidyverse)
library(sigfit)
library(ggpubr)
library(svglite)
data("cosmic_signatures_v3.2")
setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github/AMSD_cancer_mutation_spectra/scripts")
source("amsd_functions.R")

group1 <- simulate_spectra(n_samples = 10,
                           n_mutations = 200,
                           sig_probs = c(SBS1 = 0.2, SBS5 = 0.8),
                           signatures = cosmic_signatures_v3.2)
group2 <- simulate_spectra(n_samples = 10,
                           n_mutations = 200,
                           sig_probs = c(SBS1 = 0.1,
                                         SBS5 = 0.6,
                                         SBS14 = 0.2,
                                         SBS12 = 0.2),
                           signatures = cosmic_signatures_v3.2)
group3 <- simulate_spectra(n_samples = 10,
                           n_mutations = 200,
                           sig_probs = c(SBS1 = 0.1,
                                         SBS5 = 0.6,
                                         SBS18 = 0.2,
                                         SBS13 = 0.2),
                           signatures = cosmic_signatures_v3.2)
group4 <- simulate_spectra(n_samples = 10,
                           n_mutations = 200,
                           sig_probs = c(SBS1 = 0.05,
                                         SBS5 = 0.6,
                                         SBS3 = 0.25,
                                         SBS13 = 0.1),
                           signatures = cosmic_signatures_v3.2)

simple_spectra <- function(input){
  default_df <- readRDS("../inputs/default_spectrum_df.rds")
  default_df$spectra1 <- input
  
  COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
  default_df %>%
    ggplot(aes(x=trinuc, y = spectra1, fill = mut))+
    geom_col()+
    # geom_errorbar(aes(ymin=spectra1-spectra1sd, ymax=spectra1+spectra1sd),
    #               width=.2,
    #               position=position_dodge(.9))+
    scale_fill_manual(values = COLORS)+
    scale_y_continuous(expand = c(0,0)) +
    #ylim(0,0.05)+
    facet_grid(cols = vars(mut), scales = 'free')+
    theme_classic()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6),
      panel.spacing = unit(0,'lines'),
      #aspect.ratio = 1.5,
      legend.position="none",
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(face = "bold", vjust = 0.5))+
    theme(strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank()) %>%
    return()
}
simple_spectra_sd <- function(input1, input2, input3){
  default_df <- readRDS("../inputs/default_spectrum_df.rds")
  group1_summary <- cbind(default_df,
                          data.frame(a = input1,
                                     b = input2,
                                     c = input3)
  )  %>%
    rowwise() %>%
    mutate(sum = a+b+c,
           sd = sd(c(a,b,c))) %>%
    select(-spectra1sd, -a, -b, -c)
  group1_summary$spectra1 <- group1_summary$sum/sum(group1_summary$sum)
  group1_summary$spectra1sd <- group1_summary$sd/sum(group1_summary$sum)
  
  COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
  group1_summary %>%
    ggplot(aes(x=trinuc, y = spectra1, fill = mut))+
    geom_col()+
    geom_errorbar(aes(ymin=spectra1-spectra1sd, ymax=spectra1+spectra1sd),
                  width=.2,
                  position=position_dodge(.9))+
    scale_fill_manual(values = COLORS)+
    scale_y_continuous(expand = c(0,0)) +
    #ylim(0,0.05)+
    facet_grid(cols = vars(mut), scales = 'free')+
    theme_classic()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6),
      panel.spacing = unit(0,'lines'),
      #aspect.ratio = 1.5,
      legend.position="none",
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(face = "bold", vjust = 0.5))+
    theme(strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank()) %>%
    return()
}

U1 <- simple_spectra(group1[[1]])+ ggtitle("Unexposed tumor #1")
U2 <- simple_spectra(group3[[2]])+ ggtitle("Unexposed tumor #2")
U3 <- simple_spectra(group1[[3]])+ ggtitle("Unexposed tumor #3")
E1 <- simple_spectra(group2[[1]])+ ggtitle("Exposed tumor #1")
E2 <- simple_spectra(group2[[2]])+ ggtitle("Exposed tumor #2")
E3 <- simple_spectra(group4[[3]])+ ggtitle("Exposed tumor #3")
Uavg <- simple_spectra_sd(group1[[1]],group3[[2]],group1[[3]]) + ggtitle("Unexposed tumors")
Eavg <- simple_spectra_sd(group2[[1]],group2[[2]],group4[[3]]) + ggtitle("Exposed tumors")
R1 <- simple_spectra_sd(group1[[1]],group2[[2]],group3[[3]]) + ggtitle("Randomly sampled tumors")
R2 <- simple_spectra_sd(group4[[1]],group1[[2]],group1[[3]]) + ggtitle("Randomly sampled tumors")
amsd(as.data.frame(do.call(rbind, append(append(group1,group3),group4))),
     as.data.frame(do.call(rbind, append(append(append(group2,group3),group4),group2))),
     n_sim = 1000,
     seed=1234) %>%
  plot_amsd_histogram()
amsd_hist <- amsd(as.data.frame(do.call(rbind, append(append(group1,group3),group4))),
                  as.data.frame(do.call(rbind, append(append(append(group2,group3),group4),group2))),
                  n_sim = 1000,
                  seed=1234) %>%
  plot_amsd_histogram()

col1 <- ggarrange(U1,U2,U3,NULL,E1,E2,E3,
          nrow=7, ncol = 1)
col2 <- ggarrange(Uavg,NULL,Eavg,NULL, R1, NULL, R2,
          nrow=7, ncol = 1)
col3 <- ggarrange(NULL,amsd_hist, NULL,
                 nrow=3, ncol = 1)
fig1 <- ggarrange(col1,col2,col3,
          nrow=1, ncol = 3)
fig1

ggsave("../outputs/Figure1.png",
       plot = fig1,
       width = 14,
       height = 8,
       units = "in")
ggsave("../outputs/Figure1.svg",
       plot = fig1,
       width = 14,
       height = 8,
       units = "in")
