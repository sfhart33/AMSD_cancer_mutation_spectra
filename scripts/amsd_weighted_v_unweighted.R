library(tidyverse)
library(ggrepel)
library(sigfit)

# load data
# setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github/AMSD_cancer_mutation_spectra/scripts")
mouse_amsd_output <- readRDS("../outputs/mouse_amsd_output.rds")
mouse_amsd_output2 <- readRDS("../outputs/mouse_amsd_output_unweighted.rds")
mouse_amsd_output$pvalues2 <- mouse_amsd_output2$pvalues
# perms <- readRDS("../outputs/mouse_amsd_perms.rds")
# mexposuresig <- readRDS("../inputs/mouse_exposuresig.rds")
# mexposure <- readRDS("../inputs/mouse_exposure.rds")
# mouse_carcinogen_counts <- readRDS("../inputs/mouse_carcinogen_spectra.rds") # counts
# mouse_carcinogen_spectra <- mouse_carcinogen_counts/rowSums(mouse_carcinogen_counts) # spectra


# load data
ancestry_amsd_output <- readRDS("../outputs/ancestry_amsd_output.rds")
ancestry_amsd_output2 <- readRDS("../outputs/ancestry_amsd_output_unweighted.rds")
ancestry_amsd_output$pvalues2 <- ancestry_amsd_output2$pvalues
anc_spectra <- readRDS("../outputs/ancestry_spectra.rds")
perms <- readRDS("../outputs/ancestry_amsd_perms.rds")
perms2 <- readRDS("../outputs/ancestry_amsd_perms_unweighted.rds")

anc_plot <- ancestry_amsd_output %>%
  ggplot(aes(-log10(pvalues),-log10(pvalues2), color = comparison, label = tumor_type)) +
  geom_point() +
  geom_smooth(method = "lm",
              inherit.aes = FALSE,
              aes(-log10(pvalues),
                  -log10(pvalues2)),
              color = "black")+
  geom_label_repel()+
  geom_vline(xintercept = -log10(0.05))+
  geom_vline(xintercept = -log10(0.05/67))+
  geom_hline(yintercept = -log10(0.05))+
  geom_hline(yintercept = -log10(0.05/67))+
  labs(x = "-log10(AMSD p-value): all samples weighted equally",
       y = "-log10(AMSD p-value):\nsamples weighted by mutation count")+
  theme_classic()

mouse_plot <- mouse_amsd_output %>%
  ggplot(aes(-log10(pvalues),-log10(pvalues2), color = tissue, label = exposure)) +
  geom_point() +
  geom_smooth(method = "lm",
              inherit.aes = FALSE,
              aes(-log10(pvalues),
                  -log10(pvalues2)),
              color = "black")+
  geom_label_repel()+
  geom_vline(xintercept = -log10(0.05))+
  geom_vline(xintercept = -log10(0.05/67))+
  geom_hline(yintercept = -log10(0.05))+
  geom_hline(yintercept = -log10(0.05/67))+
  labs(x = "-log10(AMSD p-value): all samples weighted equally",
       y = "-log10(AMSD p-value):\nsamples weighted by mutation count")+
  theme_classic()

bothplots <- ggarrange(mouse_plot, anc_plot, nrow=2, ncol=1,  labels = c("A","B"))
bothplots

ggsave("../outputs/amsd_weighted_v_unweighted.png",
       plot = bothplots,
       width = 7,
       height = 8,
       units = "in"
)


# anc_spectra %>%
#   filter(tumor_type == "KIRP", anc3 %in% c("eur", "eas"))  %>%
#   ggplot(aes(anc3, mut_counts))+
#     geom_boxplot()+
#     geom_point()
# 
# anc_spectra %>%
#   filter(tumor_type == "BLCA") %>%
#   ggplot(aes(anc3, mut_counts))+
#   geom_boxplot()+
#   geom_point()
# 
# anc_spectra %>%
#   filter(tumor_type == "UCEC", anc3 %in% c("afr", "eas", "eur")) %>%
#   ggplot(aes(anc3, mut_counts))+
#   geom_boxplot()+
#   geom_point()
# 
# anc_spectra %>%
#   filter(tumor_type == "LUAD", anc3 %in% c("afr", "eas", "eur")) %>%
#   ggplot(aes(anc3, mut_counts))+
#   geom_boxplot()+
#   geom_point()
# 
# anc_spectra %>%
#   filter(tumor_type == "KIRP", anc3 == "eas") %>%
#   select(-IID, -tumor_type, -consensus_ancestry, -anc2, -anc3, -mut_counts) %>% 
#   plot_spectrum()
# 
# perms %>%
#   ggplot(aes(KIRP.eas_eur))+
#   geom_histogram()+
#   geom_vline(xintercept = ancestry_amsd_output[26,"cosines"])
# perms2 %>%
#   ggplot(aes(KIRP.eas_eur))+
#   geom_histogram()+
#   geom_vline(xintercept = ancestry_amsd_output2[26,"cosines"])  
# 
# 
# perms %>%
#   ggplot(aes(UCEC.eas_eur))+
#   geom_histogram()+
#   geom_vline(xintercept = ancestry_amsd_output[66,"cosines"])
# perms2 %>%
#   ggplot(aes(UCEC.eas_eur))+
#   geom_histogram()+
#   geom_vline(xintercept = ancestry_amsd_output2[66,"cosines"])  
# 
# perms %>%
#   ggplot(aes(BLCA.eas_eur))+
#   geom_histogram()+
#   geom_vline(xintercept = ancestry_amsd_output[3,"cosines"])
# perms2 %>%
#   ggplot(aes(BLCA.eas_eur))+
#   geom_histogram()+
#   geom_vline(xintercept = ancestry_amsd_output2[3,"cosines"]) 
