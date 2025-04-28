library(tidyverse)
library(ggpubr)
library(svglite)
library(sigfit)
library(ggrepel)
#library(ggbreak)
data("cosmic_signatures_v3.2")
setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github/AMSD_cancer_mutation_spectra/scripts")
source("amsd_functions.R")

# load data
  ancestry_amsd_output <- readRDS("../outputs/ancestry_amsd_output.rds")
  anc_spectra <- readRDS("../outputs/ancestry_spectra.rds")
  perms <- readRDS("../outputs/ancestry_amsd_perms.rds")
  
# difference FDR types
  ancestry_amsd_output$padj_BH <- p.adjust(ancestry_amsd_output$pvalues, method="BH")
  # ancestry_amsd_output$padj_BY <- p.adjust(ancestry_amsd_output$pvalues, method="BY")
  ancestry_amsd_output$padj_Bonf <- p.adjust(ancestry_amsd_output$pvalues, method="bonferroni")
  ggplot(ancestry_amsd_output, aes(-log10(pvalues),-log10(padj_BH)))+
    geom_point()+
    geom_hline(yintercept = -log10(0.05))+
    stat_smooth(method="lm",se=F)+
    ggtitle("Benjamini-Hochberg vs unadjusted pvalues from ancestry comparison")
  x <- -log10(ancestry_amsd_output$padj_BH)
  y <- -log10(ancestry_amsd_output$pvalues)
  BH_regression <- lm(y ~ x)
  BH_threshold <- -log10(0.05) * coef(BH_regression)[2]
  
  
# volcano plot summary of everything together
  ancestry_volcano <- ancestry_amsd_output %>%
    mutate(log10pval = -log10(pvalues),
           comparison2 = case_when(comparison == "afr_eas" ~ "AFR vs EAS",
                                   comparison == "afr_eur" ~ "AFR vs EUR",
                                   comparison == "eas_eur" ~ "EAS vs EUR")) %>%
    ggplot()+
    geom_point(aes(x=cosines,
                   y = log10pval,
                   # shape = comparison,
                   # color = tumor_type,
                   color = comparison2,
                   size = min_anc_n))+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
    geom_hline(yintercept = -log10(0.05/nrow(ancestry_amsd_output)), linetype = "dashed")+
    geom_hline(yintercept = BH_threshold, linetype = "dashed")+
    geom_text(aes(x=0.225, y = (-log10(0.05)+0.1)), label = "P = 0.05")+
    geom_text(aes(x=0.225, y = (-log10(0.05/nrow(ancestry_amsd_output))+0.1)), label = "Bonferroni FWER = 0.05")+
    geom_text(aes(x=0.225, y = BH_threshold + 0.1), label = "Benjamini-Hochberg FDR = 0.05")+
    xlim(0,0.25)+
    scale_size_continuous(
      range = c(1, 4),
      breaks = c(5,10,20,40,80,160)
    )+
    labs(x="Cosine distance",
         y="-log10(p-value)", 
         # color="Tumor type",
         # shape="Ancestry comparison",
         color="Ancestry\ncomparison",
         size="Tumor count\n(smaller group)")+
    theme_classic()+ 
    theme(legend.title.align = 0.5,
          legend.direction = "vertical",
          legend.box.just = "center")
  ancestry_volcano
  ggsave("../outputs/ancestry_amsd_output.png",
         plot = ancestry_volcano,
         width = 7, height = 7.5, units = "in")

# violin plots of AMSD null distributions
  perms2 <- perms %>%
    column_to_rownames(var = "rep") %>%
    pivot_longer(cols = everything()) %>%
    separate(name, into = c("tissue", "comparison"), sep = "\\.", remove = FALSE)
  quantiles <- group_by(perms2, tissue, comparison) %>%
    summarise(p0.05 = quantile(value, probs = 0.95),
              FDR = quantile(value, probs = 1-(0.05/nrow(ancestry_amsd_output))))
  ancestry_violin <- perms2 %>%
    ggplot(aes(x = tissue, y = value))+
    geom_violin(adjust =0.5, scale = "width")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab("")+
    ylab("Cosine distance")+
    geom_point(data = ancestry_amsd_output,
               aes(x = tumor_type, y = cosines))+
    geom_point(data = quantiles,
               aes(x = tissue, y = p0.05), shape = 95, size = 5)+
    geom_point(data = quantiles,
               aes(x = tissue, y = FDR), shape = 95, size = 5)+
    facet_grid(rows = vars(comparison))
  ancestry_violin
  ggsave("../outputs/ancestry_amsd_output_violin.png",
         plot = ancestry_violin)
  

# run for each significant ancestry
    sig_outputs <- filter(ancestry_amsd_output, pvalues < 0.05)
    signature_results <- data.frame(row.names = rownames(cosmic_signatures_v3.2))
    
  pdf("../outputs/ancestry_signature_comparisons.pdf")
  for(count in 1:nrow(sig_outputs)){
    t <- sig_outputs[count, "tumor_type"]
    a1 <- sig_outputs[count, "ancestry1"]
    a2 <- sig_outputs[count, "ancestry2"]
    set1 <- filter(anc_spectra, tumor_type == t, anc3 == a1)
    fit_result1 <- fit_signatures(counts =  t(as.matrix(colMeans(set1[,6:101])))*sum(set1$mut_counts), 
                                  signatures = cosmic_signatures_v3.2,
                                  iter = 2000, 
                                  warmup = 1000, 
                                  chains = 1, 
                                  opportunities = "human-exome",
                                  seed = 1756) %>%
      retrieve_pars(par = "exposures")
    set2 <- filter(anc_spectra, tumor_type == t, anc3 == a2)
    fit_result2 <- fit_signatures(counts =  t(as.matrix(colMeans(set2[,6:101])))*sum(set2$mut_counts), 
                                  signatures = cosmic_signatures_v3.2,
                                  iter = 2000, 
                                  warmup = 1000, 
                                  chains = 1, 
                                  opportunities = "human-exome",
                                  seed = 1756) %>%
      retrieve_pars(par = "exposures")
    sig_comp <- rbind(fit_result1$mean, 
                      fit_result1$lower_95,
                      fit_result1$upper_95,
                      fit_result2$mean,
                      fit_result2$lower_95,
                      fit_result2$upper_95) %>%
      t() %>%
      as.data.frame()
    colnames(sig_comp) <- c("a1_mean", "a1_lower95", "a1_upper95", "a2_mean", "a2_lower95", "a2_upper95")
    sig_comp2 <- sig_comp  %>%
      mutate(dif = a1_mean - a2_mean) %>%
      select(a1_mean, a2_mean, dif)
    colnames(sig_comp2) <- c(paste0(a1,"__",t, ".exp"), paste0(a2,"__",t, ".exp"), paste0(a1,"_",a2,"__",t, ".dif"))
    signature_results <- cbind(signature_results, sig_comp2)
    
    threshold = 0.01
    plot1 <- sig_comp %>%
      ggplot(aes(x = a1_mean, y = a2_mean))+
      geom_point(size = 0.25)+
      geom_pointrange(aes(xmin = a1_lower95, xmax = a1_upper95), size = 0.25)+
      geom_pointrange(aes(ymin = a2_lower95, ymax = a2_upper95), size = 0.25)+
      geom_abline(intercept = 0, slope = 1)+
      # geom_abline(intercept = threshold, slope = 1, linetype = "dashed")+
      # geom_abline(intercept = -threshold, slope = 1, linetype = "dashed")+
      # labs(title = paste(t, a1, a2,"signature exposure comparison"),
      #      x=paste(a1, "ancestry \n(signature exposure)"),
      #      y=paste(a2, "ancestry \n(signature exposure)"))+
      labs(title = paste(t, a1, a2,"signature exposure comparison"),
           x=a1,
           y=a2)+
      # lims(x=c(0,0.25),
      #      y=c(0,0.25))+
      theme_classic()
    plot2 <- plot1 +
      geom_text_repel(data = filter(sig_comp, abs(a1_mean - a2_mean) > threshold),
                      aes(label = rownames(filter(sig_comp, abs(a1_mean - a2_mean) > threshold))))
    print(plot2)
    assign(paste0(t, ".", a1, "_v_", a2, "_sigcomp"), plot1)
    # ggsave(paste0("../outputs/",t, ".", a1, "_v_", a2, "_sigcomp",".png"),
    #        plot = plot1)
    print(paste(t, a1, "signature exposure comparison done"))
  }
  dev.off()

  
# Signature comparison plots
  liver_plot <- LIHC.eas_v_eur_sigcomp +
    labs(title = "Liver",
         x = "EAS",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  esophageal_plot <- ESCA.eas_v_eur_sigcomp +
    labs(title = "Esophageal",
         x = "EAS",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  bladder_plot <- BLCA.eas_v_eur_sigcomp +
    labs(title = "Bladder",
         x = "EAS",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  uterine_plot1 <- UCEC.eas_v_eur_sigcomp +
    labs(title = "Uterine",
         x = "EAS",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  lung_plot1 <- LUAD.afr_v_eas_sigcomp +
    labs(title = "Lung Adeno.",
         x = "AFR",
         y = "EAS") +
    theme(plot.title = element_text(hjust = 0.5))
  lung_plot2 <- LUAD.afr_v_eur_sigcomp + 
    labs(title = "Lung Adeno.",
         x = "AFR",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  
  colorectal_plot <- COAD.eas_v_eur_sigcomp + 
    labs(title = "Colorectal",
         x = "EAS",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  uterine_plot2 <- UCEC.afr_v_eas_sigcomp + 
    labs(title = "Uterine",
         x = "AFR",
         y = "EAS") +
    theme(plot.title = element_text(hjust = 0.5))
  
  skin_plot <- SKCM.eas_v_eur_sigcomp +
    labs(title = "Melanoma",
         x = "EAS",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  ovarian_plot1 <- OV.afr_v_eas_sigcomp + 
    labs(title = "Ovarian",
         x = "AFR",
         y = "EAS") +
    theme(plot.title = element_text(hjust = 0.5))
  ovarian_plot2 <- OV.eas_v_eur_sigcomp + 
    labs(title = "Ovarian",
         x = "EAS",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_flip()
  head_plot1 <- HNSC.afr_v_eas_sigcomp + 
    labs(title = "Head/neck",
         x = "AFR",
         y = "EAS") +
    theme(plot.title = element_text(hjust = 0.5))
  head_plot2 <- HNSC.afr_v_eur_sigcomp + 
    labs(title = "Head/neck",
         x = "AFR",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  head_plot3 <- HNSC.eas_v_eur_sigcomp + 
    labs(title = "Head/neck",
         x = "EAS",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  uterine_plot3 <- UCEC.afr_v_eur_sigcomp + 
    labs(title = "Uterine",
         x = "AFR",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  bladder_plot2 <- BLCA.afr_v_eur_sigcomp + 
    labs(title = "Bladder",
         x = "AFR",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  colorectal_plot2 <- COAD.afr_v_eas_sigcomp + 
    labs(title = "Colorectal",
         x = "AFR",
         y = "EAS") +
    theme(plot.title = element_text(hjust = 0.5))
  breast_plot1 <- BRCA.afr_v_eas_sigcomp + 
    labs(title = "Breast",
         x = "AFR",
         y = "EAS") +
    theme(plot.title = element_text(hjust = 0.5))
  breast_plot2 <- BRCA.eas_v_eur_sigcomp + 
    labs(title = "Breast",
         x = "EAS",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  lung_plot3 <- LUSC.eas_v_eur_sigcomp + 
    labs(title = "Lung sq. cell",
         x = "EAS",
         y = "EUR") +
    theme(plot.title = element_text(hjust = 0.5))
  
  fig3 <- ggarrange(ancestry_volcano,
                    ggarrange(NULL,
                              ggarrange(esophageal_plot,
                                        liver_plot,
                                        lung_plot1,
                                        uterine_plot1,
                                        bladder_plot,
                                        lung_plot2,
                                        nrow = 2,
                                        ncol = 3),
                              NULL,
                              ncol = 3,
                              widths = c(0.1,0.8,0.1)),
                    nrow = 2,
                    labels = c("A", "B")) 
  fig3_new <- ggarrange(ancestry_volcano,
                        ggarrange(esophageal_plot + coord_flip(),
                                  liver_plot + coord_flip(),
                                  bladder_plot + coord_flip(),
                                  lung_plot1 + coord_flip(),
                                  lung_plot2 + coord_flip(),
                                  uterine_plot1 + coord_flip(),
                                  uterine_plot2,
                                  colorectal_plot + coord_flip(),
                                  colorectal_plot2,
                                  skin_plot + coord_flip(),
                                  breast_plot1,
                                  breast_plot2 + coord_flip(),
                                  ovarian_plot2,
                                  head_plot1,
                                  bladder_plot2,
                                  nrow = 3,
                                  ncol = 5),
                    nrow = 2,
                    labels = c("A", "B")) 
  fig3_new
  ggsave("../outputs/Figure3.png",
         plot = fig3_new,
         width = 10,
         height = 10,
         units = "in")
  ggsave("../outputs/Figure3.svg",
         plot = fig3_new,
         width = 10,
         height = 10,
         units = "in")

  fig3_v3 <- ggarrange(ancestry_volcano,
                        ggarrange(esophageal_plot + coord_flip() + xlim(0,0.15) +ylim(0,0.15),
                                  liver_plot + coord_flip(),
                                  bladder_plot + coord_flip(),
                                  lung_plot1 + coord_flip(),
                                  lung_plot2 + coord_flip(),
                                  uterine_plot1 + coord_flip() + xlim(0,0.15) +ylim(0,0.15),
                                  uterine_plot2 + xlim(0,0.15) +ylim(0,0.15),
                                  colorectal_plot + coord_flip() + xlim(0,0.15) +ylim(0,0.15),
                                  colorectal_plot2 + xlim(0,0.15) +ylim(0,0.15),
                                  skin_plot + coord_flip(),
                                  breast_plot1 + xlim(0,0.15) +ylim(0,0.15),
                                  breast_plot2 + coord_flip() + xlim(0,0.15) +ylim(0,0.15),
                                  ovarian_plot2,
                                  head_plot1,
                                  bladder_plot2,
                                  nrow = 3,
                                  ncol = 5),
                        nrow = 2,
                        labels = c("A", "B")) 
  #fig3_new
  fig3_v3
  ggsave("../outputs/Figure3_v3.png",
         plot = fig3_v3,
         width = 10,
         height = 10,
         units = "in")
  ggsave("../outputs/Figure3_v3.svg",
         plot = fig3_v3,
         width = 10,
         height = 10,
         units = "in")
  supp_fig <- ggarrange(esophageal_plot,
                        liver_plot,
                        bladder_plot,
                        lung_plot1,
                        uterine_plot1,
                        lung_plot2,
                        
                        skin_plot,
                        ovarian_plot2,
                        uterine_plot2,
                        breast_plot1,
                        head_plot1,
                        bladder_plot2,
                        colorectal_plot,
                        colorectal_plot2,
                        uterine_plot3,
                        breast_plot2,
                        lung_plot3,
                        head_plot2,
                        head_plot3,
                        ovarian_plot2,
                        nrow = 5,
                        ncol = 4)
  supp_fig
  ggsave("../outputs/Figure3_supp.png",
         plot = supp_fig,
         width = 7,
         height = 8,
         units = "in")
  ggsave("../outputs/Figure3_supp.svg",
         plot = supp_fig,
         width = 7,
         height = 8,
         units = "in")
  
# Individual spectra plots
  plot_tcga_spectra <- function(tissue, ancestry){
    samples <- anc_spectra %>%
      filter(anc3 == ancestry, tumor_type == tissue) %>%
      select(-IID, -tumor_type, -consensus_ancestry, -anc2, -anc3, -mut_counts)
    sample_n <- nrow(samples)
    spectra1 <- colMeans(samples, na.rm=TRUE)
    spectra1sd <- apply(samples, 2, sd, na.rm = TRUE)
    
    default_df <- readRDS("../inputs/default_spectrum_df.rds")
    
    default_df$spectra1 <- spectra1
    default_df$spectra1sd <- spectra1sd
    default_df$spectra1sem <-2 * spectra1sd / sqrt(sample_n)
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    
    default_df %>%
      ggplot(aes(x=trinuc, y = spectra1, fill = mut))+
      geom_col()+
      geom_errorbar(aes(ymin=spectra1-spectra1sem, ymax=spectra1+spectra1sem),
                    width=.2,
                    position=position_dodge(.9))+
      scale_fill_manual(values = COLORS)+
      scale_y_continuous(expand = c(0,0)) +
      #ylim(0,0.05)+
      facet_grid(cols = vars(mut), scales = 'free')+
      theme_bw()+
      theme(
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.spacing = unit(0,'lines'),
        #aspect.ratio = 1.5,
        legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(face = "bold", vjust = 0.5))+
      labs(title = paste0(tissue,", ", ancestry," (N = ", sample_n, ")"), 
           x = "Trinucleotide context",
           y = "Mutation fraction") %>%
      return()
  }
  P_lung <- ggarrange(plot_tcga_spectra("LUAD","afr"),
            plot_tcga_spectra("LUAD","eur")+ theme(strip.text = element_blank()),
            plot_tcga_spectra("LUAD","eas")+ theme(strip.text = element_blank()),
            plot_tcga_spectra("LIHC","eas")+ theme(strip.text = element_blank()),
            plot_tcga_spectra("LIHC","eur")+ theme(strip.text = element_blank()),
            plot_tcga_spectra("ESCA","eas")+ theme(strip.text = element_blank()),
            plot_tcga_spectra("ESCA","eur")+theme(strip.text = element_blank(),
                                                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)),
            nrow = 7,
            ncol = 1)
  P_lung          
            
  
  P_ut <- ggarrange(plot_tcga_spectra("UCEC","afr"),
            plot_tcga_spectra("UCEC","eas")+ theme(strip.text = element_blank()),
            plot_tcga_spectra("UCEC","eur")+ theme(strip.text = element_blank()),
            plot_tcga_spectra("COAD","afr")+ theme(strip.text = element_blank()),
            plot_tcga_spectra("COAD","eas")+ theme(strip.text = element_blank()),
            plot_tcga_spectra("COAD","eur")+theme(strip.text = element_blank(),
                                                  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6)),
            nrow = 6,
            ncol = 1)

  
  ggsave("../outputs/Figure3_lung_spectra.png",
         P_lung,
         width = 7,
         height = 11,
         units = "in")
  ggsave("../outputs/Figure3_uterine_spectra.png",
         P_ut,
         width = 7,
         height = 11,
         units = "in")
  