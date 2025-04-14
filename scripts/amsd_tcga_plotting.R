library(tidyverse)
library(ggpubr)
library(svglite)
library(sigfit)
library(ggrepel)
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
  ancestry_volcano <- mutate(ancestry_amsd_output, log10pval = -log10(pvalues)) %>%
    ggplot()+
    geom_point(aes(x=cosines,
                   y = log10pval,
                   # shape = comparison,
                   # color = tumor_type,
                   color = comparison,
                   size = min_anc_n))+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
    geom_hline(yintercept = -log10(0.05/nrow(ancestry_amsd_output)), linetype = "dashed")+
    geom_hline(yintercept = BH_threshold, linetype = "dashed")+
    geom_text(aes(x=0.225, y = (-log10(0.05)+0.1)), label = "p=0.05")+
    geom_text(aes(x=0.225, y = (-log10(0.05/nrow(ancestry_amsd_output))+0.1)), label = "Bonferroni FDR=0.05")+
    geom_text(aes(x=0.225, y = BH_threshold + 0.1), label = "Benjamini-Hochberg FDR=0.05")+
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
         size="Tumor count\n(lower count)")+
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
  
# mutational signature analysis
  # t <- "LUAD"
  # a1 <- "afr"
  # a2 <- "eas"
  # set1 <- filter(anc_spectra, tumor_type == t, anc3 == a1)
  # fit_result1 <- fit_signatures(counts =  t(as.matrix(colMeans(set1[,6:101])))*sum(set1$mut_counts),
  #                               signatures = cosmic_signatures_v3.2,
  #                               iter = 2000,
  #                               warmup = 1000,
  #                               chains = 1,
  #                               seed = 1756) %>%
  #   retrieve_pars(par = "exposures")
  # set2 <- filter(anc_spectra, tumor_type == t, anc3 == a2)
  # fit_result2 <- fit_signatures(counts =  t(as.matrix(colMeans(set2[,6:101])))*sum(set2$mut_counts),
  #                               signatures = cosmic_signatures_v3.2,
  #                               iter = 2000,
  #                               warmup = 1000,
  #                               chains = 1,
  #                               seed = 1756) %>%
  #   retrieve_pars(par = "exposures")
  # sig_comp <- rbind(fit_result1$mean,
  #                   fit_result1$lower_95,
  #                   fit_result1$upper_95,
  #                   fit_result2$mean,
  #                   fit_result2$lower_95,
  #                   fit_result2$upper_95) %>%
  #   t() %>%
  #   as.data.frame()
  # colnames(sig_comp) <- c("a1_mean", "a1_lower95", "a1_upper95", "a2_mean", "a2_lower95", "a2_upper95")
  # # sig_comp %>%
  # #   mutate(dif = a1_mean - a2_mean) %>%
  # #   ggplot(aes(x = rownames(sig_comp), y = dif))+
  # #   geom_col()+
  # #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # sig_comp %>%
  #   ggplot(aes(x = a1_mean, y = a2_mean))+
  #     geom_point(size = 4)+
  #   geom_pointrange(aes(xmin = a1_lower95, xmax = a1_upper95))+
  #   geom_pointrange(aes(ymin = a2_lower95, ymax = a2_upper95))+
  #   geom_abline(intercept = 0, slope = 1)+
  #   # geom_abline(intercept = 0.05, slope = 1, linetype = "dashed")+
  #   # geom_abline(intercept = -0.05, slope = 1, linetype = "dashed")+
  #   labs(x=paste(t, a1, "(signature exposure)"),
  #        y=paste(t, a2, "(signature exposure)"))+
  #   theme_classic()+
  #   geom_text_repel(data = filter(sig_comp, abs(a1_mean - a2_mean) > 0.05),
  #                   aes(label = rownames(filter(sig_comp, abs(a1_mean - a2_mean) > 0.05))))
  # 


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
    ggsave(paste0("../outputs/",t, ".", a1, "_v_", a2, "_sigcomp",".png"),
           plot = plot1)
    print(paste(t, a1, "signature exposure comparison done"))
  }
  dev.off()

  signature_output <- t(signature_results) %>%
    as.data.frame() %>% 
    rownames_to_column(var = "description") %>%
    separate(description, into = c("description","type"), sep = "\\.") %>%
    separate(description, into = c("anc","tumor"), sep = "__") %>%
    pivot_longer(cols = rownames(signature_results))
  filter(signature_output, anc %in% c("afr", "eas", "eur")) %>%
    separate(name, into = c(NA, "sig"), sep = "BS", remove = FALSE) %>%
    ggplot(aes(x = name, y = value, color = anc))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  filter(signature_output, anc %in% c("afr", "eas", "eur"), name %in% c("SBS1","SBS5")) %>%
    separate(name, into = c(NA, "sig"), sep = "BS", remove = FALSE) %>%
    ggplot(aes(x = name, y = value, color = anc))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(title = "Signature fitting, correcting for exons",
         x = "signature",
         y = "Fraction of mutations fit to signature")
  
  filter(signature_output, !(anc %in% c("afr", "eas", "eur"))) %>%
    ggplot(aes(x = name, y=value, color = anc)) +
      geom_boxplot(outliers = FALSE)+
      geom_point(position=position_jitterdodge(jitter.width = 0.1))+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  filter(signature_output, !(anc %in% c("afr", "eas", "eur"))) %>%
    group_by(name, anc) %>%
    summarise(mean = mean(value), sd = sd(value))#, ttest = t.test(value), wilcoxtest = wilcox.test(value))
  signature_tests <- filter(signature_output, !(anc %in% c("afr", "eas", "eur"))) %>%
    group_by(name, anc) %>%
    summarise(mean = mean(value), sd = sd(value), n = n(), ttest = t.test(value)$p.value, wilcoxtest = wilcox.test(value)$p.value)
  signature_tests
  signature_tests %>%
    arrange(wilcoxtest)
  signature_tests %>%
    arrange(ttest)
  signature_tests %>%
    ggplot(aes(x=mean, y = -log10(ttest), color = anc, size = n, label = name))+
    geom_point()+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
    scale_size_continuous(range = c(3,4), breaks = c(4,6,10))+
    geom_text(data = filter(signature_tests, ttest < 0.05),
                    aes(x=mean, y = -log10(ttest), size = 8, label = name))+
    theme_classic()+
    labs(title = "Are any signatures dispropotiantly high/low across significant comparisons in AMSD?",
         x = "Average difference in exposure for each comparison between the two ancestries",
         y = "-log10(pvalue from t-test)")
  
  
  ancestry_volcano
  
  liver_plot <- LIHC.eas_v_eur_sigcomp + ggtitle("Liver") + theme(plot.title = element_text(hjust = 0.5))
  esophageal_plot <- ESCA.eas_v_eur_sigcomp + ggtitle("Esophageal") + theme(plot.title = element_text(hjust = 0.5))
  bladder_plot <- BLCA.eas_v_eur_sigcomp + ggtitle("Bladder") + theme(plot.title = element_text(hjust = 0.5))
  uterine_plot1 <- UCEC.eas_v_eur_sigcomp + ggtitle("Uterine") + theme(plot.title = element_text(hjust = 0.5))
  lung_plot1 <- LUAD.afr_v_eas_sigcomp +
    ggtitle("Lung Adeno.") +
    theme(plot.title = element_text(hjust = 0.5))
  lung_plot2 <- LUAD.afr_v_eur_sigcomp + ggtitle("Lung Adeno.") + theme(plot.title = element_text(hjust = 0.5))
  

  colorectal_plot <- COAD.eas_v_eur_sigcomp + ggtitle("Colorectal") + theme(plot.title = element_text(hjust = 0.5))
  uterine_plot2 <- UCEC.afr_v_eas_sigcomp + ggtitle("Uterine") + theme(plot.title = element_text(hjust = 0.5))
  
  skin_plot <- SKCM.eas_v_eur_sigcomp + ggtitle("Melanoma") + theme(plot.title = element_text(hjust = 0.5))
  ovarian_plot1 <- OV.afr_v_eas_sigcomp + ggtitle("Ovarian") + theme(plot.title = element_text(hjust = 0.5))
  ovarian_plot2 <- OV.eas_v_eur_sigcomp + ggtitle("Ovarian") + theme(plot.title = element_text(hjust = 0.5)) +
    coord_flip()
  head_plot1 <- HNSC.afr_v_eas_sigcomp + ggtitle("Head/neck") + theme(plot.title = element_text(hjust = 0.5))
  head_plot2 <- HNSC.afr_v_eur_sigcomp + ggtitle("Head/neck") + theme(plot.title = element_text(hjust = 0.5))
  head_plot3 <- HNSC.eas_v_eur_sigcomp + ggtitle("Head/neck") + theme(plot.title = element_text(hjust = 0.5))
  uterine_plot3 <- UCEC.afr_v_eur_sigcomp + ggtitle("Uterine") + theme(plot.title = element_text(hjust = 0.5))
  bladder_plot2 <- BLCA.afr_v_eur_sigcomp + ggtitle("Bladder") + theme(plot.title = element_text(hjust = 0.5))
  colorectal_plot2 <- COAD.afr_v_eas_sigcomp + ggtitle("Colorectal") + theme(plot.title = element_text(hjust = 0.5))
  breast_plot1 <- BRCA.afr_v_eas_sigcomp + ggtitle("Breast") + theme(plot.title = element_text(hjust = 0.5))
  breast_plot2 <- BRCA.eas_v_eur_sigcomp + ggtitle("Breast") + theme(plot.title = element_text(hjust = 0.5))
  lung_plot3 <- LUSC.eas_v_eur_sigcomp + ggtitle("Lung sq. cell") + theme(plot.title = element_text(hjust = 0.5))
  
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
         height = 12,
         units = "in")
  ggsave("../outputs/Figure3.svg",
         plot = fig3_new,
         width = 10,
         height = 12,
         units = "in")
  # supp_fig <- ggarrange(skin_plot,
  #                       ovarian_plot1,
  #                       ovarian_plot2,
  #                       uterine_plot2,
  #                       uterine_plot3,
  #                       bladder_plot2,
  #                       head_plot1,
  #                       head_plot2,
  #                       head_plot3,
  #                       breast_plot1,
  #                       breast_plot2,
  #                       lung_plot3,
  #                       colorectal_plot,
  #                       colorectal_plot2,
  #                       nrow = 5,
  #                       ncol = 3)
  # supp_fig
  # ggsave("../outputs/Figure3_supp.png",
  #        plot = supp_fig,
  #        width = 7,
  #        height = 11,
  #        units = "in")
  # ggsave("../outputs/Figure3_supp.svg",
  #        plot = supp_fig,
  #        width = 7,
  #        height = 11,
  #        units = "in")
  