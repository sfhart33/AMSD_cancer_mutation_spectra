library(tidyverse)
library(sigfit)
library(ggrepel)
data("cosmic_signatures_v3.2")
source("amsd_functions.R")

# load data
  ancestry_amsd_output <- readRDS("../outputs/ancestry_amsd_output.rds")
  anc_spectra <- readRDS("../outputs/ancestry_spectra.rds")
  perms <- readRDS("../outputs/ancestry_amsd_perms.rds")

  
# volcano plot summary of everything together
  ancestry_volcano <- mutate(ancestry_amsd_output, log10pval = -log10(pvalues)) %>%
    ggplot()+
    geom_point(aes(x=cosines,
                   y = log10pval,
                   shape = comparison,
                   color = tumor_type,
                   size = min_anc_n))+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
    geom_hline(yintercept = -log10(0.05/nrow(ancestry_amsd_output)), linetype = "dashed")+
    #geom_hline(yintercept = -log10(1/reps))+
    geom_text(aes(x=0.225, y = (-log10(0.05/nrow(ancestry_amsd_output))+0.1)), label = "FDR=0.05")+
    geom_text(aes(x=0.225, y = (-log10(0.05)+0.1)), label = "p=0.05")+
    #geom_text(aes(x=0.225, y = (-log10(1/reps)+0.1)), label = "theoretical max")+
    theme_classic()+
    xlim(0,0.25)+
    scale_size_continuous(
      range = c(1, 6),
      breaks = c(5,10,20,40,80,160)
    )+
    labs(x="Cosine distance",
         y="-log10(p-value)", 
         color="Tumor type",
         shape="Ancestry comparison",
         size="Tumor count\n(lower anc count)")
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
  # t <- "UCEC"
  # a1 <- "eas"
  # a2 <- "eur"
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
  #   geom_abline(intercept = 0.05, slope = 1, linetype = "dashed")+
  #   geom_abline(intercept = -0.05, slope = 1, linetype = "dashed")+
  #   labs(x=paste(t, a1, "(signature exposure)"),
  #        y=paste(t, a2, "(signature exposure)"))+
  #   lims(x=c(0,0.25), y=c(0,0.25))+
  #   theme_classic()+
  #   geom_text_repel(data = filter(sig_comp, abs(a1_mean - a2_mean) > 0.05),
  #                   aes(label = rownames(filter(sig_comp, abs(a1_mean - a2_mean) > 0.05))))
  # 
  # 
  # 
  # run for each significant ancestry
  sig_outputs <- filter(ancestry_amsd_output, pvalues < 0.05)
  
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
                                  seed = 1756) %>%
      retrieve_pars(par = "exposures")
    set2 <- filter(anc_spectra, tumor_type == t, anc3 == a2)
    fit_result2 <- fit_signatures(counts =  t(as.matrix(colMeans(set2[,6:101])))*sum(set2$mut_counts), 
                                  signatures = cosmic_signatures_v3.2,
                                  iter = 2000, 
                                  warmup = 1000, 
                                  chains = 1, 
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
    # sig_comp %>%
    #   mutate(dif = a1_mean - a2_mean) %>%
    #   ggplot(aes(x = rownames(sig_comp), y = dif))+
    #   geom_col()+
    #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    threshold = 0.025
    plot1 <- sig_comp %>%
      ggplot(aes(x = a1_mean, y = a2_mean))+
      geom_point(size = 4)+
      geom_pointrange(aes(xmin = a1_lower95, xmax = a1_upper95))+
      geom_pointrange(aes(ymin = a2_lower95, ymax = a2_upper95))+
      geom_abline(intercept = 0, slope = 1)+
      # geom_abline(intercept = threshold, slope = 1, linetype = "dashed")+
      # geom_abline(intercept = -threshold, slope = 1, linetype = "dashed")+
      labs(title = paste(t, a1, a2,"signature exposure comparison"),
           x=paste(t, a1, "(signature exposure)"),
           y=paste(t, a2, "(signature exposure)"))+
      # lims(x=c(0,0.25),
      #      y=c(0,0.25))+
      theme_classic()+
      geom_text_repel(data = filter(sig_comp, abs(a1_mean - a2_mean) > threshold),
                      aes(label = rownames(filter(sig_comp, abs(a1_mean - a2_mean) > threshold))))
      
    print(plot1)
    assign(paste0(t, ".", a1, "_v_", a2, "_sigcomp"), plot1)
    ggsave(paste0("../outputs/",t, ".", a1, "_v_", a2, "_sigcomp",".png"),
           plot = plot1)
    print(paste(t, a1, "signature exposure comparison done"))
  }
  dev.off()

    