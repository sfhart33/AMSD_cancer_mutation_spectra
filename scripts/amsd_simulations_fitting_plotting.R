library(tidyverse)
library(ggpubr)

# read in results from two fitting runs and plot:
  results <- read.delim("../outputs/top_results_sigsX-40-2.tsv", header = FALSE) %>%
    rbind(read.delim("../outputs/top_results_CP-NPYR.tsv", header = FALSE)) %>%
    rbind(read.delim("../outputs/top_results_SBS2-40-CP-NPYR.tsv", header = FALSE))
  results
  colnames(results) <- c("signature", "mean_A", "mean_B", "p_ttest", "p_wilcox",
                         "present","n_signatures_present", "p_ttest_Bonf",
                         "p_wilcox_Bonf", "n_samples", "n_mutations", "additional_sig",
                         "frac_extra", "seed", "amsd_p", "replicate",  "param_set", "run_id")
  results_sigprof <- filter(results, additional_sig != "SBSX")
  
  results_sigprof2 <- results_sigprof %>%
    group_by(n_samples, n_mutations, additional_sig, frac_extra) %>%
    summarize(n_tests = n(),
              success_amsd = sum(amsd_p <= 0.05)/n(),
              success_amsd20 = sum(amsd_p <= 0.05/20)/n(),
              success_wilcox = sum(p_wilcox <= 0.05)/n(),
              success_wilcox20 = sum(p_wilcox <= 0.05/20)/n(),
              success_wilcoxBonf = sum(p_wilcox_Bonf <= 0.05)/n(),
              success_wilcoxBonf20 = sum(p_wilcox_Bonf <= 0.05/20)/n()
    )
  
  comparison_plot <- results_sigprof2 %>%
    filter(frac_extra > 0, additional_sig != "NPYR") %>%
    ggplot(aes(success_wilcoxBonf,
               success_amsd,
               color = factor(frac_extra, levels = c("0.02","0.05","0.1","0.2")), 
               size = n_samples))+
    geom_point(alpha = 0.75)+
    facet_grid(n_mutations ~ additional_sig) +
    scale_size_continuous(
      range = c(2, 5),
      breaks = c(5,25,125)
    )+
    geom_abline(slope = 1, linetype = "dashed")+
    guides(color = guide_legend(title = "Fraction extra mutations \nper exposure sample"),
           size = guide_legend(title = "Sample count \n(same # exposed\nand non-exposed)"))+
    ylab("Difference in spectra detected by AMSD \n(p<0.05, fraction of 50 simulations)")+
    xlab("Difference in signatures detected by wilcoxen rank sum test \n(p<0.05, bonferroni corrected for # signatures present, fraction of 50 simulations)")+
    # ggtitle("AMSD strength of detection in \nWES (50 mu/sample) and WGS (2500 mu/sample)")+
    theme_bw()
  comparison_plot
  ggsave("../outputs/amsd_simulations_comparison.png",
         plot = comparison_plot,
         width = 10,
         height = 5,
         units = "in"
  )

# previous plot of simulation result just with AMSD
  output <- readRDS("../outputs/amsd_simulation_output.rds")
  
  output
  output2 <- output %>%
    group_by(n_samples, n_mutations, exposure, extra_muts) %>%
    summarize(mean_p = mean(as.numeric(pvalue)),
              sd_p = sd(as.numeric(pvalue)),
              success_frac = sum(pvalue <= 0.05)/n())
  output2
  
  simulation_plot <- output2 %>%
    mutate(n_mutations = factor(n_mutations, levels = c("50", "2500"))) %>%
    ggplot(aes(x = factor(n_samples, levels = c("5","25","125","625")),
               y = success_frac,
               color = factor(extra_muts, levels = c("0.02","0.05","0.1","0.2")),
               group = factor(extra_muts, levels = c("0.02","0.05","0.1","0.2"))))+
    geom_point()+
    geom_line() +
    facet_grid(n_mutations ~ exposure) +
    guides(color = guide_legend(title = "Fraction extra mutations \nper exposure sample"))+
    xlab("Sample count (same # exposed and non-exposed)")+
    ylab("Difference detected \n(p<0.05, fraction of 100 simulations)")+
    # ggtitle("AMSD strength of detection in \nWES (50 mu/sample) and WGS (2500 mu/sample)")+
    theme_bw()
  simulation_plot

# merge and final plot
  merged_plot <- ggarrange(simulation_plot, comparison_plot, nrow=2, ncol=1,  labels = c("A","B"))
  merged_plot
  
  ggsave("../outputs/amsd_simulations_supp_figure.png",
         plot = merged_plot,
         width = 7,
         height = 8,
         units = "in"
  )
