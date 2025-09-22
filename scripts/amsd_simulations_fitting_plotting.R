library(tidyverse)
library(sigfit)

results <- read.delim("top_results_sigsX-40-2.tsv", header = FALSE) %>%
  rbind(read.delim("top_results_CP-NPYR.tsv", header = FALSE)) %>%
  rbind(read.delim("top_results_SBS2-40-CP-NPYR.tsv", header = FALSE))
results
colnames(results) <- c("signature", "mean_A", "mean_B", "p_ttest", "p_wilcox",
                       "present","n_signatures_present", "p_ttest_Bonf",
                       "p_wilcox_Bonf", "n_samples", "n_mutations", "additional_sig",
                       "frac_extra", "seed", "amsd_p", "replicate",  "param_set", "run_id")
results_sigprof <- filter(results, additional_sig != "SBSX")
results_sigprof

  output2 <- results_sigprof %>%
    group_by(n_samples, n_mutations, additional_sig, frac_extra) %>%
    summarize(n_tests = n(),
              success_amsd = sum(amsd_p <= 0.05)/n(),
              success_amsd20 = sum(amsd_p <= 0.05/20)/n(),
              success_wilcox = sum(p_wilcox <= 0.05)/n(),
              success_wilcox20 = sum(p_wilcox <= 0.05/20)/n(),
              success_wilcoxBonf = sum(p_wilcox_Bonf <= 0.05)/n(),
              success_wilcoxBonf20 = sum(p_wilcox_Bonf <= 0.05/20)/n()
    )
output2

comparison_plot <- output2 %>%
  filter(frac_extra > 0) %>%
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

######################################
# plot signatures
# load novel signatures
exp_signatures <- read.delim("../inputs/human_sbs96_filtered_v1_0.txt", 
                             sep = "\t", 
                             header = TRUE, 
                             check.names = FALSE) %>%
  separate(MutationType, into = c("a","b"), sep = "\\[", remove = FALSE) %>%
  separate(b, into = c("b","c"), sep = "\\]") %>%
  arrange(b,a,c) %>%
  select(-a,-b,-c) %>%
  column_to_rownames(var = "MutationType") %>%
  t() %>%
  as.data.frame()
colnames(exp_signatures) <- colnames(cosmic_signatures_v3.2)
cosmic_signatures_v3.2_new <- rbind(cosmic_signatures_v3.2, 
                                    CP = exp_signatures["cyclophosphamide_557117b73fe2",],
                                    NPYR = exp_signatures["n_nitrosopyrrolidine_ff630a9dde6d",])
tail(cosmic_signatures_v3.2_new)

cosmic_signatures_v3.2_new[c(c("CP", "NPYR", "SBS2", "SBS40")),] %>%
  plot_spectrum(pdf_path = "../outputs/spike-in_signatures.pdf")

