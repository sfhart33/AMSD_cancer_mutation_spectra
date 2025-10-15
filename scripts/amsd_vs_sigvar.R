library(tidyverse)
library(ggrepel)
library(ggpubr)

# load data
  comparison <- read.delim("../inputs/amsd_vs_sigvar.txt")
  comparison[1,"sigvar_p"] <- 0.001 # replace zero with limit of detection
  comparison
  sigvar_p <- readRDS("../inputs/mutsig_carcinogens_mice_bootstrap_p_vals.RDS")
  sigvar <- readRDS("../inputs/mutsig_carcinogens_mice_sigvar.RDS")
  sigvar_grouped <- sigvar %>%
    dplyr::group_by(chemical, Tissue) %>%
    dplyr::summarise(n = n(),
              mean_within = mean(within)) %>%
    filter(chemical != "Spontaneous")
  mouse_amsd_output <- readRDS("../outputs/mouse_amsd_output.rds")

# as.data.frame(sigvar_grouped)
# mouse_amsd_output
# 
# sigvar_clean <- sigvar_grouped %>%
#   mutate(
#     exposure_std = chemical %>%
#       toupper() %>%
#       str_replace_all("\\s+", "_") %>%
#       str_replace_all("\\n", "_") %>%
#       str_replace_all("_+", "_"),
#     tissue_std = toupper(Tissue)
#   )
# sigvar_clean
# 
# mouse_clean <- mouse_amsd_output %>%
#   mutate(
#     exposure_std = exposure %>%
#       toupper() %>%
#       str_replace_all("\\s+", "_") %>%
#       str_replace_all("\\n", "_") %>%
#       str_replace_all("_+", "_"),
#     tissue_std = toupper(tissue)
#   ) %>%
#   arrange(exposure_std, tissue)
# mouse_clean
# 
# merged <- sigvar_clean %>%
#   inner_join(mouse_clean,
#              by = c("exposure_std", "tissue_std"))
# 
# # clean exposures and tissue in sigvar_p
# sigvar_p_clean <- sigvar_p %>%
#   mutate(
#     exposure_std = chemical %>%
#       toupper() %>%
#       str_replace_all("\\s+", "_") %>%
#       str_replace_all("\\n", "_") %>%
#       str_replace_all("_+", "_"),
#     tissue_std = toupper(Tissue)
#   ) %>%
#   select(exposure_std, tissue_std, sigvar_p = mean_within_sample_diversity)
# 
#   
  
  
  fit_stats <- comparison %>%
    group_by(tissue) %>%
    do({
      fit <- lm(mean_within ~ cosine_dist, data = .)
      tibble(
        slope = coef(fit)[2],
        r2 = summary(fit)$r.squared,
        pval = summary(fit)$coefficients[2,4]
      )
    }) %>%
    ungroup() %>%
    mutate(
      label = sprintf("slope = %.3f\nR² = %.2f\np = %.3g", slope, r2, pval)
    )
  
  fit_stats
  
plot1 <- comparison %>%
  ggplot(aes(cosine_dist, mean_within, color = tissue, label = exposure)) +
  geom_point()+
  geom_label_repel()+
  #geom_smooth(method = "lm", color = "black") +
  geom_smooth(data = filter(comparison, tissue == "LIVER"),method = "lm", se = TRUE, color = "red") +
  geom_smooth(data = filter(comparison, tissue == "LUNG"),method = "lm", se = TRUE, color = "blue") +
  xlab("mean cosine distance")+
  ylab("mean within-sample diversity")+
  theme_classic()

plot2 <- comparison %>%
  ggplot(aes(-log10(amsd_p), -log10(sigvar_p), color = tissue, label = exposure)) +
  geom_point()+
  geom_label_repel()+
  #geom_smooth(method = "lm", color = "black") +
  geom_smooth(data = filter(comparison, tissue == "LIVER"),method = "lm", se = TRUE, color = "red") +
  geom_smooth(data = filter(comparison, tissue == "LUNG"),method = "lm", se = TRUE, color = "blue") +
  xlab("p-value from AMSD comparison (-log10)")+
  ylab("p-value from sigvar mean \nwithin-sample diversity bootstraps (-log10)")+
  theme_classic()+
  coord_cartesian(ylim = c(0, 3.5)) 

plot12 <- ggarrange(plot1, plot2, nrow=2, ncol=1,  labels = c("A","B"))
plot12

ggsave("../outputs/amsd_vs_sigvar_supp_figures.png",
       plot = plot12,
       width = 7,
       height = 8,
       units = "in"
)

  