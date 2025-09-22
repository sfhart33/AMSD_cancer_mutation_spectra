library(tidyverse)

sigvar_p <- readRDS("../inputs/mutsig_carcinogens_mice_bootstrap_p_vals.RDS")
sigvar <- readRDS("../inputs/mutsig_carcinogens_mice_sigvar.RDS")
sigvar_grouped <- sigvar %>%
  dplyr::group_by(chemical, Tissue) %>%
  dplyr::summarise(n = n(),
            mean_within = mean(within)) %>%
  filter(chemical != "Spontaneous")
mouse_amsd_output <- readRDS("../outputs/mouse_amsd_output.rds")


as.data.frame(sigvar_grouped)
mouse_amsd_output

sigvar_clean <- sigvar_grouped %>%
  mutate(
    exposure_std = chemical %>%
      toupper() %>%
      str_replace_all("\\s+", "_") %>%
      str_replace_all("\\n", "_") %>%
      str_replace_all("_+", "_"),
    tissue_std = toupper(Tissue)
  )
sigvar_clean

mouse_clean <- mouse_amsd_output %>%
  mutate(
    exposure_std = exposure %>%
      toupper() %>%
      str_replace_all("\\s+", "_") %>%
      str_replace_all("\\n", "_") %>%
      str_replace_all("_+", "_"),
    tissue_std = toupper(tissue)
  ) %>%
  arrange(exposure_std, tissue)
mouse_clean

merged <- sigvar_clean %>%
  inner_join(mouse_clean,
             by = c("exposure_std", "tissue_std"))

# clean exposures and tissue in sigvar_p
sigvar_p_clean <- sigvar_p %>%
  mutate(
    exposure_std = chemical %>%
      toupper() %>%
      str_replace_all("\\s+", "_") %>%
      str_replace_all("\\n", "_") %>%
      str_replace_all("_+", "_"),
    tissue_std = toupper(Tissue)
  ) %>%
  select(exposure_std, tissue_std, sigvar_p = mean_within_sample_diversity)

# now merge with the big merged table
merged_final <- merged %>%
  left_join(sigvar_p_clean, by = c("exposure_std", "tissue_std"))

data.frame(amsd = mouse_clean$exposure_std, sigvar = sigvar_clean$exposure_std)
  