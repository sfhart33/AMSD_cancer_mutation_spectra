library(tidyverse)

# load data
  mouse_amsd_output <- readRDS("../outputs/mouse_amsd_output.rds")
  perms <- readRDS("../outputs/mouse_amsd_perms.rds")
  mexposuresig <- readRDS("../inputs/mouse_exposuresig.rds")
  mexposure <- readRDS("../inputs/mouse_exposure.rds")
  mouse_carcinogen_counts <- readRDS("../inputs/mouse_carcinogen_spectra.rds") # counts
  mouse_carcinogen_spectra <- mouse_carcinogen_counts/rowSums(mouse_carcinogen_counts) # spectra
  
  # list of samples and attributes
  sample_table <- mouse_carcinogen_spectra %>%
    rownames_to_column(var = "label") %>%
    select(label) %>%
    separate(label, c("tissue","rest"), sep = "_", extra = "merge", remove = FALSE) %>%
    separate(rest, 
             into = c("exposure", "rep"), 
             sep = "_(?=\\d+$)", # Regex to split based on underscores near the end
             extra = "merge",            # Handle extra columns by merging into one
             fill = "right")             # Handle missing parts gracefully
  
  
# volcano plot summary of everything together
  mouse_volcano <- mutate(mouse_amsd_output, log10pval = -log10(pvalues)) %>%
    ggplot()+
    geom_point(aes(x=cosines, y = log10pval, color = tissue, size = n))+
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
    geom_hline(yintercept = -log10(0.05/nrow(mouse_amsd_output)), linetype = "dashed")+
    #geom_hline(yintercept = -log10(1/reps))+
    geom_text(aes(x=0.225, y = (-log10(0.05/nrow(mouse_amsd_output))+0.1)), label = "FDR=0.05")+
    geom_text(aes(x=0.225, y = (-log10(0.05)+0.1)), label = "p=0.05")+
    #geom_text(aes(x=0.225, y = (-log10(1/reps)+0.1)), label = "theoretical max")+
    theme_classic()+
    xlim(0,0.25)+
    scale_size_continuous(range = c(2, 5))+
    labs(x="Cosine distance",
         y="-log10(p-value)", 
         color="Tumor type",
         size="Exposed\ntumor N")
  mouse_volcano
  ggsave("../outputs/mouse_amsd_output.png",
         plot = mouse_volcano,
         width = 5, height = 5, units = "in")

  mouse_amsd_output %>%
    arrange(pvalues)

# violin plots of AMSD null distributions
  perms2 <- perms %>%
    column_to_rownames(var = "rep") %>%
    pivot_longer(cols = everything()) %>%
    separate(name, into = c("tissue", "exposure"), sep = "\\.", remove = FALSE)
  quantiles <- group_by(perms2, tissue, exposure) %>%
    summarise(p0.05 = quantile(value, probs = 0.95),
              FDR = quantile(value, probs = 1-(0.05/nrow(mouse_amsd_output))))
  mouse_violin <- perms2 %>%
    ggplot(aes(x = exposure, y = value))+
    geom_violin(adjust =0.5, scale = "width")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    xlab("")+
    ylab("Cosine distance")+
    geom_point(data = mouse_amsd_output,
               aes(x = exposure, y = cosines))+
    geom_point(data = quantiles,
               aes(x = exposure, y = p0.05), shape = 95, size = 5)+
    geom_point(data = quantiles,
               aes(x = exposure, y = FDR), shape = 95, size = 5)+
    facet_grid(rows = vars(tissue))
  mouse_violin
  ggsave("../outputs/mouse_amsd_output_violin.png",
         plot = mouse_violin)
  
# signature plots
  mexposuresigX <- mexposuresig/rowSums(mexposuresig) 
  mexposuresig2 <- mexposuresigX %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    pivot_longer(cols = -sample) %>%
    mutate(label = str_replace_all(sample, " ", "_")) %>%
    left_join(sample_table) %>%
    mutate(exposure = (replace(exposure, exposure == "1,2,3_TRICHLOROPROPANE", "TCP")))

  
  sigs <- mexposuresig2 %>%
    filter(tissue == "LIVER",
           exposure %in% c("SPONTANEOUS","TCP","OXAZEPAM")) %>%
    group_by(name) %>%
    summarize(mean = mean(value)) %>%
    filter(mean > 0) %>%
    pull(name)
  
  mexposuresig2 %>%
  filter(tissue == "LIVER",
           exposure %in% c("SPONTANEOUS","TCP","OXAZEPAM"),
           name %in% sigs) %>%
    ggplot(aes(x = as.numeric(rep), y = value, fill = name))+
    geom_col()+
    facet_grid(~ exposure, scales = "free_x", space = "free_x") +
    scale_x_continuous(breaks = seq(0,20,1)) +
    theme_classic()+
    xlab("Tumor sample")+
    ylab("Signature fraction")
  mexposuresig2 %>%
    filter(tissue == "LIVER",
           exposure %in% c("SPONTANEOUS","TCP"),
           name %in% sigs) %>%
    ggplot(aes(x = name, y = value, color = exposure))+
    geom_boxplot(outliers = FALSE)+
    geom_point(position=position_jitterdodge(jitter.width = 0.1))+
    scale_color_manual(values = c("grey", "black"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme_classic()+
    xlab("Mutational signature")+
    ylab("Signature fraction")
  mexposuresig2 %>%
    filter(tissue == "LIVER",
           exposure %in% c("SPONTANEOUS","OXAZEPAM"),
           name %in% sigs) %>%
    ggplot(aes(x = name, y = value, color = exposure))+
    geom_boxplot(outliers = FALSE)+
    geom_point(position=position_jitterdodge(jitter.width = 0.1))+
    scale_color_manual(values = c("black", "grey"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme_classic()+
    xlab("Mutational signature")+
    ylab("Signature fraction")
  