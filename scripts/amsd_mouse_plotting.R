library(tidyverse)

# load data
  mouse_amsd_output <- readRDS("../outputs/mouse_amsd_output.rds")
  perms <- readRDS("../outputs/mouse_amsd_perms.rds")
  mexposuresig <- readRDS("../inputs/mouse_exposuresig.rds")
  mexposure <- readRDS("../inputs/mouse_exposure.rds")
  mouse_carcinogen_counts <- readRDS("../inputs/mouse_carcinogen_spectra.rds") # counts
  mouse_carcinogen_spectra <- mouse_carcinogen_counts/rowSums(mouse_carcinogen_counts) # spectra
  

  mouse_amsd_output$padj_BH <- p.adjust(mouse_amsd_output$pvalues, method="BH")
  mouse_amsd_output$padj_BY <- p.adjust(mouse_amsd_output$pvalues, method="BY")
  ggplot(mouse_amsd_output, aes(-log10(pvalues),-log10(padj_BH)))+
    geom_point()+
    geom_hline(yintercept = -log10(0.05))+
    stat_smooth(method="lm",se=F)+
    ggtitle("Benjamini-Hochberg vs unadjusted pvalues from mouse comparison")
  x <- -log10(mouse_amsd_output$padj_BH)
  y <- -log10(mouse_amsd_output$pvalues)
  BH_regression <- lm(y ~ x)
  BH_threshold <- -log10(0.05) * coef(BH_regression)[2]
  
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
    geom_hline(yintercept = BH_threshold, linetype = "dashed")+
    geom_text(aes(x=0.2, y = -log10(0.05)+0.1), label = "p=0.05")+
    geom_text(aes(x=0.2, y = BH_threshold + 0.1), label = "Benjamini-Hochberg FDR=0.05")+
    geom_text(aes(x=0.2, y = -log10(0.05/nrow(mouse_amsd_output))+0.1), label = "Bonferroni FDR=0.05")+
    theme_classic()+
    theme(legend.position="top")+
    guides(
      color = guide_legend(title.position = "top", title.hjust = 0.5, direction = "horizontal"),
      size = guide_legend(title.position = "top", title.hjust = 0.5, direction = "horizontal")
    )+
    xlim(0,0.25)+
    scale_size_continuous(range = c(2, 5))+
    labs(x="Cosine distance",
         y="-log10(p-value)", 
         color="Tumor type",
         size="Exposed tumor count")
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
  
  exp_plot <- mexposuresig2 %>%
    filter(tissue == "LIVER",
           exposure %in% c("SPONTANEOUS","TCP","OXAZEPAM"),
           name %in% sigs) %>%
      ggplot(aes(x = as.numeric(rep), y = value, fill = name))+
      geom_col()+
      facet_grid(~ exposure, scales = "free_x", space = "free_x") +
      scale_x_continuous(breaks = seq(0,20,1)) +
      theme_classic()+
    theme(legend.position="top",
          legend.title = element_text(hjust = 0.5))+
    guides(fill = guide_legend(title.position = "top"))+
    labs(x = "Tumor sample",
         y = "Signature fraction",
         fill = "Signature")
  TCP_plot <- mexposuresig2 %>%
    filter(tissue == "LIVER",
           exposure %in% c("SPONTANEOUS","TCP"),
           name %in% sigs) %>%
    ggplot(aes(x = name, y = value, color = exposure))+
    geom_boxplot(outliers = FALSE)+
    geom_point(position=position_jitterdodge(jitter.width = 0.1))+
    scale_color_manual(values = c("grey", "black"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="top",
          legend.title = element_text(hjust = 0.5))+
    guides(fill = guide_legend(title.position = "top"))+
    labs(x = "Mutational signature",
         y = "Signature fraction",
         color = "Exposure")
  oxa_plot <- mexposuresig2 %>%
    filter(tissue == "LIVER",
           exposure %in% c("SPONTANEOUS","OXAZEPAM"),
           name %in% sigs) %>%
    ggplot(aes(x = name, y = value, color = exposure))+
    geom_boxplot(outliers = FALSE)+
    geom_point(position=position_jitterdodge(jitter.width = 0.1))+
    scale_color_manual(values = c("black", "grey"))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="top",
          legend.title = element_text(hjust = 0.5))+
    guides(fill = guide_legend(title.position = "top"))+
    labs(x = "Mutational signature",
         y = "Signature fraction",
         color = "Exposure")
  
# spectra plots

  samples <- pull(filter(sample_table, tissue == "LIVER", exposure == "SPONTANEOUS"), label)
  spectra1 <- colMeans(mouse_carcinogen_spectra[samples,], na.rm=TRUE)
  spectra1sd <- apply(mouse_carcinogen_spectra[samples,], 2, sd)

  default_df <- readRDS("../inputs/default_spectrum_df.rds")
  
  default_df$spectra1 <- spectra1
  default_df$spectra1sd <- spectra1sd
  COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
  
  default_df %>%
    ggplot(aes(x=trinuc, y = spectra1, fill = mut))+
    geom_col()+
    geom_errorbar(aes(ymin=spectra1-spectra1sd, ymax=spectra1+spectra1sd),
                  width=.2,
                  position=position_dodge(.9))+
    scale_fill_manual(values = COLORS)+
    scale_y_continuous(expand = c(0,0)) +
    ylim(0,0.05)+
    facet_grid(cols = vars(mut), scales = 'free')+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
      panel.spacing = unit(0,'lines'),
      strip.text = element_blank(),
      aspect.ratio = 1.5,
      panel.grid.major.x = element_blank()
    )+
    xlab("Trinucleotide context")+
    ylab("Mutation fraction")
  
  
  plot_mouse_spectra <- function(tis, exp){
    samples <- pull(filter(sample_table, tissue == tis, exposure == exp), label)
    spectra1 <- colMeans(mouse_carcinogen_spectra[samples,], na.rm=TRUE)
    spectra1sd <- apply(mouse_carcinogen_spectra[samples,], 2, sd, na.rm = TRUE)
    
    default_df <- readRDS("../inputs/default_spectrum_df.rds")
    
    default_df$spectra1 <- spectra1
    default_df$spectra1sd <- spectra1sd
    COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
    
    default_df %>%
      ggplot(aes(x=trinuc, y = spectra1, fill = mut))+
      geom_col()+
      geom_errorbar(aes(ymin=spectra1-spectra1sd, ymax=spectra1+spectra1sd),
                    width=.2,
                    position=position_dodge(.9))+
      scale_fill_manual(values = COLORS)+
      scale_y_continuous(expand = c(0,0)) +
      #ylim(0,0.05)+
      facet_grid(cols = vars(mut), scales = 'free')+
      theme_bw()+
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6),
        panel.spacing = unit(0,'lines'),
        #aspect.ratio = 1.5,
        legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(face = "bold", vjust = 0.5))+
      # theme(
      #   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6),
      #   panel.spacing = unit(0,'lines'),
      #   strip.text = element_blank(),
      #   aspect.ratio = 1.5,
      #   panel.grid.major.x = element_blank(),
      #   panel.border = element_blank(),
      #   legend.position="none"
      # )+
      xlab("Trinucleotide context")+
      ylab("Mutation fraction") %>%
      return()
  }
  p1 <- plot_mouse_spectra("LIVER","SPONTANEOUS")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank())
  p2 <- plot_mouse_spectra("LIVER","1,2,3_TRICHLOROPROPANE")+
    theme(strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank())
  p3 <- plot_mouse_spectra("LIVER","OXAZEPAM")+
    theme(strip.text = element_blank())+
    xlab("Trinucleotide context")
  spectra_plots <- ggarrange(p1, p2, p3, nrow = 3, ncol = 1, heights=c(0.4,0.3,0.4))
  top_row = ggarrange(mouse_volcano,NULL, spectra_plots, nrow=3, ncol=1, heights=c(0.55,0.02,0.4), labels = c("A","","B"))
  bottom_row <- ggarrange(exp_plot, NULL, TCP_plot,NULL, oxa_plot, nrow=5, ncol=1, heights=c(0.3,0.01,0.3,0.01,0.3), labels = c("C","","D","","E"))
  mouse_plot <- ggarrange(top_row, bottom_row, nrow=1, ncol=2, widths = c(0.6,0.4))
  mouse_plot
  ggsave("../outputs/Figure2.png",
         plot = mouse_plot,
         width = 10,
         height = 10,
         units = "in")
  ggsave("../outputs/Figure2.svg",
         plot = mouse_plot,
         width = 10,
         height = 10,
         units = "in")
  
  
  # Suppfigs
  filter(mouse_amsd_output, padj_BH < 0.05, tissue == "LUNG") 
  sigs <- mexposuresig2 %>%
    filter(tissue == "LUNG",
           exposure %in% c("SPONTANEOUS","COBALT_METAL","ISOBUTYL_NITRITE","VINYLIDENE_CHLORIDE")) %>%
    group_by(name) %>%
    summarize(mean = mean(value)) %>%
    filter(mean > 0) %>%
    pull(name)
  lung1 <- plot_mouse_spectra("LUNG","SPONTANEOUS")+
    ggtitle("SPONTANEOUS")
  lung2 <- plot_mouse_spectra("LUNG","COBALT_METAL")+
    ggtitle("COBALT_METAL")
  lung3 <- plot_mouse_spectra("LUNG","ISOBUTYL_NITRITE")+
    ggtitle("ISOBUTYL_NITRITE")
  lung4 <- plot_mouse_spectra("LUNG","VINYLIDENE_CHLORIDE")+
    ggtitle("VINYLIDENE_CHLORIDE")
  lung_plot <- mexposuresig2 %>%
    filter(tissue == "LUNG",
           exposure %in% c("SPONTANEOUS","COBALT_METAL","ISOBUTYL_NITRITE","VINYLIDENE_CHLORIDE"),
           name %in% sigs) %>%
    ggplot(aes(x = name, y = value, color = exposure))+
    geom_boxplot(outliers = FALSE)+
    geom_point(position=position_jitterdodge(jitter.width = 0.1))+
    theme_classic()+
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    #       legend.position="top",
    #       legend.title = element_text(hjust = 0.5))+
    guides(fill = guide_legend(title.position = "top"))+
    labs(x = "Mutational signature",
         y = "Signature fraction",
         color = "Exposure")
  lung_supp <- ggarrange(lung1,lung2,lung3,lung4, lung_plot, nrow=5, ncol=1, heights = c(0.15,0.15,0.15,0.15,0.4))
  lung_supp
  
  
  filter(mouse_amsd_output, padj_BH < 0.05, tissue == "LIVER") 
  sigs <- mexposuresig2 %>%
    filter(tissue == "LIVER",
           exposure %in% c("SPONTANEOUS",
                           "VINYLIDENE_CHLORIDE",
                           "ANTHRAQUINONE",
                           "BROMOCHLOROACETIC_ACID",
                           "CUMENE",
                           "DE-71",
                           "FURAN",
                           "PRIMACLONE")) %>%
    group_by(name) %>%
    summarize(mean = mean(value)) %>%
    filter(mean > 0) %>%
    pull(name)
  liver1 <- plot_mouse_spectra("LIVER","SPONTANEOUS")+
    ggtitle("SPONTANEOUS")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank())
  liver2 <- plot_mouse_spectra("LIVER","VINYLIDENE_CHLORIDE")+
    ggtitle("VINYLIDENE_CHLORIDE")+
    theme(strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank())
  liver3 <- plot_mouse_spectra("LIVER","ANTHRAQUINONE")+
    ggtitle("ANTHRAQUINONE")+
    theme(strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank())
  liver4 <- plot_mouse_spectra("LIVER","BROMOCHLOROACETIC_ACID")+
    ggtitle("BROMOCHLOROACETIC_ACID")+
    theme(strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank())
  liver5 <- plot_mouse_spectra("LIVER","CUMENE")+
    ggtitle("CUMENE")+
    theme(strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank())
  liver6 <- plot_mouse_spectra("LIVER","DE-71")+
    ggtitle("DE-71")+
    theme(strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank())
  liver7 <- plot_mouse_spectra("LIVER","FURAN")+
    ggtitle("FURAN")+
    theme(strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank())
  liver8 <- plot_mouse_spectra("LIVER","PRIMACLONE")+
    ggtitle("PRIMACLONE")+
    theme(strip.text = element_blank())+
    xlab("Trinucleotide context")
  liver_plot <- mexposuresig2 %>%
    filter(tissue == "LIVER",
           exposure %in% c("SPONTANEOUS",
                           "VINYLIDENE_CHLORIDE",
                           "ANTHRAQUINONE",
                           "BROMOCHLOROACETIC_ACID",
                           "CUMENE",
                           "DE-71",
                           "FURAN",
                           "PRIMACLONE"),
           name %in% sigs) %>%
    ggplot(aes(x = name, y = value, color = exposure))+
    geom_boxplot(outliers = FALSE)+
    geom_point(position=position_jitterdodge(jitter.width = 0.1))+
    theme_classic()+
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    #       legend.position="top",
    #       legend.title = element_text(hjust = 0.5))+
    guides(fill = guide_legend(title.position = "top"))+
    labs(x = "Mutational signature",
         y = "Signature fraction",
         color = "Exposure")
  liver_supp <- ggarrange(liver1,liver2,liver3,liver4,liver5,liver6,liver7,liver8,
                          liver_plot, nrow=9, ncol=1, heights = c(0.125,0.1,0.1,0.1,0.1,0.1,0.1,0.125,0.45))
  liver_supp
  ggsave("../outputs/mouse_liver_supp.png",
         plot = liver_supp,
         width = 12,
         height = 16,
         units = "in")
  ggsave("../outputs/mouse_liver_supp.svg",
         plot = liver_supp,
         width = 12,
         height = 16,
         units = "in")
  ggsave("../outputs/mouse_lung_supp.png",
         plot = lung_supp,
         width = 12,
         height = 16,
         units = "in")
  ggsave("../outputs/mouse_lung_supp.svg",
         plot = lung_supp,
         width = 12,
         height = 16,
         units = "in")