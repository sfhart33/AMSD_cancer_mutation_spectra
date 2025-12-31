# New Supp Figure 5 in response to reviews

# First run amsd_v_sigs_comparisons.R at least through line 275
  # source("amsd_v_sigs_comparisons.R")
  # generates these:
    mouse_amsdsig_plot
    plot1
    plot2
    ox_plot

# then run to get gene/nongene counts:
  source("amsd_mouse_snvs_oxazepam.R")

# absolute values of mouse carcinogen exposures
    exp_withcounts <- as.data.frame(rowSums(mouse_carcinogen_counts)) %>%
      rownames_to_column(var = "label") %>%
      rename('rowSums(mouse_carcinogen_counts)' = "mut_count") %>%
      right_join(mexposuresig3) %>%
      mutate(absolute = value * mut_count)
    
    sigs <- exp_withcounts %>%
      filter(tissue == "LIVER",
             exposure %in% c("SPONTANEOUS","TCP","OXAZEPAM")) %>%
      group_by(name) %>%
      summarize(mean = mean(value)) %>%
      filter(mean > 0) %>%
      pull(name)
    
    exp_withcounts_plot <- exp_withcounts %>%
      filter(tissue == "LIVER",
             exposure %in% c("SPONTANEOUS","TCP","OXAZEPAM"),
             name %in% sigs) %>%
      mutate(exposure = recode(exposure,
                               "SPONTANEOUS" = "Spontaneous",
                               "OXAZEPAM" = "Oxazepam")) %>%
      ggplot(aes(x = as.numeric(rep), y = absolute, fill = name))+
      scale_fill_brewer(palette = "Set2")+
      geom_col()+
      facet_grid(~ exposure, scales = "free_x", space = "free_x") +
      scale_x_continuous(breaks = seq(0,20,1)) +
      theme_classic()+
      theme(legend.position="top",
            legend.title = element_text(hjust = 0.5))+
      guides(fill = guide_legend(title.position = "top"))+
      labs(x = "Tumor sample",
           y = "Signature fraction",
           fill = "Mutational signature")
    
    exp_withcounts_plot
    
# Spectra of gene vs nongene for spontaneous and oxazapam

    plot_mouse_spectra <- function(tis, exp, input_spectra){
      samples <- pull(filter(sample_table, tissue == tis, exposure == exp), label)
      spectra1 <- colMeans(input_spectra[samples,], na.rm=TRUE)
      spectra1sd <- apply(input_spectra[samples,], 2, sd, na.rm = TRUE)
      
      default_df <- readRDS("../inputs/default_spectrum_df.rds")
      
      default_df$spectra1 <- spectra1
      default_df$spectra1sd <- spectra1sd
      COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")
      COLORS <- c("steelblue", "gray40", "darkred", "gray60", "darkolivegreen4", "rosybrown")
      
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
        theme_classic()+
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6),
          panel.spacing = unit(0,'lines'),
          #aspect.ratio = 1.5,
          legend.position="none",
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(face = "bold", vjust = 0.5))+
        theme(
          strip.text = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          axis.title.x = element_blank()
        )+
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
    pms1 <- plot_mouse_spectra("LIVER", "SPONTANEOUS", spon_gene) +
      ylab("     ")
    pms2 <- plot_mouse_spectra("LIVER", "SPONTANEOUS", spon_nongene) +
      ylab("     ")
    pms3 <- plot_mouse_spectra("LIVER", "OXAZEPAM", ox_gene) +
      ylab("          Mutation Fraction")
    pms4 <- plot_mouse_spectra("LIVER", "OXAZEPAM", ox_nongene)  +
      ylab("     ")
    pms_all4 <- ggarrange(pms1, pms2, pms3, pms4, nrow = 4, ncol = 1,
                          labels = c("Spontaneous, genic",
                                     "Spontaneous, intergenic",
                                     "Oxazepam, genic",
                                     "Oxazepam, intergenic"),
                          label.x = 0.5,
                          label.y = 1,
                          hjust = 0.5)
    pms_all4
# full plot
    supp5 <- ggarrange(
      ggarrange(exp_withcounts_plot,pms_all4, nrow = 1,ncol = 2, widths = c(0.4,0.6), labels = c("A","C")),
      ggarrange(ox_plot,mouse_amsdsig_plot,nrow = 1,ncol = 2, widths = c(0.35,0.65), labels = c("B","D")),
      ggarrange(plot1,plot2,nrow = 1,ncol = 2, labels = c("E","F")),
      nrow = 3,
      ncol = 1)
    
    ggsave("../outputs/mouse_amsd_supp5_revision.png",
           plot = supp5,
           width = 10,
           height = 12,
           units = "in"
    )
    