input <- tcga_3mer[,2:97] %>%
  colSums()
input
plot_spectrum(input)
input+input/2
COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", "darkolivegreen3", "rosybrown2")


as.data.frame(input) %>%
  rownames_to_column(var = "m") %>%
  mutate(m = str_replace(m, "_", ">")) %>%
  separate(m, into = c("mut", "trinuc"), sep = "\\.", remove = FALSE) %>%
  separate(mut, into = c("from", "to"), sep = ">", remove = FALSE) %>%
  ggplot(aes(x=m, y = input, fill = mut))+
    geom_col()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values = COLORS)

plot1 <- as.data.frame(input) %>%
  rownames_to_column(var = "m") %>%
  mutate(m = str_replace(m, "_", ">")) %>%
  separate(m, into = c("mut", "trinuc"), sep = "\\.", remove = FALSE) %>%
  separate(mut, into = c("from", "to"), sep = ">", remove = FALSE) %>%
  ggplot(aes(x=trinuc, y = input, fill = mut))+
    geom_col()+
    geom_errorbar(aes(ymin=input*1.5, ymax=input*0.75),
                  width=.2,
                  position=position_dodge(.9))+
    scale_fill_manual(values = COLORS)+
    scale_y_continuous(expand = c(0,0)) +
    facet_grid(cols = vars(mut), scales = 'free')+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      panel.spacing = unit(0,'lines'),
      strip.text = element_blank(),
      aspect.ratio = 1.5,
      panel.grid.major.x = element_blank()
      )+
    xlab("Trinucleotide context")+
    ylab("Mutation count")
plot1
ggsave("../outputs/test_spectrum.png", plot1)
  getwd()
  
  
  
  mouse_carcinogen_spectra <- mouse_carcinogen_counts/rowSums(mouse_carcinogen_counts) # spectra
  mouse_carcinogen_spectra 
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
  samples <- pull(filter(sample_table, tissue == "LIVER", exposure == "SPONTANEOUS"), label)
  spectra1 <- colMeans(mouse_carcinogen_spectra[samples,], na.rm=TRUE)
  spectra1sd <- apply(mouse_carcinogen_spectra[samples,], 2, sd)
  spectra1
  spectra1sd
  
  default_df <- as.data.frame(input) %>%
    rownames_to_column(var = "m") %>%
    mutate(m = str_replace(m, "_", ">")) %>%
    separate(m, into = c("mut", "trinuc"), sep = "\\.", remove = FALSE) %>%
    separate(mut, into = c("from", "to"), sep = ">", remove = FALSE) %>%
    select(-input)
  default_df$spectra <- names(spectra1)
  default_df$spectra1 <- spectra1
  default_df$spectra1sd <- spectra1sd
  default_df$spectra <- names(spectra1)
  default_df %>%
    ggplot(aes(x=trinuc, y = spectra1, fill = mut))+
    geom_col()+
    geom_errorbar(aes(ymin=spectra1-spectra1sd, ymax=spectra1+spectra1sd),
                  width=.2,
                  position=position_dodge(.9))+
    scale_fill_manual(values = COLORS)+
    scale_y_continuous(expand = c(0,0)) +
    facet_grid(cols = vars(mut), scales = 'free')+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      panel.spacing = unit(0,'lines'),
      strip.text = element_blank(),
      aspect.ratio = 1.5,
      panel.grid.major.x = element_blank()
    )+
    xlab("Trinucleotide context")+
    ylab("Mutation count")
