library(tidyverse)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicRanges)
library(Biostrings)
library(sigfit)
library(mutspecdist)
library(ggrepel)

# setwd("\\\\gs-ddn2/gs-vol1/home/sfhart/github/AMSD_cancer_mutation_spectra/scripts")
snvs <- readRDS("../inputs/snvs.rds")
head(snvs)


#################################
#### gene annotations
# Load annotation
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_gr <- genes(txdb)     # gene bodies
exons_gr <- exons(txdb)     # exons
cds_gr <- cds(txdb)

# Make GRanges from your SNVs
snvs_gr <- snvs %>%
  makeGRangesFromDataFrame(seqnames.field = "chrom",
                           start.field = "pos",
                           end.field = "pos",
                           keep.extra.columns = TRUE)

## ---- Gene overlaps ----
hits_gene <- findOverlaps(snvs_gr, genes_gr)
snvs$in_gene <- FALSE
snvs$in_gene[queryHits(hits_gene)] <- TRUE

## ---- Exon overlaps ----
hits_exon <- findOverlaps(snvs_gr, exons_gr)
snvs$in_exon <- FALSE
snvs$in_exon[queryHits(hits_exon)] <- TRUE

## ---- CDS overlaps ----
hits_cds <- findOverlaps(snvs_gr, cds_gr)
snvs$in_cds <- FALSE
snvs$in_cds[queryHits(hits_cds)] <- TRUE

# Quick check
table(snvs$in_gene)
table(snvs$in_exon)
table(snvs$in_cds)
table(snvs$in_gene,snvs$in_exon)
table(snvs$in_gene,snvs$in_cds)
head(snvs)

##############################
# Trinucleotides
snvs_all <- snvs %>%
  mutate(trinuc = as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10,
                                      chrom,
                                      start = as.numeric(pos)-1,
                                      end =  as.numeric(pos) + 1))) %>%
  mutate(
    is_purine = ref %in% c("A", "G"),
    trinuc_std = if_else(is_purine,
                         as.character(reverseComplement(DNAStringSet(trinuc))),
                         trinuc),
    ref_std = if_else(is_purine,
                      as.character(reverseComplement(DNAStringSet(ref))),
                      ref),
    alt_std = if_else(is_purine,
                      as.character(reverseComplement(DNAStringSet(alt))),
                      alt),
    trinuc_mut = paste0(trinuc_std, ">", substr(trinuc_std, 1, 1), alt_std, substr(trinuc_std, 3, 3))
  )

##############################
# Tally up and put in order
make_snv_matrix <- function(snvs,
                            sample_prefix = NULL,
                            filter_col = NULL,
                            filter_val = TRUE,
                            rowsum = TRUE) {
  df <- snvs
  
  # optional: filter by sample name prefix
  if (!is.null(sample_prefix)) {
    df <- df %>%
      filter(startsWith(sample, sample_prefix))
  }
  
  # optional: filter by annotation column (e.g. in_exon, in_gene, in_cds)
  if (!is.null(filter_col)) {
    df <- df %>%
      filter(.data[[filter_col]] == filter_val)
  }
  
  # tally and order trinucleotide mutations
  output <- df %>%
    count(sample, trinuc_mut, name = "n_snvs") %>%
    separate(trinuc_mut, into = c("trinuc", "mut"), sep = ">") %>%
    mutate(
      first = substr(trinuc, 1, 1),
      middle = substr(trinuc, 2, 2),
      middle2 = substr(mut, 2, 2),
      third = substr(trinuc, 3, 3)
    ) %>%
    arrange(middle, middle2, first, third) %>%
    mutate(trinuc_mut = paste0(trinuc, ">", mut)) %>%
    dplyr::select(-first, -middle, -third, -trinuc, -mut, -middle2) %>%
    pivot_wider(
      names_from = trinuc_mut,
      values_from = n_snvs,
      values_fill = 0
    ) %>%
    column_to_rownames("sample")
  output/rowSums(output)
  if (rowsum == TRUE) {
    output/rowSums(output)
  } else {
    output
  }
}

spon_all <- make_snv_matrix(snvs_all, sample_prefix = "LIVER_SPON")
spon_gene <- make_snv_matrix(snvs_all, sample_prefix = "LIVER_SPON", filter_col = "in_gene")
spon_nongene <- make_snv_matrix(snvs_all, sample_prefix = "LIVER_SPON", filter_col = "in_gene",filter_val = FALSE)
spon_nonexon <- make_snv_matrix(snvs_all, sample_prefix = "LIVER_SPON", filter_col = "in_exon",filter_val = FALSE)

ox_all <- make_snv_matrix(snvs_all, sample_prefix = "LIVER_OX")
ox_gene <- make_snv_matrix(snvs_all, sample_prefix = "LIVER_OX", filter_col = "in_gene")
ox_nongene <- make_snv_matrix(snvs_all, sample_prefix = "LIVER_OX", filter_col = "in_gene",filter_val = FALSE)
ox_nonexon <- make_snv_matrix(snvs_all, sample_prefix = "LIVER_OX", filter_col = "in_exon",filter_val = FALSE)

plot_spectrum(colSums(spon_nongene))
plot_spectrum(colSums(ox_nongene))

amsd_all <- amsd(spon_all,
                    ox_all,
                    mean_or_sum = "mean",
                    n_sim = 10000,
                    seed = 1234)
amsd_gene <- amsd(spon_gene,
                  ox_gene,
                  mean_or_sum = "mean",
                  n_sim = 10000,
                  seed = 1234)

amsd_nongene <- amsd(spon_nongene,
                     ox_nongene,
                     mean_or_sum = "mean",
                     n_sim = 10000,
                     seed = 1234)
amsd_nonexon <- amsd(spon_nonexon,
                     ox_nonexon,
                     mean_or_sum = "mean",
                     n_sim = 10000,
                     seed = 1234)
plot_amsd_histogram(amsd_all)
plot_amsd_histogram(amsd_nongene)
plot_amsd_histogram(amsd_gene)
plot_amsd_histogram(amsd_nonexon)

make_snv_matrix(snvs_all, sample_prefix = "LIVER_SPON", filter_col = "in_gene", rowsum = FALSE) %>% rowSums() %>% mean()
make_snv_matrix(snvs_all, sample_prefix = "LIVER_SPON", filter_col = "in_gene",filter_val = FALSE, rowsum = FALSE) %>% rowSums() %>% mean()

amsd(make_snv_matrix(snvs_all, sample_prefix = "LIVER_SPON", rowsum = FALSE),
     make_snv_matrix(snvs_all, sample_prefix = "LIVER_OX", rowsum = FALSE),
     mean_or_sum = "sum",
     n_sim = 10000,
     seed = 1234) %>%
  plot_amsd_histogram()

################
spon_all
ox_all %>%
  mutate(sample = paste0("ox_", row_number()), group = "ox") %>%
  pivot_longer(-c(sample, group), names_to = "trinuc", values_to = "freq")

################
# spectra comparison

ox_long <- ox_all %>%
  as.data.frame() %>%
  mutate(sample = paste0("ox_", row_number()), group = "ox") %>%
  pivot_longer(-c(sample, group), names_to = "trinuc", values_to = "freq")

spon_long <- spon_all %>%
  as.data.frame() %>%
  mutate(sample = paste0("spon_", row_number()), group = "spon") %>%
  pivot_longer(-c(sample, group), names_to = "trinuc", values_to = "freq")

dat_long <- bind_rows(ox_long, spon_long)

# Run Wilcoxon per trinucleotide
results_ox <- dat_long %>%
  group_by(trinuc) %>%
  summarise(
    pval = t.test(freq[group == "ox"], freq[group == "spon"])$p.value,
    mean_ox = mean(freq[group == "ox"], na.rm = TRUE),
    mean_spon = mean(freq[group == "spon"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    effect = mean_ox - mean_spon,            # raw difference
    log2fc = log2((mean_ox + 1e-6) / (mean_spon + 1e-6)),  # fold change
    padj = p.adjust(pval, method = "BH")     # FDR correction
  )
results_ox
# Volcano plot (effect size vs -log10 p-value)
ggplot(results_ox, aes(x = effect, y = -log10(pval), label = trinuc)) +
  geom_point(aes(color = padj < 0.1)) +
  geom_label_repel()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Mean difference (oxazempam - spontaneous)",
    y = "-log10(p-value)",
    color = "FDR < 0.1"
  )
