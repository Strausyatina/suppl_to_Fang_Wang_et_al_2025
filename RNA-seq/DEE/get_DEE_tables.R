## =============================================================== 
## Supplementary script for Fang et al. paper.
## Differential exon expression (DEE) analysis.
## Author: A. Mendelevich
## ===============================================================

library(tidyverse)
library(edgeR)

source("../DEG/bias_correcting_support_functions.R")

## ===============================================================
## Load merged-exon annotation and counts:

exons_merged_info = read.delim(
  file = "./data/genes.exons_merged.bed",
  header = FALSE,
  col.names = c(
    "chr", "st", "en", "strand", "gene_names", "gene_ids", "tr_names",
    "exon_ids", "genes_num", "trs_num", "min_exon_i"
  )
)
exons_merged_info$tss = exons_merged_info$min_exon_i == 1

merged_exons_tab = read.delim(
  file = "./data/RNAseq_GV_pooled_seqs_2WT_2KO.merged_exons_counts.bed.tsv",
  header = FALSE,
  col.names = c(
    "chr", "st", "en", "WT_rep1", "WT_rep2", "KO_rep1", "KO_rep2",
    "strand", "gene_names", "gene_ids", "tr_names", "exon_ids",
    "n_genes", "n_trs"
  )
)

merged_exons_tab = merge(
  merged_exons_tab,
  exons_merged_info[, c("chr", "st", "en", "tss")],
  by = c("chr", "st", "en"),
  all.x = TRUE
)

focus_chromosomes = c(as.character(1:19), "X")
merged_exons_tab = merged_exons_tab[merged_exons_tab$chr %in% focus_chromosomes, ]

libsizes = colSums(merged_exons_tab[, c("WT_rep1", "WT_rep2", "KO_rep1", "KO_rep2")])

## ===============================================================
## Run DEE:

counts_dee = data.frame(
  id = paste(merged_exons_tab$chr, merged_exons_tab$st, merged_exons_tab$en, sep = "_"),
  merged_exons_tab[, c("WT_rep1", "WT_rep2", "KO_rep1", "KO_rep2")]
)

res05 = perform_de_grouped_norm(
  list(counts_dee),
  fdr_thr = 0.05,
  thr_cpm = 0.25,
  n_cpm = 2
)

## ===============================================================
## Export significant DEE table:

dee_sig = res05[[1]]$tab_plt[, c(1, 3, 4, 5, 6, 8, 9, 10, 11, 12)]

dee_sig_coords = dee_sig$id %>%
  strsplit(split = "_") %>%
  do.call(rbind, .) %>%
  as.data.frame(stringsAsFactors = FALSE)

colnames(dee_sig_coords) = c("chr", "st", "en")
dee_sig = cbind(dee_sig_coords, dee_sig)
dee_sig$st = as.integer(dee_sig$st)
dee_sig$en = as.integer(dee_sig$en)

write_delim(
  dee_sig,
  file = "./data/DEE_edgeR.GV_pooled.WT_cKO.FC2_fdr05.tsv",
  delim = "\t",
  col_names = TRUE
)

## =============================================================== ##
## Export full annotated table with CPM and RPKM:

dee_full = merge(
  merged_exons_tab,
  dee_sig,
  by = c("chr", "st", "en"),
  all.x = TRUE
)

dee_full[, c("WT_rep1_cpm", "WT_rep2_cpm", "KO_rep1_cpm", "KO_rep2_cpm")] =
  sweep(
    dee_full[, c("WT_rep1", "WT_rep2", "KO_rep1", "KO_rep2")],
    2,
    libsizes,
    "/"
  ) * 1e6

dee_full[, c("WT_rep1_rpkm", "WT_rep2_rpkm", "KO_rep1_rpkm", "KO_rep2_rpkm")] =
  sweep(
    dee_full[, c("WT_rep1", "WT_rep2", "KO_rep1", "KO_rep2")],
    2,
    libsizes,
    "/"
  ) * 1e9 / (dee_full$en - dee_full$st)

write_delim(
  dee_full,
  file = "./data/RNAseq_GV_pooled_seqs_2WT_2KO.merged_exons_counts.DEE.tsv",
  delim = "\t",
  col_names = TRUE
)