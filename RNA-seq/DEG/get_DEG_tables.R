library(tidyverse)
library(edgeR)

source("./bias_correcting_support_functions.R")

# GO + GV (WT vs cKO):

counts_GV = read.delim(file = "./data/GV_cKO_MLL34.counts.tsv", header = T)
counts_GO = read.delim(file = "./data/GOp14_cKO_MLL34.counts.tsv", header = T)

names(counts_GV) = c("gene_id",
                     paste(rep(x = c("WT", "KO"), each = 4),
                           rep(x = c("rep1", "rep2"), each = 2),
                           c("init", "deep"), sep = "_"))
names(counts_GO) = c("gene_id",
                     paste(rep(x = c("WT", "KO"), each = 6),
                           rep(x = c("rep1", "rep2", "rep3"), each = 2),
                           c("init", "deep"), sep = "_"))

res05 = perform_de_grouped_norm(list(counts_GV, counts_GO), fdr_thr=0.05, pair_pool = T)

res05[[1]]$tab_plt[, c(1,3,4,5,6,8,9,10,11,12)] %>% 
  write_delim(file = "./data/DE_edgeR.GV_pooled.WT_cKO.FC2_fdr05.tsv", delim = "\t", col_names = T)
res05[[2]]$tab_plt[, c(1,3,4,5,6,8,9,10,11,12)] %>% 
  write_delim(file = "./data/DE_edgeR.GO_p14_pooled.WT_cKO.FC2_fdr05.tsv", delim = "\t", col_names = T)

res05[[1]]$cpm_df %>% 
  write_delim(file = "./data/cpm.GV_pooled.WT_cKO.tsv", delim = "\t", col_names = T)
res05[[2]]$cpm_df %>% 
  write_delim(file = "./data/cpm.GO_p14_pooled.WT_cKO.tsv", delim = "\t", col_names = T)


# Mll3 + Mll4 + Mll34 cKO (WT vs cKO):

comp_cKO_Mll = c("cKO_MLL3", "cKO_MLL4", "cKO_MLL34")
dRNAseq_Mll34cKO = "./data/"
counts_Mll34cKO = lapply(comp_cKO_Mll, function(x)
  read.delim(file = file.path(dRNAseq_Mll34cKO, paste0("GV_", x, ".cKO_MLL_3_4_3.counts.tsv")), header = T,
             col.names = c("gene_id",
                           paste(rep(x = c("WT", "KO"), each = 3),
                                 c("rep1", "rep2", "rep3"), sep = "_")))
)

res05_Mll34cKO = perform_de_grouped_norm(counts_Mll34cKO, fdr_thr=0.05, log2FC_thr = log2(2))

for (i in 1:3){
  res05_Mll34cKO[[i]]$tab_plt[, c(1,3,4,5,6,8,9,10,11,12)] %>% 
    write_delim(file = paste0("./data/DE_edgeR.GV.Mll_3_4_34.WT_", comp_cKO_Mll[i],".FC2_fdr05.tsv"), delim = "\t", col_names = T)
  res05_Mll34cKO[[i]]$cpm_df %>%
    write_delim(file = paste0("./data/cpm.GV_", comp_cKO_Mll[i],".WT_cKO.tsv"), delim = "\t", col_names = T)
}

