# Differential exon expression (DEE)

This folder contains the scripts and intermediate files used for differential exon expression analysis in Fang Wang et al.

## Pre-processing: merged exon annotation and count table creation

Merged exon intervals were generated from the Ensembl v94 / M19 annotation and then quantified with `featureCounts`.

Scripts:

- `get_exons.sh` -- creates merged exon annotation
- `featurecouns_on_merged_exons.sh` -- generates exon-level count table on merged exons

Outputs are written to `./data/` and used in the next step.

Main input files:

- `genes.exons_merged.bed`
- `RNAseq_GV_pooled_seqs_2WT_2KO.merged_exons_counts.bed.tsv`

## DEE analysis

Differential exon expression analysis is performed in:

- `get_DEE_tables.R`

The script:

- reads merged exon annotation and exon count tables from `./data/`
- sources helper functions from `../DEG/bias_correcting_support_functions.R`
- runs exon-level differential analysis with `edgeR`
- writes output tables to `./data/`

## Outputs

The DEE script produces:

- `DEE_edgeR.GV_pooled.WT_cKO.FC2_fdr05.tsv` -- significant differential exon expression table
- `cpm.GV_pooled.WT_cKO.merged_exons.tsv` -- CPM table for merged exons
- `RNAseq_GV_pooled_seqs_2WT_2KO.merged_exons_counts.DEE.tsv` -- full annotated output table with DEE results, CPM, and RPKM values

## Required packages:

With version used:
```
# R version 4.2.3 (2023-03-15)
tidyverse (tidyverse_1.3.2)
edgeR (edgeR_3.40.2)
```

