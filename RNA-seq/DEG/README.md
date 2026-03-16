# DEG analysis

This folder contains the scripts and input/output tables used for RNA-seq differential expression analysis in Fang Wang et al.

## Contents

- `get_DEG_tables.R` – main script to run DEG analysis and export result tables
- `bias_correcting_support_functions.R` – helper functions for filtering, normalization, model fitting, and plotting
- `data/` – input count tables and exported DEG / CPM output tables

## Method summary

Differential expression analysis was performed in `edgeR`.  
To reduce technical bias, normalization factors were estimated using genes with relatively low variance across samples ("housekeeping-like" genes).  
For datasets with repeated sequencing for the same replicate, counts were pooled before testing.

## Comparisons included

- GV WT vs MLL3/4 cKO
- GO p14 WT vs MLL3/4 cKO
- GV WT vs cKO_MLL3
- GV WT vs cKO_MLL4
- GV WT vs cKO_MLL34

## Outputs

The script exports:

- DEG tables filtered at `FDR < 0.05` and `|FC| >= 2`
- CPM tables for downstream visualization and comparison

## Required packages:

With version used:
```
# R version 4.2.3 (2023-03-15)
tidyverse (tidyverse_1.3.2)
edgeR (edgeR_3.40.2)
```
