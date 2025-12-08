ChromHMM version 1.27, build hdfd78af_0 from bioconda was used.

Input description:
------------------

- input_noH3K9me3 : folder containing all bam files listed in cellmarkfiletable.tsv
- cellmarkfiletable.tsv: file specyfing biological grouping (WT, KO GV oocytes) and chromatin marks (H3K4me1, H3K4me3 et al.) per input bam file
- chrom_sizes.tsv: file specyfying length in bases for each chromosome. Created using the bash cut command on samtools fasta index of GRCm38 genome.
- TSS0_corrected.bed: file specyfying the TSS positions for gencode m19 genes.
- input_bed: folder with bed files with functional genomic elements used for overlap enrichment calculation.
- input_eRNAs: folder with bed files with differentially expressed enhancer RNAs in KO and WT GV oocytes.


Parameter description:
----------------------
- k: number of chromatin states to emit with the LearnModel command. After swiping through a range for this parameter, value 6 was used i.e. 6 chromatin states were emitted.

Output description:
-------------------

- output_noH3K9me3: folder with results of the BinarizeBam command
- model$k'_output_noH3K9me3': folder with results of the LearnModel command for given value of parameter k.
- model$k'_output_noH3K9me3'/GV_Oocytes_WT_$k'_segments.bed': file with state assignments per genomic range for WT GV oocytes
- model$k'_output_noH3K9me3'/GV_Oocytes_KO_$k'_segments.bed': file with state assignments per genomic rage for KO GV oocytes
- model$k'_output_noH3K9me3'/state_enrichments_WT.txt: state enrichment values calculated for input bed files in WT oocytes with OverlapEnrichment command 
- model$k'_output_noH3K9me3'/state_enrichments_KO.txt: state enrichment values calculated for input bed files in KO oocytes with OverlapEnrichment command
- model$k'_output_noH3K9me3'/state_enrichments_WT.png: plot of state enrichment values calculated for input bed files in WT oocytes with OverlapEnrichment command
- model$k'_output_noH3K9me3'/state_enrichments_KO.png: plot of state enrichment values calculated for input bed files in KO oocytes with OverlapEnrichment command
- model$k'_output_noH3K9me3'/eRNA_state_enrichments_WT.txt: state enrichment values calculated for eRNA bed files in WT oocytes with OverlapEnrichment command
- model$k'_output_noH3K9me3'/eRNA_state_enrichments_KO.txt: state enrichment values calculated for eRNA bed files in KO oocytes with OverlapEnrichment command
- model$k'_output_noH3K9me3'/eRNA_state_enrichments_WT.png: plot of state enrichment values calculated for eRNA bed files in WT oocytes with OverlapEnrichment command
- model$k'_output_noH3K9me3'/eRNA_state_enrichments_KO.png: plot of state enrichment values calculated for eRNA bed files in KO oocytes with OverlapEnrichment command
- model$k'_output_noH3K9me3'/TSS0_neighbourhood_WT.txt: neighbourhood enrichment values calculated for TSS positions in WT oocytes with NeighborhoodEnrichment command
- model$k'_output_noH3K9me3'/TSS0_neighbourhood_KO.txt: neighbourhood enrichment values calculated for TSS positions in KO oocytes with NeighborhoodEnrichment command
- model$k'_output_noH3K9me3'/TSS0_neighbourhood_WT.png: plot of neighbourhood enrichment values calculated for TSS positions in WT oocytes with NeighborhoodEnrichment command
- model$k'_output_noH3K9me3'/TSS0_neighbourhood_KO.png: plot of neighbourhood enrichment values calculated for TSS positions in KO oocytes with NeighborhoodEnrichment command


