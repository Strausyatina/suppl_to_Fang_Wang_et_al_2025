#!/bin/bash


ChromHMM.sh BinarizeBam -paired chrom_sizes.tsv input_noH3K9me3 cellmarkfiletable.tsv output_noH3K9me3

let k=6
ChromHMM.sh LearnModel output_noH3K9me3 model$k'_output_noH3K9me3' $k mm10

ChromHMM.sh OverlapEnrichment model$k'_output_noH3K9me3'/GV_Oocytes_WT_$k'_segments.bed' input_bed/ model$k'_output_noH3K9me3'/state_enrichments_WT

ChromHMM.sh OverlapEnrichment model$k'_output_noH3K9me3'/GV_Oocytes_KO_$k'_segments.bed' input_bed/ model$k'_output_noH3K9me3'/state_enrichments_KO


ChromHMM.sh OverlapEnrichment model$k'_output_noH3K9me3'/GV_Oocytes_WT_$k'_segments.bed' input_eRNAs/ model$k'_output_noH3K9me3'/eRNA_state_enrichments_WT

ChromHMM.sh OverlapEnrichment model$k'_output_noH3K9me3'/GV_Oocytes_KO_$k'_segments.bed' input_eRNAs/ model$k'_output_noH3K9me3'/eRNA_state_enrichments_KO


ChromHMM.sh NeighborhoodEnrichment model$k'_output_noH3K9me3'/GV_Oocytes_WT_$k'_segments.bed' TSS0_corrected.bed model$k'_output_noH3K9me3'/TSS0_neighbourhood_WT

ChromHMM.sh NeighborhoodEnrichment model$k'_output_noH3K9me3'/GV_Oocytes_KO_$k'_segments.bed' TSS0_corrected.bed model$k'_output_noH3K9me3'/TSS0_neighbourhood_KO



echo 'done all'
