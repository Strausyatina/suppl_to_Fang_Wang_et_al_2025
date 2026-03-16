module load subread/2.0.6

bams=(`ls ../filtered_bam/*bam`)
exons=~/group3_mendelevich/references/ensembl94_M19/genes.exons_merged.saf 
odir=.
base_total="RNAseq_GV_pooled_seqs_2WT_2KO"

# note: see ~/group3_mendelevich/references/ensembl94_M19/get_exons.sh for ${exons} production.

for bam in ${bams[*]}; do
  base=$(basename $bam .bam)
  featureCounts -F SAF -a ${exons} -o ${base}.merged_exons_counts.txt -T 8 -s 0 -p -B -C -Q 10 --primary ${bam}
done

featureCounts -F SAF -a ${exons} -o ${base_total}.merged_exons_counts.txt -T 8 -s 0 -p -B -C -Q 10 --primary ${bams[*]}

module load bedtools
tail -n +3 ${base_total}.merged_exons_counts.txt | awk -F'\t' -v OFS='\t' '{$3=$3-1; print}' | cut --complement -f1,5,6 \
  | bedtools intersect -wa -wb -a - -b ${exons::-3}bed | cut --complement -f8-10  > ${base_total}.merged_exons_counts.bed.tsv
# important: if samples != 4, then it will be n=8+(N-4),m=n+3 and `cut --complement -fn-m` instead of `cut --complement -f8-10`


