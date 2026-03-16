gtf=$1

exons=${gtf::-3}exons.bed
exons_merged=${gtf::-3}exons_merged.bed

awk 'BEGIN{FS=OFS="\t"}
$0 !~ /^#/ && $3=="exon" {
  chr=$1; start=$4-1; end=$5; strand=$7; attr=$9
  gname="."; gid="."; tname="."; exname="."; exnum="."
  if (match(attr, /gene_name "([^"]+)"/, a))                     gname=a[1]
  if (match(attr, /gene_id "([^"]+)"/,   b))                     gid=b[1]
  if (match(attr, /transcript_name "([^"]+)"/, c))               tname=c[1]
  else if (match(attr, /transcript_id "([^"]+)"/, c))            tname=c[1]
  if (match(attr, /exon_id "([^"]+)"/, d))                       exname=d[1]
  else if (match(attr, /exon_name "([^"]+)"/, d))                exname=d[1]
  if (match(attr, /exon_number[[:space:]]+([^;[:space:]]+)/, e)) exnum=e[1]
  if (exname=="." && exnum!=".")                                 exname=exnum  
  print chr, start, end, strand, gname, gid, tname, exname, exnum
}' $gtf | sort -k1,1 -k2,2n  > $exons

module load bedtools 

bedtools merge -i $exons -c 4,5,6,7,8,5,7,9 -o distinct,distinct,distinct,distinct,distinct,count_distinct,count_distinct,min > $exons_merged

awk -F'\t' -v OFS='\t' '{ $2=$2+1; x = "ME_"$1"_"$2"_"$3 ; print x, $1, $2, $3, "+" }' $exons_merged > ${exons_merged::-3}saf

