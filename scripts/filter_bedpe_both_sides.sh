#!/usr/bin/bash

# Filter bedpe file by a set of bed features
# MB 20200821

BED_R1="R1.bed"
BED_R2="R2.bed"
META="metadata.txt"

awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,"loop_"NR}' $1 > $BED_R1
awk 'BEGIN{FS=OFS="\t"} {print $4,$5,$6,"loop_"NR}' $1 > $BED_R2
awk 'BEGIN{FS=OFS="\t"} {print $7,$8,$9,"loop_"NR}' $1 > $META

# Not a great idea, since it removes regions with a little bit of H3K4me3 and a lot of H3K27ac - could be active enhancers still
bedtools intersect -a $BED_R1 -b $2 > ${BED_R1/.bed/_filter.bed}

# Filter promoters by the H3K4me3 peaks
bedtools intersect -a $BED_R2 -b $2 > ${BED_R2/.bed/_filter.bed}

Rscript /data/proj/GCB_MB/CT/git_test/scCut-Tag_2020/scripts/merge_bedpe.R  ${BED_R1/.bed/_filter.bed} 4 ${BED_R2/.bed/_filter.bed} 4 $META 4 $3

rm R1*
rm R2*
rm metadata.txt