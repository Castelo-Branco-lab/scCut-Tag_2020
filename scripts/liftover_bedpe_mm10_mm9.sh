#!/usr/bin/bash

# Filter bedpe file by a set of bed features
# MB 20200821

BED_R1="R1_for_liftover.bed"
BED_R2="R2_for_liftover.bed"
META="metadata.txt"
#BEDPE_OUT="Enhancer_predictions_H3K4me3_filter_mm9.bedpe"

awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,"loop_"NR}' $1 > $BED_R1
awk 'BEGIN{FS=OFS="\t"} {print $4,$5,$6 + 500 ,"loop_"NR}' $1 > $BED_R2
awk 'BEGIN{FS=OFS="\t"} {print $7,$8,$9,"loop_"NR}' $1 > $META

[[ -f 'mm10ToMm9.over.chain.gz' ]] || wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/mm10ToMm9.over.chain.gz
liftOver $BED_R1 mm10ToMm9.over.chain.gz ${BED_R1/.bed/_lifted.bed} ${BED_R1/.bed/_unlifted.bed}
liftOver $BED_R2 mm10ToMm9.over.chain.gz ${BED_R2/.bed/_lifted.bed} ${BED_R2/.bed/_unlifted.bed}

awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3 - 500 ,$4}' ${BED_R2/.bed/_lifted.bed} > ${BED_R2/.bed/_lifted_fixed.bed} 

 
Rscript /data/proj/GCB_MB/CT/git_test/scCut-Tag_2020/scripts/merge_bedpe.R  ${BED_R1/.bed/_lifted.bed} 4 ${BED_R2/.bed/_lifted_fixed.bed} 4 $META 4 $2

rm R1*
rm R2*
rm metadata.txt