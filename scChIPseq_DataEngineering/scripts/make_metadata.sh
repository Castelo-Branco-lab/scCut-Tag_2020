#! /bin/bash

## Pacome Prompsy

## Create metadata info
when=$(date +[%Y-%m-%d]@[%H:%M:%S])
curdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
version=$(${curdir}/../bin/schip_processing.sh -v)

odir=$1
prefix=$2
cmd_line=$3
conf=$4
logdir=$5
bam_genome=$6
bam_genome_rmdup=$7
bam_barcode=$8

mkdir -p ${odir}
meta=$odir/metadata.txt
echo "######## SAMPLE DESCRIPTION ########

SAMPLE NAME = ${prefix}
RUN DATE = ${when}
VERSION = ${version} " > ${meta}

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
#echo "RUN COMMAND= ${DIR}/${logdir}
#-f ${bam_genome}
#-r ${bam_genome_rmdup}
#-c ${bam_barcode}
#-o ${odir} 
#"	 >> ${meta}

echo -e "RUN COMMAND= $cmd_line \n" >> ${meta}


echo "######## RESULTS METADATA ########
" >> ${meta}
initialReads=$(samtools view -c ${bam_barcode})

echo "NUMBER OF INITITAL READS = $initialReads" >> ${meta}

echo "
##MAPPING R2 TO BARCODE" >> ${meta}
cat ${logdir}/barcode/mapping_bowtie.log |  grep -e '^#' >> ${meta} 


echo "
##MAPPING R1 TO GENOME (paired_end) " >> ${meta}
cat ${logdir}/genome/mapping_bowtie.log |  grep -e '^#' >> ${meta} 

echo "## NUMBER OF READS MAPPING BOTH BARCODE AND GENOME =
"  >> ${meta}



reads_both=$(samtools view -c ${bam_genome})

echo "NUMBER OF READS MAPPING BOTH BARCODE AND GENOME = ${reads_both} ( ~ $(( 100*$reads_both / $initialReads )) % of intitial reads ) " >> ${meta}

echo "
## BARCODE TAGGING AND DUPLICATE REMOVAL

#BEFORE RM-DUP
" >> ${meta} 
echo "#READS BEFORE RM-DUP =  $(cat ${logdir}/addBC_removeDuplicate.log |  grep -e '^## Number of reads:' | grep -o -e '[0-9]*')" >> ${meta}

echo "
#AFTER RM-DUP
" >> ${meta} 

readsBeforeRmDup=$(cat ${logdir}/addBC_removeDuplicate.log |  grep -e '^## Number of reads: ' | grep -o -e '[0-9]*')
readsAfterRmDup=$(cat ${logdir}/addBC_removeDuplicate.log |  grep -e '^## Number of reads after duplicates removal:' | grep -o -e '[0-9]*')

echo "$(cat ${logdir}/addBC_removeDuplicate.log |  grep -e '^## Number of duplicates: ')  ( ~ $(( 100*($readsBeforeRmDup - $readsAfterRmDup)/$readsBeforeRmDup )) % of duplicates )" >> ${meta}
echo "$(cat ${logdir}/addBC_removeDuplicate.log |  grep -e '^## Number of reads after duplicates removal:')  ( ~ $(( 100*$readsAfterRmDup / $readsBeforeRmDup )) % of unique reads after rmDup )" >> ${meta}

##Create Histograms
samtools view -F 4 ${bam_genome} | cut -f 15  > ${odir}/${prefix}_flagged.ReadName_Barcode
samtools view -F 4 ${bam_genome_rmdup} | cut -f 15  > ${odir}/${prefix}_flagged_rmDup.ReadName_Barcode

##Create R report
Rscript --vanilla ${curdir}/read_stats.R ${odir}/${prefix}_flagged.ReadName_Barcode ${odir} ${prefix}_flagged $MIN_COUNT_PER_BARCODE_BEFORE_RMDUP > ${logdir}/read_stats.Rout 2>&1
Rscript --vanilla ${curdir}/read_stats.R ${odir}/${prefix}_flagged_rmDup.ReadName_Barcode ${odir} ${prefix}_flagged_rmDup $MIN_COUNT_PER_BARCODE_AFTER_RMDUP > ${logdir}/read_stats_rmdup.Rout 2>&1
Rscript --vanilla ${curdir}/plot_before_after.R ${odir} ${prefix} > ${logdir}/plot_before_after.Rout 2>&1


n500=$(awk -v limit=500 '$2>=limit && NR>1{c++} END{print c+0}' ${odir}/${prefix}_flagged_rmDup.bc_file)
n1000=$(awk -v limit=1000 '$2>=limit && NR>1{c++} END{print c+0}' ${odir}/${prefix}_flagged_rmDup.bc_file)
n1500=$(awk -v limit=1500 '$2>=limit && NR>1{c++} END{print c+0}' ${odir}/${prefix}_flagged_rmDup.bc_file)

echo "Number of barcodes with more than 500 unique reads = $n500
Number of barcodes with more than 1000 unique reads = $n1000
Number of barcodes with more than 1500 unique reads = $n1500
" >> ${meta}


echo "
######## CONFIG PARAMETERS ########
	" >> ${meta}	

cat ${conf} | awk '/PARAMETERS/{y=1;next}y' >> ${meta} 
rm ${odir}/*.ReadName_Barcode


