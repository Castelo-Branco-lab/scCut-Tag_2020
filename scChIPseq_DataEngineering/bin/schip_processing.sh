#!/bin/bash

## Mapping Research Pipeline
## Copyleft 2017 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the CECILL License
## See the LICENCE file for details


## Input : R1 fastq + R2 fastq
## Outpout : R1 fastq with a barcode flag

SOFT="schip_processing"
VERSION="0.0.3"
ARGUMENTS=$@
function usage {
    echo -e "usage : $SOFT -f FORWARD -r REVERSE -o OUTPUT -c CONFIG [-d] [-h] [-v]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo
    echo "$SOFT $VERSION"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -f|--forward R1_READ: forward fastq file"
    echo "   -r|--reverse R2_READ: forward fastq file"
    echo "   -c|--conf CONFIG: configuration file for ChIP processing"
    echo "   -o|--output OUTPUT: output folder"
    echo "   [-d|--dryrun]: dry run mode"
    echo "   [-h|--help]: help"
    echo "   [-v|--version]: version"
    exit;
}

function version {
    echo -e "$SOFT version $VERSION"
    exit
}

function set_dry_run {
    DRY_RUN=1
}

function opts_error {
    echo -e "Error : invalid parameters !" >&2
    echo -e "Use $SOFT -h for help"
    exit
}

if [ $# -lt 1 ]
then
    usage
    exit
fi

for arg in "$@"; do
  shift
  case "$arg" in
      "--reverse") set -- "$@" "-r" ;;
      "--forward") set -- "$@" "-f" ;;
      "--output") set -- "$@" "-o" ;;
      "--conf")   set -- "$@" "-c" ;;
      "--dryrun")   set -- "$@" "-d" ;;
      "--help")   set -- "$@" "-h" ;;
      "--version")   set -- "$@" "-v" ;;
      *)        set -- "$@" "$arg"
  esac
done

while getopts "f:r:o:c:dvh" OPT
do
    case $OPT in
        f) FORWARD=$OPTARG;;
        r) REVERSE=$OPTARG;;
        o) ODIR=$OPTARG;;
        c) CONF=$OPTARG;;
        d) set_dry_run ;;
        v) version ;;
        h) help ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z $FORWARD || -z $REVERSE || -z $CONF || -z $ODIR ]]; then
    usage
    exit
fi


BIN_PATH=`dirname "$0"`
BIN_NAME=`basename "$0"`
ABS_BIN_PATH=`cd "$BIN_PATH"; pwd`
SCRIPTS_PATH="$ABS_BIN_PATH/../scripts/"

. $SCRIPTS_PATH/utils.inc.sh
. $SCRIPTS_PATH/func.inc.sh

#####################
## Check Config file
#####################

if [ ! -z "$CONF" ]; then
    CONF=`abspath $CONF`
    if [ -e "$CONF" ]; then
        read_config $CONF
    else
        echo "Error - config file '$CONF' not found"
        exit
    fi
fi

CMD_LINE="$@"
LOGDIR=${ODIR}/logs
mkdir -p ${LOGDIR}
mkdir -p ${ODIR}

PREFIX=$(get_fastq_prefix ${REVERSE})

## 1- Align R2 reads on barcode indexes
barcode_index_mapping_func ${REVERSE} ${ODIR}/mapping/barcode ${PREFIX} ${LOGDIR}/barcode
BARCODE_READS=${ODIR}/mapping/barcode/${PREFIX}_read_barcodes.txt

## 2 Trim R2 reads for genome aligning	
MODE_TRIMMING='genome'
fastx_trimmer_func ${REVERSE} ${BARCODE_LINKER_LENGTH} ${ODIR}/trimming ${LOGDIR} ${MODE_TRIMMING}
REVERSE_TRIMMED_G=${ODIR}/trimming/${PREFIX}_trimmed_G.R2.fastq.gz

## 3-bis Align R2 reads on genome indexes - paired end with R1 - (STAR)
MAPPING_INDEX_STAR=${GENOME_IDX_PATH_STAR}
MAPPING_OPTS_STAR=${GENOME_MAPPING_OPTS_STAR}
PAIRED_END=(${FORWARD} ${REVERSE_TRIMMED_G})
star_func "$(echo ${PAIRED_END[@]})" ${ODIR}/mapping/genome ${LOGDIR}/genome
GENOME_BAM=${ODIR}/mapping/genome/${PREFIX}.bam

## 4- #addBC_removeDuplicate  - (STAR)
addBC_removeDuplicate_star_func ${GENOME_BAM} ${BARCODE_READS} ${ODIR}/mapping ${LOGDIR}
GENOME_BAM_FLAGGED=${ODIR}/mapping/${PREFIX}_flagged.bam
GENOME_BAM_FLAGGED_RMDUP=${ODIR}/mapping/${PREFIX}_flagged_rmPCR_RT.bam

## 4-Remove duplicates with R2 unmapped - prime (STAR)
remove_duplicates ${GENOME_BAM_FLAGGED_RMDUP} ${ODIR}/mapping/ ${LOGDIR}
GENOME_BAM_FLAGGED_RMDUP=${ODIR}/mapping/${PREFIX}_flagged_rmPCR_RT_rmDup.bam

## 5- Generate bigwig files
bw_func ${GENOME_BAM_FLAGGED_RMDUP} ${ODIR}/tracks/ ${LOGDIR} 

## 6- Use the R1 bam with the barcode flag to generate the count table (sc2counts.py)
time make_counts ${GENOME_BAM_FLAGGED_RMDUP} ${ODIR}/counts/ ${ODIR}/mapping/ ${LOGDIR}

## 7- Write Metadata
add_info_to_log ${ODIR}/mapping/genome ${ODIR}/logs ${ODIR} ${PREFIX} ${BIN_PATH} ${ARGUMENTS} ${BIN_NAME}
#write_metadata ${ODIR} ${PREFIX} "${CMD_LINE}" ${CONF} ${LOGDIR} ${GENOME_BAM_FLAGGED} ${GENOME_BAM_FLAGGED_RMDUP} ${BARCODE_BAM}

## 8- Run Downstream analysis with default parameters
DATASET_NAME=${PREFIX}_${BIN_SIZE}
COUNT_MAT=${PREFIX}_flagged_rmPCR_RT_rmDup_counts_${BIN_SIZE}.tsv

export PATH=/bioinfo/local/build/Centos/bedtools/bedtools-2.25.0/bin/:$PATH:/bioinfo/local/build/MACS2_2.0.10/bin/

Rscript ${R_DOWNSTREAM} "$(dirname ${R_DOWNSTREAM})" ${ODIR}/ ${DATASET_NAME} ${ANNOT} -1 ${ODIR}/counts/${COUNT_MAT} -n ${N_CLUSTER} -p ${MIN_PERCENT_COR} -log ${LOGDIR}/${DATASET_NAME}_R_downstream.log
/bioinfo/local/build/Centos/python/python-3.6.1/bin/multiqc -f --no-data-dir -i ${PREFIX} -o ${ODIR} -n ${DATASET_NAME}_report.html  -c ~/scChIPseq/config_defaults.yaml -f ${ODIR}/mapping/${PREFIX}_flagged.count ${ODIR}/mapping/${PREFIX}_flagged_rmPCR_RT_rmDup.count ${ODIR}/mapping/${PREFIX}_flagged_rmPCR.count ${ODIR}/mapping/${PREFIX}_flagged_rmPCR_RT.count ${ODIR}/${PREFIX}_scChIPseq_logs.txt ${ODIR}/counts/${COUNT_MAT}

for bed in ${BED_FEATURES}
do
	bed=$(basename $bed | sed 's/.bed//')
	DATASET_NAME=${PREFIX}_${bed}
	COUNT_MAT=${PREFIX}_flagged_rmPCR_RT_rmDup_counts_${bed}.tsv
	Rscript ${R_DOWNSTREAM} `dirname ${R_DOWNSTREAM}` ${ODIR}/ ${DATASET_NAME} ${ANNOT} -1 ${ODIR}/counts/${COUNT_MAT} -n ${N_CLUSTER} -p ${MIN_PERCENT_COR} -log ${LOGDIR}/${DATASET_NAME}_R_downstream.log
/bioinfo/local/build/Centos/python/python-3.6.1/bin/multiqc -f --no-data-dir -i ${PREFIX} -o ${ODIR} -n ${DATASET_NAME}_report.html  -c ~/scChIPseq/config_defaults.yaml -f ${ODIR}/mapping/${PREFIX}_flagged.count ${ODIR}/mapping/${PREFIX}_flagged_rmPCR_RT_rmDup.count ${ODIR}/mapping/${PREFIX}_flagged_rmPCR.count ${ODIR}/mapping/${PREFIX}_flagged_rmPCR_RT.count ${ODIR}/${PREFIX}_scChIPseq_logs.txt ${ODIR}/counts/${COUNT_MAT}

done

echo
echo -e "Completed on ${when}! Results are available in ${ODIR}"
echo

