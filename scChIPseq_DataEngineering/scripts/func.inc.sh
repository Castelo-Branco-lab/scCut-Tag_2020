## Mapping Research Pipeline
## Copyleft 2018 Institut Curie
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the CECILL License
## See the LICENCE file for details


## Trim reads
## -m <mode>
## $1 = fastq
## $2 = length from 5' to keep
## $3 = out dire
## $4 = log dir
fastx_trimmer_func()
{
  	       	       		  ##Logs                                                                                                                                                                                                                               
    local log=$4/trim.log
    echo -e "Running Fastx trimmer ..."
    echo -e "Logs: $log"
    echo
    
    ##Mode
    mode=$5

    ## -Q 33
    ## Centos version ?
    local out=$3
    mkdir -p ${out}
    
	
	if [ $mode == 'barcode' ] 
		then
		local ofile=${out}/$(get_fastq_prefix $1)_trimmed_BC.R2.fastq
    		cmd="${FASTX_PATH}fastx_trimmer -Q 33 -l $2 -i <(gzip -cd $1) -o ${ofile}"
    		exec_cmd ${cmd} > ${log} 2>&1	
		else			#Change option -l (last) to -f (first) to remove the barcode and linker from R2 reads in preparation for genome alignment
		local ofile=${out}/$(get_fastq_prefix $1)_trimmed_G.R2.fastq
		cmd="${FASTX_PATH}fastx_trimmer -Q 33 -f $2 -i <(gzip -cd $1) -o ${ofile}" 
    	exec_cmd ${cmd} > ${log} 2>&1	
	fi

    cmd="gzip ${ofile}"
    exec_cmd ${cmd} > ${log} 2>&1
}


## bowtie
bowtie_func()
{
   ## Logs
   mkdir -p $3
   local log=$3/mapping_bowtie.log
   echo -e "Running Bowtie mapping ..."
   echo -e "Logs: $log"
   echo

   ## input type 
   inputs=($1)
   if [[ ${#inputs[@]} -eq 1 ]]; then
       if [[ ${inputs[0]} =~ \.gz ]]; then
           cmd_in=" <(gzip -cd ${inputs[0]})"
       else
           cmd_in="${inputs[0]}"
       fi
   elif [[ ${#inputs[@]} -eq 2 ]]; then
       if [[ ${inputs[0]} =~ \.gz ]]; then
           cmd_in="-1 <(gzip -cd ${inputs[0]}) -2 <(gzip -cd ${inputs[1]})"
       else
           cmd_in="-1 ${inputs[0]} -2 ${inputs[1]}"
       fi
   fi

   ## sample_id
   if [ ! -z ${SAMPLE_ID} ]; then
       cmd_id="--sam-RG ID:${SAMPLE_ID} --sam-RG SM:${SAMPLE_ID} --sam-RG LB:${SAMPLE_ID} --sam-RG PU:${SAMPLE_ID} --sam-RG PL:ILLUMINA"
   else
       cmd_id=""
   fi

   ## Run
   local out=$2
   mkdir -p ${out}
   local out_prefix=${out}/$(get_fastq_prefix ${inputs[0]})

   cmd="bowtie -p ${NB_PROC} ${MAPPING_OPTS} $cmd_id --sam ${MAPPING_INDEX} ${cmd_in} > ${out_prefix}.sam"
   exec_cmd ${cmd} > ${log} 2>&1

   cmd="samtools view -@ ${NB_PROC} -bS ${out_prefix}.sam > ${out_prefix}.bam"
   exec_cmd ${cmd} >> ${log} 2>&1

   cmd="samtools sort -n -@ ${NB_PROC} ${out_prefix}.bam -o ${out_prefix}_nsorted.bam && mv ${out_prefix}_nsorted.bam ${out_prefix}.bam"
   exec_cmd ${cmd} >> ${log} 2>&1

   ## Clean sam file
   cmd="rm ${out_prefix}.sam"
   exec_cmd ${cmd} >> ${log} 2>&1
}

## star
star_func()
{
   ## Logs
   mkdir -p $3
   local log=$3/mapping_star.log
   echo -e "Running STAR mapping ..."
   echo -e "Logs: $log"
   echo

   ## input type 
   inputs=($1)
   if [[ ${#inputs[@]} -eq 1 ]]; then
       if [[ ${inputs[0]} =~ \.gz ]]; then
           cmd_in=" <(gzip -cd ${inputs[0]})"
       else
           cmd_in="${inputs[0]}"
       fi
   elif [[ ${#inputs[@]} -eq 2 ]]; then
       if [[ ${inputs[0]} =~ \.gz ]]; then
           cmd_in="<(gzip -cd ${inputs[0]}) <(gzip -cd ${inputs[1]})"
       else
           cmd_in="${inputs[0]} ${inputs[1]}"
       fi
   fi

   ## sample_id
   if [ ! -z ${SAMPLE_ID} ]; then
       cmd_id=""
   else
       cmd_id=""
   fi

   ## Run
   local out=$2
   mkdir -p ${out}
   local out_prefix=${out}/$(get_fastq_prefix ${inputs[0]})

   cmd="STAR --runMode alignReads --runThreadN ${NB_PROC} ${MAPPING_OPTS_STAR} --readFilesIn ${cmd_in} --genomeDir ${MAPPING_INDEX_STAR} --outFileNamePrefix ${out}/"
   exec_cmd ${cmd} > ${log} 2>&1


   cmd="samtools view -@ ${NB_PROC} -bS ${out}/Aligned.out.sam > ${out_prefix}.bam"
   exec_cmd ${cmd} >> ${log} 2>&1
	
   cmd="samtools sort -n -@ ${NB_PROC} ${out_prefix}.bam -o ${out_prefix}_nsorted.bam && mv ${out_prefix}_nsorted.bam ${out_prefix}.bam"
   exec_cmd ${cmd} >> ${log} 2>&1

   ## Clean sam file
   cmd="rm ${out}/Aligned.out.sam"
   exec_cmd ${cmd} >> ${log} 2>&1
}

barcode_index_mapping_func()
{
##Function mappin barcodes from single-cell ChIP-seq data, index by index (3 different indexes) to the reference (forward) indexes 
#Input : Read 2 (.fastq.gz) , 		out,  		${PREFIX} ,log_dir
#Input :     $1             ,              $2      ,       $3        , 	$4 
read2=$1
out=$2
prefix=$3
## logs
    mkdir -p $4 
    mkdir -p $out
    local log=$4/index_mapping_bowtie2.log
    echo -e "Running Bowtie2 index (barcode) mapping ..."
    echo -e "Logs: $log"
    echo

if [[ ${BARCODE_LENGTH} -eq 56 ]]
	then
	echo -e "Short barcode :LBC"
        start_index_1=1
	start_index_2=21
	start_index_3=41
	size_index=16
       else
	echo -e "Long barcode : Hifibio "
        start_index_1=1
	start_index_2=25
	start_index_3=49
	size_index=20
fi

#Create indexes from reads : 1 - 16 = index 1 ; 21 - 36 = index 2; 41 - 56 = index 3
cmd="gzip -cd  $read2 | awk -v start_index_1=$start_index_1 -v size_index=$size_index  'NR%4==1{print \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_1,size_index)}' > ${out}/read_indexes_1.fasta"
exec_cmd ${cmd} > ${log} 2>&1

cmd="gzip -cd  $read2 | awk -v start_index_2=$start_index_2 -v size_index=$size_index 'NR%4==1{print  \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_2,size_index)}' > ${out}/read_indexes_2.fasta"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="gzip -cd  $read2 | awk -v start_index_3=$start_index_3 -v size_index=$size_index 'NR%4==1{print \">\"substr(\$0,2)}; NR%4==2{print substr(\$0,start_index_3,size_index)}' > ${out}/read_indexes_3.fasta"
exec_cmd ${cmd} >> ${log} 2>&1

## MATCH INDEXES 1 : allowing 
cmd="bowtie2 -x ${BARCODE_BOWTIE_IDX_PATH}1 -f ${out}/read_indexes_1.fasta ${BARCODE_MAPPING_OPTS} -p ${NB_PROC} > ${out}/index_1_bowtie2.sam"
exec_cmd ${cmd} >> ${log} 2>&1

#Keep only reads that were matched by a unique index 1: 
cmd="awk -v out=\"${out}/\" '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > out\"/count_index_1\"}' ${out}/index_1_bowtie2.sam > ${out}/reads_matching_index_1.txt"
exec_cmd ${cmd} >> ${log} 2>&1

## MATCH INDEXES 2 : 
cmd="bowtie2 -x ${BARCODE_BOWTIE_IDX_PATH}2 -f ${out}/read_indexes_2.fasta ${BARCODE_MAPPING_OPTS} -p ${NB_PROC} > ${out}/index_2_bowtie2.sam"
exec_cmd ${cmd} >> ${log} 2>&1

#Keep only reads that were matched by a unique index 2: 
cmd="awk -v out=\"${out}/\" '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > out\"count_index_2\"}' ${out}/index_2_bowtie2.sam > ${out}/reads_matching_index_2.txt"
exec_cmd ${cmd} >> ${log} 2>&1

## MATCH INDEXES 3 : 
cmd="bowtie2 -x ${BARCODE_BOWTIE_IDX_PATH}3 -f ${out}/read_indexes_3.fasta ${BARCODE_MAPPING_OPTS} -p ${NB_PROC} > ${out}/index_3_bowtie2.sam"
exec_cmd ${cmd} >> ${log} 2>&1

#Keep only reads that were matched by a unique index 3: 
cmd="awk -v out=\"${out}/\" '/XS/{next} \$2!=4{print \$1,\$3;count++} ;END{print count > out\"count_index_3\"}' ${out}/index_3_bowtie2.sam > ${out}/reads_matching_index_3.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="sort -T ${TMP_DIR} --parallel=${NB_PROC} -k1,1 ${out}/reads_matching_index_1.txt > ${out}/reads_matching_index_1_sorted.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="rm ${out}/reads_matching_index_1.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="sort -T ${TMP_DIR} --parallel=${NB_PROC} -k1,1 ${out}/reads_matching_index_2.txt > ${out}/reads_matching_index_2_sorted.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="rm ${out}/reads_matching_index_2.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="sort -T ${TMP_DIR} --parallel=${NB_PROC} -k1,1 ${out}/reads_matching_index_3.txt > ${out}/reads_matching_index_3_sorted.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="rm ${out}/reads_matching_index_3.txt"
exec_cmd ${cmd} >> ${log} 2>&1

#Join the 3 together to recompose the barcode
cmd="join -t$' ' -1 1 -2 1 ${out}/reads_matching_index_1_sorted.txt ${out}/reads_matching_index_2_sorted.txt > ${out}/tmp"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="echo \$(wc -l ${out}/tmp) | cut -d' ' -f1 > ${out}/count_index_1_2"
exec_cmd ${cmd} >> ${log} 2>&1	

cmd="join -t$' ' -1 1 -2 1 ${out}/tmp ${out}/reads_matching_index_3_sorted.txt > ${out}/final"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="awk -v out=\"${out}/\" '{print substr(\$1,1)\"\tBC\"substr(\$2,2)substr(\$3,2)substr(\$4,2);count++} ;END{print count > out\"count_index_1_2_3\"}' ${out}/final > ${out}/${prefix}_read_barcodes.txt"
exec_cmd ${cmd} >> ${log} 2>&1

#Write logs
cmd="n_index_1=\$(cat ${out}/count_index_1)"
exec_cmd ${cmd} >> ${log} 2>&1
cmd="n_index_2=\$(cat ${out}/count_index_2)"
exec_cmd ${cmd} >> ${log} 2>&1
cmd="n_index_3=\$(cat ${out}/count_index_3)"
exec_cmd ${cmd} >> ${log} 2>&1
cmd="n_index_1_2=\$(cat ${out}/count_index_1_2)"
exec_cmd ${cmd} >> ${log} 2>&1
cmd="n_index_1_2_3=\$(cat ${out}/count_index_1_2_3)"
exec_cmd ${cmd} >> ${log} 2>&1

echo "## Number of matched indexes 1: $n_index_1" >> ${log}
echo "## Number of matched indexes 2: $n_index_2" >> ${log}
echo "## Number of matched indexes 1 and 2: $n_index_1_2" >> ${log}
echo "## Number of matched indexes 3: $n_index_3" >> ${log}
echo "## Number of matched barcodes: $n_index_1_2_3" >> ${log}

## Remove all non used files
cmd="rm -f ${out}/*.sam ${out}/read* ${out}/tmp ${out}/final ${out}/count"
exec_cmd ${cmd} >> ${log} 2>&1
}

## bowtie2
bowtie2_func()
{

    ## logs
    local log=$3/mapping_bowtie2.log
    echo -e "Running Bowtie2 mapping ..."
    echo -e "Logs: $log"
    echo

    ## input type
    inputs=($1)
    if [[ ${#inputs[@]} -eq 1 ]]; then
        cmd_in="-U $1"
    elif [[ ${#inputs[@]} -eq 2 ]]; then
        cmd_in="-1 ${inputs[0]} -2 ${inputs[1]}"
    fi

    ## sample_id
    if [ ! -z ${SAMPLE_ID} ]; then
        cmd_id="--rg-id ${SAMPLE_ID} --rg SM:${SAMPLE_ID} --rg PL:ILLUMINA --rg PU:${SAMPLE_ID} --rg LB:${SAMPLE_ID}"
    else
        cmd_id=""
    fi

    ## Run
    local out=$2
    mkdir -p ${out}
    local out_prefix=${out}/$(get_fastq_prefix ${inputs[0]})

    cmd="bowtie2 -t -p ${NB_PROC} -x ${MAPPING_INDEX} ${MAPPING_OPTS} ${cmd_id} ${cmd_in} > ${out_prefix}.bam"
    exec_cmd ${cmd} > ${log} 2>&1

    cmd="samtools sort -n -@ ${NB_PROC} ${out_prefix}.bam -o ${out_prefix}_sorted.bam && mv ${out_prefix}_sorted.bam ${out_prefix}.bam"
    exec_cmd ${cmd} >> ${log} 2>&1
}

## $1 = input files	
## $2 = output dir
## $3 = log dir
bw_func()
{

    source /bioinfo/local/build/Centos/miniconda/miniconda3-4.4.6/bin/activate /bioinfo/local/build/Centos/envs_conda/deeptools_3.1.0

    local out=$2
    mkdir -p ${out}
    local log=$3
    mkdir -p ${log}
    log=$log/bamCoverage.log

    echo -e "Generate bigwig file(s) ..."
    echo -e "Logs: $log"
    echo
 	
    name=$(basename $1 ".bam")
    if [[ ! -z ${ENCODE_BLACKLIST} && -e ${ENCODE_BLACKLIST} ]]; then
        local cmd="bamCoverage --bam $1 --outFileName $out/${name}.bw --numberOfProcessors ${NB_PROC} --normalizeUsing RPKM --blackListFileName ${ENCODE_BLACKLIST}"
    else
        local cmd="bamCoverage --bam $1 --outFileName $out/${name}.bw --numberOfProcessors ${NB_PROC} --normalizeUsing RPKM "
    fi
    exec_cmd ${cmd} > ${log} 2>&1

    source /bioinfo/local/build/Centos/miniconda/miniconda3-4.4.6/bin/deactivate
}



## addBC_removeDuplicate_star
#addBC_removeDuplicate_star ${GENOME_BAM} ${BARCODE_BAM} ${ODIR}/mapping ${LOGDIR}

addBC_removeDuplicate_star_func () {
## logs
    local log=$4/addBC_removeDuplicate.log 
    echo -e "Running addBC_removeDuplicate_star ..."
    echo -e "Logs: $log"
    echo
    local out=$3
    local out_prefix=${out}/$(basename $1 | sed -e 's/.bam$//')

#Tag multimapped as unmapped 
cmd="samtools view -F 256 $1 | awk -v OFS='\t' 'NR%2==1{if(\$12==\"NH:i:1\"){mapped=1;print \$0} else{mapped=0;\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0}} NR%2==0{if(mapped==1){print \$0} else{\$2=4;\$3=\"*\";\$4=0;\$6=\"*\";print \$0} }' > ${out_prefix}.sam"
exec_cmd ${cmd} > ${log} 2>&1

#If read in mapped only on R1, set R2 position as '2147483647' (identify single end reads)
cmd="cat ${out_prefix}.sam | awk -v OFS='\t' 'NR%2==1{print \$0} NR%2==0{if(\$3==\"*\"){\$4=2147483647;print \$0} else{print \$0} }' > ${out_prefix}_2.sam"

exec_cmd ${cmd} >> ${log} 2>&1

#Remove comments from the header that produce bugs in the count phase
cmd="samtools view -H $1 | sed '/^@CO/ d' > ${out_prefix}_header.sam"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="cat ${out_prefix}_2.sam >> ${out_prefix}_header.sam && mv ${out_prefix}_header.sam ${out_prefix}.sam && samtools view -b -@ ${NB_PROC} ${out_prefix}.sam > ${out_prefix}_unique.bam"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="rm -f ${out_prefix}_2.sam ${out_prefix}.sam"
exec_cmd ${cmd} >> ${log} 2>&1

#Keeping R1 aligned + R2 start as tag 'XS'
cmd="samtools view ${out_prefix}_unique.bam | awk '{OFS = \"\t\" ; if(NR%2==1 && !(\$3==\"*\")) {R1=\$0} else if(NR%2==1){R1=0}; if(NR%2==0 && !(R1==0)){tagR2Seq=\"XD:Z:\"\$10; tagR2Pos=\"XS:i:\"\$4;print R1,tagR2Pos,tagR2Seq}}' > ${out_prefix}_unique.sam"
exec_cmd ${cmd} >> ${log} 2>&1


#Sort and join on read names reads barcoded and reads mapped to genome (barcode as tag 'XB') --> filter out unbarcoded OR unmapped reads
cmd="sort -T ${TMP_DIR} --parallel=${NB_PROC} -k1,1 ${out_prefix}_unique.sam > ${out_prefix}_unique_sorted.sam"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="join -1 1  -2 1  ${out_prefix}_unique_sorted.sam <(awk -v OFS=\"\t\" '{print \$1,\"XB:Z:\"\$2}' $2) > ${out_prefix}_merged_filtered.sam"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="sed -i 's/ /\t/g' ${out_prefix}_merged_filtered.sam"
exec_cmd ${cmd} >> ${log} 2>&1

#Remove comments from the header that produce bugs in the count phase
cmd="samtools view -H ${out_prefix}_unique.bam | sed '/^@CO/ d' > ${out_prefix}_header.sam"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="cat ${out_prefix}_merged_filtered.sam >> ${out_prefix}_header.sam && mv ${out_prefix}_header.sam ${out_prefix}_merged_filtered.sam && samtools view -@ ${NB_PROC} -b ${out_prefix}_merged_filtered.sam > ${out_prefix}_merged_filtered.bam"
exec_cmd ${cmd} >> ${log} 2>&1

##Sort by barcode then chromosome then position R2
#Find the column containing the barcode tag XB
cmd="barcode_field=\$(samtools view ${out_prefix}_merged_filtered.bam  | sed -n \"1 s/XB.*//p\" | sed 's/[^\t]//g' | wc -c)"
exec_cmd ${cmd} >> ${log} 2>&1

echo "The barcode field is $barcode_field" >> ${log}

#Find the column containing the position R2 tag XS
cmd="posR2_field=\$(samtools view ${out_prefix}_merged_filtered.bam  | sed -n \"1 s/XS.*//p\" | sed 's/[^\t]//g' | wc -c)"
exec_cmd ${cmd} >> ${log} 2>&1

echo "The position R2 field is $posR2_field" >> ${log}

cmd="printf '@HD\tVN:1.4\tSO:unsorted\n' > ${out_prefix}_header.sam"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="samtools view -H ${out_prefix}_merged_filtered.bam | sed '/^@HD/ d' >> ${out_prefix}_header.sam"
exec_cmd ${cmd} >> ${log} 2>&1

#Sort by barcode then chromosome then position R2 then Position R1 (for the PCR removal) 
#It is important to sort by R1 pos also	because the removal is done by comparing consecutive lines ! 
cmd="samtools view ${out_prefix}_merged_filtered.bam | LC_ALL=C sort -T ${TMP_DIR} --parallel=${NB_PROC} -t $'\t' -k \"$barcode_field.8,$barcode_field\"n -k 3.4,3g -k \"\$posR2_field.6,\$posR2_field\"n -k 4,4n >> ${out_prefix}_header.sam && samtools view -@ ${NB_PROC} -b ${out_prefix}_header.sam > ${out_prefix}_flagged.bam"
exec_cmd ${cmd} >> ${log} 2>&1

#Create count Table from flagged (already sorted by barcode)
cmd="samtools view ${out_prefix}_flagged.bam | awk -v bc_field=$barcode_field '{print substr(\$bc_field,6)}' |  uniq -c > ${out_prefix}_flagged.count"
exec_cmd ${cmd} >> ${log} 2>&1

if [ ${UNBOUND} == 'FALSE' ] 
then
echo -e "Is NOT an unbound, removing PCR duplicates by the position"
#Remove PCR duplicates = read having same barcode, R1 position, same R2 position, same chr ("exactly equal")
cmd="samtools view ${out_prefix}_flagged.bam | awk -v bc_field=$barcode_field -v R2_field=$posR2_field -v out=${out}/ 'BEGIN{countR1unmappedR2=0;countPCR=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\"); lastR1Pos=\$4} ; NR>=2{split( \$R2_field,R2Pos,\":\");R1Pos=\$4; if(R2Pos[3]==2147483647){print \$0;countR1unmappedR2++; next}; if( (R1Pos==lastR1Pos) && (R2Pos[3]==lastR2Pos[3]) && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){countPCR++;next} {print \$0;lastR1Pos=\$4;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\") }} END {print countPCR > out\"count_PCR_duplicates\";print countR1unmappedR2 > out\"countR1unmappedR2\"}' > ${out_prefix}_flagged_rmPCR.sam"
exec_cmd ${cmd} >> $log 2>&1
else
echo -e "Is an unbound, not removing PCR duplicates by the sequence"
cmd="seqR2_field=\$(samtools view ${out_prefix}_merged_filtered.bam  | sed -n \"1 s/XD.*//p\" | sed 's/[^\t]//g' | wc -c)"
exec_cmd ${cmd} >> ${log} 2>&1


cmd="samtools view ${out_prefix}_flagged.bam | awk -v bc_field=$barcode_field -v R2_field=$seqR2_field -v out=${out}/ 'BEGIN{countR1unmappedR2=0;countPCR=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Seq,\":\"); lastR1Pos=\$4} ; NR>=2{split( \$R2_field,R2Seq,\":\");R1Pos=\$4; if( (R1Pos==lastR1Pos) && (R2Seq[3]==lastR2Seq[3]) && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){countPCR++;next} {print \$0;lastR1Pos=\$4;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Seq,\":\") }} END {print countPCR > out\"count_PCR_duplicates\";print countR1unmappedR2 > out\"countR1unmappedR2\"}' > ${out_prefix}_flagged_rmPCR.sam"
exec_cmd ${cmd} >> $log 2>&1
fi

cmd="samtools view -H ${out_prefix}_flagged.bam  | sed '/^@CO/ d' > ${out_prefix}_header.sam"
exec_cmd ${cmd} >> ${log} 2>&1
cmd="cat ${out_prefix}_flagged_rmPCR.sam >> ${out_prefix}_header.sam && samtools view -@ ${NB_PROC} -b ${out_prefix}_header.sam > ${out_prefix}_flagged_rmPCR.bam"
exec_cmd ${cmd} >> ${log} 2>&1

#Create count Table from flagged - PCR dups (already sorted by barcode)
cmd="samtools view ${out_prefix}_flagged_rmPCR.bam | awk -v bc_field=$barcode_field '{print substr(\$bc_field,6)}' |  uniq -c > ${out_prefix}_flagged_rmPCR.count"
exec_cmd ${cmd} >> $log 2>&1

## Sort flagged_rmPCR file
cmd="samtools sort -@ ${NB_PROC} ${out_prefix}_flagged_rmPCR.bam > ${out_prefix}_flagged_rmPCR_sorted.bam"
exec_cmd ${cmd} >> $log 2>&1

## Move flagged_rmPCR file
cmd="mv ${out_prefix}_flagged_rmPCR_sorted.bam ${out_prefix}_flagged_rmPCR.bam"
exec_cmd ${cmd} >> $log 2>&1

## Index flagged_rmPCR file
cmd="samtools index ${out_prefix}_flagged_rmPCR.bam"
exec_cmd ${cmd} >> $log 2>&1

if [ ${UNBOUND} == 'FALSE' ] 
then
echo -e "Not an unbound, removing RT duplicates"
#Remove RT duplicates (if two consecutive reads have the same barcode and same R2 chr&start) but not same R1 
cmd="cat ${out_prefix}_flagged_rmPCR.sam | awk -v bc_field=$barcode_field -v R2_field=$posR2_field -v out=${out}/ 'BEGIN{count=0};NR==1{print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\")} ; NR>=2{split( \$R2_field,R2Pos,\":\");if((R2Pos[3]==lastR2Pos[3]) && (R2Pos[3]!=2147483647) && (lastR2Pos[3]!=2147483647)  && ( \$3==lastChrom ) && (\$bc_field==lastBarcode) ){count++;next} {print \$0;lastChrom=\$3;lastBarcode=\$bc_field; split( \$R2_field,lastR2Pos,\":\") }} END {print count > out\"count_RT_duplicates\"}' > ${out_prefix}_flagged_rmPCR_RT.sam"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="samtools view -H ${out_prefix}_flagged.bam  | sed '/^@CO/ d' > ${out_prefix}_header.sam"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="cat ${out_prefix}_flagged_rmPCR_RT.sam >> ${out_prefix}_header.sam && samtools view -@ ${NB_PROC} -b ${out_prefix}_header.sam > ${out_prefix}_flagged_rmPCR_RT.bam"
exec_cmd ${cmd} >> ${log} 2>&1

#Create count Table from flagged - PCR dups - RT dups  (already sorted by barcode)
cmd="samtools view ${out_prefix}_flagged_rmPCR_RT.bam | awk -v bc_field=$barcode_field '{print substr(\$bc_field,6)}' | uniq -c > ${out_prefix}_flagged_rmPCR_RT.count"
exec_cmd ${cmd} >> ${log} 2>&1

## Sort flagged_rmPCR_RT file
cmd="samtools sort -@ ${NB_PROC} ${out_prefix}_flagged_rmPCR_RT.bam > ${out_prefix}_flagged_rmPCR_RT_sorted.bam"
exec_cmd ${cmd} >> $log 2>&1

## Move flagged_rmPCR_RT file
cmd="mv ${out_prefix}_flagged_rmPCR_RT_sorted.bam ${out_prefix}_flagged_rmPCR_RT.bam"
exec_cmd ${cmd} >> $log 2>&1


else
echo -e "Is an unbound, not removing RT duplicates"
## Copy flagged_rmPCR to flagged_rmPCR_RT
cmd="cp ${out_prefix}_flagged_rmPCR.bam ${out_prefix}_flagged_rmPCR_RT.bam"
exec_cmd ${cmd} >> $log 2>&1
cmd="cp ${out_prefix}_flagged_rmPCR.count ${out_prefix}_flagged_rmPCR_RT.count"
exec_cmd ${cmd} >> $log 2>&1
## Set RT duplicate count to 0
cmd="echo 0 > ${out}/count_RT_duplicates"
exec_cmd ${cmd} >> $log 2>&1
fi

## Index flagged_rmPCR_RT file
cmd="samtools index ${out_prefix}_flagged_rmPCR_RT.bam"
exec_cmd ${cmd} >> $log 2>&1


#Write logs
cmd="n_mapped_barcoded=\$(samtools view -c ${out_prefix}_merged_filtered.bam)"
exec_cmd ${cmd} >> $log 2>&1
cmd="n_pcr_duplicates=\$(cat ${out}/count_PCR_duplicates)"
exec_cmd ${cmd} >> $log 2>&1
cmd="n_rt_duplicates=\$(cat ${out}/count_RT_duplicates)"
exec_cmd ${cmd} >> $log 2>&1
cmd="n_R1_mapped_R2_unmapped=\$(cat ${out}/countR1unmappedR2)"
exec_cmd ${cmd} >> $log 2>&1

cmd="n_unique_except_R1_unmapped_R2=\$(($n_mapped_barcoded - $n_pcr_duplicates - $n_rt_duplicates))"
exec_cmd ${cmd} >> $log 2>&1

echo "## Number of reads mapped and barcoded: $n_mapped_barcoded" >> ${log}
echo "## Number of pcr duplicates: $n_pcr_duplicates" >> ${log}
echo "## Number of rt duplicates: $n_rt_duplicates" >> ${log}
echo "## Number of R1 mapped but R2 unmapped: $n_R1_mapped_R2_unmapped" >> ${log}
echo "## Number of reads after PCR and RT removal (not R1 unmapped R2): $n_unique_except_R1_unmapped_R2" >> ${log}
## Remove all non used files
cmd="rm -f ${out_prefix}_unique.bam ${out_prefix}_genome_nsorted.bam ${out_prefix}_barcode_nsorted.bam  ${out}/count* ${out_prefix}_merged_filtered.bam ${out}/*.sam"
exec_cmd ${cmd} >> $log 2>&1
}

## Rm dup for single cell ChIP-seq
remove_duplicates()
{
    ## logs
    local log=$3/rmDup.log
    echo -e "Removing duplicates ..."
    echo -e "Logs: $log"
    echo

    local odir=$2
    mkdir -p ${odir}
    local prefix=${odir}/$(basename $1 | sed -e 's/.bam$//')

    if [ ! -z ${DUPLICATES_WINDOW} ]; then
	cmd="${PYTHON_PATH}/python ${SCRIPTS_PATH}/rmDup.py -i ${prefix}.bam -o ${prefix}_rmDup.bam -d ${DUPLICATES_WINDOW} -v "
    else
	cmd="${PYTHON_PATH}/python ${SCRIPTS_PATH}/rmDup.py -i ${prefix}.bam -o ${prefix}_rmDup.bam -v "
    fi
    exec_cmd ${cmd} >> ${log} 2>&1
    #Create count Table from flagged - PCR dups - RT dups and window-based rmDup (need to sort by b arcode)
    cmd="barcode_field=\$(samtools view ${prefix}_rmDup.bam  | sed -n \"1 s/XB.*//p\" | sed 's/[^\t]//g' | wc -c)"
    exec_cmd ${cmd} >> $log 2>&1
    cmd="samtools view ${prefix}_rmDup.bam | awk -v bc_field=$barcode_field '{print substr(\$bc_field,6)}' | sort | uniq -c > ${prefix}_rmDup.count"		
    exec_cmd ${cmd} >> $log 2>&1
 
    ## Index BAM file
    cmd="samtools index ${prefix}_rmDup.bam"
    exec_cmd ${cmd} >> $log 2>&1
}

## Generate genomic count table
make_counts(){

    ## logs
    local log=$4/make_counts.log
    echo -e "Generating counts table ..."
    echo -e "Logs: $log"
    echo
	
    local odir=$2
    mkdir -p ${odir}
    local prefix=${odir}/$(basename $1 | sed -e 's/.bam$//')
    local mapping_dir=$3

    #Counting unique BCs
    local bc_prefix=$(basename $1 | sed -e 's/.bam$//')
    cmd="barcodes=\$(wc -l ${mapping_dir}/${bc_prefix}.count | awk '{print \$1}')"
    exec_cmd ${cmd} >> ${log} 2>&1

    echo "Barcodes found = $barcodes" >> ${log} 
    for bsize in ${BIN_SIZE}
    do
	echo -e "at ${bsize} resolution ..."
	opts="-b ${bsize} "
        if [ ! -z ${MIN_COUNT_PER_BARCODE_AFTER_RMDUP} ]; then
	    opts="${opts} -f ${MIN_COUNT_PER_BARCODE_AFTER_RMDUP} "
	fi
        cmd="${PYTHON_PATH}/python ${SCRIPTS_PATH}/sc2counts.py -i $1 -o ${prefix}_counts_${bsize}.tsv ${opts} -s $barcodes -v"
	exec_cmd ${cmd} >> ${log} 2>&1
    done

    for bed in ${BED_FEATURES}
    do
        echo -e "at ${bed} resolution ..."
        opts="-B ${bed} "
        if [ ! -z ${MIN_COUNT_PER_BARCODE_AFTER_RMDUP} ]; then
            opts="${opts} -f ${MIN_COUNT_PER_BARCODE_AFTER_RMDUP} "
        fi
	osuff=$(basename ${bed} | sed -e 's/.bed//')
        cmd="${PYTHON_PATH}/python ${SCRIPTS_PATH}/sc2counts.py -i $1 -o ${prefix}_counts_${osuff}.tsv ${opts} -s $barcodes -v"
        exec_cmd ${cmd} >> ${log} 2>&1
    done
    echo

}

add_info_to_log(){

    
genomedir=$1
logdir=$2
odir=$3
prefix=$4
bin_path=$5
arguments=$6
bin_name=$7
R_downstream_dir="$(dirname ${R_DOWNSTREAM})"
echo ${R_downstream_dir}

 ## logs
    local log=$2/add_info_to_log.log
    echo -e "add_info_to_log ..."
    echo -e "Logs: $log"
    echo

cmd="cat ${genomedir}/Log.final.out > ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1
cmd="echo \"Barcode mapping\" >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="grep '##' ${logdir}/barcode/index_mapping_bowtie2.log >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="echo \"Barcode Flagging And Exact Duplicate Removal\"  >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="grep '##' ${logdir}/addBC_removeDuplicate.log >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="echo \"Duplicate Removal Window (rmDup)\"  >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="grep '## Number' ${logdir}/rmDup.log >> ${odir}/${prefix}_scChIPseq_logs.txt"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="cp -r ${bin_path}/../scripts/ ${odir}/.data_engineering_run/"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="cp -r ${bin_path}/../bin/ ${odir}/.data_engineering_run/"
exec_cmd ${cmd} >> ${log} 2>&1	

cmd="cp -r ${CONF} ${odir}/.data_engineering_run/"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="mkdir -p ${odir}/.data_analysis_run/"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="cp ${R_downstream_dir}/R_scChIP_seq_analysis.R ${odir}/.data_analysis_run/"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="cp -r ${R_downstream_dir}/Modules/ ${odir}/.data_analysis_run/"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="cp ${R_downstream_dir}/report_scChIPseq_analysis.Rmd ${odir}/.data_analysis_run/"
exec_cmd ${cmd} >> ${log} 2>&1

cmd="echo '${bin_name} ${arguments}' > ${odir}/.data_engineering_run/call.sh"
exec_cmd ${cmd} >> ${log} 2>&1

}


##Write metadata

#write_metadata ${ODIR} ${PREFIX} "${CMD_LINE}" ${CONF} ${LOGDIR} ${GENOME_BAM_FLAGGED} ${GENOME_BAM_FLAGGED_RMDUP} ${BARCODE_BAM}
write_metadata(){

    ## logs
    local log=$5/metadata.log
    echo -e "Generating metadata table ..."
    echo -e "Logs: $log"
    echo
 
    cmd="bash ${SCRIPTS_PATH}/make_metadata.sh $1/metadata $2 \"$3\" $4 $5 $6 $7 $8"
    exec_cmd ${cmd} > ${log} 2>&1
}

get_command_line(){
    line=$(history 1)
    line=${line#*[0-9]  }
    echo "$line"
}
