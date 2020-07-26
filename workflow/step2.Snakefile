import os
import sys

# Define paths and sample names in config file
configfile: "config.yaml"

# Linking with a rule did not work for me
# So here is python code to linkg the files

# Create folder ./data
if not os.path.isdir("data"):
    os.mkdir("data")
    
# Link the source 10x files defined in config to ./data
for key in config["samples"].keys():
    
    data_dir = config["samples"][key]
    link_dir = " data/" + key

    if not os.path.exists(link_dir):
        shell("ln -sf {0} {1}".format(data_dir, link_dir))

# print(config['clustering_params'])
##################################

rule all:
     input:
         # bigwig files
         expand("bigwig/{sample}_all_reads.bw", sample=config["samples"]),

         # snap files
#         expand("logs/{sample}_bmat_out.log", sample=config["samples"]),

         # macs peak files
         expand("macs2_peaks/narrow/{sample}/{sample}_peaks.narrowPeak",sample=config["samples"]),
         expand("macs2_peaks/broad/{sample}/{sample}_peaks.broadPeak",sample=config["samples"]),
         "macs2_peaks.combined/narrow/combined_peaks.narrowPeak",
         "macs2_peaks.combined/broad/combined_peaks.broadPeak",
         
         # Merged bam file
         "data/combined/outs/possorted_bam.bam",

         # Remove duplicates from bam file
         expand("barcode_statistics/{sample}/possorted_bam_nodup.bam",sample=config["samples"]),

         # Barcode statistics log file
         expand("barcode_statistics/{sample}/barcode_statistics_all.log",sample=config["samples"]),
         expand("barcode_statistics/{sample}/barcode_statistics_peaks.log",sample=config["samples"]),
         
         # Modified fragments files:
         expand("data/combined/outs/fragments_{sample}.tsv",sample=config["samples"]),
         "data/combined/outs/fragments.tsv.gz",
         
         #Clustering
         expand("clustering_snakemake/01.clustering/Seurat_{binwidth}_{umi}.Rds",binwidth = config['binwidth'] + ["peaks"],umi = ['UMI']),
           
         # bigwig
         expand("clustering_snakemake/01.clustering/bigwig/{binwidth}_{umi}/{clusterMethod}_all_clusters.bw",binwidth=config['binwidth'] + ["peaks"], umi=['UMI'],clusterMethod=['seurat']),

         # MEME motif search
#         expand("MEME/out_{npeaks}/meme-chip.html",npeaks = config['MEME_npeaks'])
 
############ Clustering
rule clustering_final:
  input:
      peaks           = "macs2_peaks.combined/broad/combined_peaks.broadPeak",
      fragments_UMI   = "data/combined/outs/fragments.tsv.gz",
#      fragments_noUMI = "data/combined/outs/noUMI/fragments.tsv.gz",
  output:
      "clustering_snakemake/01.clustering/Seurat_{binwidth}_{umi}.Rds",
  params:
      min_reads      =  "--reads_min " + str(config["clustering_params"]['min_reads_log10']),
      max_reads      =  "--reads_max " + str(config["clustering_params"]['max_reads_log10']),
      min_peaks      =  "--peaks_min "  + str(config["clustering_params"]['min_peaks_ratio']),
      max_peaks      =  "--peaks_max "  + str(config["clustering_params"]['max_peaks_ratio']),
      peaks_file     = lambda wildcards, input: "--peaks_file " + input.peaks if wildcards.binwidth == "peaks" else "",
      fragments_file = lambda wildcards, input: "--fragments_file " + input.fragments_UMI if wildcards.umi == "UMI" else "--fragments_file " + input.fragments_noUMI,
      feature        = lambda wildcards: "--feature peaks" if wildcards.binwidth == "peaks" else "--feature bins --window " + str(wildcards.binwidth),
      binwidth       = "{binwidth}",
      umi            = "{umi}"

  shell:
      "Rscript /home/marek/bin/snaptools_scripts/snakemake/R/01.clustering_peaks.R {params.min_reads} "
                                                                                  "{params.max_reads} "
                                                                                  "{params.min_peaks} "
                                                                                  "{params.max_peaks} "
                                                                                  "{params.feature}   "
                                                                                  "{params.peaks_file} "
                                                                                  "{params.fragments_file} "
                                                                                  "--out_prefix {params.binwidth}_{params.umi}"


rule clustering_export_bigwig:
  input:
    "clustering_snakemake/01.clustering/Seurat_{binwidth}_{umi}.Rds" 
  output:   
    "clustering_snakemake/01.clustering/bigwig/{binwidth}_{umi}/{clusterMethod}_all_clusters.bw"
  params:
    clusterMethod = "{clusterMethod}"
  shell:
    "Rscript /home/marek/bin/snaptools_scripts/snakemake/R/02.export_bw.R {input} {params.clusterMethod}"

######################
rule merge_bam:
    input:
        expand("data/{sample}/outs/possorted_bam.bam",sample=config["samples"])
    output:
        "data/combined/outs/possorted_bam.bam"
    shell:
        "samtools merge {output} {input}"

rule modify_fragments_file_part1:
    input:
        "data/{sample}/outs/fragments.tsv.gz"
    output:
        "data/combined/outs/fragments_{sample}.tsv"
    params:
        sample="{sample}_"
    shell:
        """gzip -dc {input} | awk -v SAMPLE={params.sample} 'BEGIN {{FS=OFS="\t"}} {{print $1,$2,$3,SAMPLE$4,$5}}' -  > {output} """

rule modify_fragments_file_part2:
    input:
        expand("data/combined/outs/fragments_{sample}.tsv",sample=config["samples"])
    output:
        "data/combined/outs/fragments.tsv.gz"
    threads: 8
    shell:
        """sort -m -k 1,1V -k2,2n {input} > fragments.tsv &&
    bgzip -c -@ {threads} fragments.tsv  > {output} && 
    tabix -p bed {output} """


rule bam_to_bw:
    input:
        "data/{sample}/outs/possorted_bam.bam"
    output:
        "bigwig/{sample}_all_reads.bw"
    threads: 8
    shell:
        "bamCoverage -b {input} -o {output} -p {threads} --minMappingQuality 5 "
        " --binSize 100 --centerReads --smoothLength 500 --normalizeUsing RPKM --ignoreDuplicates"
        
        
rule run_macs_narrow:
    input:
        "data/{sample}/outs/possorted_bam.bam"
    output:
        "macs2_peaks/narrow/{sample}/{sample}_peaks.narrowPeak"
    params:
        macs_outdir="macs2_peaks/narrow/{sample}"
    shell:
        "macs2 callpeak -t {input} -g mm -f BAMPE -n {wildcards.sample} "
        "--outdir {params.macs_outdir} -q 0.05 -B --SPMR --keep-dup=1 2>&1 "
        

rule run_macs_broad:
    input:
        "data/{sample}/outs/possorted_bam.bam"
    output:
        "macs2_peaks/broad/{sample}/{sample}_peaks.broadPeak"
    params:
        macs_outdir="macs2_peaks/broad/{sample}"
    shell:
        "macs2 callpeak -t {input} -g mm -f BAMPE -n {wildcards.sample} "
        "--outdir {params.macs_outdir} -q 0.05 -B --SPMR --keep-dup=1 --broad-cutoff=0.1 --broad 2>&1 "

rule run_macs_merged_narrow:
    input:
        "data/combined/outs/possorted_bam.bam"
    output:
        "macs2_peaks.combined/narrow/combined_peaks.narrowPeak"
#        "macs2_peaks.combined/narrow/combined_summits.bed" 
    params:
        macs_outdir="macs2_peaks.combined/narrow"
    shell:
        "macs2 callpeak -t {input} -g mm -f BAMPE -n combined "
        "--outdir {params.macs_outdir} -q 0.05 -B --SPMR --keep-dup=1 2>&1 "
        

rule run_macs_merged_broad:
    input:
        "data/combined/outs/possorted_bam.bam"
    output:
        "macs2_peaks.combined/broad/combined_peaks.broadPeak"
    params:
        macs_outdir="macs2_peaks.combined/broad"
    shell:
        "macs2 callpeak -t {input} -g mm -f BAMPE -n combined "
        "--outdir {params.macs_outdir} -q 0.05 -B --SPMR --keep-dup=1 --broad-cutoff=0.1 --broad 2>&1 "

rule barcode_statistics_peaks:
  input:
    bam          = "data/{sample}/outs/possorted_bam.bam",
    peaks_broad  = "macs2_peaks/broad/{sample}/{sample}_peaks.broadPeak",
    peaks_narrow = "macs2_peaks/narrow/{sample}/{sample}_peaks.narrowPeak"
  output:
    "barcode_statistics/{sample}/barcode_statistics_peaks.log"
  params:
    out_narrow = "barcode_statistics/{sample}/peaks_barcodes_narrow.txt",
    out_broad  = "barcode_statistics/{sample}/peaks_barcodes_broad.txt"
  shell:
    "(bedtools intersect -abam {input.bam} -b {input.peaks_broad} -u | samtools view -f2 | "
    " awk -f ~/bin/snaptools_scripts/get_cell_barcode.awk | sed 's/CB:Z://g' |sort | uniq -c > {params.out_broad} && "
    " bedtools intersect -abam {input.bam} -b {input.peaks_narrow} -u | samtools view -f2 | "
    " awk -f ~/bin/snaptools_scripts/get_cell_barcode.awk | sed 's/CB:Z://g' |sort | uniq -c > {params.out_narrow} ) "
    " 2> {output} &&"
    " echo 'DONE' >> {output} "

rule remove_duplicates:
    input:
        "data/{sample}/outs/possorted_bam.bam"
    output:
        "barcode_statistics/{sample}/possorted_bam_nodup.bam"
    params:
        metrics = "barcode_statistics/{sample}/picard_markduplicates.log"
    shell:
        "java -Xmx4g -jar ~/bin/picard/picard.jar MarkDuplicates REMOVE_DUPLICATES=true INPUT={input} OUTPUT={output} METRICS_FILE={params.metrics}"


rule barcode_statistics_all:
  input:
    bam          = "data/{sample}/outs/possorted_bam.bam",
  output:
    "barcode_statistics/{sample}/barcode_statistics_all.log"
  params:
    out_narrow = "barcode_statistics/{sample}/all_barcodes_narrow.txt",
    out_broad  = "barcode_statistics/{sample}/all_barcodes_broad.txt"
  shell:
    "(samtools view -f2 {input.bam} | "
    " awk -f ~/bin/snaptools_scripts/get_cell_barcode.awk | sed 's/CB:Z://g' |sort | uniq -c > {params.out_broad} && "
    " samtools view -f2 {input.bam} | "
    " awk -f ~/bin/snaptools_scripts/get_cell_barcode.awk | sed 's/CB:Z://g' |sort | uniq -c > {params.out_narrow} ) "
    " 2> {output} &&"
    " echo 'DONE' >> {output} "

        
rule barcode_into_read_name:
    input:
        "data/{sample}/outs/possorted_bam.bam"
    output:
        "snaptools_scripts/{sample}_for_snaptools.bam"
    shell:
        "samtools view -h {input} | "
        "awk -f ~/bin/snaptools_scripts/barcode_read_name.awk | "
        "samtools view -bS - > {output}"

rule sort_bam:
    input:
        "snaptools_pre/{sample}_for_snaptools.bam"
    output:
        "snaptools_pre/{sample}_sort.bam"
    threads: 16
    shell:
        "samtools sort -n -@ {threads} {input} -o {output}"

rule download_mm10:
    output:
        os.getenv("HOME") + "/bin/snaptools_scripts/snaptools_resources/mm10/mm10.chrom.sizes"
    shell:
        "wget -O {output} http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"


rule download_blacklist:
  output:
    "mm10.blacklist.bed.gz"
  shell:
    "wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz"

        
# rule snap_pre:
#     input:
#         snap="snaptools_pre/{sample}_sort.bam",
#         genome=os.getenv("HOME") + "/bin/snaptools_scripts/snaptools_resources/mm10/mm10.chrom.sizes"
#     output:
#         "snaptools_pre/{sample}.snap"
#     shell:
#         "python2 ~/.local/bin/snaptools snap-pre --input-file={input.snap}\
#                             --output-snap={output}\
#                             --genome-name=mm10 \
#                             --genome-size={input.genome} \
#                             --min-mapq=30  \
#                             --min-flen=50  \
#                             --max-flen=1000  \
#                             --keep-chrm=TRUE  \
#                             --keep-single=FALSE  \
#                             --keep-secondary=False  \
#                             --overwrite=True  \
#                             --max-num=20000  \
#                             --min-cov=500  \
#                             --verbose=True"

# rule snap_pre:
#     input:
#         "snaptools_pre/{sample}.snap"
#     output:
#         "logs/{sample}_sort_out.log"
#     shell:
#         "(python2 ~/.local/bin/snaptools snap-add-bmat --snap-file={input} "
#                                 "--bin-size 1000 5000 10000 25000 "
#                                 "--verbose True && echo 'done') 2>&1 > {output}"



# rule motif_search:
#   input:
#       peaks     = "macs2_peaks.combined/narrow/combined_summits.bed",
#       blacklist = "mm10.blacklist.bed.gz"
#   output:
#       meme_out           = "MEME/out_{npeaks}/meme-chip.html"
#   params:
#       summits_filtered   = "MEME/summits_filtered_{npeaks}.bed",
#       top_summits        = "MEME/top_summits_{npeaks}.bed",
#       top_summits_padded = "MEME/top_summits_padded_{npeaks}.bed",
#       top_summits_fa     = "MEME/top_summits_{npeaks}.fa",
#       genome_fa          = config['genome_fa'],
#       npeaks             = "{npeaks}",
#       out                = "MEME/out_{npeaks}"
#   shell:
#       """
#       set +o pipefail;
#       cat {input.peaks} | grep -v -e 'chrM' | sort-bed - | bedops -n 1 - {input.blacklist} > {params.summits_filtered};
#       sort -k5gr {params.summits_filtered} | head -{params.npeaks} | sort-bed - > {params.top_summits};
#       bedops --range 150 -u {params.top_summits} > {params.top_summits_padded};
#       bedtools getfasta -fi {params.genome_fa} -bed {params.top_summits_padded} -fo {params.top_summits_fa};
#       /home/marek/anaconda3/envs/CT/bin/meme-chip -oc {params.out} -dreme-m 10 -meme-nmotifs 10 {params.top_summits_fa};
#       """
#
