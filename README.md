# single cell in the brain Cut-Tag 2020
### Marek Bartosovic
### Goncalo Castelo-Branco lab


This repo contains code needed to generate figures for the paper XY (link doi)

Code is organized in snakemake pipelines

1. Fastq files are processed using standard cellranger-atac count
2. First preprocess snakemake pipeline is run to generate bulk tracks, cell barcode statistics and do cells identification
3. Second snakemake pipeline contains R markdown notebooks used for merging of the samples, analysis and figure generation

## Step 1. 

Each sequencing run is processed separately with cellranger:

`cellranger-atac count --fastqs=./fastq/ --reference=PATH_TO_CELLRANGER_REFERENCE_MM10 --sample=SEQUNCING_ID --id=SAMPLE_ID`


## Step 2.

Clone git repo
`git clone https://github.com/mardzix/scCut-Tag_2020/master`

Create conda environment with all necessary tools installed by:
`conda env create -f scCut-Tag_2020/envs/CT_snakemake.yaml `

Provide cell ranger oudput folder to snakemake pipeline step2.Snakemake







