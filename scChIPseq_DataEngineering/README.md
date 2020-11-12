---
author: "Nicolas Servant; Pacome Prompsy"
date: "01st Octobre 2018"
output: html_document
---

# scChIPseq_DataEngineering
This pipeline is dedicated to create single-cell count matrix from paired-end raw fastQ files coming from single-cell ChIP-seq experiments.

# single-cell ChIP-seq pipeline
# v0.0.2

```
usage : schip_processing -f FORWARD -r REVERSE -o OUTPUT -c CONFIG [-d] [-h] [-v]
Use option -h|--help for more information

schip_processing 0.0.2
---------------
OPTIONS

   -f|--forward R1_READ: forward fastq file
   -r|--reverse R2_READ: forward fastq file
   -c|--conf CONFIG: configuration file for ChIP processing
   -o|--output OUTPUT: output folder
   [-d|--dryrun]: dry run mode
   [-h|--help]: help
   [-v|--version]: version
```


=======
