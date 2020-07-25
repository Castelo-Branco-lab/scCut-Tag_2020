#!/usr/bin/env python3

import sys
import yaml
import pysam

# Required are 4 arguments
bam_path        = sys.argv[1]    # Original possorted_bam.bam file
config_path     = sys.argv[2]    # Path to snakemake config file
pass_bcd        = sys.argv[3]    # Csv of passed barcodes exported from R script
out_path        = sys.argv[4]    # Path to the output bam file
out_path_nopass = sys.argv[5]    # Path to the output bam with reads that do not pass
# Read the yaml config file
config         = yaml.safe_load(open(config_path,'r'))
config_samples = config["samples"]
config_names   = {config_samples[sample].split("/")[-1]: sample for sample in config_samples}

# Read the barcodes csv file
bcd = []
for line in open(pass_bcd):
    bcd.append(line.rstrip())
    

# Open the possorted_bam.bam file
samfile = pysam.AlignmentFile(bam_path, "rb")


# Iterate over the bam file and filter it
with pysam.AlignmentFile(out_path, "wb",header=samfile.header) as f, pysam.AlignmentFile(out_path_nopass, "wb",header=samfile.header) as g:
    for read in samfile:
        try:
            sample  = read.get_tag("RG:Z")
            barcode = read.get_tag("CB:Z")
        except(KeyError):
            g.write(read)
            continue      # If any of the tags is missing skip
        
        sample     = sample.split(":")[0]                     # sample name is first part in the RG:Z: tag in the bam file
        barcode_id = config_names[sample] + "_" + barcode     # Create sample_barcode string
        
        if barcode_id in bcd:                                 # Compare sample_barcode string with barcodes in the csv file
            f.write(read)
        else:
            g.write(read)
            



