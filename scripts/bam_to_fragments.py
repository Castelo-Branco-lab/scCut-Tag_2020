#!/usr/bin/env python3

import sys
import pysam
import yaml
import gzip

COUNTER_ALL       = 0
COUNTER_PROPER    = 0
COUNTER_NOMATCH   = 0
COUNTER_NOBARCODE = 0

verbose = False

bam_path = sys.argv[1]
out_path = sys.argv[2]

sys.stderr.write("\t".join(sys.argv) + "\n")

config_path = "config.yaml"

# Read config yaml file into dictionary
config         = yaml.safe_load(open(config_path,'r'))
config_samples = config["samples"]
config_names   = {config_samples[sample].split("/")[-1]: sample for sample in config_samples}

# Open bam file using pysam
samfile = pysam.AlignmentFile(bam_path, "rb")

with open(out_path,"w") as f:
  for read in samfile:
      COUNTER_ALL += 1
    
      # First check if read is a proper pair
      if not read.is_proper_pair:
          if verbose:
              sys.stderr.write("***Warning: Read is not part of a proper pair*** %s \n" % read.query_name )
          COUNTER_PROPER += 1
          continue
    
      # Match RG:Z tag to sample names in config file: (Returns index of the sample)
      sample_index = [i for i,x in enumerate(config_names.keys()) if x in read.get_tag("RG:Z")]    
    
      # Should return only one match, if not continue
      if len(sample_index) != 1:
          if verbose:
              sys.stderr.write("***Warning: could not match read to single sample: Skipping*** %s \n" % read.query_name )
          COUNTER_NOMATCH +=1 
          continue
        
      # List to integer - index of sample in config_names.keys()
      sample_index = int(sample_index[0])
    
      # Get sample id from index
      sample_id = [config_names[x] for x in config_names.keys()][sample_index]
    
      # Get barcode from bam file
      try:
          barcode = read.get_tag("CB:Z")
      except KeyError:
          if verbose:
              sys.stderr.write("*** WARNING: Not cell barcode tag present for read %s, skipping *** \n" % read.query_name)
          COUNTER_NOBARCODE += 1
          continue
    
      BED = [str(x) for x in [read.reference_name,read.pos,read.reference_end,sample_id + "_" + barcode,1]]
      f.write("%s\n" % "\t".join(BED))
    
    
    

sys.stderr.write("Reads processed: {0}\nSkipped reads:\n\tnot proper:\t\t{1}\t{2:.2%}\t\n\tno library match:\t{3}\t{4:.2%}\n\tno cell barcode:\t{5}\t{6:.2%}\n".format(COUNTER_ALL,COUNTER_PROPER,COUNTER_PROPER/COUNTER_ALL,COUNTER_NOMATCH,COUNTER_NOMATCH/COUNTER_ALL,COUNTER_NOBARCODE,COUNTER_NOBARCODE/COUNTER_ALL))
pysam.tabix_index(out_path,preset="bed",force="True")



    
    
    