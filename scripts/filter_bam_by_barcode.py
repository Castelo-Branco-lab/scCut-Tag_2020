#!/usr/bin/env python3

import sys,os
import pysam
from contextlib import ExitStack


bamFile             = sys.argv[1]      # Possorted bam file from 10x with barcodes           e.g. AAACGAAAGTTGCTTG-1
barcode_annotations = sys.argv[2]      # cluster-barcode csv.                                e.g. mOL,H3K27ac_N1_AAACGAAAGTTGCTTG-1
barcode_prefix      = sys.argv[3]      # Sample name to prefix to bam barcode                e.g. mOL,H3K27ac_N1  (Use "NA" if no prefix is needed)
out_prefix          = sys.argv[4]      # folder to write the output to                       e.g. /path/to/output

if barcode_prefix == "NA":
  barcode_prefix = ""

sys.stderr.write("*** Reading cluster - barcode csv file ***\n\n")

# Parse cluster file into dictionary
clusters_dic = {}
for line in open(barcode_annotations,'r'):
  if line.startswith("#"):
    continue
  line = line.rstrip().split(',')
  clusters_dic[line[1]] = line[0]

# Iterate over clusters and generate list of paths for output files
clusters = list(set(clusters_dic.values()))
clusters_outfiles = {x: out_prefix + "/" + x + "_out.bam" for x in clusters}

sys.stderr.write("*** Found following clusters in cluster - barcode file ***\n")
sys.stderr.write("\n".join(clusters) + "\n\n")

sys.stderr.write("*** Creating following output files ***\n")
print("\n".join(clusters_outfiles.values()) + "\n\n")

# Open bam file and save the header 
bamfile = pysam.AlignmentFile(bamFile, "rb")
header  = bamfile.header


# Iterate over bamfile
# Write as bam

with ExitStack() as stack:
    files = {fname: stack.enter_context(pysam.AlignmentFile(fname,'wb',header=header)) for fname in list(clusters_outfiles.values())}
    # Iterate over the bam file
    n = 0
    for line in bamfile:
      n+=1
      if n % 5000000 == 0:
        sys.stderr.write("*** {} lines processed ***\n".format(n))
      try:
        barcode = line.get_tag("CB")
      except KeyError:
        continue
      if barcode_prefix != "":
        barcode = barcode_prefix + "_" + barcode
      if barcode in clusters_dic:
        cluster = clusters_dic[barcode]
        files[clusters_outfiles[cluster]].write(line)