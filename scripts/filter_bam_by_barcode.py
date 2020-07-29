#!/usr/bin/env python3

import sys
import pysam
from contextlib import ExitStack


bamFile             = sys.argv[1]
barcode_annotations = sys.argv[2]
sample_id           = sys.argv[3]

sys.stderr.write("*** Reading cluster - barcode csv file ***\n\n")

# Parse cluster file into dictionary
clusters_dic = {}
for line in open(barcode_annotations,'r'):
  if line.startswith("#"):
    continue
  line = line.rstrip().split(',')
  clusters_dic[line[1]] = line[0]

clusters = list(set(clusters_dic.values()))
clusters_outfiles = {x: x + "_out.sam" for x in clusters}

sys.stderr.write("*** Found following clusters in cluster - barcode file ***\n")
sys.stderr.write("\n".join(clusters) + "\n\n")

sys.stderr.write("*** Creating following output files ***\n")
print("\n".join(clusters_outfiles.values()) + "\n\n")

# Open bam file
bamfile = pysam.AlignmentFile(bamFile, "rb")
header  = bamfile.header


    # Create sam headers in output files
for key in clusters_outfiles:
  f = pysam.AlignmentFile(clusters_outfiles[key],'w',header = header)
  f.close()
    

# for line in bamfile:
#   n+=1
#   if n % 10000 == 0:
#     sys.stderr.write("*** {} lines processed ***\n".format(n))
#   try:
#     barcode = line.get_tag("CB")
#   except KeyError:
#     continue
#
#   barcode = sample_id + "_" + barcode
#   if barcode in clusters_dic:
#     with open(clusters_outfiles[clusters_dic[barcode]],'a+') as f:
#       f.write("{}\n".format(line.tostring(bamFile)))
      
with ExitStack() as stack:
    files = {fname: stack.enter_context(open(fname,'w')) for fname in list(clusters_outfiles.values())}
    # Iterate over the bam file
    n = 0
    for line in bamfile:
      n+=1
      if n % 100000 == 0:
        sys.stderr.write("*** {} lines processed ***\n".format(n))
      try:
        barcode = line.get_tag("CB")
      except KeyError:
        continue
    
      barcode = sample_id + "_" + barcode
      if barcode in clusters_dic:
        cluster = clusters_dic[barcode]
        files[clusters_outfiles[cluster]].write("{}\n".format(line.tostring(bamFile)))