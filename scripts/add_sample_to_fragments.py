#!/usr/bin/env python3

import sys
import gzip

fragments_file = sys.argv[1]
sample_name    = sys.argv[2]

with gzip.open(fragments_file,'rb') as f:
   for line in f:
      line = line.decode('UTF-8').rstrip().split('\t')
      line[3] = sample_name + "_" + line[3]

      sys.stdout.write("\t".join(line) + "\n")