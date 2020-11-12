#!/usr/bin/env python3

import argparse
import pysam
import sys
import gzip

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("--bam_file","-b",type=str,action="store")
parser.add_argument("--out","-o",type=str,action="store")

args = parser.parse_args()

def main(args):
  N=1
  with pysam.AlignmentFile(args.bam_file,'r') as p:
    with pysam.AlignmentFile(args.out,'wb',template = p) as out:
      for line in p:
        line.tags += [("RG",line.query_name)]
        line.query_name = line.query_name + "_read_" + str(N)
        out.write(line)
        N += 1

main(args)