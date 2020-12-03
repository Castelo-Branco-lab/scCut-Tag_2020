#!/usr/bin/env python3

import argparse
import pysam
import csv
import random
import gzip

parser = argparse.ArgumentParser(description='foo')
parser.add_argument('-f','--fragments',action='store')
parser.add_argument('-b','--barcodes',action='store')
parser.add_argument('-o','--output',action='store')
args = parser.parse_args()

def parse_barcodes(args):
  bcd = {}
  with open(args.barcodes) as f:
    barcodes = csv.reader(f)
    for line in barcodes:
      bcd[line[1]] = line[0]
  return(bcd)

def filter_bam(args,bcd):
  reads     = {}
  replicate = {}
  with gzip.open(args.fragments) as f:
    tbx = pysam.tabix_iterator(f,pysam.asBed())
    for line in tbx:
      if line.name in bcd:
        try:
          reads[bcd[line.name] + "_rep" + line.name.split("_")[-2]].append(str(line))
        except KeyError:
          reads[bcd[line.name] + "_rep" + line.name.split("_")[-2]]     = [str(line)]
          
  return(reads)
  
  
def main(args):
  barcodes  = parse_barcodes(args)
  reads_dic = filter_bam(args,barcodes)
  for key in reads_dic:
    with open(args.output + key + ".bed",'w') as f:
      f.write("\n".join(reads_dic[key]))
  
main(args)