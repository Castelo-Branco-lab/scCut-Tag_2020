#!/usr/bin/env python3

import argparse
import pysam
import csv
import random
import gzip
import sys, os

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
  sys.stderr.write("{} barcodes read\n".format(len(bcd)))
  return(bcd)

def filter_bam(args,bcd):
  reads     = {}
  replicate = {}
  i = 0
  with gzip.open(args.fragments) as f:
    tbx = pysam.tabix_iterator(f,pysam.asBed())
    for line in tbx:
      i += 1
      if line.name in bcd:
        try:
          reads[bcd[line.name] + "_rep" + line.name.split("_")[-2]].append(str(line))
        except KeyError:
          reads[bcd[line.name] + "_rep" + line.name.split("_")[-2]]     = [str(line)]
  sys.stderr.write("{} reads processed \n".format(i))
  return(reads)

def add_forward_slash(s):
  if not s.endswith("/"):
    s = s + "/"
    return s
    
def create_directory(d):
  if not os.path.isdir(d):
    os.makedirs(d)
  
def main(args):
  args.output = add_forward_slash(args.output)
  create_directory(args.output)
  
  barcodes  = parse_barcodes(args)
  reads_dic = filter_bam(args,barcodes)
  for key in reads_dic:
    with open(args.output + key + ".bed",'w') as f:
      f.write("\n".join(reads_dic[key]))
  
main(args)