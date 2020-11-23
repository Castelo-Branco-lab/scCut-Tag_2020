#!/usr/bin/env python3

import argparse
import pysam
import csv
import random
import gzip

parser = argparse.ArgumentParser(description='foo')
parser.add_argument('-f','--fragments',action='store')
parser.add_argument('-b','--barcodes',action='store')
parser.add_argument('-n','--number', action='store',type=int)
parser.add_argument('-c','--cluster',action='store')
parser.add_argument('-o','--output',action='store')

args = parser.parse_args()



def select_barcodes(args):
  bcd = []
  with open(args.barcodes) as f:
    barcodes = csv.reader(f)
    for line in barcodes:
      if args.cluster in line[0]:
        bcd.append(line[1])
  if args.number > len(bcd):
    args.number = len(bcd)
  bcd_selected = random.sample(bcd,int(args.number))
  return(bcd_selected)

def filter_bam(args,bcd):
  with open(args.output,'w') as o:
    with gzip.open(args.fragments) as f:
      tbx = pysam.tabix_iterator(f,pysam.asBed())
      for line in tbx:
        if line.name in bcd:
          o.write("{}\n".format(str(line)))
  return 0
  
  
def main(args):
  barcodes = select_barcodes(args)
  filter_bam(args,barcodes)
  
main(args)