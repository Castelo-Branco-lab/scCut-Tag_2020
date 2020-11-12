"""
Mapping Research Pipeline                                                                                                                                                                                
Copyleft 2018 Institut Curie                                                                                                                                                                             
Author(s): Nicolas Servant                                                                                                                                                                               
Contact: nicolas.servant@curie.fr                                                                                                                                                                        
This software is distributed without any guarantee under the terms of the CECILL License                                                                                                                 
See the LICENCE file for details                                                                                                                                                                         

Remove duplicates for single-cell ChIP-seq data
"""

import getopt
import sys
import os
import re
import pysam
from itertools import izip

def usage():
    """Usage function"""
    print "Usage : python addBarcodeFlag.py"
    print "-i/--input < mapped file [BAM]>"
    print "[-d/--dist] <distance to consider reads as duplicates [INT]>"
    print "[-o/--ofile] <output file [BAM]>"
    print "[-t/--tag] <tag>"
    print "[-v/--verbose] <Verbose>"
    print "[-h/--help] <Help>"
    return


def get_args():
    """Get arguments"""
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "i:d:o:t:vh",
            ["inputFile=",
             "dist=",
             "output=",
             "tag=", 
             "verbose", "help"])
    except getopt.GetoptError:
        usage()
        sys.exit(-1)
    return opts

def get_read_tag(read, tag):
    """
    Extract a flag from a read alignment
    """
    for t in read.tags:
        if t[0] == tag:
            return t[1]
    return None


def get_read_start(read):
    """                                                                                                                                                                                        
    Return the 5' end of the read                                                                                                                                                              
    """
    if read.is_reverse:
        pos = read.pos + read.alen -1
    else:
        pos = read.pos
    return pos


if __name__ == "__main__":

    opts = get_args()
    verbose = False
    tag = "XB"
    output = "-"
    reads_counter = 0
    dup_counter = 0
    dist = 150
    ref_perbarcode = {}
    #ref_perbarcode_fwd = {}
    #ref_perbarcode_rev = {}

    if len(opts) == 0:
        usage()
        sys.exit()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            inputFile = arg
        elif opt in ("-o", "--ofile"):
            output = arg
        elif opt in ("-t", "--tag"):
            tag = arg
        elif opt in ("-d", "--dist"):
            dist = int(arg)
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    if verbose:
        print "## rmDup.py"
        print "## input =", inputFile
        print "## dist =", dist
        print "## tag =", tag
        print "## output =", output
        print "## verbose =", verbose
        print

    ## Read bam file
    if verbose:
        print "## Reading BAM file ..."
    samfile = pysam.Samfile(inputFile, "rb")

    ## output
    if output == "-":
        outfile = pysam.AlignmentFile(output, "w", template=samfile)
    else:
        outfile = pysam.AlignmentFile(output, "wb", template=samfile)

    for r1 in samfile.fetch(until_eof=True):
        reads_counter += 1
        
        ## Get Barcode
        barcode = str(get_read_tag(r1, tag))
        
        ## Distinguish reverse and forward reads
##        if not r1.is_reverse:
            ## Set new reference
        if barcode not in ref_perbarcode:
            ref_perbarcode[barcode] = r1
            outfile.write(r1)                  
        else:
            ## Compare with existing reference
            ref = ref_perbarcode[barcode]
            ref_start = get_read_start(ref)
            r1_start = get_read_start(r1)

            ## Is duplicates
            if ref.tid == r1.tid and r1_start < (ref_start + dist) :
                dup_counter += 1
            else:
                    ## update reference
                ref_perbarcode[barcode] = r1
                outfile.write(r1)
##        else:
##            if barcode not in ref_perbarcode_rev:
##                ref_perbarcode_rev[barcode] = r1
##                outfile.write(r1)
##            else:
##                ref = ref_perbarcode_rev[barcode]
##                ref_start = get_read_start(ref)
##                r1_start = get_read_start(r1)

##                if ref.tid == r1.tid and r1_start < (ref_start + dist) and r1_start > (ref_start - dist)  :
##                    dup_counter += 1
##                else:
##                    ref_perbarcode_rev[barcode] = r1
# #                   outfile.write(r1)
        
        if (reads_counter % 1000000 == 0 and verbose):
            print "##", reads_counter

samfile.close()

if verbose:
    print "## Number of reads: " + str(reads_counter)
    print "## Number of duplicates: " + str(dup_counter)
    print "## Number of reads after duplicates removal: " + str(reads_counter - dup_counter)

