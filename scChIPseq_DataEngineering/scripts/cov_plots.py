"""
Produce pngs for the main chromosomes based on bigwig file
"""

import sys
import os
import re
import pyBigWig as pyBW
import argparse
import collections
import subprocess

if __name__ == "__main__":
    # Reads args
    parser = argparse.ArgumentParser(prog='cov_plots.py', description='''
        Produce pngs for the main chromosomes based on bigwig file
        ''', epilog='''
        ''')
    parser.add_argument('-i', '--input', help="BigWig file with barcode tag (.bw/.bigwig)", required=True)
    parser.add_argument('-o', '--output', help="Output directory", required=True, type=str)

    args = parser.parse_args()

    # check args
    if not os.path.isdir(args.output):
        try:
            os.makedirs(args.output)
            print("Creating directory " + args.output)
        except:
            print("Can't create directory " + args.output +", path not found")
            sys.exit(-1)




    if args.input.endswith(".bw") or args.input.endswith(".bigwig"):
        track = pyBW.open(args.input)
        if track is not None:
            chrom_to_plot = collections.OrderedDict()

            for chrom in track.chroms():
                if "_" not in chrom:
                    chrom_to_plot[chrom] = (track.chroms()[chrom])

            if len(chrom_to_plot) > 0:
                #Run the make_track_file from pyGenomeTracks library to create the config file including the .bw file
                bashCommand = str("python3 /home/pprompsy/.local/bin/make_tracks_file --trackFiles " + args.input + " -o " + args.output + "/tracks.ini")
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()

                for chrom in chrom_to_plot:
                    #Run the pyGenomeTracks_file from pyGenomeTracks library to create the coverage plots as .png
                    bashCommand2 = str(
                        "python3 /home/pprompsy/.local/bin/pyGenomeTracks --tracks " + args.output + "/tracks.ini --region " + str(
                            chrom + ":1-" + str(chrom_to_plot[
                                                    chrom])) + " --outFileName " + args.output + "/" + chrom + "_mqc.png")
                    process = subprocess.Popen(bashCommand2.split(), stdout=subprocess.PIPE)
                    output, error = process.communicate()

    else:
        print("Input bigwig file is in the wrong format, exiting")
        sys.exit(-1)
    print("Finished creating the " +str(len(chrom_to_plot)) + " coverage plots")
