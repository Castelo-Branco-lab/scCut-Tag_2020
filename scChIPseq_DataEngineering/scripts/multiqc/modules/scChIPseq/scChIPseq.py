#!/usr/bin/env python

""" MultiQC module to parse output from scChIPseq pipeline """

from __future__ import print_function
from collections import OrderedDict
import logging
import os
import re
import pandas as pd
import subprocess
import pyBigWig as pyBW
from multiqc import config
from multiqc.plots import bargraph
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule
from itertools import chain
from multiqc.plots import linegraph
import math
# Initialise the logger
log = logging.getLogger(__name__)
# Initialise your class and so on
class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='scChIPseq', anchor='scChIPseq',
        href="https://gitlab.curie.fr/data-analysis/ChIP-seq_single-cell_LBC",
        info="is a DNA alignment pipeline dedicated to single-cell ChIP-seq experiments")

        # Find and load any scChIPseq reports
        self.scChIPseq_data = dict()
        for f in self.find_log_files('scChIPseq/all_logs'):
            log.info('Found the all_logs!')
            parsed_data = self.parse_scChIPseq_report(f['f'])
            if parsed_data is not None:
                s_name = f['s_name']
                if s_name == '':
                    s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
                if s_name in self.scChIPseq_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.add_data_source(f, section='SummaryLog')
                self.scChIPseq_data[s_name] = parsed_data

        # Read in flagged_count
        self.scChIPseq_flagged_count = dict()
        for f in self.find_log_files('scChIPseq/flagged_count'):
            log.info('Found the flagged_count !')
            colnames = ['count', 'barcode']
            if not f['root']:
                log.info("is empty")
                count = pd.read_csv("./" + f['fn'], delim_whitespace=True, names=colnames)

            else:
                log.info("is not empty")
                count = pd.read_csv(f['root'] +"/" + f['fn'], delim_whitespace=True, names=colnames)
            s_name = f['s_name']
            if count is not None:
                if s_name == '':
                    s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
                if s_name in self.scChIPseq_flagged_count:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.scChIPseq_flagged_count[s_name] = count


        #Read in flagged_PCR_count
        self.scChIPseq_flagged_PCR_count = dict()
        for f in self.find_log_files('scChIPseq/flagged_PCR_count'):
            log.info('Found the scChIPseq_flagged_PCR_count !')
            colnames = ['count', 'barcode']
            if not f['root']:
                log.info("is empty")
                count = pd.read_csv("./" + f['fn'], delim_whitespace=True, names=colnames)
            else:
                log.info("is not empty")
                count = pd.read_csv(f['root'] + "/" + f['fn'], delim_whitespace=True, names=colnames)
            s_name = f['s_name']
            if count is not None:
                if s_name == '':
                    s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
                if s_name in self.scChIPseq_flagged_PCR_count:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.scChIPseq_flagged_PCR_count[s_name] = count

        #Read in flagged_PCR_RT_rmDup_count
        self.scChIPseq_flagged_PCR_RT_count = dict()
        for f in self.find_log_files('scChIPseq/flagged_PCR_RT_count'):
            log.info('Found the scChIPseq_flagged_PCR_RT_count !')
            colnames = ['count', 'barcode']
            if not f['root']:
                log.info("is empty")
                count = pd.read_csv("./" + f['fn'], delim_whitespace=True, names=colnames)
            else:
                log.info("is not empty")
                count = pd.read_csv(f['root'] + "/" + f['fn'], delim_whitespace=True, names=colnames)
            s_name = f['s_name']
            if count is not None:
                if s_name == '':
                    s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
                if s_name in self.scChIPseq_flagged_PCR_RT_count:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.scChIPseq_flagged_PCR_RT_count[s_name] = count

        # Read in flagged_PCR_RT_rmDup_count
        self.scChIPseq_flagged_PCR_RT_rmDup_count = dict()
        for f in self.find_log_files('scChIPseq/flagged_PCR_RT_rmDup_count'):
            log.info('FOUND THE scChIPseq_flagged_PCR_RT_rmDup_count !')
            colnames = ['count', 'barcode']
            if not f['root']:
                log.info("is empty")
                count = pd.read_csv("./" + f['fn'], delim_whitespace=True, names=colnames)
            else:
                log.info("is not empty")
                count = pd.read_csv(f['root'] + "/" + f['fn'], delim_whitespace=True, names=colnames)
            s_name = f['s_name']
            if count is not None:
                if s_name == '':
                    s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
                if s_name in self.scChIPseq_flagged_PCR_RT_rmDup_count:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.scChIPseq_flagged_PCR_RT_rmDup_count[s_name] = count

        # Read in count_matrix
        self.scChIPseq_count_matrix = dict()
        for f in self.find_log_files('scChIPseq/count_matrix'):
            log.info('FOUND THE scChIPseq_count_matrix !')
            if not f['root']:
                log.info("is empty")
                count = pd.read_csv("./" + f['fn'], delim_whitespace=True)
            else:
                log.info("is not empty")
                count = pd.read_csv(f['root'] + "/" + f['fn'], delim_whitespace=True)
            s_name = f['s_name']
            if count is not None:
                if s_name == '':
                    s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
                if s_name in self.scChIPseq_count_matrix:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.scChIPseq_count_matrix[s_name] = count
        # Filter to strip out ignored sample names
        self.scChIPseq_data = self.ignore_samples(self.scChIPseq_data)
        #self.scChIPseq_flagged_count = self.ignore_samples(self.scChIPseq_flagged_count)
        #self.scChIPseq_flagged_PCR_count = self.ignore_samples(self.scChIPseq_flagged_PCR_count)
        #self.scChIPseq_flagged_PCR_RT_count = self.ignore_samples(self.scChIPseq_flagged_PCR_RT_count)

       # if len(self.scChIPseq_data) == 0 or len(self.scChIPseq_flagged_count) == 0 or len(self.scChIPseq_flagged_PCR_count) == 0 or len(self.scChIPseq_flagged_PCR_RT_count) == 0:
       #     raise UserWarning

        if len(self.scChIPseq_data) > 0:
            log.info("Found {} reports".format(len(self.scChIPseq_data)))
        if len(self.scChIPseq_flagged_count) > 0:
            log.info("Found {} reports".format(len(self.scChIPseq_flagged_count)))
        if len(self.scChIPseq_flagged_PCR_count) > 0:
            log.info("Found {} reports".format(len(self.scChIPseq_flagged_PCR_count)))
        if len(self.scChIPseq_flagged_PCR_RT_count) > 0:
            log.info("Found {} reports".format(len(self.scChIPseq_flagged_PCR_RT_count)))
        if len(self.scChIPseq_flagged_PCR_RT_rmDup_count) > 0:
            log.info("Found {} reports".format(len(self.scChIPseq_flagged_PCR_RT_rmDup_count)))
        if len(self.scChIPseq_count_matrix) > 0:
            log.info("Found {} reports".format(len(self.scChIPseq_flagged_PCR_RT_rmDup_count)))

        if len(self.scChIPseq_data) > 0:
            # Write parsed report data to a file
            self.write_data_file(self.scChIPseq_data, 'multiqc_scChIPseq')
            # Basic Stats Table
            self.scChIPseq_stats_table()


            # Barcode matching bar plot
            self.add_section(
                name='Barcode Matching',
                anchor='scChIPseq_barcode',
                plot=self.scChIPseq_barcode_chart()
            )

            # Alignment bar plot
            self.add_section (
                name = 'Alignment Scores',
                anchor = 'scChIPseq_alignments',
                plot = self.scChIPseq_alignment_chart()
            )
        if len(self.scChIPseq_flagged_count) > 0:
            self.scChIPseq_flagged_coverage_chart()

        if len(self.scChIPseq_flagged_PCR_count) > 0:
            self.scChIPseq_flagged_PCR_coverage_chart()
        if len(self.scChIPseq_flagged_PCR_RT_count) > 0:
            self.scChIPseq_flagged_PCR_RT_coverage_chart()
        if len(self.scChIPseq_flagged_PCR_RT_rmDup_count) > 0:
            self.scChIPseq_flagged_PCR_RT_rmDup_coverage_chart()
        if len(self.scChIPseq_count_matrix) > 0:
            self.scChIPseq_count_matrix_region_coverage_chart()
            self.scChIPseq_count_matrix_cell_coverage_chart()

    def parse_scChIPseq_report (self, raw_data):
        """ Parse the combined scChIPseq log file. """

        regexes = {
            'total_reads':                  r"Number of input reads \|\s+(\d+)",
            'avg_input_read_length':        r"Average input read length \|\s+([\d\.]+)",
            'uniquely_mapped':              r"Uniquely mapped reads number \|\s+(\d+)",
            'uniquely_mapped_percent':      r"Uniquely mapped reads % \|\s+([\d\.]+)",
            'avg_mapped_read_length':       r"Average mapped length \|\s+([\d\.]+)",
            'num_splices':                  r"Number of splices: Total \|\s+(\d+)",
            'num_annotated_splices':        r"Number of splices: Annotated \(sjdb\) \|\s+(\d+)",
            'num_GTAG_splices':             r"Number of splices: GT/AG \|\s+(\d+)",
            'num_GCAG_splices':             r"Number of splices: GC/AG \|\s+(\d+)",
            'num_ATAC_splices':             r"Number of splices: AT/AC \|\s+(\d+)",
            'num_noncanonical_splices':     r"Number of splices: Non-canonical \|\s+(\d+)",
            'mismatch_rate':                r"Mismatch rate per base, % \|\s+([\d\.]+)",
            'deletion_rate':                r"Deletion rate per base \|\s+([\d\.]+)",
            'deletion_length':              r"Deletion average length \|\s+([\d\.]+)",
            'insertion_rate':               r"Insertion rate per base \|\s+([\d\.]+)",
            'insertion_length':             r"Insertion average length \|\s+([\d\.]+)",
            'multimapped':                  r"Number of reads mapped to multiple loci \|\s+(\d+)",
            'multimapped_percent':          r"% of reads mapped to multiple loci \|\s+([\d\.]+)",
            'multimapped_toomany':          r"Number of reads mapped to too many loci \|\s+(\d+)",
            'multimapped_toomany_percent':  r"% of reads mapped to too many loci \|\s+([\d\.]+)",
            'unmapped_mismatches_percent':  r"% of reads unmapped: too many mismatches \|\s+([\d\.]+)",
            'unmapped_tooshort_percent':    r"% of reads unmapped: too short \|\s+([\d\.]+)",
            'unmapped_other_percent':       r"% of reads unmapped: other \|\s+([\d\.]+)",
            'match_index_1':                r"## Number of matched indexes 1:\s+([\d\.]+)",
            'match_index_2':                r"## Number of matched indexes 2:\s+([\d\.]+)",
            'match_index_1_2':              r"## Number of matched indexes 1 and 2:\s+([\d\.]+)",
            'match_index_3':                r"## Number of matched indexes 3:\s+([\d\.]+)",
            'match_barcode':                r"## Number of matched barcodes:\s+([\d\.]+)",
            'uniquely_mapped_and_barcoded': r"## Number of reads mapped and barcoded:\s+([\d\.]+)",
            'pcr_duplicates':               r"## Number of pcr duplicates:\s+([\d\.]+)",
            'rt_duplicates':                r"## Number of rt duplicates:\s+([\d\.]+)",
            'R1_mapped_R2_unmapped':        r"## Number of R1 mapped but R2 unmapped:\s+([\d\.]+)",
            'reads_after_pcr_rt_rm':        r"## Number of reads after PCR and RT removal \(not R1 unmapped R2\):\s+([\d\.]+)",
            'R2_unmapped_duplicates':       r"## Number of duplicates:\s+([\d\.]+)",
            'unique_reads':                 r"## Number of reads after duplicates removal:\s+([\d\.]+)"
        }
        parsed_data = {}
        for k, r in regexes.items():
            r_search = re.search(r, raw_data, re.MULTILINE)
            if r_search:
                parsed_data[k] = float(r_search.group(1))
        # Figure out the numbers for unmapped as for some reason only the percentages are given
        try:
            total_mapped = parsed_data['uniquely_mapped'] + parsed_data['multimapped'] + parsed_data['multimapped_toomany']
            unmapped_count = parsed_data['total_reads'] - total_mapped
            total_unmapped_percent = parsed_data['unmapped_mismatches_percent'] + parsed_data['unmapped_tooshort_percent'] + parsed_data['unmapped_other_percent']

            parsed_data['uniquely_mapped_unbarcoded'] =  int(round(parsed_data['uniquely_mapped']-parsed_data['uniquely_mapped_and_barcoded']))
            parsed_data['multimapped'] = int(round(parsed_data['multimapped'] + parsed_data['multimapped_toomany']))
            parsed_data['unmapped'] = unmapped_count
            #Data for the barcode matching graph
            parsed_data['reads_after_pcr_rt_rm']=parsed_data['reads_after_pcr_rt_rm'] - parsed_data['R1_mapped_R2_unmapped']
            parsed_data['index_1_2_not_3'] = int(round(parsed_data['match_index_1_2'] - parsed_data['match_barcode']))
            parsed_data['index_1_not_2_not_3'] = int(round(parsed_data['match_index_1'] - parsed_data['index_1_2_not_3'] - parsed_data['match_barcode']))
            parsed_data['index_2_not_1_3'] = int(round(parsed_data['match_index_2'] - parsed_data['match_index_1_2']))
            parsed_data['index_3_not_1_2'] = int(round(parsed_data['match_index_3'] - parsed_data['match_barcode']))
            parsed_data['no_index_found'] = int(round(parsed_data['total_reads'] - parsed_data['match_barcode'] - parsed_data['index_1_2_not_3'] - parsed_data['index_1_not_2_not_3'] - parsed_data['index_2_not_1_3'] - parsed_data['index_3_not_1_2']))
            parsed_data['uniquely_mapped_and_barcoded_percent'] = 100*parsed_data['uniquely_mapped_and_barcoded'] / parsed_data['total_reads']
            parsed_data['unique_reads_percent'] = 100 * parsed_data['unique_reads'] / \
                                                                  parsed_data['total_reads']
            log.info(parsed_data['uniquely_mapped_and_barcoded_percent'])
            try:
                parsed_data['unmapped_mismatches'] = int(round(unmapped_count * (parsed_data['unmapped_mismatches_percent'] / total_unmapped_percent), 0))
                parsed_data['unmapped_tooshort'] = int(round(unmapped_count * (parsed_data['unmapped_tooshort_percent'] / total_unmapped_percent), 0))
                parsed_data['unmapped_other'] = int(round(unmapped_count * (parsed_data['unmapped_other_percent'] / total_unmapped_percent), 0))
            except ZeroDivisionError:
                parsed_data['unmapped_mismatches'] = 0
                parsed_data['unmapped_tooshort'] = 0
                parsed_data['unmapped_other'] = 0
        except KeyError:
            pass

        if len(parsed_data) == 0: return None
        return parsed_data

    def scChIPseq_stats_table(self):
        """ Take the parsed stats from the STAR report and add them to the
        basic stats table at the top of the report """

        headers = OrderedDict()
        headers['uniquely_mapped_percent'] = {
            'title': '% Aligned',
            'description': '% Uniquely mapped reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }
        headers['uniquely_mapped_and_barcoded_percent'] = {
            'title': '% Aligned and Barcoded',
            'description': '% Aligned and Barcoded reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }
        headers['unique_reads_percent'] = {
            'title': '% Unique Reads',
            'description': '% Unique Reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'YlGn'
        }
        self.general_stats_addcols(self.scChIPseq_data, headers)

    def scChIPseq_alignment_chart (self):
        """ Make the plot showing alignment rates """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['unique_reads'] = {'color': '#00bf00', 'name': 'Deduplicated reads'}
        keys['R2_unmapped_duplicates'] = {'color': '#00e887', 'name': '\"Window\" duplicates'}
        keys['rt_duplicates'] = {'color': '#0c7bd1', 'name': 'RT duplicates'}
        keys['pcr_duplicates'] = {'color': '#4914e8', 'name': 'PCR duplicates'}
        keys['uniquely_mapped_unbarcoded'] = {'color': '#b5d30c', 'name': 'Uniquely mapped, not barcoded'}
        keys['multimapped'] = {'color': '#edb900', 'name': 'Mapped to multiple loci'}
        keys['unmapped'] = {'color': '#ff2c20', 'name': 'Unmapped'}

        # Config for the plot
        pconfig = {
            'id': 'scChIPseq_alignment_plot',
            'title': 'scChIPseq: Alignment Scores',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }

        return bargraph.plot(self.scChIPseq_data, keys, pconfig)

    def scChIPseq_barcode_chart (self):
        """ Make the plot showing alignment rates """

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys['match_barcode'] = {'color': '#00bf00', 'name': 'Barcoded'}
        keys['index_1_2_not_3'] =      { 'color': '#b5d30c', 'name': 'Index 1 and 2 found, not 3'}
        keys['index_1_not_2_not_3'] =          { 'color': '#edb900', 'name': 'Index 1 found, not 2 and 3' }
        keys['index_2_not_1_3'] =  { 'color': '#8922ff', 'name': 'Index 2 found, not 1 and 3' }
        keys['index_3_not_1_2'] =  { 'color': '#fb21ff', 'name': 'Index 3 found, not 1 and 2' }
        keys['no_index_found'] =    { 'color': '#ff2c20', 'name': 'No Index Found ~ genomic DNA' }

        # Config for the plot
        pconfig = {
            'id': 'scChIPseq_barcode_plot',
            'title': 'scChIPseq: Barcode Mapping',
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        return bargraph.plot(self.scChIPseq_data, keys, pconfig)

    def scChIPseq_flagged_coverage_chart (self):
        """ Make the plot showing alignment rates """

        for keys in self.scChIPseq_flagged_count.keys():
            flagged =  pd.Series(self.scChIPseq_flagged_count[keys]['count']).value_counts()
            flagged = pd.DataFrame(data=[flagged.values.tolist(), flagged.keys().to_list()])
            flagged=flagged.transpose()
            flagged.columns = ['Barcodes_Number', 'Reads_per_barcode']
            flagged = flagged[flagged.Reads_per_barcode >= 500]

            max_bins=math.ceil(flagged['Reads_per_barcode'].quantile(0.95))
            step= math.ceil((max_bins-500)/40)
            bins = list(range(500,max_bins,step))
            flagged_dict = dict()
            for index, row in flagged.iterrows():
                for i in bins:
                    if row['Reads_per_barcode'] >= i and row['Reads_per_barcode'] < (i + step):
                        if i not in flagged_dict:
                            flagged_dict[i] = int(row['Barcodes_Number'])
                        else:
                            flagged_dict[i] = flagged_dict[i] + int(row['Barcodes_Number'])

                if row['Reads_per_barcode'] >= (max_bins + step):
                    if (max_bins + step) not in flagged_dict:
                        flagged_dict[max_bins + step] = int(row['Barcodes_Number'])
                    else:
                        flagged_dict[max_bins + step] = flagged_dict[max_bins + step] + int(row['Barcodes_Number'])
            data = dict()
            data[keys] = flagged_dict
            data_color = dict()
            data_color[keys] = "#15a594"
            log.info(data)

        #log.info(dict(list(data.items())[0:2]))
        # Config for the plot
        pconfig = {
            'id': 'scChIPseq_flagged_coverage_plot',
            'title': "Read distribution across barcodes before duplicate removals",
            'ylab': '# Barcodes',
            'xlab': '# Reads per barcode',
            'cpswitch_counts_label': 'Number of Reads',
            'colors': data_color,
            'smooth_points': 100,  # Supply a number to limit number of points / smooth data
            'smooth_points_sumcounts': True,
        }
        desc = "**Number of barcodes with more than 500 reads: **" + str(sum(flagged[flagged['Reads_per_barcode']>=500].Barcodes_Number)) +"<br>"+ "**Number of barcodes with more than 1000 reads: **" + str(sum(flagged[flagged['Reads_per_barcode']>=1000].Barcodes_Number)) + "<br>"+ "**Number of barcodes with more than 1500 reads: **" + str(sum(flagged[flagged['Reads_per_barcode']>=1500].Barcodes_Number))
        self.add_section(
            name='Read distribution across barcodes before duplicate removal',
            anchor='scChIPseq_coverage_flagged',
            description=desc,
            plot=linegraph.plot(data, pconfig)
        )

    def scChIPseq_flagged_PCR_coverage_chart (self):
        """ Make the plot showing alignment rates """

        for keys in self.scChIPseq_flagged_PCR_count.keys():
            flagged_PCR = pd.Series(self.scChIPseq_flagged_PCR_count[keys]['count']).value_counts()
            flagged_PCR = pd.DataFrame(data=[flagged_PCR.values.tolist(), flagged_PCR.keys().to_list()])
            flagged_PCR = flagged_PCR.transpose()
            flagged_PCR.columns = ['Barcodes_Number', 'Reads_per_barcode']
            flagged_PCR = flagged_PCR[flagged_PCR.Reads_per_barcode >= 500]

            max_bins=math.ceil(flagged_PCR['Reads_per_barcode'].quantile(0.95))
            step= math.ceil((max_bins-500)/40)
            bins = list(range(500,max_bins,step))
            flagged_PCR_dict = dict()
            for index, row in flagged_PCR.iterrows():
                for i in bins:
                    if row['Reads_per_barcode'] >= i and row['Reads_per_barcode'] < (i + step):
                        if i not in flagged_PCR_dict:
                            flagged_PCR_dict[i] = int(row['Barcodes_Number'])
                        else:
                            flagged_PCR_dict[i] = flagged_PCR_dict[i] + int(row['Barcodes_Number'])

                if row['Reads_per_barcode'] >= (max_bins + step):
                    if (max_bins + step) not in flagged_PCR_dict:
                        flagged_PCR_dict[max_bins + step] = int(row['Barcodes_Number'])
                    else:
                        flagged_PCR_dict[max_bins + step] = flagged_PCR_dict[max_bins + step] + int(row['Barcodes_Number'])
            data = dict()
            data[keys] = flagged_PCR_dict
            data_color=dict()
            data_color[keys]="#4914e8"
            log.info(data)

        #log.info(dict(list(data.items())[0:2]))
        # Config for the plot
        pconfig = {
            'id': 'scChIPseq_flagged_PCR_coverage_plot',
            'title': "Barcodes distribution across reads after PCR duplicate removals\nNumber of barcodes: "+str(sum(flagged_PCR['Barcodes_Number'])),
            'ylab': '# Barcodes',
            'xlab': '# Reads per barcode',
            'cpswitch_counts_label': 'Number of Reads',
            'colors': data_color,
            'smooth_points': 100,  # Supply a number to limit number of points / smooth data
            'smooth_points_sumcounts': True,
        }
        desc = "**Number of barcodes with more than 500 reads: **" + str(sum(flagged_PCR[flagged_PCR['Reads_per_barcode']>=500].Barcodes_Number)) +"<br>"+ "**Number of barcodes with more than 1000 reads: **" + str(sum(flagged_PCR[flagged_PCR['Reads_per_barcode']>=1000].Barcodes_Number)) + "<br>"+ "**Number of barcodes with more than 1500 reads: **" + str(sum(flagged_PCR[flagged_PCR['Reads_per_barcode']>=1500].Barcodes_Number))
        self.add_section(
            name='Barcodes distribution across reads after PCR duplicate removal',
            anchor='scChIPseq_coverage_flagged_PCR',
            description=desc,
            plot=linegraph.plot(data, pconfig)
        )

    def scChIPseq_flagged_PCR_RT_coverage_chart (self):
        """ Make the plot showing alignment rates """

        for keys in self.scChIPseq_flagged_PCR_RT_count.keys():
            flagged_PCR_RT = pd.Series(self.scChIPseq_flagged_PCR_RT_count[keys]['count']).value_counts()
            flagged_PCR_RT = pd.DataFrame(data=[flagged_PCR_RT.values.tolist(), flagged_PCR_RT.keys().to_list()])
            flagged_PCR_RT = flagged_PCR_RT.transpose()
            flagged_PCR_RT.columns = ['Barcodes_Number', 'Reads_per_barcode']
            flagged_PCR_RT = flagged_PCR_RT[flagged_PCR_RT.Reads_per_barcode >= 500]

            max_bins=math.ceil(flagged_PCR_RT['Reads_per_barcode'].quantile(0.95))
            step= math.ceil((max_bins-500)/40)
            bins = list(range(500,max_bins,step))
            flagged_PCR_RT_dict = dict()
            for index, row in flagged_PCR_RT.iterrows():
                for i in bins:
                    if row['Reads_per_barcode'] >= i and row['Reads_per_barcode'] < (i + step):
                        if i not in flagged_PCR_RT_dict:
                            flagged_PCR_RT_dict[i] = int(row['Barcodes_Number'])
                        else:
                            flagged_PCR_RT_dict[i] = flagged_PCR_RT_dict[i] + int(row['Barcodes_Number'])

                if row['Reads_per_barcode'] >= (max_bins + step):
                    if (max_bins + step) not in flagged_PCR_RT_dict:
                        flagged_PCR_RT_dict[max_bins + step] = int(row['Barcodes_Number'])
                    else:
                        flagged_PCR_RT_dict[max_bins + step] = flagged_PCR_RT_dict[max_bins + step] + int(row['Barcodes_Number'])
            data = dict()
            data[keys]= flagged_PCR_RT_dict
            data_color = dict()
            data_color[keys] = "#0c7bd1"
            log.info(flagged_PCR_RT)
            log.info(data)

        # Config for the plot
        pconfig = {
            'id': 'scChIPseq_flagged_PCR_RT_coverage_plot',
            'title': "Barcodes distribution across reads after PCR and RT duplicate removals",
            'ylab': '# Barcodes',
            'xlab': '# Reads per barcode',
            'colors': data_color,
            'cpswitch_counts_label': 'Number of Reads',

        }
        desc = "**Number of barcodes with more than 500 reads: **" + str(sum(flagged_PCR_RT[flagged_PCR_RT['Reads_per_barcode']>=500].Barcodes_Number)) +"<br>"+ "**Number of barcodes with more than 1000 reads: **" + str(sum(flagged_PCR_RT[flagged_PCR_RT['Reads_per_barcode']>=1000].Barcodes_Number)) + "<br>"+ "**Number of barcodes with more than 1500 reads: **" + str(sum(flagged_PCR_RT[flagged_PCR_RT['Reads_per_barcode']>=1500].Barcodes_Number))

        self.add_section(
            name='Barcodes distribution across reads after PCR and RT duplicate removal',
            anchor='scChIPseq_coverage_flagged_PCR_RT',
            description=desc,
            plot=linegraph.plot(data, pconfig)
        )
    def scChIPseq_flagged_PCR_RT_rmDup_coverage_chart (self):
        """ Make the plot showing alignment rates """

        for keys in self.scChIPseq_flagged_PCR_RT_rmDup_count.keys():
            flagged_PCR_RT_rmDup = pd.Series(self.scChIPseq_flagged_PCR_RT_rmDup_count[keys]['count']).value_counts()
            flagged_PCR_RT_rmDup = pd.DataFrame(data=[flagged_PCR_RT_rmDup.values.tolist(), flagged_PCR_RT_rmDup.keys().to_list()])
            flagged_PCR_RT_rmDup = flagged_PCR_RT_rmDup.transpose()
            flagged_PCR_RT_rmDup.columns = ['Barcodes_Number', 'Reads_per_barcode']
            flagged_PCR_RT_rmDup = flagged_PCR_RT_rmDup[flagged_PCR_RT_rmDup.Reads_per_barcode >= 500]

            max_bins=math.ceil(flagged_PCR_RT_rmDup['Reads_per_barcode'].quantile(0.95))
            step= math.ceil((max_bins-500)/40)
            bins = list(range(500,max_bins,step))
            flagged_PCR_RT_rmDup_dict = dict()
            for index, row in flagged_PCR_RT_rmDup.iterrows():
                for i in bins:
                    if row['Reads_per_barcode'] >= i and row['Reads_per_barcode'] < (i + step):
                        if i not in flagged_PCR_RT_rmDup_dict:
                            flagged_PCR_RT_rmDup_dict[i] = int(row['Barcodes_Number'])
                        else:
                            flagged_PCR_RT_rmDup_dict[i] = flagged_PCR_RT_rmDup_dict[i] + int(row['Barcodes_Number'])

                if row['Reads_per_barcode'] >= (max_bins + step):
                    if (max_bins + step) not in flagged_PCR_RT_rmDup_dict:
                        flagged_PCR_RT_rmDup_dict[max_bins + step] = int(row['Barcodes_Number'])
                    else:
                        flagged_PCR_RT_rmDup_dict[max_bins + step] = flagged_PCR_RT_rmDup_dict[max_bins + step] + int(row['Barcodes_Number'])
            data = dict()
            data[keys]= flagged_PCR_RT_rmDup_dict
            data_color = dict()
            data_color[keys] = "#00bf00"
            log.info(flagged_PCR_RT_rmDup)
            log.info(data)

        # Config for the plot
        pconfig = {
            'id': 'scChIPseq_flagged_PCR_RT_rmDup_coverage_plot',
            'title': "Barcodes distribution across reads after PCR, RT and window-based duplicate removals",
            'ylab': '# Barcodes',
            'xlab': '# Reads per barcode',
            'colors': data_color,
            'cpswitch_counts_label': 'Number of Reads',

        }
        desc = "**Number of barcodes with more than 500 reads: **" + str(sum(flagged_PCR_RT_rmDup[flagged_PCR_RT_rmDup['Reads_per_barcode']>=500].Barcodes_Number)) +"<br>"+ "**Number of barcodes with more than 1000 reads: **" + str(sum(flagged_PCR_RT_rmDup[flagged_PCR_RT_rmDup['Reads_per_barcode']>=1000].Barcodes_Number)) + "<br>"+ "**Number of barcodes with more than 1500 reads: **" + str(sum(flagged_PCR_RT_rmDup[flagged_PCR_RT_rmDup['Reads_per_barcode']>=1500].Barcodes_Number))

        self.add_section(
            name='Barcode distribution across reads after PCR, RT and window-based duplicate removal',
            anchor='scChIPseq_coverage_flagged_PCR_RT_rmDup',
            description=desc,
            plot=linegraph.plot(data, pconfig)
        )
        
    def scChIPseq_count_matrix_region_coverage_chart (self):
        """ Make the plot showing alignment rates """

        for keys in self.scChIPseq_count_matrix.keys():
            count_matrix = self.scChIPseq_count_matrix[keys]
            log.info(count_matrix.head())
            count_matrix[count_matrix>=1]=1
            colsum = count_matrix.sum(axis=1)

            #Distribution by region
            max_bins=math.ceil(colsum.quantile(0.999))
            step= math.ceil((max_bins)/40)
            bins = list(range(0,max_bins,step))
            regions_dict = dict()
            for row in colsum:
                for i in bins:
                    if row >= i and row < (i + step):
                        if i not in regions_dict:
                            regions_dict[i] = int(row)
                        else:
                            regions_dict[i] = regions_dict[i] + int(row)

                if row >= (max_bins + step):
                    if (max_bins + step) not in regions_dict:
                        regions_dict[max_bins + step] = int(row)
                    else:
                        regions_dict[max_bins + step] = regions_dict[max_bins + step] + int(row)
            data = dict()
            data[keys]= regions_dict
            data_color = dict()
            data_color[keys] = "#00bf00"
            log.info(regions_dict)
            log.info(data)

        # Config for the plot
        pconfig = {
            'id': 'scChIPseq_count_matrix_region_coverage_plot',
            'title': "Read distribution across regions",
            'ylab': '# Reads',
            'xlab': '# Reads per region',
            'cpswitch_counts_label': 'Number of Reads',

        }
        desc = "Read distribution across regions (bins)"

        self.add_section(
            name='Read distribution across regions',
            anchor='scChIPseq_region_coverage_count_matrix',
            description=desc,
            plot=linegraph.plot(data, pconfig)
        )

    def scChIPseq_count_matrix_cell_coverage_chart (self):
        """ Make the plot showing alignment rates """

        for keys in self.scChIPseq_count_matrix.keys():
            count_matrix = self.scChIPseq_count_matrix[keys]
            log.info(count_matrix.head())
            count_matrix[count_matrix>=1]=1
            rowsum = count_matrix.sum(axis=0)

            #Distribution by region
            max_bins=math.ceil(rowsum.quantile(0.999))
            step= math.ceil((max_bins-500)/40)
            bins = list(range(500,max_bins,step))
            cells_dict = dict()
            for row in rowsum:
                for i in bins:
                    if row >= i and row < (i + step):
                        if i not in cells_dict:
                            cells_dict[i] = int(row)
                        else:
                            cells_dict[i] = cells_dict[i] + int(row)

                if row >= (max_bins + step):
                    if (max_bins + step) not in cells_dict:
                        cells_dict[max_bins + step] = int(row)
                    else:
                        cells_dict[max_bins + step] = cells_dict[max_bins + step] + int(row)
            data = dict()
            data[keys]= cells_dict
            data_color = dict()
            data_color[keys] = "#00bf00"
            log.info(cells_dict)
            log.info(data)

        # Config for the plot
        pconfig = {
            'id': 'scChIPseq_count_matrix_cells_coverage_plot',
            'title': "Read distribution across cells",
            'ylab': '# Reads',
            'xlab': '# Reads per cell',
            'cpswitch_counts_label': 'Number of Reads',

        }
        desc = "Read distribution across cells (bins)"

        self.add_section(
            name='Read distribution across cells',
            anchor='scChIPseq_cell_coverage_count_matrix',
            description=desc,
            plot=linegraph.plot(data, pconfig)
        )


