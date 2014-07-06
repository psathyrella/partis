""" Process observed counts of recombination parameters into a form useful for sampling. """

import csv
import subprocess
import math
import os
import json
import sys
import operator

from opener import opener
from event import RecombinationEvent
import utils

# ----------------------------------------------------------------------------------------
class VersionCounter(object):
    """ Process observed counts of recombination parameters into a form useful for sampling.

    Input: csv with counts for each vdj gene choice, cdr3_length, and erosion lengths.
    Output: probabilities for each (vdj gene choice + cdr3_length combo + erosion length), eg: (IGHV3-23*04, IGHD4-17*01, IGHJ4*02_F, 42, 0, 2, 5, 2) : 3.4e-5
    """
    def __init__(self, infnames, base_outdir, min_counts=0, index_columns=utils.index_columns):
        assert len(index_columns) > 0
        self.min_counts = min_counts  # if a combo has fewer counts than this, ignore it completely
        self.version_freqs = {}  # Dict with entries like: (IGHV3-23*04, IGHD4-17*01, IGHJ4*02_F, 42) : 360, i.e. we saw this combination of genes with cdr3_length 42 a total of 360 times
        self.infnames = infnames
        self.base_outdir = base_outdir
        self.index_columns = index_columns  # By default, use all the columns that are needed to specify a recombination event.
                                            # If you just want, say, the frequency distribution of d_gene choices you want index_columns='d_gene'.
                                            # For the d right-hand erosion length, because of correlations, you'd need (d_gene, d_5p_del, d_3p_del)
        self.outfname = self.base_outdir + '/'
        for ic in self.index_columns:
            self.outfname += ic + '-'
        self.outfname = self.outfname.rstrip('-')
        self.outfname += '-probs.csv.bz2'

        self.all_seqs = {}
        if 'vd_insertion' in self.index_columns or 'dj_insertion' in self.index_columns:  # the insertions aren't listed in the input file, so we need to calculate them
            self.all_seqs = utils.read_germlines()
            reco_data_dir = '/home/dralph/Dropbox/work/recombinator'
            with opener('r')(reco_data_dir + '/data/v-meta.json') as json_file:  # get location of <begin> cysteine in each v region
                self.cyst_positions = json.load(json_file)
            with opener('r')(reco_data_dir + '/data/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in each j region (TGG)
                tryp_reader = csv.reader(csv_file)
                self.tryp_positions = {row[0]:row[1] for row in tryp_reader}  # WARNING: this doesn't filter out the header line

    # ----------------------------------------------------------------------------------------
    def count(self):
        """ Run the whole counting procedure. """
        total = self.parse_input()
        self.write_gene_choice_probs(total)
 
    # ----------------------------------------------------------------------------------------
    def parse_input(self):
        """ Convert input csv file to dicts. Also sum totals for gene freqs. """
        print '  parsing input files '
        total = 0.0
        for infname in self.infnames:  # usually only one file, but one of the patients was split into two files
            print '    %s' % infname,
            broke = False
            with opener('r')(infname) as infile:
                in_data = csv.DictReader(infile)
                # progress_count = 0
                # n_lines = 5000000 #subprocess.check_output(['bzcat ', infname, ' | wc'])
                for line in in_data:
                    # if progress_count % 100000 == 0:
                    #     print '  %3.1e / %3.0e = %5.2f  %s' % (progress_count, n_lines, float(progress_count) / n_lines, line['count'])
                    # progress_count += 1
                    if int(line['count']) < self.min_counts:
                        print '            BREAK count < %d.' % self.min_counts
                        broke = True
                        break
                    line['vd_insertion'], line['dj_insertion'] = 0,0  # TODO store the actual bases, not just the length
                    if 'vd_insertion' in self.index_columns or 'dj_insertion' in self.index_columns:  # the insertions aren't listed in the input file, so we need to calculate them
                        event = RecombinationEvent(self.all_seqs)
                        full_event_index = tuple(line[column] for column in utils.index_columns)  # tuplized index of this entire row
                        try:
                            event.set_vdj_combo(full_event_index,
                                                self.cyst_positions[full_event_index[utils.index_keys['v_gene']]]['cysteine-position'],
                                                int(self.tryp_positions[full_event_index[utils.index_keys['j_gene']]]),
                                                self.all_seqs)
                        except AssertionError:  # conserved codons screwed up. TODO fix that shit
                            print 'hrg'
                            continue
                        total_insertion_length = event.total_deletion_length + event.net_length_change  # TODO get this written in the original freq file, don't infer it
                        if total_insertion_length < 0:  # hrg. see comments in recombinator. TODO needs to be fixed, for reals this time
                            continue
                        partition_point = int(round(float(total_insertion_length) / 2))  # TODO this is kinda totally wrong, but an ok approximation for the moment
                        line['vd_insertion'] = int(round(partition_point))
                        line['dj_insertion'] = total_insertion_length - line['vd_insertion']

                    index = tuple(line[ic] for ic in self.index_columns)
                    if index not in self.version_freqs:  # duplicates are either from non-sensical 'v_5p_del' and 'j_3p_del' erosions or from the splitting of patient A's naive cells into two batches
                        self.version_freqs[index] = 0
                    self.version_freqs[index] += int(line['count'])
                    total += int(line['count'])
                    # if progress_count < 5:
                    #     print index,int(line['count']),self.version_freqs[index]

            if not broke:
                print ''
        assert total != 0.0
        return total
    
    # ----------------------------------------------------------------------------------------
    def write_gene_choice_probs(self, total):
        """ Write vdj combo freqs to file. """
        print '      writing gene choice probabilities to %s' % self.outfname
        if os.path.isfile(self.outfname):
            os.remove(self.outfname)
        elif not os.path.exists(self.base_outdir):
            os.makedirs(self.base_outdir)
        with opener('w')(self.outfname) as outfile:
            out_fieldnames = list(self.index_columns)
            out_fieldnames.append('prob')
            out_data = csv.DictWriter(outfile, out_fieldnames)
            out_data.writeheader()
            # NOTE this will in general not be sorted
            tmp_total = 0.0  # neurotic double check
            for index,count in self.version_freqs.iteritems():
                line = {}
                for ic in range(len(index)):
                    line[self.index_columns[ic]] = index[ic]
                prob = float(count) / total
                tmp_total += prob
                line['prob'] = prob
                out_data.writerow(line)
            assert math.fabs(tmp_total - 1.0) < 1e-8
