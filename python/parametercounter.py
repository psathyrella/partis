#!/usr/bin/env python
import os
import csv
import time
import sys

import utils
import glutils
from opener import opener
import plotting
from hist import Hist
from mutefreqer import MuteFreqer

# ----------------------------------------------------------------------------------------
class ParameterCounter(object):
    """ class to keep track of how many times we've seen each gene version, erosion length,
    insertion (length and base content), and mutation """
    def __init__(self, glfo, args):
        self.glfo = glfo
        self.args = args
        self.mfreqer = MuteFreqer(self.glfo)
        self.reco_total = 0  # total number of recombination events
        self.mute_total = 0  # total number of sequences
        self.counts = {}
        self.counts['all'] = {}
        for column in utils.column_dependencies:
            self.counts[column] = {}
        for bound in utils.boundaries:
            self.counts[bound + '_insertion_content'] = {n : 0 for n in utils.nukes}  # base content of each insertion
        self.counts['seq_content'] = {n : 0 for n in utils.nukes}

    # ----------------------------------------------------------------------------------------
    def get_index(self, info, deps):
        index = []
        for ic in deps:
            if ic[2:] == '_insertion':  # insertion length
                index.append(len(info[ic]))
            else:
                assert 'insertion' not in ic
                assert 'content' not in ic
                index.append(info[ic])
        return tuple(index)

    # ----------------------------------------------------------------------------------------
    def increment_all_params(self, info):
        self.increment_per_sequence_params(info)
        self.increment_per_family_params(info)

    # ----------------------------------------------------------------------------------------
    def increment_per_sequence_params(self, info):
        """ increment parameters that differ for each sequence within the clonal family """
        self.mute_total += 1
        self.mfreqer.increment(info)
        seq = info['seq']
        for nuke in seq:
            if nuke in utils.ambiguous_bases:
                continue
            self.counts['seq_content'][nuke] += 1

    # ----------------------------------------------------------------------------------------
    def increment_per_family_params(self, info):
        """ increment parameters that are the same for the entire clonal family """
        self.reco_total += 1

        all_index = self.get_index(info, tuple(list(utils.index_columns) + ['cdr3_length', ]))
        if all_index not in self.counts['all']:
            self.counts['all'][all_index] = 0
        self.counts['all'][all_index] += 1

        for deps in utils.column_dependency_tuples:
            column = deps[0]
            index = self.get_index(info, deps)
            if index not in self.counts[column]:
                self.counts[column][index] = 0
            self.counts[column][index] += 1

        for bound in utils.boundaries:
            for nuke in info[bound + '_insertion']:
                if nuke in utils.ambiguous_bases:
                    continue
                self.counts[bound + '_insertion_content'][nuke] += 1

    # ----------------------------------------------------------------------------------------
    def clean_plots(self, plotdir, subset_by_gene):
        self.mfreqer.clean_plots(plotdir + '/mute-freqs')
        utils.prep_dir(plotdir + '/overall')  #, multilings=('*.csv', '*.svg'))
        for column in self.counts:
            if subset_by_gene and ('_del' in column or column == 'vd_insertion' or column == 'dj_insertion'):  # option to subset deletion and (real) insertion plots by gene
                thisplotdir = plotdir + '/' + column
                utils.prep_dir(thisplotdir, wildlings=['*.csv', '*.svg'])

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, subset_by_gene=False, cyst_positions=None, tryp_positions=None, only_csv=False):
        print '  plotting parameters',
        sys.stdout.flush()
        start = time.time()

        self.clean_plots(plotdir, subset_by_gene)

        self.mfreqer.plot(plotdir + '/mute-freqs', cyst_positions, tryp_positions, only_csv=only_csv)  #, mean_freq_outfname=base_outdir + '/REGION-mean-mute-freqs.csv')  # REGION is replace by each region in the three output files

        overall_plotdir = plotdir + '/overall'

        for column in self.counts:
            if column == 'all':
                continue
            values, gene_values = {}, {}
            if len(self.counts[column]) == 0:
                raise Exception('no counts in %s' % column)
            for index, count in self.counts[column].iteritems():
                gene = None
                if subset_by_gene and ('_del' in column or column == 'vd_insertion' or column == 'dj_insertion'):  # option to subset deletion and (real) insertion plots by gene
                    if '_del' in column:
                        region = column[0]
                    else:
                        region = column[1]
                    assert region in utils.regions
                    assert 'IGH' + region.upper() in index[1]  # NOTE this is hackey, but it works find now and will fail obviously
                    gene = index[1]                            #   if I ever change the correlations to be incompatible. so screw it
                    if gene not in gene_values:
                        gene_values[gene] = {}

                column_val = index[0]
                if gene is not None:
                    if column_val not in gene_values[gene]:
                        gene_values[gene][column_val] = 0.0
                    gene_values[gene][column_val] += count
                if column_val not in values:
                    values[column_val] = 0.0
                values[column_val] += count

                try:  # figure out whether this is an integer or string (only used outside this loop when we make the plots)
                    int(column_val)
                    var_type = 'int'
                except:
                    var_type = 'string'

            if subset_by_gene and ('_del' in column or column == 'vd_insertion' or column == 'dj_insertion'):  # option to subset deletion and (real) insertion plots by gene
                thisplotdir = plotdir + '/' + column
                for gene in gene_values:
                    plotname = utils.sanitize_name(gene) + '-' + column
                    hist = plotting.make_hist_from_dict_of_counts(gene_values[gene], var_type, plotname, sort=True)
                    plotting.draw_no_root(hist, plotname=plotname, plotdir=thisplotdir, errors=True, write_csv=True, only_csv=only_csv)
                if not only_csv:
                    plotting.make_html(thisplotdir)

            plotname = column
            hist = plotting.make_hist_from_dict_of_counts(values, var_type, plotname, sort=True)
            plotting.draw_no_root(hist, plotname=plotname, plotdir=overall_plotdir, errors=True, write_csv=True, only_csv=only_csv)

        if not only_csv:
            plotting.make_html(overall_plotdir)

        print '(%.1f sec)' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def write(self, base_outdir):
        print '    writing parameters',
        sys.stdout.flush()
        start = time.time()

        utils.prep_dir(base_outdir, subdirs=('hmms', 'mute-freqs', 'germline-sets'), wildlings=('*.csv', '*.yaml', '*.fasta'))  # it's kind of hackey to specify the /hmms dir here, but as soon as we write the parameters below, the previous yamels are out of date, so it's pretty much necessary

        self.mfreqer.write(base_outdir + '/mute-freqs', mean_freq_outfname=base_outdir + '/REGION-mean-mute-freqs.csv')  # REGION is replace by each region in the three output files) 
        glutils.write_glfo(base_outdir + '/' + glutils.glfo_dir, glfo=self.glfo, debug=True)

        for column in self.counts:
            index = None
            outfname = None
            if column == 'all':
                index = tuple(list(utils.index_columns) + ['cdr3_length', ])
                outfname = base_outdir + '/' + utils.get_parameter_fname(column='all')
            elif '_content' in column:
                index = [column,]
                outfname = base_outdir + '/' + column + '.csv'
            else:
                index = [column,] + utils.column_dependencies[column]
                outfname = base_outdir + '/' + utils.get_parameter_fname(column_and_deps=index)
            if os.path.isfile(outfname):
                os.remove(outfname)
            elif not os.path.exists(base_outdir):
                os.makedirs(base_outdir)
            with opener('w')(outfname) as outfile:
                out_fieldnames = list(index)
                out_fieldnames.append('count')
                out_data = csv.DictWriter(outfile, out_fieldnames)
                out_data.writeheader()
                # NOTE this will in general not be sorted
                for key, count in self.counts[column].iteritems():
                    line = {}
                    for ic in range(len(key)):
                        line[index[ic]] = key[ic]
                    line['count'] = count
                    out_data.writerow(line)

        print '(%.1f sec)' % (time.time()-start)
