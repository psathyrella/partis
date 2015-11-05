#!/usr/bin/env python
import os
import csv
import time
import sys
from subprocess import check_call

import utils
from opener import opener
import plotting
from mutefreqer import MuteFreqer
has_root = plotting.check_root()

# ----------------------------------------------------------------------------------------
class ParameterCounter(object):
    """ class to keep track of how many times we've seen each gene version, erosion length,
    insertion (length and base content), and mutation """
    def __init__(self, germline_seqs):   #, base_outdir='', plotdir='', write_parameters=True, plot_parameters=True):
        self.reco_total = 0  # total number of recombination events
        self.mute_total = 0  # total number of sequences
        self.counts = {}
        self.counts['all'] = {}
        for column in utils.column_dependencies:
            self.counts[column] = {}
        for bound in utils.boundaries:
            self.counts[bound + '_insertion_content'] = {n : 0 for n in utils.nukes}  # base content of each insertion
        self.counts['seq_content'] = {n : 0 for n in utils.nukes}
        self.mutefreqer = MuteFreqer(germline_seqs)  #, self.base_outdir, self.plotdir, write_parameters=self.write_parameters, plot_parameters=self.plot_parameters)

    # ----------------------------------------------------------------------------------------
    def clean(self):
        """ remove all the parameter files """
        self.mutefreqer.clean()
        for column in self.counts:
            if column == 'all':
                os.remove(self.base_outdir + '/' + utils.get_parameter_fname(column='all'))
            else:
                index = [column,] + utils.column_dependencies[column]
                os.remove(self.base_outdir + '/' + utils.get_parameter_fname(column_and_deps=index))

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
    def increment_per_sequence_params(self, info):
        """ increment parameters that differ for each sequence within the clonal family """
        self.mute_total += 1
        self.mutefreqer.increment(info)
        seq = info['seq']
        # if 'indels' in info:
        #     seq = info['indels']['reversed_seq']  # TODO unhackify this
        for nuke in seq:
            if nuke in utils.ambiguous_bases:
                continue
            self.counts['seq_content'][nuke] += 1

        for region in ['v', 'j']:
            self.

    # ----------------------------------------------------------------------------------------
    def increment_per_family_params(self, info):
        """ increment parameters that are the same for the entire clonal family """
        self.reco_total += 1

        all_index = self.get_index(info, utils.index_columns)
        if all_index not in self.counts['all']:
            self.counts['all'][all_index] = 0
        self.counts['all'][all_index] += 1

        for deps in utils.column_dependency_tuples:
            if '_read_truncation' in deps[0]:  # these are sequence-level
                continue
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
    def __str__(self):
        return_str = []
        print 'hm I think I was too lazy to put \'all\' in this string'
        print '  or [vdj]_insertion_content or seq_content'
        for column in self.counts:
            return_str.append('%s\n' % column)
            return_str.append('%20s' % column)
            for dep in utils.column_dependencies[column]:
                return_str.append('%20s' % dep)
            return_str.append('\n')
            for index, count in self.counts[column].iteritems():
                for val in index:
                    return_str.append('%20s' % str(val))
                return_str.append('   %d / %d = %f\n' % (count, self.reco_total, float(count) / self.reco_total))
        return ''.join(return_str)

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, subset_by_gene=False, cyst_positions=None, tryp_positions=None):
        print '  plotting parameters'
        # start = time.time()
        utils.prep_dir(plotdir + '/plots')  #, multilings=('*.csv', '*.svg'))
        for column in self.counts:
            if column == 'all':
                continue
            values, gene_values = {}, {}
            if len(self.counts[column]) == 0:
                print 'ERROR no counts in %s' % column
                assert False
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
                utils.prep_dir(thisplotdir + '/plots', multilings=['*.csv', '*.svg'])
                for gene in gene_values:
                    plotname = utils.sanitize_name(gene) + '-' + column
                    hist = plotting.make_hist_from_dict_of_counts(gene_values[gene], var_type, plotname, sort=True)
                    plotting.draw(hist, var_type, plotname=plotname, plotdir=thisplotdir, errors=True, write_csv=True)
                check_call(['./bin/makeHtml', thisplotdir, '3', 'null', 'svg'])
                check_call(['./bin/permissify-www', thisplotdir])  # NOTE this should really permissify starting a few directories higher up

            plotname = column
            hist = plotting.make_hist_from_dict_of_counts(values, var_type, plotname, sort=True)
            plotting.draw(hist, var_type, plotname=plotname, plotdir=plotdir, errors=True, write_csv=True)

        self.mutefreqer.plot(plotdir, cyst_positions, tryp_positions)  #, mean_freq_outfname=base_outdir + '/REGION-mean-mute-freqs.csv')  # REGION is replace by each region in the three output files

        if has_root:
            check_call(['./bin/makeHtml', plotdir, '3', 'null', 'svg'])
            check_call(['./bin/permissify-www', plotdir])  # NOTE this should really permissify starting a few directories higher up

        # print '    parameter plot time: %.3f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def write(self, base_outdir):
        print '    writing parameters'
        # start = time.time()

        utils.prep_dir(base_outdir, multilings=('*.csv', '*.svg'))
        # mute_start = time.time()
        self.mutefreqer.write(base_outdir, mean_freq_outfname=base_outdir + '/REGION-mean-mute-freqs.csv')  # REGION is replace by each region in the three output files) 
        # print '      mut freq write time: %.3f' % (time.time() - mute_start)
        # print ' %d / %d cached' % (self.mutefreqer.n_cached, self.mutefreqer.n_cached + self.mutefreqer.n_not_cached)
        for column in self.counts:
            index = None
            outfname = None
            if column == 'all':
                index = utils.index_columns
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

        # print '    parameter write time: %.3f' % (time.time()-start)
