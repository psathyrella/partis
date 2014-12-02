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
    def __init__(self, germline_seqs, base_outdir, plotdir=''):
        self.base_outdir = base_outdir
        self.plotdir = plotdir
        self.total = 0
        self.counts = {}
        utils.prep_dir(self.base_outdir, '*.csv')
        self.counts['all'] = {}
        for column in utils.column_dependencies:
            self.counts[column] = {}
        for bound in utils.boundaries:
            self.counts[bound + '_insertion_content'] = {'A':0, 'C':0, 'G':0, 'T':0}  # base content of each insertion TODO add correlation to previous base
        self.counts['seq_content'] = {'A':0, 'C':0, 'G':0, 'T':0}
        self.mutefreqer = MuteFreqer(self.base_outdir, self.plotdir, germline_seqs)

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
    def increment(self, info):
        self.total += 1

        all_index = self.get_index(info, utils.index_columns)
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
                self.counts[bound + '_insertion_content'][nuke] += 1
        for nuke in info['seq']:
            self.counts['seq_content'][nuke] += 1

        self.mutefreqer.increment(info)

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
                return_str.append('   %d / %d = %f\n' % (count, self.total, float(count) / self.total))
        return ''.join(return_str)

    # ----------------------------------------------------------------------------------------
    def plot(self):
        for column in self.counts:
            if column == 'all':
                continue
            values = {}
            for index, count in self.counts[column].iteritems():
                column_val = index[0]
                try:
                    int(column_val)
                    var_type = 'int'
                except:
                    var_type = 'string'
                if column_val not in values:
                    values[column_val] = 0.0
                values[column_val] += count
            hist = plotting.make_hist(values, var_type, column, sort=True)

            plotting.draw(hist, var_type, plotname=column, plotdir=self.plotdir, errors=('_content' in column))
        if has_root:
            check_call(['./permissify-www', self.plotdir])  # NOTE this should really permissify starting a few directories higher up
            check_call(['makeHtml', self.plotdir, '3', 'null', 'svg'])

    # ----------------------------------------------------------------------------------------
    def write_counts(self):
        if self.plotdir != '' and has_root:
            self.plot()
        print 'write mute freqs'
        mute_start = time.time()
        n_cached, n_not_cached = self.mutefreqer.write(self.base_outdir)
        print 'mute freq write time: %.3f' % (time.time() - mute_start)
        print ' %d / %d cached' % (n_cached, n_cached + n_not_cached)
        for column in self.counts:
            index = None
            outfname = None
            if column == 'all':
                index = utils.index_columns
                outfname = self.base_outdir + '/' + utils.get_parameter_fname(column='all')
            elif '_content' in column:
                index = [column,]
                outfname = self.base_outdir + '/' + column + '.csv'
            else:
                index = [column,] + utils.column_dependencies[column]
                outfname = self.base_outdir + '/' + utils.get_parameter_fname(column_and_deps=index)
            if os.path.isfile(outfname):
                os.remove(outfname)
            elif not os.path.exists(self.base_outdir):
                os.makedirs(self.base_outdir)
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
