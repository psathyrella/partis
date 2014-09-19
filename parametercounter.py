#!/usr/bin/env python
import os
import csv
import types
import sys
from subprocess import check_call

from utils import utils
from utils.opener import opener
from utils import plotting
from mutefreqer import MuteFreqer

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
        for column in utils.column_dependencies:
            self.counts[column] = {}
        self.mutefreqer = MuteFreqer(self.base_outdir, self.plotdir, germline_seqs)

    # ----------------------------------------------------------------------------------------
    def get_index(self, info, deps):
        index = []
        for ic in deps:
            if '_insertion' in ic:
                index.append(len(info[ic]))
            else:
                index.append(info[ic])
        return tuple(index)

    # ----------------------------------------------------------------------------------------
    def increment(self, info, gl_match):
        self.total += 1
        for deps in utils.column_dependency_tuples:
            if deps == utils.index_columns:  # don't need the fully marginalized one. Oh *snap* did I just use 'marginalized' correctly?
                continue
            column = deps[0]
            index = self.get_index(info, deps)
            if index not in self.counts[column]:
                self.counts[column][index] = 0
            self.counts[column][index] += 1
        self.mutefreqer.increment(info, gl_match)

    # ----------------------------------------------------------------------------------------
    def __str__(self):
        return_str = []
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
            values = {}
            var_type = 'int'
            for index, count in self.counts[column].iteritems():
                column_val = index[0]
                if type(column_val) == types.StringType:
                    var_type = 'string'
                if column_val not in values:
                    values[column_val] = 0.0
                values[column_val] += count
            hist = plotting.make_hist(values, var_type, column)
            plotting.draw(hist, var_type, plotname=column, plotdir=self.plotdir)
        check_call(['./permissify-www', self.plotdir])  # NOTE this should really permissify starting a few directories higher up
        check_call(['makeHtml', self.plotdir, '4', 'null', 'svg'])

    # ----------------------------------------------------------------------------------------
    def write_counts(self):
        if self.plotdir != '':
            self.plot()
        self.mutefreqer.write(self.base_outdir)
        for column in self.counts:
            index_columns = [column,] + utils.column_dependencies[column]
            outfname = self.base_outdir + '/' + utils.get_prob_fname_tuple(index_columns)
            if os.path.isfile(outfname):
                os.remove(outfname)
            elif not os.path.exists(self.base_outdir):
                os.makedirs(self.base_outdir)
            with opener('w')(outfname) as outfile:
                out_fieldnames = index_columns
                out_fieldnames.append('count')
                out_data = csv.DictWriter(outfile, out_fieldnames)
                out_data.writeheader()
                # NOTE this will in general not be sorted
                for index, count in self.counts[column].iteritems():
                    line = {}
                    for ic in range(len(index)):
                        line[index_columns[ic]] = index[ic]
                    line['count'] = count
                    out_data.writerow(line)
