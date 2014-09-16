#!/usr/bin/env python
import os
import csv
import sys
from utils import utils
from utils.opener import opener
from mutefreqer import MuteFreqer

class ParameterCounter(object):
    """ class to keep track of how many times we've seen each gene version, erosion length,
    insertion (length and base content), and mutation """
    def __init__(self, pdriver, base_outdir):
        self.pdriver = pdriver
        self.total = 0
        self.counts = {}
        self.base_outdir = base_outdir
        utils.prep_dir(self.base_outdir, '*.csv.bz2')
        for column in utils.column_dependencies:
            self.counts[column] = {}
        self.mutefreqer = MuteFreqer(self.base_outdir, os.getenv('www') + '/partis/' + self.pdriver.args.human + '/' + self.pdriver.args.naivety, self.pdriver.germline_seqs)

    # ----------------------------------------------------------------------------------------
    def increment(self, info, gl_match):
        self.total += 1
        for deps in utils.column_dependency_tuples:
            if deps == utils.index_columns:
                continue
            column = deps[0]
            index = tuple(info[ic] for ic in deps)
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
    def write_counts(self):
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
