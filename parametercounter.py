#!/usr/bin/env python
from utils import utils

class ParameterCounter(object):
    """ class to keep track of how many times we've seen each gene version, erosion length,
    insertion (length and base content), and mutation """
    def __init__(self):
        self.total = 0
        self.counts = {}
        for column in utils.column_dependencies:
            self.counts[column] = {}

    # ----------------------------------------------------------------------------------------
    def increment(self, params):
        self.total += 1
        for deps in utils.column_dependency_tuples:
            if deps == utils.index_columns:
                continue
            column = deps[0]
            index = tuple(params[ic] for ic in deps)
            if index not in self.counts[column]:  # duplicates are either from non-sensical 'v_5p_del' and 'j_3p_del' erosions or from the splitting of patient A's naive cells into two batches
                self.counts[column][index] = 0
            self.counts[column][index] += 1

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
