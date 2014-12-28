#!/usr/bin/env python
""" Read the inferred tree parameters from Connor's json files, and generate a bunch of trees to later sample from. """

import sys
import os
import random
import json
import numpy
import math
import tempfile
from subprocess import check_call
from opener import opener
import plotting
import utils

# ----------------------------------------------------------------------------------------
class Hist(object):  # NOTE I wrote this class when I couldn't get pyroot to work. Now it does, so...
    """ a simple histogram """
    def __init__(self, n_bins = 20, xmin = 0.0, xmax = 30.0):
        self.n_bins = n_bins
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.low_edges = []  # lower edge of each bin
        self.centers = []  # center of each bin
        self.bin_contents = []
        dx = (self.xmax - self.xmin) / self.n_bins
        for ib in range(self.n_bins + 2):  # ROOT conventions: zero is underflow and last bin is overflow
            self.low_edges.append(self.xmin + (ib-1)*dx)  # subtract one from ib so underflow bin has upper edge xmin
            self.centers.append(self.low_edges[-1] + 0.5*dx)
            self.bin_contents.append(0.0)

    def fill(self, value, weight=1.0):
        if value <= self.low_edges[0]:  # underflow
            self.bin_contents[0] += weight
        elif value > self.low_edges[self.n_bins + 1]:  # overflow
            self.bin_contents[self.n_bins + 1] += weight
        else:
            for ib in range(self.n_bins + 2):  # loop over the rest of the bins
                if value > self.low_edges[ib] and value <= self.low_edges[ib+1]:
                    self.bin_contents[ib] += weight

    def normalize(self):
        sum_value = 0.0
        for ib in range(1, self.n_bins + 1):  # don't include under/overflows in sum_value
            sum_value += self.bin_contents[ib]
        if sum_value == 0.0:
            print 'WARNING sum zero in Hist::normalize'
            return
        # make sure there's not too much stuff in the under/overflows
        if self.bin_contents[0]/sum_value > 1e-10 or self.bin_contents[self.n_bins+1]/sum_value > 1e-10:
            print 'WARNING under/overflows'
        for ib in range(1, self.n_bins + 1):
            self.bin_contents[ib] /= sum_value
        check_sum = 0.0
        for ib in range(1, self.n_bins + 1):  # check it
            check_sum += self.bin_contents[ib]
        assert math.fabs(check_sum - 1.0) < 1e-10

    def write(self, outfname):
        with opener('w')(outfname) as outfile:
            for ib in range(self.n_bins + 2):
                outfile.write('%.15e %.15e\n' % (self.centers[ib], self.bin_contents[ib]))

# ----------------------------------------------------------------------------------------
class TreeGenerator(object):
    def __init__(self, args, mute_freq_fname, seed):
        self.args = args
        self.tree_generator = 'TreeSim'  # ape
        self.read_mute_freqs(mute_freq_fname)
        assert self.args.outfname != None
        assert self.args.n_leaves > 1
        random.seed(seed)
        numpy.random.seed(seed)
        if self.args.debug:
            print 'generating %d trees from %s' % (self.args.n_trees, mute_freq_fname)
            if self.args.random_number_of_leaves:
                print '    with random number of leaves in [2, %d]' % self.args.n_leaves
            else:
                print '    with %d leaves' % self.args.n_leaves

    #----------------------------------------------------------------------------------------
    def read_mute_freqs(self, mute_freq_fname):
        # NOTE mute_freq_fname is mute freqs written as *percents* (it makes sense for plotting)
        #  so we need to convert to branch lengths
        # NOTE also I'm not doing it really correctly yet, but it's ok for now
        self.branch_lengths = {}
        self.branch_lengths['lengths'], self.branch_lengths['probs'] = [], []
        mutehist = plotting.make_hist_from_bin_entry_file(mute_freq_fname, 'mute-freqs')

        if mutehist.GetBinContent(0) > 0.0 or mutehist.GetBinContent(mutehist.GetNbinsX()+1) > 0.0:
            print 'WARNING nonzero under/overflow bins read from %s' % mute_freq_fname

        check_sum = 0.0
        for ibin in range(1, mutehist.GetNbinsX()+1):  # ignore under/overflow bins
            freq = mutehist.GetBinCenter(ibin)
            branch_length = float(freq) / 100
            prob = mutehist.GetBinContent(ibin)
            self.branch_lengths['lengths'].append(branch_length)
            self.branch_lengths['probs'].append(prob)
            check_sum += self.branch_lengths['probs'][-1]
        assert utils.is_normed(check_sum)

    #----------------------------------------------------------------------------------------
    def choose_branch_length(self):
        iprob = numpy.random.uniform(0,1)
        sum_prob = 0.0
        for ibin in range(len(self.branch_lengths['lengths'])):
            sum_prob += self.branch_lengths['probs'][ibin]
            if iprob < sum_prob:
                return self.branch_lengths['lengths'][ibin]
                
        assert False  # shouldn't fall through to here
    
    # ----------------------------------------------------------------------------------------
    def generate_trees(self, seed, outfname):
        if os.path.exists(outfname):
            os.remove(outfname)

        # build up the R command line
        r_command = 'R --slave'
        if self.tree_generator == 'ape':
            assert False  # needs updating
            # r_command += ' -e \"'
            # r_command += 'library(ape); '
            # r_command += 'set.seed(0); '
            # r_command += 'for (itree in 1:' + str(self.args.n_trees)+ ') { write.tree(rtree(' + str(self.args.n_leaves) + ', br = rexp, rate = 10), \'test.tre\', append=TRUE) }'
            # r_command += '\"'
            # check_call(r_command, shell=True)
        elif self.tree_generator == 'TreeSim':
            # set parameters. NOTE gee, I'm not really sure these parameters are all right
            speciation_rate = '1'
            extinction_rate = '0.5'
            n_trees_each_run = '1'
            # build command line, one (painful) tree at a time
            with tempfile.NamedTemporaryFile() as commandfile:
                commandfile.write('library(TreeSim)\n')
                commandfile.write('set.seed(' + str(seed)+ ')\n')
                for it in range(self.args.n_trees):
                    if self.args.random_number_of_leaves:
                        n_leaves = random.randint(2, self.args.n_leaves)  # NOTE interval is inclusive!
                    else:
                        n_leaves = self.args.n_leaves
                    age = self.choose_branch_length()
                    commandfile.write('trees <- sim.bd.taxa.age(' + str(n_leaves) + ', ' + n_trees_each_run + ', ' + speciation_rate + ', ' + extinction_rate + ', frac = 1, ' + str(age) + ', ' + 'mrca = FALSE)\n')
                    commandfile.write('write.tree(trees[[1]], \"' + outfname + '\", append=TRUE)\n')
                r_command += ' -f ' + commandfile.name
                commandfile.flush()
                check_call(r_command, shell=True)
        else:
            assert False
