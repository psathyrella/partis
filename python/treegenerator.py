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
from Bio import Phylo
from opener import opener
import plotting
import utils

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
            branch_length = float(freq)
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
    def check_tree_lengths(self, treefname, ages):
        trees = list(Phylo.parse(treefname, 'newick'))
        print 'checking branch lengths... ',
        assert len(trees) == len(ages)
        treetotal, agetotal = 0.0, 0.0
        for itree in range(len(ages)):
            # print '%7.4f  %7.4f' % (ages[itree], trees[itree].distance('t1'))
            treetotal += trees[itree].distance('t1')  # NOTE should be the same for t[0-9]... but I guess I should check at some point
            agetotal += ages[itree]

        print '  mean: %7.4f (asked for %7.4f)' % (agetotal / len(ages), treetotal / len(ages))

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
                ages = []
                for itree in range(self.args.n_trees):
                    if self.args.random_number_of_leaves:
                        n_leaves = random.randint(2, self.args.n_leaves)  # NOTE interval is inclusive!
                    else:
                        n_leaves = self.args.n_leaves
                    ages.append(self.choose_branch_length())
                    commandfile.write('trees <- sim.bd.taxa.age(' + str(n_leaves) + ', ' + n_trees_each_run + ', ' + speciation_rate + ', ' + extinction_rate + ', frac = 1, ' + str(ages[itree]) + ', ' + 'mrca = FALSE)\n')
                    commandfile.write('write.tree(trees[[1]], \"' + outfname + '\", append=TRUE)\n')
                r_command += ' -f ' + commandfile.name
                commandfile.flush()
                check_call(r_command, shell=True)
            self.check_tree_lengths(outfname, ages)
        else:
            assert False
