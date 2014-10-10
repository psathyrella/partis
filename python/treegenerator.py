#!/usr/bin/env python
""" Read the inferred tree parameters from Connor's json files, and generate a bunch of trees to later sample from. """

import sys
import os
import json
import numpy
import math
import tempfile
from subprocess import check_call
from opener import opener
import utils

# ----------------------------------------------------------------------------------------
class Hist(object):
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
    def __init__(self, args):
        self.args = args
        self.tree_generator = 'TreeSim'  # ape
        self.n_trees = '1000'
        self.branch_length_fname = self.args.hackey_extra_data_dir + '/branch-lengths.txt'
        self.treefname = self.args.hackey_extra_data_dir + '/trees.tre'
        self.make_branch_length_hists()
        self.generate_trees()

    #----------------------------------------------------------------------------------------
    def make_branch_length_hists(self):
        for human in utils.humans:
            with opener('r')(self.args.tree_parameter_file) as infile:
            # with opener('r')(data_dir + '/' + human + '/' + naivety + '/tree-parameters.json.gz') as infile:
                js = json.load(infile)
                hist = Hist(100, 0.0, 1.3)
                il = 0
                for length in js['branchLengths']:
                    hist.fill(length)
                    il += 1
                    if il > 100:
                        break
                hist.normalize()
                hist.write(self.branch_length_fname)
    
    #----------------------------------------------------------------------------------------
    def choose_branch_length(self, bin_centers, contents):
        iprob = numpy.random.uniform(0,1)
        sum_prob = 0.0
        # print '  %f' % iprob
        for ibin in range(len(bin_centers)):
            sum_prob += contents[ibin]
            # print '    %5f %5f' % (contents[ibin], sum_prob)
            if iprob < sum_prob:
                # print '    returning %s' % bin_centers[ibin]
                return bin_centers[ibin]
                
        assert False  # shouldn't fall through to here
    
    # ----------------------------------------------------------------------------------------
    def generate_trees(self):
        # read branch length histogram
        bin_centers, contents = [], []
        with opener('r')(self.branch_length_fname) as branchfile:
            check_sum = 0.0
            for line in branchfile:
                bin_centers.append(line.split()[0])  # keep it as a string
                contents.append(float(line.split()[1]))  # probability of each bin
                check_sum += contents[-1]
            assert utils.is_normed(check_sum)
    
        if os.path.exists(self.treefname):
            os.remove(self.treefname)

        # build up the R command line
        r_command = 'R --slave'
        if self.tree_generator == 'ape':
            r_command += ' -e \"'
            r_command += 'library(ape); '
            r_command += 'set.seed(0); '
            r_command += 'for (itree in 1:' + self.n_trees+ ') { write.tree(rtree(' + str(n_leaves) + ', br = rexp, rate = 10), \'test.tre\', append=TRUE) }'
            r_command += '\"'
            check_call(r_command, shell=True)
        elif self.tree_generator == 'TreeSim':
            # set parameters. TODO gee, I'm not really sure these parameters are all right
            n_species = '5'
            speciation_rate = '1'
            extinction_rate = '0.5'
            n_trees_each_run = '1'
            # build command line, one (painful) tree at a time
            with tempfile.NamedTemporaryFile() as commandfile:
                commandfile.write('library(TreeSim)\n')
                for it in range(int(self.n_trees)):
                    age = self.choose_branch_length(bin_centers, contents)
                    commandfile.write('trees <- sim.bd.taxa.age(' + n_species + ', ' + n_trees_each_run + ', ' + speciation_rate + ', ' + extinction_rate + ', frac = 1, ' + age + ', ' + 'mrca = FALSE)\n')
                    commandfile.write('write.tree(trees[[1]], \"' + self.treefname + '\", append=TRUE)\n')
                r_command += ' -f ' + commandfile.name
                commandfile.flush()
                check_call(r_command, shell=True)
        else:
            assert False
    
           # # plot the tree to a png
           # r_command += 'png(file=\\"foo.png\\"); '
           # r_command += 'plot.phylo(tree); '
           # r_command += 'add.scale.bar(); '
           # r_command += 'dev.off(); '
