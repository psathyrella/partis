#!/usr/bin/env python
""" Read the inferred tree parameters from Connor's json files, and generate a bunch of trees to later sample from. """

import sys
import os
import re
import random
import json
import numpy
import math
from cStringIO import StringIO
import tempfile
from subprocess import check_call
from Bio import Phylo
from opener import opener
from hist import Hist
import utils


# ----------------------------------------------------------------------------------------
def get_leaf_node_depths(treestr):
    """ return mapping from leaf node names to depths """
    tree = Phylo.read(StringIO(treestr), 'newick')
    nodes = {}
    for tclade in tree.get_terminals():
        length = tree.distance(tclade.name) + tree.root.branch_length
        nodes[tclade.name] = length
    return nodes

# ----------------------------------------------------------------------------------------
def rescale_tree(treestr, factor):
    """ 
    Rescale the branch lengths in <treestr> (newick-formatted) by <factor>
    I.e. multiply each float in <treestr> by <factor>.
    """
    for match in re.findall('[0-9]*\.[0-9][0-9]*', treestr):
        treestr = treestr.replace(match, str(factor*float(match)))
    return treestr

# ----------------------------------------------------------------------------------------
class TreeGenerator(object):
    def __init__(self, args, mute_freq_dir, seed):
        self.args = args
        self.tree_generator = 'TreeSim'  # other option: ape
        self.branch_lengths = {}
        self.read_mute_freqs(mute_freq_dir)
        assert self.args.outfname is not None
        random.seed(seed)
        numpy.random.seed(seed)
        if self.args.debug:
            print 'generating %d trees from %s' % (self.args.n_trees, mute_freq_dir),
            if self.args.constant_number_of_leaves:
                print ' with %d leaves' % self.args.n_leaves
            else:
                print ' with random number of leaves with parameter %d' % self.args.n_leaves

    #----------------------------------------------------------------------------------------
    def read_mute_freqs(self, mute_freq_dir):
        # NOTE these are mute freqs, not branch lengths, but it's ok for now
        for mtype in ['all',] + utils.regions:
            infname = mute_freq_dir + '/' + mtype + '-mean-mute-freqs.csv'
            self.branch_lengths[mtype] = {}
            self.branch_lengths[mtype]['lengths'], self.branch_lengths[mtype]['probs'] = [], []
            mutehist = Hist(fname=infname)
            self.branch_lengths[mtype]['mean'] = mutehist.get_mean()

            # if mutehist.GetBinContent(0) > 0.0 or mutehist.GetBinContent(mutehist.GetNbinsX()+1) > 0.0:
            #     print 'WARNING nonzero under/overflow bins read from %s' % infname
            mutehist.normalize(include_overflows=False)  # if it was written with overflows included, it'll need to be renormalized
            check_sum = 0.0
            for ibin in range(1, mutehist.n_bins + 1):  # ignore under/overflow bins
                freq = mutehist.get_bin_centers()[ibin]
                branch_length = float(freq)
                prob = mutehist.bin_contents[ibin]
                self.branch_lengths[mtype]['lengths'].append(branch_length)
                self.branch_lengths[mtype]['probs'].append(prob)
                check_sum += self.branch_lengths[mtype]['probs'][-1]
            if not utils.is_normed(check_sum):
                raise Exception('not normalized %f' % check_sum)

        if self.args.debug:
            print '  mean branch lengths'
            for mtype in ['all',] + utils.regions:
                print '     %4s %7.3f (ratio %7.3f)' % (mtype, self.branch_lengths[mtype]['mean'], self.branch_lengths[mtype]['mean'] / self.branch_lengths['all']['mean'])

    #----------------------------------------------------------------------------------------
    def add_branch_lengths_and_things(self, treefname, lonely_leaves, ages):
        """ 
        Each tree is written with branch length the mean branch length over the whole sequence
        So we need to add the length for each region afterward, so each line looks e.g. like
        (t2:0.003751736951,t1:0.003751736951):0.001248262937;v:0.98,d:1.8,j:0.87
        """
        # first read the newick info for each tree
        with opener('r')(treefname) as treefile:
            treestrings = treefile.readlines()
        # then add the region-specific branch info
        length_list = ['%s:%f'% (region, self.branch_lengths[region]['mean'] / self.branch_lengths['all']['mean']) for region in utils.regions]
        for itree in range(len(ages)):
            if lonely_leaves[itree]:
                treestrings.insert(itree, 't1:%f;\n' % ages[itree])
            treestrings[itree] = treestrings[itree].replace(';', ';' + ','.join(length_list))
        # and finally write out the final lines
        with opener('w')(treefname) as treefile:
            for line in treestrings:
                treefile.write(line)
                
    #----------------------------------------------------------------------------------------
    def choose_mean_branch_length(self):
        """ mean for entire sequence, i.e. weighted average over v, d, and j """
        iprob = numpy.random.uniform(0,1)
        sum_prob = 0.0
        for ibin in range(len(self.branch_lengths['all']['lengths'])):
            sum_prob += self.branch_lengths['all']['probs'][ibin]
            if iprob < sum_prob:
                return self.branch_lengths['all']['lengths'][ibin]
                
        assert False  # shouldn't fall through to here
    
    # ----------------------------------------------------------------------------------------
    def check_tree_lengths(self, treefname, ages):
        treestrs = []
        with opener('r')(treefname) as treefile:
            for line in treefile:
                treestrs.append(line.split(';')[0] + ';')  # ignore the info I added after the ';'
        if self.args.debug > 1:
            print '  checking branch lengths... '
        assert len(treestrs) == len(ages)
        total_length, total_leaves = 0.0, 0
        for itree in range(len(ages)):
            if self.args.debug > 1:
                print '    asked for', ages[itree],
            for name, depth in get_leaf_node_depths(treestrs[itree]).items():
                if self.args.debug > 1:
                    print '%s:%f' % (name, depth),
                if not utils.is_normed(depth / ages[itree], this_eps=1e-6):
                    raise Exception('asked for branch length %f but got %f\n   %s' % (ages[itree], depth, treestrs[itree]))  # ratio of <age> (requested length) and <length> (length in the tree file) should be 1 within float precision
            total_length += ages[itree]
            total_leaves += len(re.findall('t', treestrs[itree]))
            if self.args.debug > 1:
                print ''
        if self.args.debug:
            print '    mean branch length %.5f' % (total_length / len(ages))
            print '    mean n leaves %.2f' % (float(total_leaves) / len(ages))

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
            # from docs:
            #   frac: each tip is included into the final tree with probability frac
            #   age: the time since origin / most recent common ancestor
            #   mrca: if FALSE, time since the origin of the process, else time since the most recent common ancestor of the sampled species.
            speciation_rate = '1'
            extinction_rate = '0.5'
            n_trees_each_run = '1'
            # build command line, one (painful) tree at a time
            with tempfile.NamedTemporaryFile() as commandfile:
                commandfile.write('require(TreeSim, quietly=TRUE)\n')
                commandfile.write('set.seed(' + str(seed)+ ')\n')
                ages, lonely_leaves = [], []  # keep track of which trees should have one leaft, so we can go back and add them later in the proper spots
                for itree in range(self.args.n_trees):
                    if self.args.constant_number_of_leaves:
                        n_leaves = self.args.n_leaves
                    else:
                        n_leaves = numpy.random.geometric(1./self.args.n_leaves)
                        # n_leaves = max(1, int(numpy.random.exponential(scale=self.args.n_leaves)))
                        # n_leaves = random.randint(2, self.args.n_leaves)  # NOTE interval is inclusive!
                    age = self.choose_mean_branch_length()
                    ages.append(age)
                    if n_leaves == 1:  # TODO doesn't work yet
                        lonely_leaves.append(True)
                        continue
                    lonely_leaves.append(False)
                    commandfile.write('trees <- sim.bd.taxa.age(' + str(n_leaves) + ', ' + n_trees_each_run + ', ' + speciation_rate + ', ' + extinction_rate + ', frac=1, age=' + str(age) + ', mrca = FALSE)\n')
                    commandfile.write('write.tree(trees[[1]], \"' + outfname + '\", append=TRUE)\n')
                r_command += ' -f ' + commandfile.name
                commandfile.flush()
                if lonely_leaves.count(True) == len(ages):
                    print '  all lonely leaves'
                    open(outfname, 'w').close()
                else:
                    check_call(r_command, shell=True)
            self.add_branch_lengths_and_things(outfname, lonely_leaves, ages)
            self.check_tree_lengths(outfname, ages)
        else:
            assert False
