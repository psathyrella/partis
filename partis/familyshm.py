from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import collections
import math
import numpy

from . import utils
from . import treeutils

# Family-level shm parameters (issue #460): how a family's mean shm depends on its size, and how much shm varies among the sequences in a family.
# Counted from multi-sequence annotations at parameter-counting time (cache-parameters --input-partition-fname, or partition --count-parameters), used by treegenerator.
# All three quantities are relative (a percentile, a relative variance, and a fraction), so they can be pooled across samples with different shm levels (see data/recombinator/scratch-parameters/).

fname = 'family-shm.yaml'
size_bins = [(1, 1), (2, 2), (3, 5), (6, 20), (21, None)]  # inclusive bounds (None: no upper bound). These are the bins in which #460 measured family means and within-family variance
n_percentile_bins = 100  # resolution of the family mean percentile hists
nuke_bytes = numpy.array([n.encode() for n in utils.nukes], dtype='S1')
split_half_block_len = 30  # the two halves of each sequence for the split-half covariance are alternating blocks of this many positions. Alternating single positions (odd/even) give a fraction 0.02-0.03 higher on real data, since hotspots and multi-base mutations span adjacent positions (measured for #460)

# ----------------------------------------------------------------------------------------
def bin_str(lo, hi):
    return str(lo) if lo == hi else ('%d+' % lo if hi is None else '%d-%d' % (lo, hi))
bin_strs = [bin_str(lo, hi) for lo, hi in size_bins]
def get_bin(n):
    for lo, hi in size_bins:
        if n >= lo and (hi is None or n <= hi):
            return bin_str(lo, hi)
    raise Exception('family size %d not in any size bin' % n)

# ----------------------------------------------------------------------------------------
class FamilyShmCounter(object):
    def __init__(self):
        self.seq_n_muted = collections.Counter()  # n mutated (same definition as mutefreqer's all-mean-n-muted.csv) -> number of sequences, for every sequence in every family
        self.fam_mean_n_muted = {b : [] for b in bin_strs}  # size bin -> list of family mean n mutated
        self.within_sums = {b : {'sum-sq-dev' : 0., 'sum-n-mean-sq' : 0.} for b in bin_strs[1:]}  # size bin -> sums over families of: sum over seqs of (mut freq - family mean)^2, and n seqs * (family mean)^2. Their ratio is the within-family relative variance
        self.split_half_sums = {'sum-cov' : 0., 'sum-sq-dev' : 0.}  # over all multi-seq families: sum over seqs of (deviation from family mean in half a) * (in half b), and of (full-seq deviation)^2

    # ----------------------------------------------------------------------------------------
    def increment(self, info):
        n_muted, freqs, halves = [], [], []
        for iseq in range(len(info['unique_ids'])):
            naive_seq, mseq = utils.subset_iseq(info, iseq, restrict_to_region='')  # same positions as mutefreqer's all-mean-n-muted.csv
            freq, nmut = utils.hamming_fraction(naive_seq, mseq, also_return_distance=True)
            n_muted.append(nmut)
            freqs.append(freq)
            if len(info['unique_ids']) > 1:  # singletons add nothing to the split-half sums
                halves.append(split_half_freqs(naive_seq, mseq))
        self.seq_n_muted.update(n_muted)
        sbin = get_bin(len(n_muted))
        self.fam_mean_n_muted[sbin].append(numpy.mean(n_muted))
        if len(freqs) < 2:
            return
        freqs = numpy.array(freqs)
        fmean = freqs.mean()
        self.within_sums[sbin]['sum-sq-dev'] += numpy.sum((freqs - fmean)**2)
        self.within_sums[sbin]['sum-n-mean-sq'] += len(freqs) * fmean**2
        avals, bvals = [numpy.array(v) for v in zip(*halves)]
        self.split_half_sums['sum-cov'] += numpy.sum((avals - avals.mean()) * (bvals - bvals.mean()))
        self.split_half_sums['sum-sq-dev'] += numpy.sum((freqs - fmean)**2)

    # ----------------------------------------------------------------------------------------
    def percentile_counts(self):  # for each size bin, a hist of each family mean's percentile within the per-sequence n-muted distribution
        n_tot = float(sum(self.seq_n_muted.values()))
        max_n = max(self.seq_n_muted)
        midcdf, below = [], 0.  # mid-cdf at each integer n mutated: fraction of sequences below n, plus half those at n
        for nval in range(max_n + 1):
            midcdf.append((below + 0.5 * self.seq_n_muted[nval]) / n_tot)
            below += self.seq_n_muted[nval]
        counts = {}
        for sbin, fmeans in self.fam_mean_n_muted.items():
            pcts = numpy.interp(fmeans, list(range(max_n + 1)), midcdf)  # linear between integers
            counts[sbin] = [int(c) for c in numpy.histogram(pcts, bins=n_percentile_bins, range=(0., 1.))[0]]
        return counts

    # ----------------------------------------------------------------------------------------
    def has_families(self):  # true if there were any multi-sequence families (i.e. parameters were counted on a partition)
        return any(len(self.fam_mean_n_muted[b]) > 0 for b in bin_strs[1:])

    # ----------------------------------------------------------------------------------------
    def to_dict(self):
        return {'size-bins' : bin_strs, 'mean-percentile-counts' : self.percentile_counts(), 'within-family' : self.within_sums, 'split-half' : self.split_half_sums}

    # ----------------------------------------------------------------------------------------
    def write(self, outdir):
        if self.has_families():
            utils.jsdump('%s/%s' % (outdir, fname), self.to_dict())

# ----------------------------------------------------------------------------------------
def split_half_freqs(naive_seq, mseq):  # mutated fraction in each of the two halves of the sequence (alternating blocks of split_half_block_len positions), ignoring ambiguous positions
    nv, ms = [numpy.frombuffer(s.encode(), dtype='S1') for s in (naive_seq, mseq)]
    ok = numpy.isin(nv, nuke_bytes) & numpy.isin(ms, nuke_bytes)
    muted = (nv != ms) & ok
    in_a = (numpy.arange(len(nv)) // split_half_block_len) % 2 == 0
    return [(muted & h).sum() / float(max(1, (ok & h).sum())) for h in (in_a, ~in_a)]

# ----------------------------------------------------------------------------------------
def merge_counts(fshm_list, weights):  # sum the counts in several FamilyShmCounter.to_dict() results (e.g. for pooling samples), each multiplied by its weight
    merged = {'size-bins' : bin_strs, 'mean-percentile-counts' : {b : [0.] * n_percentile_bins for b in bin_strs},
              'within-family' : {b : {'sum-sq-dev' : 0., 'sum-n-mean-sq' : 0.} for b in bin_strs[1:]}, 'split-half' : {'sum-cov' : 0., 'sum-sq-dev' : 0.}}
    for fshm, wgt in zip(fshm_list, weights):
        if fshm['size-bins'] != bin_strs:
            raise Exception('size bins %s don\'t match %s' % (fshm['size-bins'], bin_strs))
        for sbin in bin_strs:
            merged['mean-percentile-counts'][sbin] = [m + wgt * c for m, c in zip(merged['mean-percentile-counts'][sbin], fshm['mean-percentile-counts'][sbin])]
        for sbin in bin_strs[1:]:
            for key in merged['within-family'][sbin]:
                merged['within-family'][sbin][key] += wgt * fshm['within-family'][sbin][key]
        for key in merged['split-half']:
            merged['split-half'][key] += wgt * fshm['split-half'][key]
    return merged

# ----------------------------------------------------------------------------------------
def pooled_fname(locus):
    return '%s/data/recombinator/scratch-parameters/%s/hmm/%s' % (utils.get_partis_dir(), locus, fname)

# ----------------------------------------------------------------------------------------
def read_family_shm(parameter_dir, locus, debug=False):  # returns a FamilyShm from the parameter dir's family-shm.yaml or the pooled one for <locus>, or None (with a warning) if neither exists, e.g. for tr loci, which have no pooled one
    own_fn = '%s/%s' % (parameter_dir, fname) if parameter_dir is not None else None
    if (own_fn is None or not os.path.exists(own_fn)) and not os.path.exists(pooled_fname(locus)):
        print('  %s no %s in parameter dir %s or pooled for %s, so each family\'s mean shm won\'t depend on its size, and leaves within each family will have the same depth' % (utils.wrnstr(), fname, parameter_dir, locus))
        return None
    return FamilyShm(own_fn, pooled_fname(locus), debug=debug)

# ----------------------------------------------------------------------------------------
class FamilyShm(object):  # reads family-shm.yaml from a parameter dir (<own_fn>), falling back to the pooled one in scratch-parameters (<pooled_fn>), for the whole file if there isn't one, or for each size bin that has no families
    def __init__(self, own_fn, pooled_fn, debug=False):
        self.srcs = {}  # size bin -> file it came from
        pooled = utils.read_json_yaml(pooled_fn) if os.path.exists(pooled_fn) else None
        own = utils.read_json_yaml(own_fn) if own_fn is not None and os.path.exists(own_fn) else None
        assert own is not None or pooled is not None  # read_family_shm() checks this
        self.pct_probs, self.relvars = {}, {}
        for sbin in bin_strs:
            fshm = own if own is not None and (pooled is None or sum(own['mean-percentile-counts'][sbin]) > 0) else pooled
            self.srcs[sbin] = own_fn if fshm is own else pooled_fn
            counts = numpy.array(fshm['mean-percentile-counts'][sbin], dtype=float)
            if counts.sum() == 0:  # no families of this size in the only file there is: use a uniform percentile (the same draw as without family shm parameters), and add no within-family variance
                counts = numpy.ones(n_percentile_bins)
                self.srcs[sbin] = 'none (no families of this size in %s)' % self.srcs[sbin]
            self.pct_probs[sbin] = counts / counts.sum()
            if sbin != bin_strs[0]:
                wsums = fshm['within-family'][sbin]
                self.relvars[sbin] = wsums['sum-sq-dev'] / wsums['sum-n-mean-sq'] if wsums['sum-n-mean-sq'] > 0 and not self.srcs[sbin].startswith('none') else 0.
        fshm = own if own is not None else pooled  # the split-half fraction is one number over all family sizes
        self.split_half_fn = own_fn if fshm is own else pooled_fn
        shsums = fshm['split-half']
        self.split_half_fraction = min(1., max(0., shsums['sum-cov'] / shsums['sum-sq-dev'])) if shsums['sum-sq-dev'] > 0 else 0.
        if debug:
            print('    family shm parameters (split-half fraction %.3f from %s):' % (self.split_half_fraction, self.split_half_fn))
            for sbin in bin_strs:
                print('      %5s  %s  %s' % (sbin, '' if sbin not in self.relvars else 'within-family relative variance %.4f' % self.relvars[sbin], self.srcs[sbin]))

    # ----------------------------------------------------------------------------------------
    def choose_percentile(self, n_leaves):  # draw a family mean percentile for a family with <n_leaves> leaves
        sbin = get_bin(n_leaves)
        ibin = numpy.random.choice(n_percentile_bins, p=self.pct_probs[sbin])
        return (ibin + numpy.random.uniform()) / n_percentile_bins

    # ----------------------------------------------------------------------------------------
    def target_relvar(self, n_leaves):  # relative variance of leaf depths to put in a tree with <n_leaves> leaves: the data's within-family relative variance, times the fraction of it that isn't from random mutation counts (bppseqgen adds those)
        return self.split_half_fraction * self.relvars[get_bin(n_leaves)]

# ----------------------------------------------------------------------------------------
def value_at_percentile(pct, values, probs):  # inverse of the mid-cdf used in FamilyShmCounter.percentile_counts(), for a binned distribution with bin centers <values>
    vals, mids, below = [], [], 0.
    for val, prob in zip(values, probs):
        if prob > 0:
            vals.append(val)
            mids.append(below + 0.5 * prob)
        below += prob
    return float(numpy.interp(pct, mids, vals))

# ----------------------------------------------------------------------------------------
def add_within_family_variance(treestr, target_relvar):  # multiply each branch by an independent lognormal factor with mean 1, with variance chosen so the expected relative variance (variance / mean^2) of leaf depths is <target_relvar>, then rescale to the original mean depth
    # For an ultrametric tree and independent mean-1 factors with variance v, the expected variance over leaves of leaf depth is v * sum over branches of length^2 * f * (1 - f), where f is the fraction of leaves below the branch
    dtree = treeutils.get_dendro_tree(treestr=treestr, suppress_internal_node_taxa=True)
    n_leaves = treeutils.get_n_leaves(dtree)
    depth = treeutils.get_mean_leaf_height(tree=dtree)
    ssum = 0.
    for node in dtree.preorder_node_iter():
        if node is dtree.seed_node or node.edge_length is None:
            continue
        frac = len(node.leaf_nodes()) / float(n_leaves)
        ssum += node.edge_length**2 * frac * (1 - frac)
    if ssum <= 0 or depth <= 0 or target_relvar <= 0:
        return treestr
    sigma = math.sqrt(math.log1p(target_relvar * depth**2 / ssum))
    for node in dtree.preorder_node_iter():
        if node is not dtree.seed_node and node.edge_length is not None and node.edge_length > 0:
            node.edge_length *= numpy.random.lognormal(-sigma**2 / 2, sigma)
    treeutils.rescale_tree(depth, dtree=dtree)
    return dtree.as_string(schema='newick').strip()
