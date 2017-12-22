import sys
import os
import re
import math
import collections
import yaml
from scipy.stats import norm
import csv
import time

import utils
import glutils
import paramutils
from hist import Hist

# ----------------------------------------------------------------------------------------
def get_bin_list(values, bin_type):
    assert bin_type == 'all' or bin_type == 'empty' or bin_type == 'full'
    lists = {}
    lists['all'] = []
    lists['empty'] = []
    lists['full'] = []
    for bin_val, bin_contents in values.iteritems():
        lists['all'].append(bin_val)
        if bin_contents < utils.eps:
            lists['empty'].append(bin_val)
        else:
            lists['full'].append(bin_val)

    return sorted(lists[bin_type])

# ----------------------------------------------------------------------------------------
def find_full_bin(bin_val, full_bins, side):
    """
    Find the member of <full_bins> which is closest to <bin_val> on the <side>.
    NOTE if it can't find it, i.e. if <bin_val> is equal to or outside the limits of <full_bins>, returns the outermost value of <full_bins>
    """
    assert full_bins == sorted(full_bins)
    assert len(full_bins) > 0
    assert side == 'lower' or side == 'upper'

    if side == 'lower':
        nearest_bin = full_bins[0]
        for ib in full_bins:
            if ib < bin_val and ib > nearest_bin:
                nearest_bin = ib
    elif side == 'upper':
        nearest_bin = full_bins[-1]
        for ib in sorted(full_bins, reverse=True):
            if ib > bin_val and ib < nearest_bin:
                nearest_bin = ib

    return nearest_bin

# ----------------------------------------------------------------------------------------
def add_empty_bins(values):
    # add an empty bin between any full ones
    all_bins = get_bin_list(values, 'all')
    for ib in range(all_bins[0], all_bins[-1]):
        if ib not in values:
            values[ib] = 0.0

# ----------------------------------------------------------------------------------------
def interpolate_bins(values, n_max_to_interpolate, bin_eps, debug=False, max_bin=-1):
    """
    Interpolate the empty (less than utils.eps) bins in <values> if the neighboring full bins have fewer than <n_max_to_interpolate> entries.
    Otherwise, fill with <bin_eps>.
    NOTE there's some shenanigans if you have empty bins on the edges
    <max_bin> specifies not to add any bins after <max_bin>
    """
    if debug:
        print '---- interpolating with %d' % n_max_to_interpolate
        for x in sorted(values.keys()):
            print '    %3d %f' % (x, values[x])
    add_empty_bins(values)
    full_bins = get_bin_list(values, 'full')
    if debug:
        print '----'
        for x in sorted(values.keys()):
            print '     %3d %f' % (x, values[x])
    for empty_bin in get_bin_list(values, 'empty'):
        lower_full_bin = find_full_bin(empty_bin, full_bins, side='lower')
        upper_full_bin = find_full_bin(empty_bin, full_bins, side='upper')
        if n_max_to_interpolate == -1 or values[lower_full_bin] + values[upper_full_bin] < n_max_to_interpolate:
            lower_weight = 1.0 / max(1, abs(empty_bin - lower_full_bin))
            upper_weight = 1.0 / max(1, abs(empty_bin - upper_full_bin))
            values[empty_bin] = lower_weight*values[lower_full_bin] + upper_weight*values[upper_full_bin]
            values[empty_bin] /= lower_weight + upper_weight
        else:
            values[empty_bin] = math.sqrt(values[lower_full_bin] + values[upper_full_bin])
    if debug:
        print '----'
        for x in sorted(values.keys()):
            print '     %3d %f' % (x, values[x])

    if values[full_bins[-1]] < n_max_to_interpolate:  # if the last full bin doesn't have enough entries, we add on a linearly-decreasing tail with slope such that it hits zero the same distance out as the last full bin is from zero
        slope = - float(values[full_bins[-1]]) / full_bins[-1]
        new_bin_val = values[full_bins[-1]]
        for new_bin in range(full_bins[-1] + 1, 2*full_bins[-1] + 1):
            new_bin_val += slope
            if new_bin_val <= 0.0 or new_bin >= max_bin:
                break
            values[new_bin] = new_bin_val

    if debug:
        print '----'
        for x in sorted(values.keys()):
            print '     %3d %f' % (x, values[x])

# ----------------------------------------------------------------------------------------
class Track(object):
    def __init__(self, name, letters):
        self.name = name
        self.letters = list(letters)  # should be a list (make a copy of it so we can modify it)
    def getdict(self):
        return {self.name : self.letters}

# ----------------------------------------------------------------------------------------
class State(object):
    def __init__(self, name):
        self.name = name
        self.transitions = {}
        self.emissions = None
        self.extras = {}  # any extra info you want to add

    def add_emission(self, track, emission_probs):
        if self.emissions == None:
            self.emissions = {}
        for letter in track.letters:
            assert letter in emission_probs
        assert 'track' not in self.emissions
        assert 'probs' not in self.emissions
        self.emissions['track'] = track.name
        self.emissions['probs'] = emission_probs

    def add_transition(self, to_name, prob):
        assert to_name not in self.transitions
        self.transitions[to_name] = prob

    def check(self):
        total = 0.0
        for _, prob in self.transitions.iteritems():
            assert prob >= 0.0
            total += prob
        if not utils.is_normed(total):
            raise Exception('transition probs not normed in %s: %s' % (self.name, self.transitions))

        if self.name == 'init':  # no emissions for 'init' state
            return

        if self.emissions is not None:
            total = 0.0
            for _, prob in self.emissions['probs'].iteritems():
                assert prob >= 0.0
                total += prob
            assert utils.is_normed(total)

# ----------------------------------------------------------------------------------------
class HMM(object):
    def __init__(self, name, tracks):
        self.name = name
        self.tracks = tracks
        self.states = []
        self.extras = {}  # any extra info you want to add
    def add_state(self, state):
        state.check()
        self.states.append(state)

# ----------------------------------------------------------------------------------------
class HmmWriter(object):
    def __init__(self, base_indir, outdir, gene_name, glfo, args, debug=False):
        self.region = utils.get_region(gene_name)
        self.raw_name = gene_name  # i.e. unsanitized
        self.germline_seqs = glfo['seqs']  # all germline alleles
        self.germline_seq = self.germline_seqs[self.region][gene_name]  # germline sequence for this hmm
        self.indir = base_indir
        self.args = args
        self.debug = debug
        self.codon_positions = {r : glfo[c + '-positions'] for r, c in utils.conserved_codons[args.locus].items()}

        print self.raw_name  # xXX
        # parameters with values that I more or less made up
        self.precision = '16'  # number of digits after the decimal for probabilities
        self.eps = 1e-6  # NOTE I also have an eps defined in utils, and they should in principle be combined
        self.n_max_to_interpolate = args.min_observations_to_write
        self.min_mean_unphysical_insertion_length = {'fv' : 1.5, 'jf' : 25}  # jf has to be quite a bit bigger, since besides account for the variation in J length from the tryp position to the end, it has to account for the difference in cdr3 lengths

        self.erosion_pseudocount_length = 10  # if we're closer to the end of the gene than this, make sure erosion probability isn't zero

        self.outdir = outdir
        self.smallest_entry_index = -1  # keeps track of the first state that has a chance of being entered from init -- we want to start writing (with add_internal_state) from there

        self.insertions = []
        if self.region == 'v':
            self.insertions.append('fv')
        elif self.region == 'd':
            self.insertions.append('vd')
        elif self.region == 'j':
            self.insertions.append('dj')
            self.insertions.append('jf')

        assert len(utils.ambiguous_bases) == 1 and utils.ambiguous_bases[0] == 'N'  # maybe need to update some stuff below if this changes

        if self.debug:
            print '%s' % utils.color_gene(gene_name)

        self.n_occurences = utils.read_single_gene_count(self.indir, gene_name, debug=self.debug)  # how many times did we observe this gene in data?
        replacement_genes = None
        # NOTE this never happens any more, since partitiondriver.cache_parameters() resets <args.min_observations_to_write> if it's arger than 10*(number of sequences)
        if self.n_occurences < self.args.min_observations_to_write:  # if we didn't see it enough, average over all the genes that find_replacement_genes() gives us
            if self.debug:
                print '      only saw it %d times (wanted %d), so use info from all other genes' % (self.n_occurences, self.args.min_observations_to_write)
            replacement_genes = utils.find_replacement_genes(self.indir, self.args.min_observations_to_write, gene_name, debug=self.debug)

        self.erosion_probs = self.read_erosion_info(gene_name, replacement_genes)
        self.insertion_probs, self.insertion_content_probs = self.read_insertion_info(gene_name, replacement_genes)
        self.mute_freqs = paramutils.read_mute_freqs(self.indir, this_gene=gene_name, locus=self.args.locus, approved_genes=replacement_genes)  # weighted averages over genes
        self.mute_counts = paramutils.read_mute_counts(self.indir, gene_name, self.args.locus)  # raw per-{ACGT} counts

        self.track = Track('nukes', utils.nukes)
        self.saniname = utils.sanitize_name(gene_name)
        self.hmm = HMM(self.saniname, self.track.getdict())  # pass the track as a dict rather than a Track object to keep the yaml file a bit more readable
        self.hmm.extras['gene_prob'] = max(self.eps, utils.read_overall_gene_probs(self.indir, only_gene=gene_name))  # if we really didn't see this gene at all, take pity on it and kick it an eps
        tmp_mean_freq_hist = Hist(fname=self.indir + '/all-mean-mute-freqs.csv')
        self.hmm.extras['overall_mute_freq'] = tmp_mean_freq_hist.get_mean()
        self.hmm.extras['per_gene_mute_freq'] = self.mute_freqs['unweighted_overall_mean']  # the other (weighted) one might be technically more accurate, depending on what you want, but it's probably not what anyone is expecting, so we write the unweighted one

    # ----------------------------------------------------------------------------------------
    def write(self):
        self.add_states()
        assert os.path.exists(self.outdir)
        with open(self.outdir + '/' + self.saniname + '.yaml', 'w') as outfile:
            yaml.dump(self.hmm, outfile, width=150)

    # ----------------------------------------------------------------------------------------
    def add_states(self):
        # NOTE it'd kinda make more sense for the fv and jf insertions to only have one state (rather than 4), but for the moment I'm just leaving with 4 'cause it's easier to leave them the same as the physical insertions
        self.add_init_state()
        # then left side insertions
        for insertion in self.insertions:
            if insertion == 'jf':
                continue
            if self.raw_name == glutils.dummy_d_genes[self.args.locus]:
                continue
            self.add_lefthand_insert_states(insertion)
        # then write internal states
        assert self.smallest_entry_index >= 0  # should have been set in add_region_entry_transitions
        for inuke in range(self.smallest_entry_index, len(self.germline_seq)):
            self.add_internal_state(inuke)
        # and finally right side insertions
        if self.region == 'j':
            self.add_righthand_insert_state(insertion='jf')

    # ----------------------------------------------------------------------------------------
    def add_init_state(self):
        init_state = State('init')
        lefthand_insertion = ''
        if len(self.insertions) > 0:
            lefthand_insertion = self.insertions[0]  # assumes it's impossible to have only a righthand deletion... which it oughta be
            assert 'jf' not in lefthand_insertion
        self.add_region_entry_transitions(init_state, lefthand_insertion)
        self.hmm.add_state(init_state)

    # ----------------------------------------------------------------------------------------
    def add_lefthand_insert_states(self, insertion):
        if not insertion in utils.boundaries:
            nukelist = ['N', ]
        else:
            nukelist = utils.nukes
        for nuke in nukelist:  # these states emit Ns, so we don't need a separate N insert state
            insert_state = State('insert_left_' + nuke)
            insert_state.extras['germline'] = nuke
            self.add_region_entry_transitions(insert_state, insertion)
            self.add_emissions(insert_state, germline_nuke=nuke)
            self.hmm.add_state(insert_state)

    # ----------------------------------------------------------------------------------------
    def add_internal_state(self, inuke):
        germline_nuke = self.germline_seq[inuke]

        # initialize
        state = State('%s_%d' % (self.saniname, inuke))
        state.extras['germline'] = germline_nuke

        # transitions
        exit_probability = self.get_exit_probability(inuke) # probability of ending this region here, i.e. excising the rest of the germline gene
        distance_to_end = len(self.germline_seq) - inuke - 1
        if distance_to_end > 0:  # if we're not at the end of this germline gene, add a transition to the next state
            state.add_transition('%s_%d' % (self.saniname, inuke+1), 1.0 - exit_probability)
        if exit_probability >= utils.eps or distance_to_end == 0:  # add transitions to end and righthand insertions if there's a decent chance of eroding to here, or if we're at the end of the germline sequence
            self.add_region_exit_transitions(state, exit_probability)

        # emissions
        self.add_emissions(state, inuke=inuke, germline_nuke=germline_nuke)

        self.hmm.add_state(state)

    # ----------------------------------------------------------------------------------------
    def add_righthand_insert_state(self, insertion):
        if not insertion in utils.boundaries:
            nukelist = ['N', ]
        else:
            nukelist = utils.nukes
        for nuke in nukelist:  # these states emit Ns, so we don't need a separate N insert state
            insert_state = State('insert_right_' + nuke)
            insert_state.extras['germline'] = nuke
            self.add_region_exit_transitions(insert_state, exit_probability=1.0)  # the 1.0 means we're not really exiting a region, but already exited it and are in the righthand insert states

            self.add_emissions(insert_state, germline_nuke=nuke)
            self.hmm.add_state(insert_state)

    # ----------------------------------------------------------------------------------------
    def add_pseudocounts(self, erosion_probs):
        for n_eroded in range(self.erosion_pseudocount_length):
            if n_eroded not in erosion_probs:
                erosion_probs[n_eroded] = 0
            if erosion_probs[n_eroded] == 0:
                erosion_probs[n_eroded] += 1

    # ----------------------------------------------------------------------------------------
    def read_erosion_info(self, this_gene, approved_genes=None):
        # NOTE that d erosion lengths depend on each other... but I don't think that's modellable with an hmm. At least for the moment we integrate over the other erosion
        if approved_genes is None:
            approved_genes = [this_gene, ]
        eprobs = {}
        genes_used = set()
        for erosion in utils.all_erosions:
            if erosion[0] != self.region:
                continue
            eprobs[erosion] = {}
            if this_gene == glutils.dummy_d_genes[self.args.locus]:
                eprobs[erosion][0] = 1.  # always erode zero bases
                continue
            deps = utils.column_dependencies[erosion + '_del']
            with open(self.indir + '/' + utils.get_parameter_fname(column=erosion + '_del', deps=deps), 'r') as infile:
                reader = csv.DictReader(infile)
                for line in reader:
                    # first see if we want to use this line (if <region>_gene isn't in the line, this erosion doesn't depend on gene version)
                    if self.region + '_gene' in line and line[self.region + '_gene'] not in approved_genes:  # NOTE you'll need to change this if you want it to depend on another region's genes
                        continue
                    # then skip nonsense erosions that're too long for this gene, but were ok for another
                    if int(line[erosion + '_del']) >= len(self.germline_seq):
                        continue

                    # then add in this erosion's counts
                    n_eroded = int(line[erosion + '_del'])
                    if n_eroded not in eprobs[erosion]:
                        eprobs[erosion][n_eroded] = 0.0
                    eprobs[erosion][n_eroded] += float(line['count'])

                    if self.region + '_gene' in line:
                        genes_used.add(line[self.region + '_gene'])

            if len(eprobs[erosion]) == 0:
                raise Exception('didn\'t read any %s erosion probs from %s' % (erosion, self.indir + '/' + utils.get_parameter_fname(column=erosion + '_del', deps=deps)))

            # do some smoothingy things NOTE that we normalize *after* interpolating
            if erosion in utils.real_erosions:  # for real erosions, don't interpolate if we lots of information about neighboring bins (i.e. we're pretty confident this bin should actually be zero)
                n_max = self.n_max_to_interpolate
            else:  # for fake erosions, always interpolate
                n_max = -1
            # print '   interpolate erosions'
            interpolate_bins(eprobs[erosion], n_max, bin_eps=self.eps, max_bin=len(self.germline_seq))
            self.add_pseudocounts(eprobs[erosion])

            # and finally, normalize
            total = 0.0
            for _, val in eprobs[erosion].iteritems():
                total += val

            test_total = 0.0
            for n_eroded in eprobs[erosion]:
                eprobs[erosion][n_eroded] /= total
                test_total += eprobs[erosion][n_eroded]
            assert utils.is_normed(test_total)

        if len(genes_used) > 1 and self.debug:  # if length is 1, we will have just used the actual gene
            print '    used erosion info from:', ' '.join(genes_used)

        return eprobs

    # ----------------------------------------------------------------------------------------
    def read_insertion_info(self, this_gene, approved_genes=None):
        if approved_genes is None:  # if we aren't explicitly passed a list of genes to use, we just use the gene for which we're actually writing the hmm
            approved_genes = [this_gene, ]
        iprobs, icontentprobs = {}, {}
        genes_used = set()
        for insertion in self.insertions:
            iprobs[insertion] = {}
            if this_gene == glutils.dummy_d_genes[self.args.locus]:
                iprobs[insertion][0] = 1.  # always insert zero bases
                icontentprobs[insertion] = {n : 0.25 for n in utils.nukes}
                continue
            deps = utils.column_dependencies[insertion + '_insertion']
            with open(self.indir + '/' + utils.get_parameter_fname(column=insertion + '_insertion', deps=deps), 'r') as infile:
                reader = csv.DictReader(infile)
                for line in reader:
                    # first see if we want to use this line (if <region>_gene isn't in the line, this erosion doesn't depend on gene version)
                    if self.region + '_gene' in line and line[self.region + '_gene'] not in approved_genes:  # NOTE you'll need to change this if you want it to depend on another region's genes
                        continue

                    # then add in this insertion's counts
                    n_inserted = 0
                    n_inserted = int(line[insertion + '_insertion'])
                    if n_inserted not in iprobs[insertion]:
                        iprobs[insertion][n_inserted] = 0.0
                    iprobs[insertion][n_inserted] += float(line['count'])

                    if self.region + '_gene' in line:
                        genes_used.add(line[self.region + '_gene'])

            if len(iprobs[insertion]) == 0:
                raise Exception('didn\'t read any %s insertion probs from %s' % (insertion, self.indir + '/' + utils.get_parameter_fname(column=insertion + '_insertion', deps=deps)))

            # print '   interpolate insertions'
            interpolate_bins(iprobs[insertion], self.n_max_to_interpolate, bin_eps=self.eps)  #, max_bin=len(self.germline_seq))  # NOTE that we normalize *after* this

            if 0 not in iprobs[insertion] or len(iprobs[insertion]) < 2:  # all hell breaks loose lower down if we haven't got shit in the way of information
                if self.debug:
                    print '    WARNING adding pseudocount to 1-bin in insertion probs'
                iprobs[insertion][0] = 1
                iprobs[insertion][1] = 1
                if self.debug:
                    print '      ', iprobs[insertion]

            assert 0 in iprobs[insertion] and len(iprobs[insertion]) >= 2  # all hell breaks loose lower down if we haven't got shit in the way of information

            # and finally, normalize
            total = 0.0
            for _, val in iprobs[insertion].iteritems():
                total += val
            test_total = 0.0
            for n_inserted in iprobs[insertion]:
                iprobs[insertion][n_inserted] /= total
                test_total += iprobs[insertion][n_inserted]
            assert utils.is_normed(test_total)

            if 0 not in iprobs[insertion] or iprobs[insertion][0] == 1.0:
                print 'ERROR cannot have all or none of the probability mass in the zero bin:', iprobs[insertion]
                assert False

            icontentprobs[insertion] = self.read_insertion_content(insertion)  # also read the base content of the insertions

        if len(genes_used) > 1:  # if length is 1, we will have just used the actual gene
            if self.debug:
                print '    insertions used:', ' '.join(genes_used)

        return iprobs, icontentprobs

    # ----------------------------------------------------------------------------------------
    def read_insertion_content(self, insertion):
        icontentprobs = {}  # NOTE this is only the probs for <insertion>, even though name is the same as in the previous function
        if insertion in utils.boundaries:  # i.e. if it's a real insertion
            with open(self.indir + '/' + insertion + '_insertion_content.csv', 'r') as icfile:
                reader = csv.DictReader(icfile)
                total = 0
                for line in reader:
                    icontentprobs[line[insertion + '_insertion_content']] = int(line['count'])
                    total += int(line['count'])

                if total == 0. and self.debug:
                    print '\n    WARNING zero insertion content probs read from %s, so setting to uniform distribution' % self.indir + '/' + insertion + '_insertion_content.csv'
                for nuke in utils.nukes:
                    if total == 0.:
                        icontentprobs[nuke] = 1. / len(utils.nukes)
                    else:
                        if nuke not in icontentprobs:
                            print '    %s not in insertion content probs, adding with zero' % nuke
                            icontentprobs[nuke] = 0
                        icontentprobs[nuke] /= float(total)
        else:  # just return uniform probs for effective (fv and jf) insertions
            icontentprobs = {n : 0.25 for n in utils.nukes}

        assert utils.is_normed(icontentprobs)

        return icontentprobs

    # ----------------------------------------------------------------------------------------
    def get_mean_insert_length(self, insertion, debug=False):
        if insertion in utils.boundaries:  # for real insertions, use the mean of the inferred histogram
            total, n_tot = 0.0, 0.0
            for length, prob in self.insertion_probs[insertion].iteritems():
                total += prob*length
                n_tot += prob
            if n_tot == 0.0:
                return -1
            else:
                return total / n_tot
        else:  # but for fv and jf, use the relative germline allele lengths
            lengths = []
            if insertion == 'fv':  # use the mean difference between other V gene cdr3 positions and ours (or zero, if its negative, i.e. if this V gene is really long)
                our_cpos = self.codon_positions['v'][self.raw_name]  # cpos of this hmm's gene
                if debug:
                    print '  %15s  cpos  delta_cpos' % ''
                for gene in self.germline_seqs['v']:
                    cpos = self.codon_positions['v'][gene]
                    delta_cpos = cpos - our_cpos
                    if delta_cpos < 0:  # count shorter genes as zeros in the averaging
                        lengths.append(0)
                    else:
                        lengths.append(delta_cpos)
                    if debug:
                        print '  %15s %3d %3d   %3d' % (gene, cpos, delta_cpos, lengths[-1])
            elif insertion == 'jf':
                our_tpos_to_end = len(self.germline_seq) - self.codon_positions['j'][self.raw_name]  # tpos of this hmm's gene
                if debug:
                    print '  our tpos to end: %d - %d = %d' % (len(self.germline_seq), self.codon_positions['j'][self.raw_name], len(self.germline_seq) - self.codon_positions['j'][self.raw_name])
                    print '  %15s  tpos-to-end  delta' % ''
                for gene in self.germline_seqs['j']:
                    tpos_to_end = len(self.germline_seqs['j'][gene]) - self.codon_positions['j'][gene]
                    delta_tpos_to_end = tpos_to_end - our_tpos_to_end
                    if delta_tpos_to_end < 0:  # count shorter genes as zeros in the averaging
                        lengths.append(0)
                    else:
                        lengths.append(delta_tpos_to_end)
                    if debug:
                        print '  %15s %3d %3d    %3d' % (gene, tpos_to_end, delta_tpos_to_end, lengths[-1])
            else:
                raise Exception('bad unphysical insertion %s' % insertion)

            mean_len = float(sum(lengths)) / len(lengths)
            return_val = mean_len
            if mean_len <= self.min_mean_unphysical_insertion_length[insertion]:  # if the mean is very small, return 1.5 instead, since we want to allow a little fuzz (i.e. some insertion) just in case (we invert this to get the prob of zero length insertion, so we need it to be bigger than 1.)
                return_val = self.min_mean_unphysical_insertion_length[insertion]

            if debug:
                print 'mean: %.2f  (return %.2f)' % (mean_len, return_val)

            return return_val

    # ----------------------------------------------------------------------------------------
    def get_inverse_insert_length(self, insertion):
        mean_length = self.get_mean_insert_length(insertion)  # it's kind of wasteful to call this in a few places, but it's actually really fast
        assert mean_length >= 0.0
        inverse_length = 0.0
        if mean_length > 0.0:
            inverse_length = 1.0 / mean_length
        if insertion != 'fv' and insertion != 'jf' and mean_length < 1.0:
            if self.debug:
                print '      small mean insert length %f' % mean_length

        return inverse_length

    # ----------------------------------------------------------------------------------------
    def get_insert_self_transition_prob(self, insertion):
        """ Probability of insert state transitioning to itself """
        inverse_length = self.get_inverse_insert_length(insertion)
        if inverse_length < 1.0:  # if (mean length) > 1, approximate insertion length as a geometric distribution
            return 1.0 - inverse_length  # i.e. the prob of remaining in the insert state is [1 - 1/mean_insert_length]
        else:  # while if (mean length) <=1, return the fraction of entries in bins greater than zero. NOTE this is a weird heuristic, *but* it captures the general features (it gets bigger if we have more insertions longer than zero)
            non_zero_sum = 0.0
            for length, prob in self.insertion_probs[insertion].iteritems():
                if length != 0:
                    non_zero_sum += prob
            self_transition_prob = non_zero_sum / float(non_zero_sum + self.insertion_probs[insertion][0])  # NOTE this otter be less than 1, since we only get here if the mean length is less than 1
            assert self_transition_prob >= 0.0 and self_transition_prob <= 1.0
            if insertion != 'fv' and insertion != 'jf':  # we pretty much expect this for unphysical insertions
                if self.debug:
                    print '    WARNING using insert self-transition probability hack for %s insertion p(>0) / p(0) = %f / %f = %f' % (insertion, non_zero_sum, non_zero_sum + self.insertion_probs[insertion][0], self_transition_prob)
                    print '      ', self.insertion_probs[insertion]
            return self_transition_prob

    # ----------------------------------------------------------------------------------------
    def get_zero_length_insertion_prob(self, insertion):
        if insertion in utils.boundaries:  # for real insertions, it's just the zero bin in the input inferred histogram
            return self.insertion_probs[insertion][0]
        else:  # but for fv and jf insertions, we need the insertion length to be determined by the relative lengths of the germline alleles
            mean_length = self.get_mean_insert_length(insertion)  # it's kind of wasteful to call this in a few places, but it's actually really fast
            assert mean_length > 0.0  # should have returned min_mean_unphysical_insertion_length if the mean was small
            zero_length_prob = 1. / mean_length  # return the prob of zero for a geometric distribution with the same mean
            if zero_length_prob < 0. or zero_length_prob > 1.:
                raise Exception('bad zero length insertion prob %f' % zero_length_prob)
            return zero_length_prob

    # ----------------------------------------------------------------------------------------
    def add_region_entry_transitions(self, state, insertion):
        """
        Add transitions *into* the v, d, or j regions. Called on either the 'init' state or the 'insert_left' state.
        For v, this is (mostly) the prob that the read doesn't extend all the way to the left side of the v gene.
        For d and j, this is (mostly) the prob to actually erode on the left side.
        The two <mostly>s are there because in both cases, we're starting from *approximate* smith-waterman alignments, so we need to add some fuzz in case the s-w is off.
        """
        assert 'jf' not in insertion  # need these to only be *left*-hand insertions
        assert state.name == 'init' or 'insert' in state.name

        # first add transitions to the insert state
        region_entry_prob = 0.0  # Prob to go to an internal germline state (i.e. not to an insert state)
        if state.name == 'init':
            if insertion == '':
                region_entry_prob = 1.0  # if no insert state on this side (i.e. we're on left side of v), we have no choice but to enter the region (the internal states)
            else:
                region_entry_prob = self.get_zero_length_insertion_prob(insertion)  # prob of entering the region from 'init' is the prob of a zero-length insertion
        elif 'insert' in state.name:
            region_entry_prob = 1.0 - self.get_insert_self_transition_prob(insertion)  # the 'insert_left' state has to either go to itself, or else enter the region
        else:
            assert False

        # If this is an 'init' state, we add a transition to 'insert' with probability the observed probability of a non-zero insertion
        # Whereas if this is an 'insert' state, we add a *self*-transition with probability 1/<mean observed insert length>
        # update: now, we also multiply by the insertion content prob, since we now have four insert states (and can thus no longer use this prob in the emissions)
        if insertion != '' and region_entry_prob < 1.0:
            if not insertion in utils.boundaries:
                nukelist = ['N', ]
            else:
                nukelist = utils.nukes
            for nuke in nukelist:
                content_prob = 1. if nuke == 'N' else self.insertion_content_probs[insertion][nuke]
                state.add_transition('insert_left_' + nuke, (1.0 - region_entry_prob) * content_prob)

        # then add transitions to the region's internal states
        total = 0.0
        if self.region == 'v':  # only add a transition to the zeroth internal state
            state.add_transition('%s_%d' % (self.saniname, 0), region_entry_prob)
            total += region_entry_prob
            self.smallest_entry_index = 0
        else:
            erosion = self.region + '_5p'
            for inuke in range(len(self.germline_seq)):
                erosion_length = inuke
                if erosion_length in self.erosion_probs[erosion]:
                    prob = self.erosion_probs[erosion][erosion_length]
                    total += prob * region_entry_prob
                    if region_entry_prob != 0.0:  # only add the line if there's a chance of entering the region from this state
                        state.add_transition('%s_%d' % (self.saniname, inuke), prob * region_entry_prob)
                        if self.smallest_entry_index == -1 or inuke < self.smallest_entry_index:  # tells us where we need to start adding internal states (the smallest internal state index we add is the first one that has nonzero transition probability here)
                            self.smallest_entry_index = inuke
                    else:
                        assert state.name == 'init' or self.raw_name == glutils.dummy_d_genes[self.args.locus]  # if there's *no* chance of entering the region, this better *not* be the 'insert_left' state (UPDATE: or, it can be the dummy d)

        if region_entry_prob != 0.0 and not utils.is_normed(total / region_entry_prob):
            raise Exception('normalization problem in add_region_entry_transitions():\n  region_entry_prob: %f   total / region_entry_prob: %f' % (region_entry_prob, total / region_entry_prob))

    # ----------------------------------------------------------------------------------------
    def add_region_exit_transitions(self, state, exit_probability):  # <exit_probability> is the prob of skipping the remainder of the internal states (of eroding up to here) [it's just 1.0 if <state> is an insert state]
        """ add transitions from <state> to righthand insert and end states (<state> can be internal or a righthand insert). """
        insertion = ''
        if self.region == 'j':
            insertion = 'jf'

        if 'insert' in state.name:
            insert_self_transition_prob = self.get_insert_self_transition_prob(insertion)
            end_prob = 1.0 - insert_self_transition_prob  # prob of going to end state from this insert state
        else:  # internal state
            if insertion == '':  # if no jf insertion, there's nowhere else to go but the end
                end_prob = 1.0
            else:
                end_prob = self.get_zero_length_insertion_prob(insertion)  # prob of skipping insertions altogether

        if insertion != '':  # add transition to righthand insert states with probability the observed probability of a non-zero insertion (times the exit_probability)
            if not insertion in utils.boundaries:
                nukelist = ['N', ]
            else:
                nukelist = utils.nukes
            for nuke in nukelist:
                content_prob = 1. if nuke == 'N' else self.insertion_content_probs[insertion][nuke]
                if 'insert' in state.name:
                    state.add_transition('insert_right_' + nuke, insert_self_transition_prob * exit_probability * content_prob)  # exit_probability should be 1.0 in this case
                else:  # internal state
                    state.add_transition('insert_right_' + nuke, (1.0 - end_prob) * exit_probability * content_prob)

        state.add_transition('end', end_prob * exit_probability)

    # ----------------------------------------------------------------------------------------
    def get_exit_probability(self, inuke):
        """
        Prob of exiting the chain of states for this region at <inuke>.
        In other words, what is the prob that we will erode all the bases to the right of <inuke>.
        """
        distance_to_end = len(self.germline_seq) - inuke - 1
        if distance_to_end == 0:  # last state has to exit region
            return 1.0

        if self.region == 'j':
            return 0.0

        erosion = self.region + '_3p'
        erosion_length = distance_to_end
        if erosion_length in self.erosion_probs[erosion]:
            prob = self.erosion_probs[erosion][erosion_length]
            if prob > utils.eps:
                return prob
            else:
                return 0.0
        else:
            return 0.0

    # ----------------------------------------------------------------------------------------
    def get_emission_prob(self, nuke1, is_insert=True, inuke=-1, germline_nuke='', insertion=''):
        if nuke1 not in utils.nukes + utils.ambiguous_bases:
            raise Exception('bad nuke (%s)' % nuke1)
        if is_insert:
            if germline_nuke == '' or germline_nuke == 'N':
                assert insertion == 'fv' or insertion == 'jf'
                return 1. / len(utils.nukes)
            else:
                assert germline_nuke in utils.nukes + ['N', ]
                mute_freq = self.mute_freqs['overall_mean']  # for now, just use the mean mute freq over the whole sequence for the insertion mute freq
                if nuke1 == germline_nuke:
                    return 1.0 - mute_freq  # I can think of some other ways to arrange this, but this seems ok
                else:
                    return mute_freq / 3.0
        else:
            assert inuke >= 0 and germline_nuke != ''

            if germline_nuke in utils.ambiguous_bases:
                return 1. / len(utils.nukes)
            else:
                mute_freq = self.mute_freqs.get(inuke, self.mute_freqs['overall_mean'])  # if it isn't there, that means we want to make an hmm state for a position that wasn't observed... which I think'll happen mostly with shorter read lengths
                assert mute_freq != 1.0 and mute_freq != 0.0

                if nuke1 == germline_nuke:
                    return 1.0 - mute_freq
                else:
                    return mute_freq / 3.0

        assert False  # shouldn't fall through to here

    # ----------------------------------------------------------------------------------------
    def get_ambiguos_emission_prob(self, inuke):
        if inuke in self.mute_counts:
            prob = 1. / len(utils.nukes)  # hm, on second thought maybe 0.25 is better?
            # obs = [c for c in self.mute_counts[inuke].values()]  # list of length four with the number of times we observed A, C, G, and T at this position
            # sqobs = [o*o for o in obs]
            # prob = float(sum(sqobs)) / (sum(obs))**2  # mean prob per base without Ns (i.e., the use of this value for N emission should ensure that sequences with lots of Ns on average have the same total probability as those without any Ns)
        else:
            prob = 1. / len(utils.nukes)  # NOTE it's debatable whether this is the best value. This will mostly (only?) get used for parts of the germline that we never saw in data (e.g. way on the left side of V when we have short reads). For these parts we use the mean mutation rate, which gives something like [0.94, 0.02, 0.02, 0.02]... but these are in practice nearly always going to be observed as ambiguous, so screw it

        return prob

    # ----------------------------------------------------------------------------------------
    def add_emissions(self, state, inuke=-1, germline_nuke=''):
        insertion = ''
        if 'insert' in state.name:
            assert len(self.insertions) == 1 or len(self.insertions) == 2
            if len(self.insertions) == 1:
                insertion = self.insertions[0]
            elif 'left' in state.name:
                insertion = self.insertions[0]
            elif 'right' in state.name:
                insertion = self.insertions[1]
            assert insertion != ''

        emission_probs = {}
        total = 0.0
        for base in utils.nukes:
            emission_probs[base] = self.get_emission_prob(base, is_insert=('insert' in state.name), inuke=inuke, germline_nuke=germline_nuke, insertion=insertion)
            total += emission_probs[base]
        if math.fabs(total - 1.0) >= self.eps:
            raise Exception('emission not normalized in state %s (total %f)   %s' % (state.name, total, emission_probs))
        state.add_emission(self.track, emission_probs)
        if len(utils.ambiguous_bases) > 0:  # and 'insert' not in state.name:
            state.extras['ambiguous_emission_prob'] = self.get_ambiguos_emission_prob(inuke)
            assert len(utils.ambiguous_bases) == 1
            state.extras['ambiguous_char'] = utils.ambiguous_bases[0]
            # print '%30s, %4d %f' % (state.name, inuke, self.get_ambiguos_emission_prob(inuke))
