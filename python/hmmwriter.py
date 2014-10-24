import sys
import os
import re
import math
import collections
import yaml
from scipy.stats import norm
import csv
import utils
from opener import opener

# ----------------------------------------------------------------------------------------
class Track(object):
    def __init__(self, name, letters):
        self.name = name
        self.letters = letters  # should be a list

# ----------------------------------------------------------------------------------------
class State(object):
    def __init__(self, name):
        self.name = name
        self.transitions = {}
        self.emissions = {}  # partially implement emission to multiple tracks (I say 'partially' because I think I haven't written it into ham yet)
        self.pair_emissions = {}
        self.extras = {}  # any extra info you want to add

    def add_emission(self, track, emission_probs):  # NOTE we only allow one single (i.e. non-pair) emission a.t.m
        for letter in track.letters:
            assert letter in emission_probs
        assert 'track' not in self.emissions
        assert 'probs' not in self.emissions
        self.emissions['track'] = track.name
        self.emissions['probs'] = emission_probs

    def add_pair_emission(self, track, pair_emission_probs):  # NOTE we only allow one pair emission a.t.m
        for letter1 in track.letters:
            assert letter1 in pair_emission_probs
            for letter2 in track.letters:
                assert letter2 in pair_emission_probs[letter1]
        assert 'tracks' not in self.pair_emissions
        assert 'probs' not in self.pair_emissions
        self.pair_emissions['tracks'] = [track.name, track.name]
        self.pair_emissions['probs'] = pair_emission_probs

    def add_transition(self, to_name, prob):
        assert to_name not in self.transitions
        self.transitions[to_name] = prob

# ----------------------------------------------------------------------------------------
class HMM(object):
    def __init__(self, name, tracks):
        self.name = name
        self.tracks = tracks
        self.states = []
        self.extras = {}  # any extra info you want to add
    def add_state(self, state):
        self.states.append(state)

# ----------------------------------------------------------------------------------------
class HmmWriter(object):
    def __init__(self, base_indir, outdir, gene_name, naivety, germline_seq):
        self.indir = base_indir
        self.precision = '16'  # number of digits after the decimal for probabilities. TODO increase this?
        self.eps = 1e-6  # TODO I also have an eps defined in utils
        self.min_occurences = 10

        self.insert_mute_prob = 0.0
        self.mean_mute_freq = 0.0

        self.outdir = outdir
        self.region = utils.get_region(gene_name)

        n_occurences = utils.read_overall_gene_probs(self.indir, only_gene=gene_name, normalize=False)  # how many times did we observe this gene in data?
        # self.gene_name = gene_name
        replacement_gene = gene_name
        if n_occurences < self.min_occurences:  # if we didn't see it enough, use <replacement_gene> for all the parameters
            replacement_gene = utils.find_replacement_gene(self.indir, gene_name, self.min_occurences)
            print '\n    only saw it %d times, use info from %s instead' % (n_occurences, replacement_gene)
        self.saniname = utils.sanitize_name(gene_name)  # TODO make this not a member variable to make absolutely sure you don't confuse gene_name and replacement_gene
        self.naivety = naivety
        self.germline_seq = germline_seq
        self.smallest_entry_index = -1  # keeps track of the first state that has a chance of being entered from init -- we want to start writing (with add_internal_state) from there
        self.insertion = ''
        if self.region == 'd':
            self.insertion = 'vd'
        elif self.region == 'j':
            self.insertion = 'dj'

        self.erosion_probs = {}
        self.read_erosion_probs(replacement_gene)  # try this exact gene, but...
        if len(self.erosion_probs) == 0:  # ...if we don't have info for it use other alleles
            print '\n    no erosion probs for', self.gene_name, 'so try other alleles'
            self.read_erosion_probs(use_other_alleles=True)
        #     if len(self.erosion_probs) == 0:  # ...if we don't have info for it use other 'primary versions'
        #         print '    couldn\'t find any alleles, either, so try other primary versions'
        #         self.read_erosion_probs(use_other_alleles=False, use_other_primary_versions=True)
        assert len(self.erosion_probs) > 0

        self.insertion_probs = {}
        if self.region != 'v':
            self.read_insertion_probs(replacement_gene)
            # if len(self.insertion_probs) == 0:  # try again using all the alleles for this gene
            #     print '\n    no insertion probs found for', self.gene_name, 'so try other alleles'
            #     self.read_insertion_probs(use_other_alleles=True)
            #     if len(self.insertion_probs) == 0:  # ...if we don't have info for it use other 'primary versions'
            #         print '    couldn\'t find any alleles, either, so try other primary versions'
            #         self.read_insertion_probs(use_other_alleles=False, use_other_primary_versions=True)
            #         if len(self.insertion_probs) == 0:  # ...totally punt and just use info for some random version
            #             assert self.region == 'j'
            #             print '     hrg, can\'t find that either. Using IGHJ4*02_F, he\'s a popular fellow'
            #             self.read_insertion_probs(use_other_alleles=False, use_other_primary_versions=False, use_specific_version='IGHJ4*02_F')
            assert len(self.insertion_probs) > 0

        self.mute_freqs = {}  # TODO make sure that the overall 'normalization' of the mute freqs here agrees with the branch lengths in the tree simulator in recombinator. I kinda think it doesn't
        if self.naivety == 'M':  # mutate if not naive
            self.read_mute_freqs(replacement_gene)

        self.track = Track('nukes', list(utils.nukes))
        self.hmm = HMM(self.saniname, {'nukes':list(utils.nukes)})  # pass the track as a dict rather than a Track object to keep the yaml file a bit more readable
        self.hmm.extras['gene_prob'] = utils.read_overall_gene_probs(self.indir, only_gene=gene_name)

    # ----------------------------------------------------------------------------------------
    def write(self):
        self.add_states()
        assert os.path.exists(self.outdir)
        with opener('w')(self.outdir + '/' + self.saniname + '.yaml') as outfile:
            yaml.dump(self.hmm, outfile, width=200)
    
    # ----------------------------------------------------------------------------------------
    def add_states(self):
        self.add_init_state()
        if self.region != 'v':  # for d and j regions add insert state on the left-hand side
            self.add_insert_state()
        # then write internal states
        assert self.smallest_entry_index >= 0  # should have been set in add_region_entry_transitions
        for inuke in range(self.smallest_entry_index, len(self.germline_seq)):
            self.add_internal_state(inuke)

    # ----------------------------------------------------------------------------------------
    def add_init_state(self):
        init_state = State('init')
        self.add_region_entry_transitions(init_state)
        self.hmm.add_state(init_state)

    # ----------------------------------------------------------------------------------------
    def add_insert_state(self):
        insert_state = State('insert')
        self.add_region_entry_transitions(insert_state)  # TODO allow d region to be entirely eroded? Um, I don't think I have any probabilities for how frequently this happens though.
        self.add_emissions(insert_state)
        self.hmm.add_state(insert_state)

    # ----------------------------------------------------------------------------------------
    def add_internal_state(self, inuke):
        # arbitrarily replace ambiguous nucleotides with A TODO figger out something better
        germline_nuke = self.germline_seq[inuke]
        if germline_nuke == 'N' or germline_nuke == 'Y':
            print '\n    WARNING replacing %s with A' % germline_nuke
            germline_nuke = 'A'

        # initialize
        state = State('%s_%d' % (self.saniname, inuke))
        state.extras['germline'] = germline_nuke

        # transitions
        exit_probability = self.get_exit_probability(inuke) # probability of ending this region here, i.e. excising the rest of the germline gene
        distance_to_end = len(self.germline_seq) - inuke - 1
        if distance_to_end > 0:  # if we're not at the end of this germline gene, add a transition to the next state
            state.add_transition('%s_%d' % (self.saniname, inuke+1), 1 - exit_probability)
        if exit_probability >= utils.eps or distance_to_end == 0:  # add transition to 'end' if there's a non-zero chance of eroding to, here or if we're at the end of the germline sequence
            state.add_transition('end', exit_probability)

        # emissions
        self.add_emissions(state, inuke=inuke, germline_nuke=germline_nuke)

        self.hmm.add_state(state)

    # ----------------------------------------------------------------------------------------
    def read_erosion_probs(self, gene_name, use_other_alleles=False, use_other_primary_versions=False):
        # NOTE reads info for <gene_name>, which is *not* necessarily the gene for which we're writing an hmm
        # TODO in cases where the bases which are eroded are the same as those inserted (i.e. cases that *suck*) I seem to *always* decide on the choice with the shorter insertion. not good!
        # TODO fill in cases where we don't have info from data with some (small) default value
        # NOTE that d erosion lengths depend on each other... but I don't think that's modellable with an hmm. At least for the moment we integrate over the other erosion
        for erosion in utils.real_erosions + utils.effective_erosions:
            if erosion[0] != self.region:
                continue
            self.erosion_probs[erosion] = {}
            deps = utils.column_dependencies[erosion + '_del']
            with opener('r')(self.indir + '/' + utils.get_parameter_fname(column=erosion + '_del', deps=deps)) as infile:
                reader = csv.DictReader(infile)
                total = 0.0
                genes_used = set()
                for line in reader:
                    # first see if we want to skip this line
                    if self.region != 'j':  # j erosion doesn't depend on gene choice, so we use everything
                        if use_other_primary_versions:
                            if not utils.are_same_primary_version(line[self.region + '_gene'], gene_name):
                                continue
                            genes_used.add(line[self.region + '_gene'])
                        elif use_other_alleles:
                            if not utils.are_alleles(line[self.region + '_gene'], gene_name):
                                continue
                            genes_used.add(line[self.region + '_gene'])
                        elif line[self.region + '_gene'] != gene_name:  # skip other genes
                            continue
                    if int(line[erosion + '_del']) >= len(self.germline_seq):  # nonsense erosions (that're too long) should only occur if we're looking at other genes
                        assert self.region == 'j' or use_other_alleles or use_other_primary_versions
                        continue

                    # then add in this erosion's counts
                    n_eroded = int(line[erosion + '_del'])
                    if n_eroded not in self.erosion_probs[erosion]:
                        self.erosion_probs[erosion][n_eroded] = 0.0
                    self.erosion_probs[erosion][n_eroded] += float(line['count'])
                    total += float(line['count'])

                if len(genes_used) > 0:
                    print '    used', ' '.join(genes_used)

                if len(self.erosion_probs[erosion]) == 0:  # if we didn't find any erosion information return without doing anything else (we'll try again, looking at other genes)
                    self.erosion_probs = {}
                    return


                # # for fake v_5p and j_3p, interpolate the bases for which we didn't find info
                # if erosion in utils.effective_erosions:
                #     for inuke in range(len(self.germline_seq)):
                #         n_eroded = inuke if '5p' in erosion else len(self.germline_seq)) - inuke - 1
                #         nearest_left_neighbor = -1
                #         ipos = inuke
                #         while ipos >= 0
                #             if 

                # and finally, normalize
                test_total = 0.0
                for n_eroded in self.erosion_probs[erosion]:
                    self.erosion_probs[erosion][n_eroded] /= total
                    test_total += self.erosion_probs[erosion][n_eroded]
                assert utils.is_normed(test_total)

    # ----------------------------------------------------------------------------------------
    def read_insertion_probs(self, gene_name, use_other_alleles=False, use_other_primary_versions=False, use_specific_version=''):
        # NOTE reads info for <gene_name>, which is *not* necessarily the gene for which we're writing an hmm
        # TODO fill in cases where we don't have info from data with some (small) default value
        self.insertion_probs[self.insertion] = {}
        deps = utils.column_dependencies[self.insertion + '_insertion']
        with opener('r')(self.indir + '/' + utils.get_parameter_fname(column=self.insertion + '_insertion', deps=deps)) as infile:
            reader = csv.DictReader(infile)
            total = 0.0
            genes_used = set()
            for line in reader:
                # first see if we should skip this line
                if self.insertion == 'dj':  # for dj insertion, skip rows for genes other than the current one (vd insertion doesn't depend on gene choice, so use everything for it)
                    if use_other_primary_versions:
                        if not utils.are_same_primary_version(line[self.region + '_gene'], gene_name):
                            continue
                        genes_used.add(line[self.region + '_gene'])
                    if use_other_alleles:
                        if not utils.are_alleles(line['j_gene'], gene_name):
                            continue
                        genes_used.add(line[self.region + '_gene'])
                    elif use_specific_version != '':
                        if line['j_gene'] != use_specific_version:
                            continue
                        genes_used.add(line[self.region + '_gene'])
                    elif line['j_gene'] != gene_name:
                        continue

                # then add in this insertion's counts
                n_inserted = 0
                n_inserted = int(line[self.insertion + '_insertion'])
                if n_inserted not in self.insertion_probs[self.insertion]:
                    self.insertion_probs[self.insertion][n_inserted] = 0.0
                self.insertion_probs[self.insertion][n_inserted] += float(line['count'])
                total += float(line['count'])

            if len(genes_used) > 0:
                print '    used', ' '.join(genes_used)

            if len(self.insertion_probs[self.insertion]) == 0:  # if we didn't find any insertion information, return without doing anything else (we'll try again, looking at other genes)
                self.insertion_probs = {}
                return

            if 0 not in self.insertion_probs[self.insertion] or len(self.insertion_probs[self.insertion]) < 2:  # all hell breaks loose lower down if we haven't got shit in the way of information. May as well just use another gene
                print '\n    resetting insertion probs', self.insertion_probs[self.insertion]
                self.insertion_probs = {}
                return

            # if 0 in self.insertion_probs[self.insertion] and len(self.insertion_probs[self.insertion]) == 1:  # if we only observed zero-length insertions, add a single 1-insertion TODO shouldn't need this if you're running on lots of sequences
            #     self.insertion_probs[self.insertion][1] = 1
            #     total += 1

            # and finally, normalize
            test_total = 0.0
            for n_inserted in self.insertion_probs[self.insertion]:
                self.insertion_probs[self.insertion][n_inserted] /= total
                test_total += self.insertion_probs[self.insertion][n_inserted]
            assert utils.is_normed(test_total)
        
    # ----------------------------------------------------------------------------------------
    def read_mute_freqs(self, gene_name):
        # NOTE reads info for <gene_name>, which is *not* necessarily the gene for which we're writing an hmm
        mutefname = self.indir + '/mute-freqs/' + utils.sanitize_name(gene_name) + '.csv'
        i_try = 0
        assert os.path.exists(mutefname)
        if not os.path.exists(mutefname):
            print '\n    no mute file for',gene_name,'so try other alleles:',
        while not os.path.exists(mutefname):  # loop through other possible alleles
            # TODO double check and unify how I look for other alleles/versions when I don't have info for a particular gene version
            # TODO if I have *no* info for a position, I should really use an average rather than a very small value
            allele_id = re.findall('_star_[0-9][0-9]*', mutefname)
            assert len(allele_id) == 1  # make sure we found only one match
            allele_id = allele_id[0]
            if i_try == 0:
                allele_number = 1
            else:
                allele_number = int(allele_id.replace('_star_', ''))
                # assert allele_number > 0 and allele_number < 100
                allele_number += 1
            mutefname = mutefname.replace(allele_id, '_star_%02d' % allele_number)
            print allele_number,
            if i_try > 99:
                assert False  # er, shouldn't need this any more
                break
            i_try += 1
        if i_try > 0:
            print ''
        
        assert os.path.exists(mutefname)

        total = 0.0
        n_positions = 0
        with opener('r')(mutefname) as mutefile:
            reader = csv.DictReader(mutefile)
            for line in reader:
                freq = float(line['mute_freq'])
                lo_err = float(line['lo_err'])
                hi_err = float(line['hi_err'])
                assert freq >= 0.0 and lo_err >= 0.0 and hi_err >= 0.0  # you just can't be too careful
                if freq < utils.eps or abs(1.0 - freq) < utils.eps:  # if <freq> too close to 0 or 1, replace it with the midpoint of its uncertainty band
                    freq = 0.5 * (lo_err + hi_err)
                self.mute_freqs[int(line['position'])] = freq
                total += freq
                n_positions += 1
        assert n_positions > 0
        self.mean_mute_freq = total / n_positions
        self.insert_mute_prob = self.mean_mute_freq

    # ----------------------------------------------------------------------------------------
    def add_region_entry_transitions(self, state):
        """
        Add transitions *into* the v, d, or j regions. Called from either the 'init' state or the 'insert' state.
        For v, this is (mostly) the prob that the read doesn't extend all the way to the left side of the v gene.
        For d and j, this is (mostly) the prob to actually erode on the left side.
        The two <mostly>s are there because in both cases, we're starting from *approximate* smith-waterman alignments, so we need to add some fuzz in case the s-w is off.
        For insert states, 
        """
        assert state.name == 'init' or state.name == 'insert'

        # if self.region == 'v':
        #     # probability to enter v at a given base is given by a normal distribution with mean <mean_start_pos> and width <fuzz_around_v_left_edge>
        #     assert self.v_right_length != -1
        #     mean_start_pos = len(self.germline_seq) - self.v_right_length
        #     probs = collections.OrderedDict()
        #     total = 0.0
        #     for inuke in range(len(self.germline_seq)):  # start at far left side of v, but only write out probs that are greater than self.eps
        #         tmp_prob = float(norm.pdf(float(inuke), float(mean_start_pos), self.fuzz_around_v_left_edge))  # NOTE not yet normalized
        #         if tmp_prob < self.eps:
        #             continue
        #         probs[inuke] = tmp_prob
        #         total += tmp_prob
        #         if self.smallest_entry_index == -1 or inuke < self.smallest_entry_index:
        #             self.smallest_entry_index = inuke
        #     test_total = 0.0
        #     for inuke in probs:  # normalize and check
        #         probs[inuke] /= total
        #         test_total += probs[inuke]
        #     assert utils.is_normed(test_total)
        #     for inuke in probs:
        #         state.add_transition('%s_%d' % (self.saniname, inuke), probs[inuke])
        # else:
        non_zero_insertion_prob = 0.0  # Prob of a non-zero-length insertion (i.e. prob to *not* go directly into the region)
                                       # The sum of the region entry probs must be (1 - non_zero_insertion_prob) for d and j
                                       # (i.e. such that [prob of transitions to insert] + [prob of transitions *not* to insert] is 1.0)

        # first add transitions to the insert state
        if self.region != 'v':
            if state.name == 'init':
                if 0 not in self.insertion_probs[self.insertion] or self.insertion_probs[self.insertion][0] == 1.0:  # if there is a non-zero prob of a zero-length insertion, subtract that prob from 1.0
                    print 'arg',self.insertion_probs[self.insertion]
                    assert False
                non_zero_insertion_prob = 1.0 - self.insertion_probs[self.insertion][0]
            elif state.name == 'insert':  # we want the prob of *leaving* the insert state to be 1/insertion_length, so multiply all the region entry probs (below) by this
                mean_length = self.get_mean_insert_length()
                inverse_length = 0.0
                if mean_length > 0.0:
                    inverse_length = 1.0 / mean_length
                non_zero_insertion_prob = 1.0 - inverse_length  # set the prob of *remaining* in the insert state to [1 - 1/mean_insert_length]
            state.add_transition('insert', non_zero_insertion_prob)

        # then add transitions to the region's internal states
        total = 0.0
        for inuke in range(len(self.germline_seq)):
            erosion = self.region + '_5p'
            erosion_length = inuke
            if erosion_length in self.erosion_probs[erosion]:  # TODO this disallows erosions that we didn't see in data
                prob = self.erosion_probs[erosion][erosion_length]
                # if prob < utils.eps:
                #     continue
                total += prob * (1.0 - non_zero_insertion_prob)
                if non_zero_insertion_prob != 1.0:  # only add the line if there's a chance of zero-length insertion
                    state.add_transition('%s_%d' % (self.saniname, inuke), prob * (1.0 - non_zero_insertion_prob))
                    if self.smallest_entry_index == -1 or inuke < self.smallest_entry_index:
                        self.smallest_entry_index = inuke
        assert non_zero_insertion_prob == 1.0 or utils.is_normed(total / (1.0 - non_zero_insertion_prob))

    # ----------------------------------------------------------------------------------------
    def get_exit_probability(self, inuke):
        """
        Prob of exiting the chain of states for this region at <inuke>.
        In other words, what is the prob that we will erode all the bases to the right of <inuke>.
        """
        # TODO note that taking these numbers straight from data, with no smoothing, means that we are *forbidding* erosion lengths that we do not see in the training sample. Good? Bad? t.b.d.
        distance_to_end = len(self.germline_seq) - inuke - 1
        if distance_to_end == 0:  # last state has to go to END
            return 1.0
        # if self.region == 'j':  # always go to end of germline j region
        #     return 0.0
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
    def get_mean_insert_length(self):
        total, n_tot = 0.0, 0.0
        for length, count in self.insertion_probs[self.insertion].iteritems():
            total += count*length
            n_tot += count
        if n_tot == 0.0:
            return -1
        else:
            return total / n_tot

    # ----------------------------------------------------------------------------------------
    def get_emission_prob(self, nuke1, nuke2='', is_insert=True, inuke=-1, germline_nuke=''):
        assert nuke1 in utils.nukes
        assert nuke2 == '' or nuke2 in utils.nukes
        prob = 1.0
        if is_insert:
            assert self.insert_mute_prob != 0.0
            if nuke2 == '':  # single (non-pair) emission
                prob = 1./len(utils.nukes)
            else:
                if nuke1 == nuke2:
                    prob = (1. - self.insert_mute_prob) / 4
                else:
                    prob = self.insert_mute_prob / 12
        else:
            assert inuke >= 0
            assert germline_nuke != ''

            # first figure out the mutation frequency we're going to use
            mute_freq = self.mean_mute_freq
            if inuke in self.mute_freqs:  # if we found this base in this gene version in the data parameter file
                mute_freq = self.mute_freqs[inuke]

            # then calculate the probability
            if nuke2 == '':
                assert mute_freq != 1.0 and mute_freq != 0.0
                if nuke1 == germline_nuke:  # TODO note that if mute_freq is 1.0 this gives zero
                    prob = 1.0 - mute_freq
                else:
                    prob = mute_freq / 3.0  # TODO take into account different frequency of going to different bases
            else:  # TODO change this back to the commented block
                for nuke in (nuke1, nuke2):
                    if nuke == germline_nuke:
                        prob *= 1.0 - mute_freq
                    else:
                        prob *= mute_freq / 3.0
                # cryptic_factor_from_normalization = (math.sqrt(3.)*math.sqrt(-8.*mute_freq**2 + 16.*mute_freq + 27.) - 9.) / 12.
                # if nuke1 == germline_nuke and nuke2 == germline_nuke:
                #     prob = (1.0 - mute_freq)**2
                # elif nuke1 == nuke2 and nuke1 != germline_nuke:  # assume this requires *one* mutation event (i.e. ignore higher-order terms, I think)
                #     prob = cryptic_factor_from_normalization
                # elif nuke1 == germline_nuke or nuke2 == germline_nuke:
                #     prob = cryptic_factor_from_normalization
                # else:
                #     prob = cryptic_factor_from_normalization**2
        # print '  prob',prob,'for',nuke1,nuke2,is_insert,inuke,germline_nuke
        return prob
            
    # ----------------------------------------------------------------------------------------
    def add_emissions(self, state, inuke=-1, germline_nuke=''):
        # first add single emission
        emission_probs = {}
        total = 0.0
        for nuke in utils.nukes:
            emission_probs[nuke] = self.get_emission_prob(nuke, is_insert=(state.name=='insert'), inuke=inuke, germline_nuke=germline_nuke)
            total += emission_probs[nuke]
        if math.fabs(total - 1.0) >= self.eps:
            print 'ERROR emission not normalized in state %s in %s (%f)' % (state.name, 'X', total)  #utils.color_gene(gene_name), total)
            assert False
        state.add_emission(self.track, emission_probs)

        # then the pair emission
        pair_emission_probs = {}
        total = 0.0
        for nuke1 in utils.nukes:
            pair_emission_probs[nuke1] = {}
            for nuke2 in utils.nukes:
                pair_emission_probs[nuke1][nuke2] = self.get_emission_prob(nuke1, nuke2, is_insert=(state.name=='insert'), inuke=inuke, germline_nuke=germline_nuke)
                total += pair_emission_probs[nuke1][nuke2]
        if math.fabs(total - 1.0) >= self.eps:
            print 'ERROR pair emission not normalized in state %s in %s (%f)' % (state.name, 'X', total)  #utils.color_gene(gene_name), total)
            assert False
        state.add_pair_emission(self.track, pair_emission_probs)
