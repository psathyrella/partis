#!/usr/bin/env python

import sys
import os
import math
import collections
from scipy.stats import norm
import csv
import utils
from opener import opener

# define this up here so the multi line string doesn't mess up the indentation below
header_base_text = """#STOCHHMM MODEL FILE
MODEL INFORMATION
======================================================
MODEL_NAME:	bcell
MODEL_DESCRIPTION:  maturation model
MODEL_CREATION_DATE:	today

TRACK SYMBOL DEFINITIONS
======================================================
NUKES: """

class HmmWriter(object):
    def __init__(self, base_indir, base_outdir, gene_name, naivety, germline_seq, v_right_length=-1):
        self.indir = base_indir
        self.precision = '16'  # number of digits after the decimal for probabilities. TODO increase this?
        self.v_right_length = v_right_length  # only take *this* much of the v gene, starting from the *right* end. mimics the fact that our reads don't extend all the way through v
        self.fuzz_around_v_left_edge = 5.0  # width of the normal distribution I'm using to account for uncertainty about where we jump into the v on the left side. TODO maybe change this?
        self.outdir = base_outdir  # + '/' + region
        self.region = utils.get_region(gene_name)
        self.gene_name = gene_name
        self.naivety = naivety
        self.germline_seq = germline_seq
        self.insertion = ''
        self.smallest_entry_index = -1
        if self.region == 'd':
            self.insertion = 'vd'
        elif self.region == 'j':
            self.insertion = 'dj'
        self.text = ''  # text of the hmm description, to be written to file on completion
        self.erosion_probs = {}
        self.read_erosion_probs()
        self.insertion_probs = {}
        if self.region != 'v':
            try:
                self.read_insertion_probs(False)
            except:
                self.read_insertion_probs(True)  # try again using all the alleles for this gene. TODO do something better
        self.mute_freqs = {}
        if self.naivety == 'M':  # mutate if not naive
            self.read_mute_freqs()

    # # ----------------------------------------------------------------------------------------
    # def get_erosion_key(erosion, line):
    #     key = [line[erosion + '_del']]
    #     for dep in deps:
    #         key.append(line[dep])
    #     return tuple(key)

    # ----------------------------------------------------------------------------------------
    def read_erosion_probs(self):
        for erosion in utils.erosions:
            if erosion[0] != self.region:
                continue
            self.erosion_probs[erosion] = {}
            deps = utils.column_dependencies[erosion + '_del']
            with opener('r')(self.indir + '/' + utils.get_prob_fname_key_val(erosion + '_del', deps)) as infile:
                reader = csv.DictReader(infile)
                total = 0.0
                for line in reader:
                    if self.region != 'j' and line[self.region + '_gene'] != self.gene_name:  # skip other genes (j erosion doesn't depend on gene choice)
                        continue
                    if int(line[erosion + '_del']) >= len(self.germline_seq):  # j erosion lengths don't depend on j gene, so we have to skip the ones that're too long
                        assert self.region == 'j'
                        continue
                    n_eroded = int(line[erosion + '_del'])
                    if n_eroded not in self.erosion_probs[erosion]:  # d erosion lengths depend on each other, but with the current hmm structure we cannot model this, so for the time being we integrate over the erosion on the other side. TODO fix that (well... maybe. it's kinda probably really not a big deal)
                        self.erosion_probs[erosion][n_eroded] = 0.0
                    self.erosion_probs[erosion][n_eroded] += float(line['count'])
                    total += float(line['count'])

                assert len(self.erosion_probs[erosion]) != 0
                test_total = 0.0
                for n_eroded in self.erosion_probs[erosion]:  # then normalize
                    self.erosion_probs[erosion][n_eroded] /= total
                    test_total += self.erosion_probs[erosion][n_eroded]
                assert utils.is_normed(test_total)

    # ----------------------------------------------------------------------------------------
    def read_insertion_probs(self, use_other_alleles=False):
        self.insertion_probs[self.insertion] = {}
        deps = utils.column_dependencies[self.insertion + '_insertion']
        with opener('r')(self.indir + '/' + utils.get_prob_fname_key_val(self.insertion + '_insertion', deps)) as infile:
            reader = csv.DictReader(infile)
            total = 0.0
            for line in reader:
                if self.region == 'j':  # for dj insertion, skip rows for genes other than the current one
                    assert self.insertion == 'dj'  # neuroticism is a positive personality trait
                    if use_other_alleles and not utils.are_alleles(line['j_gene'], self.gene_name):
                        continue
                    if not use_other_alleles and line['j_gene'] != self.gene_name:
                        continue
                n_inserted = int(line[self.insertion + '_insertion'])
                if n_inserted not in self.insertion_probs[self.insertion]:
                    self.insertion_probs[self.insertion][n_inserted] = 0.0
                self.insertion_probs[self.insertion][n_inserted] += float(line['count'])
                total += float(line['count'])

            assert len(self.insertion_probs[self.insertion]) != 0  # designed to fail

            test_total = 0.0
            for n_inserted in self.insertion_probs[self.insertion]:  # then normalize
                self.insertion_probs[self.insertion][n_inserted] /= total
                test_total += self.insertion_probs[self.insertion][n_inserted]
            assert utils.is_normed(test_total)
        
    # ----------------------------------------------------------------------------------------
    def read_mute_freqs(self):
        mutefname = self.indir + '/mute-freqs/' + utils.sanitize_name(self.gene_name) + '.csv'
        with opener('r')(mutefname) as mutefile:
            reader = csv.DictReader(mutefile)
            for line in reader:
                self.mute_freqs[int(line['position'])] = float(line['mute_freq'])  # TODO is there some way to incorporate the uncertainty on this?

    # ----------------------------------------------------------------------------------------
    def write(self, outfname='', remove=False):
        self.add_header()
        self.add_states()
        if outfname == '':
            outfname = self.outdir + '/' + utils.sanitize_name(self.gene_name) + '.hmm'
        with opener('w')(outfname) as outfile:
            outfile.write(self.text)
        self.text = ''
        if remove:  # if using mkfifo, remove afterwards
            os.remove(outfname)

    # # ----------------------------------------------------------------------------------------
    # def get_erosion_prob(self):  
    #     for inuke in range(len(self.germline_seq)):
    #         erosion_length = 0
    #         if '5p' in erosion:  
    #             erosion_length = inuke
    #         else:  # right erosions: prob to *leave* the region *after* inuke, i.e. to erode the bases to the rigth of inuke
    #             assert '3p' in erosion
    #             erosion_length = len(self.germline_seq) - inuke -1
    #         if str(erosion_length) in self.erosion_probs[erosion]:
    #             print '  %4d%4d%10f' % (inuke, erosion_length, self.erosion_probs[erosion][str(erosion_length)])
    #         else:
    #             continue
    #             print '  %4d%4d  hrg' % (inuke, erosion_length)
                        
    # ----------------------------------------------------------------------------------------
    def add_region_entry_probs(self, for_insert_state=False):
        """
        Probabilities to enter germline gene at point <inuke>.
        In the hmm file they appear as lines with the prob to go from INIT to the <inuke> state.
        For v, this is (mostly) the prob that our reads stop after <inuke> (going from right to left).
        For d and j, this is (mostly) the prob to *erode* up to (but not including) inuke.
        The two <mostly>s are there because in both cases, we're starting from *approximate* smith-waterman alignments, so we need to add some fuzz in case the s-w is off.
        <for_insert_state>: if we're adding the region entry probs for the insert state, we don't want t
        """

        # prob for non-zero-length insertion (i.e. prob to *not* go directly into the region
        non_zero_insertion_prob = 1.0
        if self.region != 'v':  # no insertion state in v hmm
            if for_insert_state:  # for the insert state, we want the prob of *leaving* the insert state to be 1./insertion_length, so multiply all the region entry probs by this
                insertion_length = self.get_mean_insert_length()
                self.text += (' insert: %.' + self.precision + 'f\n') % (1.0 - 1./insertion_length)
                non_zero_insertion_prob -= 1./insertion_length  # TODO I have no real idea if this is right. well, it's probably ok, but I really need to read through it again
            else:
                if 0 in self.insertion_probs[self.insertion]:  # if there is a non-zero prob of a zero-length insertion, subtract that prob from 1.0       *giggle*
                    non_zero_insertion_prob -= self.insertion_probs[self.insertion][0]
                else:
                    pass  # TODO think of something to check here
                self.text += (' insert: %.' + self.precision + 'f\n') % non_zero_insertion_prob

        # assert False  # arg, just realized I can't use the insertion length probs like this. It's an *hmm*, after all.
        #               # Well, looking at the distributions I made before... it's not a *horrible* approximation to use a geometric distribution or whatever the hell I get if I just use 1./<mean length>
        #               # TBD

        # prob to actually enter region
        # NOTE must be normalized to 1.0 - non_zero_insertion_prob for d and j
        if self.region == 'v':
            istart = 0
            if self.v_right_length != -1:
                istart = len(self.germline_seq) - self.v_right_length
            probs = collections.OrderedDict()
            total = 0.0
            for inuke in range(len(self.germline_seq)):  # start at far left side of v, but only write out probs that are greater than utils.eps (so will only write out probs near to 'location of' v_right_length
                tmp_prob = norm.pdf(float(inuke), float(istart), self.fuzz_around_v_left_edge)  # NOTE not yet normalized
                if tmp_prob < utils.eps:  # TODO this really probably goes further out into the tails than we need to, so probably slows us down quite a bit
                    continue
                probs[inuke] = tmp_prob  # normal distribution centered around istart
                total += tmp_prob
                if self.smallest_entry_index == -1 or inuke < self.smallest_entry_index:  # keep track of the first state that has a chance of being entered from INIT -- we want to start writing (with add_internal_state) from there
                    self.smallest_entry_index = inuke
            test_total = 0.0
            for inuke in probs:  # normalize and check
                probs[inuke] /= total
                test_total += probs[inuke]
            assert utils.is_normed(test_total)
            for inuke in probs:  # add to text
                self.text += ('  %18s_%s: %.' + self.precision + 'f\n') % (utils.sanitize_name(self.gene_name), inuke, probs[inuke])  # see gene probs in recombinator/data/human-beings/A/M/ighv-probs.txt
        else:  # TODO note that taking these numbers straight from data, with no smoothing, means that we are *forbidding* erosion lengths that we do not see in the training sample. Good? Bad? t.b.d.
            self.smallest_entry_index = 0
            total = 0.0
            for inuke in range(len(self.germline_seq)):
                erosion = self.region + '_5p'
                erosion_length = inuke
                if erosion_length in self.erosion_probs[erosion]:
                    prob = self.erosion_probs[erosion][erosion_length]
                    if prob < utils.eps:
                        continue
                    total += prob * (1.0 - non_zero_insertion_prob)
                    if non_zero_insertion_prob != 1.0:  # only add the line if there's a chance of zero-length insertion
                        self.text += ('  %18s_%d: %.' + self.precision + 'f\n') % (utils.sanitize_name(self.gene_name), inuke, prob * (1.0 - non_zero_insertion_prob))  # see gene probs in recombinator/data/human-beings/A/M/ighv-probs.txt
            assert non_zero_insertion_prob == 1.0 or utils.is_normed(total / (1.0 - non_zero_insertion_prob))

    # ----------------------------------------------------------------------------------------
    def get_exit_probability(self, seq, inuke):
        """
        Prob of exiting the chain of states for this region.
        In other words, what is the prob that we will erode all the bases to the right of <inuke>.
        """
        # TODO note that taking these numbers straight from data, with no smoothing, means that we are *forbidding* erosion lengths that we do not see in the training sample. Good? Bad? t.b.d.
        if inuke == len(seq) - 1:  # last state has to go to END
            return 1.0
        if self.region == 'j':  # always go to end of germline j region
            return 0.0
        erosion = self.region + '_3p'
        erosion_length = len(self.germline_seq) - inuke -1
        if erosion_length in self.erosion_probs[erosion]:
            prob = self.erosion_probs[erosion][erosion_length]
            if prob > utils.eps:
                return prob
            else:
                return 0.0
        else:
            return 0.0
    
    # ----------------------------------------------------------------------------------------
    def add_header(self):
        header_string = header_base_text
        for nuke in utils.nukes:
            header_string += nuke + ','
        header_string = header_string.rstrip(',')
        self.text += header_string + '\n'

    # ----------------------------------------------------------------------------------------
    def add_state_header(self, name, label=''):
        self.text += '#############################################\n'
        self.text += 'STATE:\n'
        self.text += '  NAME:        %s\n' % name
        if label != '':
            self.text += '  PATH_LABEL:  %s\n' % label
            self.text += '  GFF_DESC:    %s\n' % name
    
    # ----------------------------------------------------------------------------------------
    def add_transition_header(self):
        self.text += 'TRANSITION: STANDARD: P(X)\n'

    # ----------------------------------------------------------------------------------------
    def add_emission_header(self):
        self.text += 'EMISSION: NUKES: P(X)\n'
        self.text += '  ORDER: 0\n'

    # ----------------------------------------------------------------------------------------
    def add_init_state(self):
        self.add_state_header('INIT')
        self.add_transition_header()
        self.add_region_entry_probs()

    # ----------------------------------------------------------------------------------------
    def get_mean_insert_length(self):
        total, n_tot = 0.0, 0.0
        for length,count in self.insertion_probs[self.insertion].iteritems():
            total += count*length
            n_tot += count
        return total / n_tot

    # ----------------------------------------------------------------------------------------
    def add_insert_state(self):
        self.add_state_header('insert', 'i')
        self.add_transition_header()
        self.add_region_entry_probs(True)
        # TODO allow d region to be entirely eroded. Um, I don't think I have any probabilities for how frequentlyt his happens
        # if self.region == 'd':  # allow erosion of entire d region
        #     self.text += '  END: 1\n'
        self.add_emission_header()
        emission_probability_string = ''
        for nuke in utils.nukes:
            emission_probability_string += (' %18.' + self.precision + 'f') % (1./len(utils.nukes))  # TODO use insertion base composition from data
        emission_probability_string = emission_probability_string.rstrip()
        self.text += emission_probability_string + '\n'

    # ----------------------------------------------------------------------------------------
    def add_internal_state(self, seq, inuke, germline_nuke):
        # TODO unify utils.eps and self.precision
        # TODO the transition probs out of a state should add to one. But should this sum include the transitions to END, which stochhm requires
        #   to be 1 by themselves? AAAGGGGHGHGHHGHG I don't know. What the hell does stochhmm do?

        saniname = utils.sanitize_name(self.gene_name)
        self.add_state_header('%s_%d' % (saniname, inuke), '%s_%d' % (saniname, inuke))

        # transitions
        self.add_transition_header()
        exit_probability = self.get_exit_probability(seq, inuke) # probability of ending this region here, i.e. excising the rest of the germline gene
        if inuke < len(seq) - 1:  # if we're not at the end of this germline gene, add a transition to the next state
            self.text += ('  %s_%d:  %.' + self.precision + 'f\n') % (saniname, inuke+1, 1 - exit_probability)
        # if exit_probability >= utils.eps:  #10**(-int(self.precision)+1):  # don't write transitions that have zero probability
        #     self.text += ('  insert:  %.' + self.precision + 'f\n') % exit_probability  # and one to the insert state
        distance_to_end = len(seq) - inuke - 1
        if exit_probability >= utils.eps or distance_to_end == 0:
            self.text += ('  END: %.' + self.precision + 'f\n') % exit_probability

        # emissions
        self.add_emission_header()
        emission_label_string = '@'  # add line with human-readable headers
        for nuke in utils.nukes:
            emission_label_string += ' %18s' % nuke
        self.text += emission_label_string + '\n'
        emission_probability_string = ' '  # then add the line with actual probabilities
        for nuke in utils.nukes:
            mute_freq = 0.0
            if inuke in self.mute_freqs:
                mute_freq = self.mute_freqs[inuke]
            prob = 0.0
            if nuke == germline_nuke:
                prob = 1.0 - mute_freq
            else:
                prob = mute_freq / 3.0
            emission_probability_string += (' %18.' + self.precision + 'f') % prob
        emission_probability_string = emission_probability_string.rstrip()
        self.text += emission_probability_string + '\n'
    
    # ----------------------------------------------------------------------------------------
    def add_states(self):
        self.text += '\nSTATE DEFINITIONS\n'
        self.add_init_state()

        # for d and j regions add insert state to left-hand side of hmm
        if self.region != 'v':
            self.add_insert_state()

        # write internal states
        # istart = 0
        # if self.region == 'v' and self.v_right_length != -1:  # chop off the left side of the v
        #     istart = len(self.germline_seq) - self.v_right_length
        for inuke in range(self.smallest_entry_index, len(self.germline_seq)):
            nuke = self.germline_seq[inuke]
            self.add_internal_state(self.germline_seq, inuke, nuke)
    
        # finish up
        self.text += '#############################################\n'
        self.text += '//END\n'
