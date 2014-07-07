#!/usr/bin/env python

import sys
import os
import math
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
    def __init__(self, base_indir, base_outdir, region, gene_name, germline_seq):
        self.indir = base_indir
        self.precision = '16'  # number of digits after the decimal for probabilities. TODO increase this?
        self.v_right_length = 100  # only take *this* much of the v gene, starting from the *right* end. mimics the fact that our reads don't extend all the way through v
        self.fuzz_around_v_left_edge = 5.0  # width of the normal distribution I'm using to account for uncertainty about where we jump into the v on the left side
        self.outdir = base_outdir + '/' + region
        self.region = region
        self.gene_name = gene_name
        self.germline_seq = germline_seq
        self.insertion = ''
        if self.region == 'd':
            self.insertion = 'vd'
        elif self.region == 'j':
            self.insertion = 'dj'
        self.text = ''  # text of the hmm description, to be written to file on completion
        self.erosion_probs = {}
        self.read_erosion_probs()
        self.insertion_probs = {}
        if self.region != 'v':
            self.read_insertion_probs()

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
                    n_eroded = line[erosion + '_del']  # NOTE n_eroded is a string
                    if n_eroded not in self.erosion_probs[erosion]:  # d erosion lengths depend on each other, but with the current hmm structure we cannot model this, so for the time being we integrate over the erosion on the other side. TODO fix that (well... maybe. it's kinda probably really not a big deal)
                        self.erosion_probs[erosion][n_eroded] = 0.0
                    self.erosion_probs[erosion][n_eroded] += float(line['count'])
                    total += float(line['count'])

                test_total = 0.0
                for n_eroded in self.erosion_probs[erosion]:  # then normalize
                    self.erosion_probs[erosion][n_eroded] /= total
                    test_total += self.erosion_probs[erosion][n_eroded]
                assert utils.is_normed(test_total)

    # ----------------------------------------------------------------------------------------
    def read_insertion_probs(self):

        self.insertion_probs[self.insertion] = {}
        deps = utils.column_dependencies[self.insertion + '_insertion']
        with opener('r')(self.indir + '/' + utils.get_prob_fname_key_val(self.insertion + '_insertion', deps)) as infile:
            reader = csv.DictReader(infile)
            total = 0.0
            for line in reader:
                if self.region == 'j' and line['j_gene'] != self.gene_name:  # for dj insertion, skip rows for genes other than the current one
                    assert self.insertion == 'dj'  # neuroticism is a positive personality trait
                    continue
                n_inserted = line[self.insertion + '_insertion']  # NOTE n_inserted is a string
                if n_inserted not in self.insertion_probs[self.insertion]:
                    self.insertion_probs[self.insertion][n_inserted] = 0.0
                self.insertion_probs[self.insertion][n_inserted] += float(line['count'])
                total += float(line['count'])

            test_total = 0.0
            for n_inserted in self.insertion_probs[self.insertion]:  # then normalize
                self.insertion_probs[self.insertion][n_inserted] /= total
                test_total += self.insertion_probs[self.insertion][n_inserted]
            assert utils.is_normed(test_total)
        
    # ----------------------------------------------------------------------------------------
    def write(self):
        self.add_header()
        self.add_states()
        outfname = self.outdir + '/' + utils.sanitize_name(self.gene_name) + '.hmm'
        with opener('w')(outfname) as outfile:
            outfile.write(self.text)
        self.text = ''

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
    def add_region_entry_probs(self):
        """
        Probabilities to enter germline gene at point <inuke>.
        In the hmm file they appear as lines with the prob to go from INIT to the <inuke> state.
        For v, this is (mostly) the prob that our reads stop after <inuke> (going from right to left).
        For d and j, this is (mostly) the prob to *erode* up to (but not including) inuke.
        The two <mostly>s are there because in both cases, we're starting from *approximate* smith-waterman alignments, so we need to add some fuzz in case the s-w is off.
        """

        # prob for non-zero-length insertion (i.e. prob to *not* go directly into the region
        if self.region != 'v':  # no insertion state in v hmm
            prob = 1.0
            if '0' in self.insertion_probs[self.insertion]:  # if there is a non-zero prob of a zero-length insertion, subtract that prob from 1.0       *giggle*
                assert '1' in self.insertion_probs[self.insertion] or '2' in self.insertion_probs[self.insertion] or '3' in self.insertion_probs[self.insertion] # this is kind of weird, but it's just to make sure I don't switch to indexing by integers instead of strings
                prob -= self.insertion_probs[self.insertion]['0']
            self.text += (' insert: %.' + self.precision + 'f\n') % prob

        assert False  # arg, just realized I can't use the insertion length probs like this. It's an *hmm*, after all.
                      # Well, looking at the distributions I made before... it's not a *horrible* approximation to use a geometric distribution or whatever the hell I get if I just use 1./<mean length>
                      # TBD

        # prob to actually enter region
        if self.region == 'v':
            istart = 0
            if self.v_right_length != -1:
                istart = len(self.germline_seq) - self.v_right_length
            probs = {}
            total = 0.0
            for inuke in range(len(self.germline_seq)):  # start at far left side of v, but only write out probs that are greater than utils.eps (so will only write out probs near to 'location of' v_right_length
                tmp_prob = norm.pdf(float(inuke), float(istart), self.fuzz_around_v_left_edge)  # NOTE not yet normalized
                if tmp_prob < utils.eps:  # TODO this really probably goes further out into the tails than we need to, so probably slows us down quite a bit
                    continue
                probs[inuke] = tmp_prob  # normal distribution centered around istart
                total += tmp_prob
            test_total = 0.0
            for inuke in probs:  # normalize and check
                probs[inuke] /= total
                test_total += probs[inuke]
            assert utils.is_normed(test_total)
            for inuke in probs:  # add to text
                self.text += ('  %35s_%s: %.' + self.precision + 'f\n') % (utils.sanitize_name(self.gene_name), inuke, probs[inuke])  # see gene probs in recombinator/data/human-beings/A/M/ighv-probs.txt
        else:  # TODO note that taking these numbers straight from data, with no smoothing, means that we are *forbidding* erosion lengths that we do not see in the training sample. Good? Bad? t.b.d.
            total = 0.0
            for inuke in range(len(self.germline_seq)):
                erosion = self.region + '_5p'
                erosion_length = inuke
                if str(erosion_length) in self.erosion_probs[erosion]:
                    prob = self.erosion_probs[erosion][str(erosion_length)]
                    if prob < utils.eps:
                        continue
                    total += prob
                    self.text += ('  %35s_%d: %.' + self.precision + 'f\n') % (utils.sanitize_name(self.gene_name), inuke, prob)  # see gene probs in recombinator/data/human-beings/A/M/ighv-probs.txt
            assert utils.is_normed(total)

    # ----------------------------------------------------------------------------------------
    def get_exit_probability(self, seq, inuke):
        """
        Prob of exiting the chain of states for this region.
        In other words, what is the prob that we will erode all the bases to the right of <inuke>.
        """
        # TODO note that taking these numbers straight from data, with no smoothing, means that we are *forbidding* erosion lengths that we do not see in the training sample. Good? Bad? t.b.d.
        if self.region == 'j':  # always go to end of germline j region
            return 0.0
        erosion = self.region + '_3p'
        erosion_length = len(self.germline_seq) - inuke -1
        if str(erosion_length) in self.erosion_probs[erosion]:
            prob = self.erosion_probs[erosion][str(erosion_length)]
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
    def add_insert_state(self):
        self.add_state_header('insert', 'i')
        self.add_transition_header()
        insertion_length = 7  # TODO use lengths from data
        self.text += (' insert: %.' + self.precision + 'f\n') % (1./insertion_length)
        self.text += '  END: 1\n'
        self.add_emission_header()
        emission_probability_string = ''
        for nuke in utils.nukes:
            emission_probability_string += (' %18.' + self.precision + 'f') % (1./len(utils.nukes))  # TODO use emission probs from data
        emission_probability_string = emission_probability_string.rstrip()
        self.text += emission_probability_string + '\n'

    # ----------------------------------------------------------------------------------------
    def add_internal_state(self, seq, inuke, germline_nuke):
        # TODO unify utils.eps and self.precision
        # TODO the transition probs out of a state should add to one. But should this sum include the transitions to END, which stochhm requires to be 1 by themselves? AAAGGGGHGHGHHGHG I don't know. What the hell does stochhmm do?

        saniname = utils.sanitize_name(self.gene_name)
        self.add_state_header('%s_%d' % (saniname, inuke), '%s_%d' % (saniname, inuke))

        # transitions
        self.add_transition_header()
        exit_probability = self.get_exit_probability(seq, inuke) # probability of ending this region here, i.e. excising the rest of the germline gene
        if inuke < len(seq) - 1:  # if we're not at the end of this germline gene, add a transition to the next state
            self.text += ('  %s_%d:  %.' + self.precision + 'f\n') % (saniname, inuke+1, 1 - exit_probability)
        if exit_probability >= utils.eps:  #10**(-int(self.precision)+1):  # don't write transitions that have zero probability
            self.text += ('  insert:  %.' + self.precision + 'f\n') % exit_probability  # and one to the insert state
        distance_to_end = len(seq) - inuke - 1
        if exit_probability >= utils.eps or distance_to_end == 0:
            self.text += '  END: 1\n'

        # emissions
        self.add_emission_header()
        emission_label_string = '@'  # add line with human-readable headers
        for nuke in utils.nukes:
            emission_label_string += ' %18s' % nuke
        self.text += emission_label_string + '\n'
        emission_probability_string = ' '  # then add the line with actual probabilities
        for nuke in utils.nukes:
            prob = 0.0
            if nuke == germline_nuke:  # TODO use data mute probs
                prob = 0.94
            else:
                prob = 0.02
            emission_probability_string += (' %18.' + self.precision + 'f') % prob
        emission_probability_string = emission_probability_string.rstrip()
        self.text += emission_probability_string + '\n'
    
    # ----------------------------------------------------------------------------------------
    def add_states(self):
        self.text += '\nSTATE DEFINITIONS\n'
        self.add_init_state()

        # write internal states
        istart = 0
        if self.region == 'v' and self.v_right_length != -1:  # chop off the left side of the v
            istart = len(self.germline_seq) - self.v_right_length
        for inuke in range(istart, len(self.germline_seq)):
            nuke = self.germline_seq[inuke]
            self.add_internal_state(self.germline_seq, inuke, nuke)
    
        # for v and d regions add insert state to right-hand side of hmm
        if self.region == 'v' or self.region == 'd':
            self.add_insert_state()

        # finish up
        self.text += '#############################################\n'
        self.text += '//END\n'
