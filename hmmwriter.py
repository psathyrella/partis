#!/usr/bin/env python

import sys
import math
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
    def __init__(self, outfile, region, data_dir):
        self.precision = '16'  # number of digits after the decimal for probabilities. TODO increase this. I just have it at two for human readibility while debugging
        self.region = region
        self.germline_seqs = utils.read_germlines(data_dir)
        self.text = ''  # text of the hmm description, to be written to file on completion

    def write(self):
        self.add_header()
        self.add_states()
        outfile.write(self.text)

    # ----------------------------------------------------------------------------------------
    def add_region_entry_probs(self):
        """ Probabilities to enter <region> at point <inuke> in version <gene_name> """
        # TODO use probs from data
        max_erosion = 7
        n_genes = len(self.germline_seqs[self.region])
        total = 0.0
        for gene,seq in self.germline_seqs[self.region].iteritems():
            for inuke in range(max_erosion):  #len(seq)):
                gene_prob = float(1./n_genes)  # prob of choosing this gene
                erosion_prob = float(1./max_erosion)  # prob of entering this germline gene at this position, i.e. of eroding until here
                probability = gene_prob * erosion_prob
                total += probability
                self.text += ('  %35s_%d: %.' + self.precision + 'f\n') % (utils.sanitize_name(gene), inuke, probability)  # see gene probs in recombinator/data/human-beings/A/M/ighv-probs.txt
        assert math.fabs(total - 1.0) < 1e-10  # make sure probs sum to one

    # ----------------------------------------------------------------------------------------
    def get_exit_probability(self, seq, inuke):
        """
        Prob of exiting the chain of states for this region and entering the insert state.
        In other words, how much do we erode?
        """
        if self.region == 'j':  # no insert state to right of j -- we go until sequence ends, i.e. we get a transition to END
            return 0.0
        distance_to_end = len(seq) - inuke - 1
        decay_length = 7  # number of bases over which probability decays to 1/e of initial value (using a e^(-x) a.t.m.
        return math.exp(-distance_to_end / decay_length)
    
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

    # ----------------------------------------------------------------------------------------
    def add_init_state(self):
        self.add_state_header('INIT')
        self.add_transition_header()
        self.add_region_entry_probs()

    # ----------------------------------------------------------------------------------------
    def add_insert_state(self):
        self.add_state_header('insert', 'i')
        self.add_transition_header()
        insertion_length = 7
        self.text += (' insert: %.' + self.precision + 'f\n') % (1./insertion_length)
        self.text += '  END: 1\n'
        self.add_emission_header()
        self.text += '  ORDER: 0\n'
        emission_probability_string = ''
        for nuke in utils.nukes:
            emission_probability_string += (' %18.' + self.precision + 'f') % (1./len(utils.nukes))
        emission_probability_string = emission_probability_string.rstrip()
        self.text += emission_probability_string + '\n'

    # ----------------------------------------------------------------------------------------
    def add_internal_state(self, gene_name, seq, inuke, germline_nuke):
        saniname = utils.sanitize_name(gene_name)
        self.add_state_header('%s_%d' % (saniname, inuke), '%s_%d' % (saniname, inuke))

        self.add_transition_header()
        exit_probability = self.get_exit_probability(seq, inuke) # probability of ending this region here, i.e. excising the rest of the germline gene. TODO use data
        if inuke < len(seq) - 1:
            self.text += ('  %s_%d:  %.' + self.precision + 'f\n') % (saniname, inuke+1, 1 - exit_probability)  # add a transition to next state
        if exit_probability > 0.0:
            self.text += ('  insert:  %.' + self.precision + 'f\n') % exit_probability  # and one to the insert state
        self.text += '  END: 1\n'
        self.add_emission_header()
        self.text += '  ORDER: 0\n'
        emission_label_string = '@'
        for nuke in utils.nukes:
            emission_label_string += ' %18s' % nuke
        self.text += emission_label_string + '\n'
        emission_probability_string = ' '
        for nuke in utils.nukes:
            prob = 0.0
            if nuke == germline_nuke:  # TODO use real mute probs
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
        for gene_name in self.germline_seqs[self.region]:
            for inuke in range(len(self.germline_seqs[self.region][gene_name])):
                nuke = self.germline_seqs[self.region][gene_name][inuke]
                self.add_internal_state(gene_name, self.germline_seqs[self.region][gene_name], inuke, nuke)
    
        # for v and d regions add insert state to right-hand side of hmm
        if self.region == 'v' or self.region == 'd':
            self.add_insert_state()

        # finish up
        self.text += '#############################################\n'
        self.text += '//END\n'
    
# ----------------------------------------------------------------------------------------
region = 'd'
outfname = 'bcell/' + region + '.hmm'
with opener('w')(outfname) as outfile:
    writer = HmmWriter(outfile, region, '.')  #'/home/dralph/Dropbox/work/recombinator')
    writer.write()
