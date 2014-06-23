#!/usr/bin/env python

import sys
import os
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
    def __init__(self, base_outdir, region, gene_name, germline_seq):
        self.precision = '16'  # number of digits after the decimal for probabilities. TODO increase this. I just have it at two for human readibility while debugging
        self.v_right_length = 50  # only take *this* much of the v gene, starting from the *right* end. mimics the fact that our reads don't extend all the way through v
        self.outdir = outdir
        self.region = region
        self.gene_name = gene_name
        self.germline_seq = germline_seq
        self.text = ''  # text of the hmm description, to be written to file on completion

    # ----------------------------------------------------------------------------------------
    def write(self):
        self.add_header()
        self.add_states()
        outfname = self.outdir + '/' + utils.sanitize_name(self.gene_name) + '.hmm'
        with opener('w')(outfname) as outfile:
            outfile.write(self.text)
        self.text = ''

    # ----------------------------------------------------------------------------------------
    def add_region_entry_probs(self):
        """ Probabilities to enter germline gene at point <inuke> """
        # TODO use probs from data
        max_erosion = 7
        total = 0.0
        istart = 0
        if region == 'v' and self.v_right_length != -1:
            istart = len(self.germline_seq) - self.v_right_length
        for inuke in range(istart, istart+max_erosion):  #len(seq)):
            probability = float(1./max_erosion)  # prob of entering the germline gene at this position, i.e. of eroding until here (*left* side erosion)
            total += probability
            self.text += ('  %35s_%d: %.' + self.precision + 'f\n') % (utils.sanitize_name(self.gene_name), inuke, probability)  # see gene probs in recombinator/data/human-beings/A/M/ighv-probs.txt
        assert math.fabs(total - 1.0) < 1e-10  # make sure probs sum to one

    # ----------------------------------------------------------------------------------------
    def get_exit_probability(self, seq, inuke):
        """
        Prob of exiting the chain of states for this region and entering the insert state.
        In other words, how much do we erode? (*right* side erosion)
        """
        if self.region == 'j':  # no insert state to right of j -- we go until query sequence ends, i.e. until we get a transition to END
            return 0.0
        distance_to_end = len(seq) - inuke - 1
        decay_length = 7  # number of bases over which probability decays to 1/e of initial value (using a e^(-x) a.t.m.
        return math.exp(-float(distance_to_end + 5) / decay_length)  # TODO use probs from data
    
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
        saniname = utils.sanitize_name(self.gene_name)
        self.add_state_header('%s_%d' % (saniname, inuke), '%s_%d' % (saniname, inuke))

        # transitions
        self.add_transition_header()
        exit_probability = self.get_exit_probability(seq, inuke) # probability of ending this region here, i.e. excising the rest of the germline gene
        if inuke < len(seq) - 1:  # if we're not at the end of this germline gene, add a transition to the next state
            self.text += ('  %s_%d:  %.' + self.precision + 'f\n') % (saniname, inuke+1, 1 - exit_probability)
        if exit_probability > 10**(-int(self.precision)+1):  # don't write transitions that have zero probability
            self.text += ('  insert:  %.' + self.precision + 'f\n') % exit_probability  # and one to the insert state
        self.text += '  END: 1\n'  # TODO don't write this transition to END for states (e.g. the start of v) that have no real chance to transition to END)

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
    
# ----------------------------------------------------------------------------------------

n_max_versions = 5  # only look at the first n gene versions (speeds things up for testing)
germline_seqs = utils.read_germlines('/home/dralph/Dropbox/work/recombinator')

for region in utils.regions:
    outdir = 'bcell/' + region
    if os.path.exists(outdir):
        for hmmfile in os.listdir(outdir):
            if hmmfile.endswith(".hmm"):
                os.remove(outdir + "/" + hmmfile)
    else:
        os.makedirs(outdir)
    
    igene = 0
    for gene_name in germline_seqs[region]:
        if igene >= n_max_versions:
            print 'breaking after %d gene versions' % n_max_versions
            break
        print '  %d / %d (%s)' % (igene, len(germline_seqs[region]), gene_name)
        igene += 1
        writer = HmmWriter('bcell', region, gene_name, germline_seqs[region][gene_name])
        writer.write()
