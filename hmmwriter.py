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
        self.precision = '2'  # number of digits after the decimal for probabilities
        self.germline_seqs = utils.read_germlines(data_dir)  #'/home/dralph/Dropbox/work/recombinator')
        self.text = ''  # text of the hmm description, to be written to file on completion

    def write(self):
        self.add_header()
        self.add_states()
        outfile.write(self.text)

    def get_end_probability(self, gene_name, seq, inuke):
        """ In other words, how much do we erode? """
        if inuke == len(seq) - 1:  # end of gene
            return 1.0
        elif len(seq) - 1 - inuke < 3:  # could erode to here
            return 0.2
        else:
            return 0.0
    
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
    def add_state(self, region, gene_name, seq, inuke, germline_nuke):
        saniname = utils.sanitize_name(gene_name)
        self.add_state_header('%s_%d' % (saniname, inuke), '%s_%d' % (saniname, inuke))

        self.add_transition_header()
        end_probability = self.get_end_probability(gene_name, seq, inuke) # probability of ending this region here, i.e. excising the rest of the germline gene. TODO use data
        if inuke < len(seq) - 1:
            self.text += ('  %s_%d:  %.' + self.precision + 'f\n') % (gene_name, inuke+1, 1 - end_probability)  # add a transition to next state
        if end_probability > 0.0:
            self.text += ('  insert:  %.' + self.precision + 'f\n') % end_probability  # and one to the insert state
        self.text += '  END: 1\n'  # TODO figure out how to swap back to <region>_end in previous line
        self.add_emission_header()
        self.text += '  ORDER: 0\n'
        emission_label_string = '@'
        for nuke in utils.nukes:
            emission_label_string += '%18s' % nuke
        self.text += emission_label_string + '\n'
        emission_probability_string = ''
        for nuke in utils.nukes:
            prob = 0.0
            if nuke == germline_nuke:
                prob = 0.94
            else:
                prob = 0.02
            emission_probability_string += ('%18.' + self.precision + 'f') % prob
        emission_probability_string = emission_probability_string.rstrip()
        self.text += emission_probability_string + '\n'
    
    # ----------------------------------------------------------------------------------------
    def add_states(self):
        self.text += '\nSTATE DEFINITIONS\n'
    
        # write init state
        self.add_state_header('INIT')
        self.add_transition_header()
        n_v_genes = len(self.germline_seqs['v'])
        total = 0.0
        for v_gene in self.germline_seqs['v']:
            probability = float(1./n_v_genes)  # TODO use the actual probabilities
            total += probability
            self.text += ('  %35s_0: %.' + self.precision + 'f\n') % (utils.sanitize_name(v_gene), probability)  # see gene probs in recombinator/data/human-beings/A/M/ighv-probs.txt
        assert math.fabs(total - 1.0) < 1e-10
    
        # write internal states
        for region in 'v':  #utils.regions:
            for gene_name in self.germline_seqs[region]:
                for inuke in range(len(self.germline_seqs[region][gene_name])):
                    nuke = self.germline_seqs[region][gene_name][inuke]
                    self.add_state(region, gene_name, self.germline_seqs[region][gene_name], inuke, nuke)
    
        # write insert state
        self.add_state_header('insert', 'i')
        self.add_transition_header()
        insertion_length = 7
        self.text += (' insert: %.' + self.precision + 'f\n') % (1./insertion_length)
        self.text += '  END: 1\n'
        self.add_emission_header()
        self.text += '  ORDER: 0\n'
        emission_probability_string = ''
        for nuke in utils.nukes:
            emission_probability_string += ('%18.' + self.precision + 'f') % (1./len(utils.nukes))
        emission_probability_string = emission_probability_string.rstrip()
        self.text += emission_probability_string + '\n'
    
        # TODO figure out how to swap back to <region>_end in previous line
        # # write region end state
        # self.text += '#############################################\n'
        # self.text += 'STATE:\n'
        # self.text += '  NAME:        %s_end\n' % region
        # self.text += '  PATH_LABEL:  %s_end\n' % region
        # self.text += '  GFF_DESC:    %s_end\n' % region
        # self.add_transition_header()
        # self.text += '   END: 1\n'
    
        # self.text += '  END:   1\n'  # TODO need to write the v_end state, and the END state. how to make v_end silent? did the way you tried *work*?
        self.text += '#############################################\n'
        self.text += '//END\n'
    
# ----------------------------------------------------------------------------------------
outfname = 'bcell/bcell.hmm'
with opener('w')(outfname) as outfile:
    writer = HmmWriter(outfile, 'v', '.')
    writer.write()
