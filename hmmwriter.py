#!/usr/bin/env python

import sys
import math
import utils
from opener import opener

germline_seqs = utils.read_germlines('.')  #'/home/dralph/Dropbox/work/recombinator')
precision = '2'  # number of digits after the decimal for probabilities

def get_end_probability(gene_name, seq, inuke):
    """ In other words, how much do we erode? """
    if inuke == len(seq) - 1:  # end of gene
        return 1.0
    elif len(seq) - 1 - inuke < 3:  # could erode to here
        return 0.2
    else:
        return 0.0

def write_header(outfile):
    header_string = """#STOCHHMM MODEL FILE
MODEL INFORMATION
======================================================
MODEL_NAME:	bcell
MODEL_DESCRIPTION:  maturation model
MODEL_CREATION_DATE:	today

TRACK SYMBOL DEFINITIONS
======================================================
NUKES: """
    for nuke in utils.nukes:
        header_string += nuke + ','
    header_string = header_string.rstrip(',')
    outfile.write(header_string + '\n')

def write_state(outfile, region, gene_name, seq, inuke, germline_nuke):
    saniname = utils.sanitize_name(gene_name)
    outfile.write('#############################################\n')
    outfile.write('STATE:\n')
    outfile.write('  NAME:        %s_%d\n' % (saniname, inuke))
    outfile.write('  PATH_LABEL:  %s_%d\n' % (saniname, inuke))
    outfile.write('  GFF_DESC:    %s_%d\n' % (saniname, inuke))
    outfile.write('TRANSITION: STANDARD: P(X)\n')
    end_probability = get_end_probability(gene_name, seq, inuke) # probability of ending this region here, i.e. excising the rest of the germline gene. TODO use data
    if inuke < len(seq) - 1:
        outfile.write(('   %s_%d:  %.' + precision + 'f\n') % (gene_name, inuke+1, 1 - end_probability))  # add a transition to next state
    outfile.write(('  %s_end:  %.' + precision + 'f\n') % (region, end_probability))  # and one the end_region state
    outfile.write('EMISSION: NUKES: P(X)\n')
    outfile.write('  ORDER: 0\n')
    emission_label_string = '@'
    for nuke in utils.nukes:
        emission_label_string += '%18s' % nuke
    outfile.write(emission_label_string + '\n')
    emission_probability_string = ''
    for nuke in utils.nukes:
        prob = 0.0
        if nuke == germline_nuke:
            prob = 0.94
        else:
            prob = 0.02
        emission_probability_string += ('%18.' + precision + 'f') % prob
    emission_probability_string = emission_probability_string.rstrip()
    outfile.write(emission_probability_string + '\n')

def write_states(outfile):
    outfile.write('\nSTATE DEFINITIONS\n')
    # write init state
    outfile.write('#############################################\n')
    outfile.write('STATE:\n')
    outfile.write('  NAME: INIT\n')
    outfile.write('TRANSITION: STANDARD: P(X)\n')
    n_v_genes = len(germline_seqs['v'])
    total = 0.0
    for v_gene in germline_seqs['v']:
        probability = float(1./n_v_genes)  # TODO use the actual probabilities
        total += probability
	outfile.write(('  %35s_0: %.' + precision + 'f\n') % (utils.sanitize_name(v_gene), probability))  # see gene probs in recombinator/data/human-beings/A/M/ighv-probs.txt
    assert math.fabs(total - 1.0) < 1e-10

    # write internal states
    for region in 'v':  #utils.regions:
        for gene_name in germline_seqs[region]:
            for inuke in range(len(germline_seqs[region][gene_name])):
                nuke = germline_seqs[region][gene_name][inuke]
                write_state(outfile, region, gene_name, germline_seqs[region][gene_name], inuke, nuke)

    # write region end state
    outfile.write('#############################################\n')
    outfile.write('STATE:\n')
    outfile.write('  NAME:        %s_end\n' % region)
    outfile.write('  PATH_LABEL:  %s_end\n' % region)
    outfile.write('  GFF_DESC:    %s_end\n' % region)
    outfile.write('TRANSITION: STANDARD: P(X)\n')
    outfile.write('   END: 1\n')

    # outfile.write('  END:   1\n')  # TODO need to write the v_end state, and the END state. how to make v_end silent? did the way you tried *work*?
    outfile.write('#############################################\n')
    outfile.write('//END\n')

outfname = 'bcell/bcell.hmm'
with opener('w')(outfname) as outfile:
    write_header(outfile)
    write_states(outfile)
