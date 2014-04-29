""" A few utility functions. At the moment simply functions used in recombinator which do not
require member variables. """

import sys

#----------------------------------------------------------------------------------------
regions = ['v', 'd', 'j']
erosions = ['v_3p', 'd_5p', 'd_3p', 'j_5p']
boundaries = ('vd', 'dj')
# Infrastrucure to allow hashing all the columns together into a dict key.
# Uses a tuple with the variables that are used to index selection frequencies
index_columns = ('v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del')
index_keys = {}
for i in range(len(index_columns)):  # dict so we can access them by name instead of by index number
    index_keys[index_columns[i]] = i

#----------------------------------------------------------------------------------------
def int_to_nucleotide(number):
    """ Convert between (0,1,2,3) and (A,C,G,T) """
    if number == 0:
        return 'A'
    elif number == 1:
        return 'C'
    elif number == 2:
        return 'G'
    elif number == 3:
        return 'T'
    else:
        print 'ERROR nucleotide number not in [0,3]'
        sys.exit()

#----------------------------------------------------------------------------------------                    
def check_conserved_cysteine(seq, cyst_position):
    """ Ensure there's a cysteine at <cyst_position> in <seq>. """
    cyst_word = str(seq[cyst_position:cyst_position+3])
    if cyst_word != 'TGT' and cyst_word != 'TGC':
        print 'ERROR cysteine in V is messed up: %s' % cyst_word
        assert False
def check_conserved_tryptophan(seq, tryp_position):
    """ Ensure there's a tryptophan at <tryp_position> in <seq>. """
    tryp_word = str(seq[tryp_position:tryp_position+3])
    if tryp_word != 'TGG':
        print 'ERROR tryptophan in J is messed up: %s' % tryp_word
        assert False
def check_conserved_codons(seq, cyst_position, tryp_position):
    """ Double check that we conserved the cysteine and the tryptophan. """
    check_conserved_cysteine(seq, cyst_position)
    check_conserved_tryptophan(seq, tryp_position)
def are_conserved_codons_screwed_up(reco_event):
    """ Version that checks all the final seqs in reco_event.

    Returns True if codons are screwed up, or if no sequences have been added.
    """
    if len(reco_event.final_seqs) == 0:
        return True
    for seq in reco_event.final_seqs:
        try:
            check_conserved_codons(seq, reco_event.cyst_position, reco_event.final_tryp_position)
        except:
            return True

    return False

#----------------------------------------------------------------------------------------
def is_position_protected(protected_positions, prospective_position):
    """ Would a mutation at <prospective_position> screw up a protected codon? """
    for position in protected_positions:
        if (prospective_position == position or
            prospective_position == (position + 1) or
            prospective_position == (position + 2)):
            return True
    return False

#----------------------------------------------------------------------------------------
def would_erode_conserved_codon(reco_event):
    """ Would any of the erosion <lengths> delete a conserved codon? """
    lengths = reco_event.erosions
    # check conserved cysteine
    if len(reco_event.seqs['v']) - lengths['v_3p'] <= reco_event.cyst_position + 2:
        print '      about to erode cysteine (%d), try again' % lengths['v_3p']
        return True  # i.e. it *would* screw it up
    # check conserved tryptophan
    if lengths['j_5p'] - 1 >= reco_event.tryp_position:
        print '      about to erode tryptophan (%d), try again' % lengths['j_5p']
        return True

    return False  # *whew*, it won't erode either of 'em

#----------------------------------------------------------------------------------------
def is_erosion_longer_than_seq(reco_event):
    """ Are any of the proposed erosion <lengths> longer than the seq to be eroded? """
    lengths = reco_event.erosions
    if lengths['v_3p'] > len(reco_event.seqs['v']):  # NOTE not actually possible since we already know we didn't erode the cysteine
        print '      v_3p erosion too long (%d)' % lengths['v_3p']
        return True
    if lengths['d_5p'] + lengths['d_3p'] > len(reco_event.seqs['d']):
        print '      d erosions too long (%d)' % (lengths['d_5p'] + lengths['d_3p'])
        return True
    if lengths['j_5p'] > len(reco_event.seqs['j']):  # NOTE also not possible for the same reason
        print '      j_5p erosion too long (%d)' % lengths['j_5p']
        return True
    return False

#----------------------------------------------------------------------------------------
#def find_tryp_in_joined_seq(reco_event):
    # NOTE remember this is the tryp position, in j alone, *before* erosion
#                                                      tryp_position,  
#                                                      reco_event.seqs['v'],
#                                                      reco_info['vd_insertion'],
#                                                      reco_event.seqs['d'],
#                                                      reco_info['dj_insertion'],
#                                                      reco_event.seqs['j'],
#                                                      reco_info['j_5p_del'])
#
def find_tryp_in_joined_seq(tryp_position_in_j, v_seq, vd_insertion, d_seq, dj_insertion, j_seq, j_erosion):
    """ Find the <end> tryptophan in a joined sequence.

    Given local tryptophan position in the j region, figure
    out what position it's at in the final sequence.
    NOTE tryp_position_in_j is the position *before* the j was eroded,
    but this fcn assumes that the j *has* been eroded.
    """
#    print 'checking tryp with: %s, %d -%d = %d' % (j_seq, tryp_position_in_j, j_erosion, tryp_position_in_j - j_erosion)
    check_conserved_tryptophan(j_seq, tryp_position_in_j - j_erosion)  # make sure tryp is where it's supposed to be
    length_to_left_of_j = len(v_seq + vd_insertion + d_seq + dj_insertion)
#    print '  finding tryp position as'
#    print '    length_to_left_of_j = len(v_seq + vd_insertion + d_seq + dj_insertion) = %d + %d + %d + %d' % (len(v_seq), len(vd_insertion), len(d_seq), len(dj_insertion))
#    print '    result = tryp_position_in_j - j_erosion + length_to_left_of_j = %d - %d + %d = %d' % (tryp_position_in_j, j_erosion, length_to_left_of_j, tryp_position_in_j - j_erosion + length_to_left_of_j)
    return tryp_position_in_j - j_erosion + length_to_left_of_j

class Colors:
    head = '\033[95m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    end = '\033[0m'

def red(seq):
    return Colors.red + seq + Colors.end

def is_mutated(original, final):
    if original == final:
        return final
    else:
        return red(final)

def print_reco_event(germlines, line):
    """ Print ascii summary of recombination event and mutation. """
    original_seqs = {}
    for region in regions:
        original_seqs[region] = germlines[region][line[region+'_gene']]
    v_start = 0
    v_length = len(original_seqs['v']) - int(line['v_3p_del'])
    d_start = v_length + len(line['vd_insertion'])
    d_length = len(original_seqs['d']) - int(line['d_5p_del']) - int(line['d_3p_del'])
    j_start = v_length + len(line['vd_insertion']) + d_length + len(line['dj_insertion'])
    j_length = len(original_seqs['j']) - int(line['j_5p_del'])
    eroded_seqs = {}
    eroded_seqs['v'] = original_seqs['v'][:len(original_seqs['v'])-int(line['v_3p_del'])]
    eroded_seqs['d'] = original_seqs['d'][int(line['d_5p_del']) : len(original_seqs['d'])-int(line['d_3p_del'])]
    eroded_seqs['j'] = original_seqs['j'][int(line['j_5p_del']) :]

    germline_v_end = len(original_seqs['v']) - 1
    germline_d_start = len(original_seqs['v']) - int(line['v_3p_del']) + len(line['vd_insertion']) - int(line['d_5p_del'])
    germline_d_end = germline_d_start + len(original_seqs['d'])
    germline_j_start = germline_d_end + 1 - int(line['d_3p_del']) + len(line['dj_insertion']) - int(line['j_5p_del'])

    final_seq = ''
    for inuke in range(len(line['seq'])):
        ilocal = inuke
        if ilocal < v_length:
            final_seq += is_mutated(eroded_seqs['v'][ilocal], line['seq'][inuke])
        else:
            ilocal -= v_length
            if ilocal < len(line['vd_insertion']):
                final_seq += is_mutated(line['vd_insertion'][ilocal], line['seq'][inuke])
            else:
                ilocal -= len(line['vd_insertion'])
                if ilocal < d_length:
                    final_seq += is_mutated(eroded_seqs['d'][ilocal], line['seq'][inuke])
                else:
                    ilocal -= d_length
                    if ilocal < len(line['dj_insertion']):
                        final_seq += is_mutated(line['dj_insertion'][ilocal], line['seq'][inuke])
                    else:
                        ilocal -= len(line['dj_insertion'])
                        final_seq += is_mutated(eroded_seqs['j'][ilocal], line['seq'][inuke])

    # pad with dots
    eroded_seqs['v'] = eroded_seqs['v'] + int(line['v_3p_del']) * '.'
    eroded_seqs['d'] = int(line['d_5p_del']) * '.' + eroded_seqs['d'] + int(line['d_3p_del']) * '.'
    eroded_seqs['j'] = int(line['j_5p_del']) * '.' + eroded_seqs['j']

    insertions = v_length * ' ' + line['vd_insertion'] + d_length * ' ' + line['dj_insertion'] + j_length * ' '
    vj = germline_d_start * ' ' + eroded_seqs['d'] + (len(original_seqs['j']) - int(line['j_5p_del']) + len(line['dj_insertion']) - int(line['d_3p_del'])) * ' '
    d = eroded_seqs['v'] + (germline_j_start - germline_v_end - 2) * ' ' + eroded_seqs['j']

    print '    ',insertions
    print '    ',final_seq
    print '    ',vj
    print '    ',d
    if 'ack' in line and line['ack']:
        sys.exit()
#    assert len(line['seq']) == line['v_5p_del'] + len(hmms['v']) + len(outline['vd_insertion']) + len(hmms['d']) + len(outline['dj_insertion']) + len(hmms['j']) + outline['j_3p_del']
