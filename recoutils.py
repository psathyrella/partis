""" A few utility functions. At the moment simply functions used in recombinator which do not
require member variables. """
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
        sys.exit()
def check_conserved_tryptophan(seq, tryp_position):
    """ Ensure there's a tryptophan at <tryp_position> in <seq>. """
    tryp_word = str(seq[tryp_position:tryp_position+3])
    if tryp_word != 'TGG':
        print 'ERROR tryptophan in J is messed up: %s' % tryp_word
        sys.exit()
def check_conserved_codons(seq, cyst_position, tryp_position):
    """ Double check that we conserved the cysteine and the tryptophan. """
    check_conserved_cysteine(seq, cyst_position)
    check_conserved_tryptophan(seq, tryp_position)

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
def would_erode_conserved_codon(lengths, seqs, cyst_position, tryp_position):
    """ Would any of the erosion <lengths> delete a conserved codon? """
    # check conserved cysteine
    if len(seqs['v']) - lengths['v_3p_del'] <= cyst_position + 2:
        print '      about to erode cysteine (%d), try again' % lengths['v_3p_del']
        return True  # i.e. it *would* screw it up
    # check conserved tryptophan
    if lengths['j_5p_del'] - 1 >= tryp_position:
        print '      about to erode tryptophan (%d), try again' % lengths['j_5p_del']
        return True

    return False  # *whew*, it won't erode either of 'em

#----------------------------------------------------------------------------------------
def is_erosion_longer_than_seq(lengths, seqs):
    """ Are any of the proposed erosion <lengths> longer than the seq to be eroded? """
    if lengths['v_3p_del'] > len(seqs['v']):  # NOTE not actually possible since we already know we didn't erode the cysteine
        print '      v_3p_del erosion too long (%d)' % lengths['v_3p_del']
        return True
    if lengths['d_5p_del'] + lengths['d_3p_del'] > len(seqs['d']):
        print '      d erosions too long (%d)' % (lengths['d_5p_del'] + lengths['d_3p_del'])
        return True
    if lengths['j_5p_del'] > len(seqs['j']):  # NOTE also not possible for the same reason
        print '      j_5p_del erosion too long (%d)' % lengths['j_5p_del']
        return True
    return False

#----------------------------------------------------------------------------------------
def find_tryp_in_joined_seq(tryp_position_in_j, v_seq, vd_insertion, d_seq, dj_insertion, j_seq, j_erosion):
    """ Find the <end> tryptophan in a joined sequence.

    Given local tryptophan position in the j region, figure
    out what position it's at in the final sequence.
    NOTE tryp_position_in_j is the position *before* the j was eroded,
    but this fcn assumes that the j *has* been eroded.
    """
    check_conserved_tryptophan(j_seq, tryp_position_in_j - j_erosion)  # make sure tryp is where it's supposed to be
    length_to_left_of_j = len(v_seq + vd_insertion + d_seq + dj_insertion)
#    print '  finding tryp position as'
#    print '    length_to_left_of_j = len(v_seq + vd_insertion + d_seq + dj_insertion) = %d + %d + %d + %d' % (len(v_seq), len(vd_insertion), len(d_seq), len(dj_insertion))
#    print '    result = tryp_position_in_j - j_erosion + length_to_left_of_j = %d - %d + %d = %d' % (tryp_position_in_j, j_erosion, length_to_left_of_j, tryp_position_in_j - j_erosion + length_to_left_of_j)
    return tryp_position_in_j - j_erosion + length_to_left_of_j

