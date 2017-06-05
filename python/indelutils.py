import copy

# NOTE that even the uncommented functions are kinda not completely right for multiple indels (they're fine as long as we're asking for cyst/tryp/phen positions, and indels are in v or j, which should always be true a.t.m.)

# # ----------------------------------------------------------------------------------------
# ick, needs germline info, and, screw it, I'll just go back to writing the indel-reversed sequence to file
# def reverse_indels(input_seq, indelfo):  # reverse the action of indel reversion
#     if indelfo['reversed_seq'] == '':
#         return input_seq

#     for indel in indelfo['indels']:
#         if indel['type'] == 'insertion':
#             return_seq = return_seq[ : indel['pos']] + indel['seqstr'] + return_seq[indel['pos'] : ]
#         elif indel['type'] == 'deletion':
#             return_seq = return_seq[ : indel['pos']] + return_seq[indel['pos'] + indel['len'] : ]
#         else:
#             assert False

#     return return_seq

# # ----------------------------------------------------------------------------------------
#
# This is *not* correct (the 'pos' for each indel is needs to be different -- see waterer.py).
# Nevertheless, leaving this here to remind you not to reimplement it.
#
# def get_seq_with_indels_reinstated(line, iseq=0):  # reverse the action of indel reversion
#     indelfo = line['indelfos'][iseq]
#     return_seq = line['seqs'][iseq]
#     if indelfo['reversed_seq'] == '':
#         return return_seq
#
#     # for indel in reversed(indelfo['indels']):
#     for indel in indelfo['indels']:
#         if indel['type'] == 'insertion':
#             return_seq = return_seq[ : indel['pos']] + indel['seqstr'] + return_seq[indel['pos'] : ]
#         elif indel['type'] == 'deletion':
#             excision = return_seq[indel['pos'] : indel['pos'] + indel['len']]
#             if excision != indel['seqstr']:
#                 # raise Exception('ack %s %s' % (excision, indel['seqstr']))
#                 print '%s %s %s %s' % (color('red', 'ack'), ' '.join(line['unique_ids']), excision, indel['seqstr'])
#             return_seq = return_seq[ : indel['pos']] + return_seq[indel['pos'] + indel['len'] : ]
#         else:
#             assert False
#
#     return return_seq

# ----------------------------------------------------------------------------------------
def adjust_single_position_for_reinstated_indels(indel, position):
    if indel['pos'] > position:  # NOTE this just ignores the case where the indel's in the middle of the codon, because, like, screw that I don't want to think about it
        return position
    if indel['type'] == 'insertion':
        return position + indel['len']
    elif indel['type'] == 'deletion':
        return position - indel['len']
    else:
        assert False

# ----------------------------------------------------------------------------------------
def get_codon_positions_with_indels_reinstated(line, iseq, codon_positions):
    # NOTE as long as the indels are reversed, all the sequences have the same codon positions. But as soon as we reinstate the indels, all heck breaks loose.
    indelfo = line['indelfos'][iseq]
    reinstated_codon_positions = copy.deepcopy(codon_positions)
    if indelfo['reversed_seq'] == '':
        return reinstated_codon_positions

    for indel in indelfo['indels']:
        for region in reinstated_codon_positions:
            reinstated_codon_positions[region] = adjust_single_position_for_reinstated_indels(indel, reinstated_codon_positions[region])
    return reinstated_codon_positions

# ----------------------------------------------------------------------------------------
def get_regional_bounds_with_indels_reinstated(line, iseq):
    indelfo = line['indelfos'][iseq]
    regional_bounds = copy.deepcopy(line['regional_bounds'])
    if indelfo['reversed_seq'] == '':
        return regional_bounds

    for indel in indelfo['indels']:
        for region in regional_bounds:
            regional_bounds[region] = (adjust_single_position_for_reinstated_indels(indel, regional_bounds[region][0]),
                                       adjust_single_position_for_reinstated_indels(indel, regional_bounds[region][1]))
    return regional_bounds

# ----------------------------------------------------------------------------------------
def get_qr_seqs_with_indels_reinstated(line, iseq):
    rbounds = get_regional_bounds_with_indels_reinstated(line, iseq)
    # assert line['input_seqs'][iseq] == get_seq_with_indels_reinstated(line, iseq)
    inseq = line['input_seqs'][iseq]
    qr_seqs = {r : inseq[rbounds[r][0] : rbounds[r][1]] for r in regions}
    return qr_seqs
