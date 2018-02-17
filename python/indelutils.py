import random
import numpy
import copy

import utils

# NOTE that even the uncommented functions are kinda not completely right for multiple indels (they're fine as long as we're asking for cyst/tryp/phen positions, and indels are in v or j, which should always be true a.t.m.)

# # ----------------------------------------------------------------------------------------
# ick, needs germline info, and, screw it, I'll just go back to writing the indel-reversed sequence to file
# def reverse_indels(input_seq, indelfo):  # reverse the action of indel reversion
#     if not has_indels(indelfo):
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
#     if not has_indels(indelfo):
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
def get_empty_indel():
    return {'reversed_seq' : '', 'indels' : []}

# ----------------------------------------------------------------------------------------
def has_indels(indelfo):
    return len(indelfo['indels']) > 0

# ----------------------------------------------------------------------------------------
def sign(ifo):
    if ifo['type'] == 'insertion':
        return 1
    elif ifo['type'] == 'deletion':
        return -1
    else:
        assert False

# ----------------------------------------------------------------------------------------
def net_length(indelfo):
    return sum([sign(ifo) * ifo['len'] for ifo in indelfo['indels']])

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
    if not has_indels(indelfo):
        return reinstated_codon_positions

    for indel in indelfo['indels']:
        for region in reinstated_codon_positions:
            reinstated_codon_positions[region] = adjust_single_position_for_reinstated_indels(indel, reinstated_codon_positions[region])
    return reinstated_codon_positions

# ----------------------------------------------------------------------------------------
def get_regional_bounds_with_indels_reinstated(line, iseq):
    indelfo = line['indelfos'][iseq]
    regional_bounds = copy.deepcopy(line['regional_bounds'])
    if not has_indels(indelfo):
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
    qr_seqs = {r : inseq[rbounds[r][0] : rbounds[r][1]] for r in utils.regions}
    return qr_seqs

# ----------------------------------------------------------------------------------------
def add_insertion(indelfo, seq, pos, length, debug=False):
    # should probably be moved inside of add_single_indel()
    """ insert a random sequence with <length> beginning at <pos> """
    inserted_sequence = ''
    for ipos in range(length):
        inuke = random.randint(0, len(utils.nukes) - 1)  # inclusive
        inserted_sequence += utils.nukes[inuke]
    return_seq = seq[ : pos] + inserted_sequence + seq[pos : ]
    indelfo['indels'].append({'type' : 'insertion', 'pos' : pos, 'len' : length, 'seqstr' : inserted_sequence})
    if debug:
        print '          inserting %s at %d' % (inserted_sequence, pos)
    return return_seq

# ----------------------------------------------------------------------------------------
def add_single_indel(seq, indelfo, mean_length, codon_positions, indel_location=None, pos=None, keep_in_frame=False, debug=False):  # NOTE modifies <indelfo> and <codon_positions>
    # if <pos> is specified we use that, otherwise we use <indel_location> to decide the region of the sequence from which to choose a position
    if pos is None:
        if indel_location is None:  # uniform over entire sequence
            pos = random.randint(5, len(seq) - 6)  # this will actually exclude either before the first index or after the last index. No, I don't care.
        elif indel_location == 'v':  # within the meat of the v
            pos = random.randint(5, codon_positions['v'])
        elif indel_location == 'cdr3':  # inside cdr3
            pos = random.randint(codon_positions['v'], codon_positions['j'])
        else:
            assert False

    length = numpy.random.geometric(1. / mean_length)
    if keep_in_frame:
        itry = 0
        while length % 3 != 0:
            length = numpy.random.geometric(1. / mean_length)
            itry += 1
            if itry > 99:
                raise Exception('tried too many times to get in-frame indel length')

    if numpy.random.uniform(0, 1) < 0.5:  # fifty-fifty chance of insertion and deletion
        new_seq = add_insertion(indelfo, seq, pos, length, debug=debug)
    else:
        deleted_seq = seq[ : pos] + seq[pos + length : ]  # delete <length> bases beginning with <pos>
        indelfo['indels'].append({'type' : 'deletion', 'pos' : pos, 'len' : length, 'seqstr' : seq[pos : pos + length]})
        if debug:
            print '          deleting %d bases at %d' % (length, pos)
        new_seq = deleted_seq

    for region in codon_positions:
        if pos < codon_positions[region]:  # not sure that this is right if the indel is actually in the codon, but who fucking cares, right?
            codon_positions[region] += sign(indelfo['indels'][-1]) * length
    assert utils.codon_unmutated('cyst', new_seq, codon_positions['v'], debug=True)

    return new_seq

# ----------------------------------------------------------------------------------------
def get_indelfo_from_cigar(cigarstr, qrseq, glseq, gene):
    cigars = re.findall('[0-9][0-9]*[A-Z]', cigarstr)  # split cigar string into its parts
    cigars = [(cstr[-1], int(cstr[:-1])) for cstr in cigars]  # split each part into the code and the length

    codestr = ''
    qpos = 0  # position within query sequence
    indelfo = indelutils.get_empty_indel()  # replacement_seq: query seq with insertions removed and germline bases inserted at the position of deletions
    tmp_indices = []
    for code, length in cigars:
        codestr += length * code
        if code == 'I':  # advance qr seq but not gl seq
            indelfo['indels'].append({'type' : 'insertion', 'pos' : qpos, 'len' : length, 'seqstr' : ''})  # insertion begins at <pos>
            tmp_indices += [len(indelfo['indels']) - 1  for _ in range(length)]  # indel index corresponding to this position in the alignment
        elif code == 'D':  # advance qr seq but not gl seq
            indelfo['indels'].append({'type' : 'deletion', 'pos' : qpos, 'len' : length, 'seqstr' : ''})  # first deleted base is <pos> (well, first base which is in the position of the first deleted base)
            tmp_indices += [len(indelfo['indels']) - 1  for _ in range(length)]  # indel index corresponding to this position in the alignment
        else:
            tmp_indices += [None  for _ in range(length)]  # indel index corresponding to this position in the alignment
        qpos += length

    qrprintstr, glprintstr = '', ''
    iqr, igl = 0, 0
    for icode in range(len(codestr)):
        code = codestr[icode]
        if code == 'M':
            qrbase = qrseq[iqr]
            if qrbase != glseq[igl]:
                qrbase = utils.color('red', qrbase)
            qrprintstr += qrbase
            glprintstr += glseq[igl]
            indelfo['reversed_seq'] += qrseq[iqr]  # add the base to the overall sequence with all indels reversed
        elif code == 'S':
            continue
        elif code == 'I':
            qrprintstr += utils.color('light_blue', qrseq[iqr])
            glprintstr += utils.color('light_blue', '*')
            indelfo['indels'][tmp_indices[icode]]['seqstr'] += qrseq[iqr]  # and to the sequence of just this indel
            igl -= 1
        elif code == 'D':
            qrprintstr += utils.color('light_blue', '*')
            glprintstr += utils.color('light_blue', glseq[igl])
            indelfo['reversed_seq'] += glseq[igl]  # add the base to the overall sequence with all indels reversed
            indelfo['indels'][tmp_indices[icode]]['seqstr'] += glseq[igl]  # and to the sequence of just this indel
            iqr -= 1
        else:
            raise Exception('unhandled code %s' % code)

        iqr += 1
        igl += 1

    dbg_str_list = ['          %20s %s' % (gene, glprintstr),
                    '          %20s %s' % ('query', qrprintstr)]
    for idl in indelfo['indels']:
        dbg_str_list.append('          %10s: %d bases at %d (%s)' % (idl['type'], idl['len'], idl['pos'], idl['seqstr']))
    indelfo['dbg_str'] = '\n'.join(dbg_str_list)

    return indelfo
