import sys
import string
import re
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
        if pos < codon_positions[region]:  # this isn\'t right if the indel is actually in the codon, but in that case we just let the messed up codon through below
            codon_positions[region] += sign(indelfo['indels'][-1]) * length
    if not utils.codon_unmutated('cyst', new_seq, codon_positions['v']):
        print '  adding indel within %s codon' % 'cyst'

    return new_seq

# ----------------------------------------------------------------------------------------
def color_cigar(cigarstr):
    return ''.join([utils.color('bold', utils.color('blue', c)) if c in 'MIDS' else c for c in cigarstr])

# ----------------------------------------------------------------------------------------
def split_cigarstr(cstr):
    assert len(cstr) > 0
    code = cstr[-1]
    if code not in 'MIDS':
        raise Exception('unhandled cigar code %s' % code)
    lstr = cstr[:-1] if len(cstr) > 1 else '1'  # stupid vsearch doesn't write the 1 (e.g. 'D' instead of '1D')
    return code, int(lstr)

# ----------------------------------------------------------------------------------------
def get_indelfo_from_cigar(cigarstr, full_qrseq, qrbounds, full_glseq, glbounds, gene, vsearch_conventions=False, debug=False):
    if debug:
        print '  initial:'
        print '    %s' % color_cigar(cigarstr)
        print '    qr %3d %3d %s' % (qrbounds[0], qrbounds[1], full_qrseq)
        print '    gl %3d %3d %s' % (glbounds[0], glbounds[1], full_glseq)

    cigars = [split_cigarstr(cstr) for cstr in re.findall('[0-9]*[A-Z]', cigarstr)]  # split cigar string into its parts, then split each part into the code and the length 
    if vsearch_conventions:
        assert utils.get_region(gene) == 'v'  # would need to be generalized
        cigars = [(code.translate(string.maketrans('ID', 'DI')), length) for code, length in cigars]  # vsearch reverses what's the query and what's the target/gene/whathaveyou compared to what ig-sw does
        for iend in [0, -1]:
            if cigars[iend][0] == 'I':  # qr extends beyond gl: ig-sw calls these soft-clips, vsearch calls them insertions
                cigars[iend] = ('S', cigars[iend][1])
            elif cigars[iend][0] == 'D':  # gl goes past qr: ig-sw just calls them not part of the alignment, vsearch calls them deletions
                cigars.pop(iend)
    cigars = [(code, length) for code, length in cigars if code != 'S']  # remove soft-clipping
    cigarstr = ''.join(['%d%s' % (l, c) for c, l in cigars])
    qrseq = full_qrseq[qrbounds[0] : qrbounds[1]]  # ...and trim qrseq and glseq
    glseq = full_glseq[glbounds[0] : glbounds[1]]

    if debug:
        print '  parsed:'
        print '    %s' % color_cigar(cigarstr)
        print '    %s' % '   '.join(['%s %d' % (c, l) for c, l in cigars])
        print '    qr %s' % qrseq
        print '    gl %s' % glseq

    # check consistency between cigar and qr/gl seqs
    for seqtype, tmpseq, tmpcode in (('qr', qrseq, 'D'), ('gl', glseq, 'I')):
        cigar_len = sum([length for code, length in cigars if code != tmpcode])
        if cigar_len != len(tmpseq):
            raise Exception('cigar length %d doesn\'t match %s seq length %d' % (cigar_len, seqtype, len(tmpseq)))

    indelfo = get_empty_indel()  # replacement_seq: query seq with insertions removed and germline bases inserted at the position of deletions
    # TODO should probably also ignore indels on either end (I think only relevant for vsearch)
    if 'I' not in cigarstr and 'D' not in cigarstr:  # has to happen after we've changed from vsearch conventions
        if debug:
            print '  no indels'
        return indelfo

    # add each indel to <indelfo['indels']>, and build <codestr> and <tmp_indices> to keep track of what's going on at each position
    codestr = ''.join([length * code for code, length in cigars])  # each position is cigar code corresponding to that position in the alignment
    qpos = 0  # position within query sequence
    tmp_indices = []  # integer for each position in the alignment, giving the index of the indel that we're within (None if we're not in an indel)
    if debug:
        print '      code  length'
    for code, length in cigars:
        if debug:
            print '        %s     %3d' % (code, length)
        if code == 'I':  # advance qr seq but not gl seq
            indelfo['indels'].append({'type' : 'insertion', 'pos' : qpos, 'len' : length, 'seqstr' : []})  # insertion begins at <pos> (note that 'seqstr' later on gets converted from a list to a string)
            tmp_indices += [len(indelfo['indels']) - 1  for _ in range(length)]  # indel index corresponding to this position in the alignment
        elif code == 'D':  # advance qr seq but not gl seq
            indelfo['indels'].append({'type' : 'deletion', 'pos' : qpos, 'len' : length, 'seqstr' : []})  # first deleted base is <pos> (well, first base which is in the position of the first deleted base)
            tmp_indices += [len(indelfo['indels']) - 1  for _ in range(length)]  # indel index corresponding to this position in the alignment
        else:
            tmp_indices += [None  for _ in range(length)]  # indel index corresponding to this position in the alignment
        qpos += length

    if debug:
        print '      %s  codestr' % ''.join([c if c not in 'ID' else utils.color('blue', c) for c in codestr])
        print '      %s  indel index' % ''.join([str(ti if ti is not None else ' ') for ti in tmp_indices])

    # then construct the dbg strings, indel-reversed input sequence, and 'seqstr' entries in indelfo
    qrprintstr, glprintstr, reversed_seq = [], [], []
    iqr, igl = 0, 0
    for icode in range(len(codestr)):
        code = codestr[icode]
        if code == 'M':
            qrbase = qrseq[iqr]
            if qrbase != glseq[igl]:
                qrbase = utils.color('red', qrbase)
            qrprintstr.append(qrbase)
            glprintstr.append(glseq[igl])
            reversed_seq.append(qrseq[iqr])  # add the base to the overall sequence with all indels reversed
        elif code == 'S':
            continue
        elif code == 'I':
            qrprintstr.append(utils.color('light_blue', qrseq[iqr]))
            glprintstr.append(utils.color('light_blue', '*'))
            indelfo['indels'][tmp_indices[icode]]['seqstr'].append(qrseq[iqr])  # and to the sequence of just this indel
            igl -= 1
        elif code == 'D':
            qrprintstr.append(utils.color('light_blue', '*'))
            glprintstr.append(utils.color('light_blue', glseq[igl]))
            reversed_seq.append(glseq[igl])  # add the base to the overall sequence with all indels reversed
            indelfo['indels'][tmp_indices[icode]]['seqstr'].append(glseq[igl])  # and to the sequence of just this indel
            iqr -= 1
        else:
            raise Exception('unhandled cigar code %s' % code)

        iqr += 1
        igl += 1

    # convert character lists to strings (indels are rare enough that this probably isn't that much faster, but it just feels wrong not to)
    qrprintstr = ''.join(qrprintstr)
    glprintstr = ''.join(glprintstr)
    indelfo['reversed_seq'] = ''.join(reversed_seq)
    indelfo['reversed_seq'] = full_qrseq[ : qrbounds[0]] + indelfo['reversed_seq'] + full_qrseq[qrbounds[1] : ]  # add to reversed seq the bits to left and right of the aligned region TODO
    for ifo in indelfo['indels']:
        ifo['seqstr'] = ''.join(ifo['seqstr'])

    # make the dbg str for indelfo
    gwidth = str(len(gene))  # doesn't account for color abbreviation, but oh well
    dbg_str_list = [('%' + gwidth + 's  %s') % (utils.color_gene(gene, width=int(gwidth), leftpad=True), glprintstr),
                    ('%' + gwidth + 's  %s') % ('query', qrprintstr)]
    for idl in indelfo['indels']:
        dbg_str_list.append('%10s: %d base%s at %d (%s)' % (idl['type'], idl['len'], utils.plural(idl['len']), idl['pos'], idl['seqstr']))
    indelfo['dbg_str'] = '\n'.join(dbg_str_list)

    if debug:
        print utils.pad_lines(indelfo['dbg_str'], 0)

    return indelfo

# ----------------------------------------------------------------------------------------
def pad_indel_info(indelfo, leftstr, rightstr):
    indelfo['reversed_seq'] = leftstr + indelfo['reversed_seq'] + rightstr

    # dbg_indel_lines = indelfo['dbg_str'].split('\n')[2:]  # i.e. just the lines describing each indel, not the first two lines that have the qrprintstr and glprintstr
    # if len(dbg_indel_lines) != len(indelfo['indels']):
    #     print '  %s wrong number of dbg indel lines %d %d' % (utils.color('red', 'error'), len(dbg_indel_lines), len(indelfo['indels']))

    # dbg_indel_lines = indelfo['dbg_str'].split('\n')
    # if len(dbg_indel_lines) != len(indelfo['indels']) - 2:
    #     print '  %s wrong number of dbg indel lines %d %d' % (utils.color('red', 'error'), len(dbg_indel_lines) - 2, len(indelfo['indels']))
    # for iline in range(2, len(dbg_indel_lines)):
    #     def len_str( = ' %d ' % indelfo['indels'][iline]['len']
    #     assert dbg_indel_lines[iline].count(len_str) == 1
    #     dbg_indel_lines[iline] = dbg_indel_lines[iline].replace(len_str, )

    # aw, screw it, that ^ is hard, just remove the dbg_str
    if 'dbg_str' in indelfo:
        del indelfo['dbg_str']

    for indel in indelfo['indels']:
        indel['pos'] += len(leftstr)
