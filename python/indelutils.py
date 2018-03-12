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
def get_empty_indel():  # TODO update, probably should just be None or something? it at least should be simpler
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
def get_dbg_str(indelfo):
    qrprintstr, glprintstr = [], []
    for ich in range(len(indelfo['qr_gap_seq'])):
        qrb, glb = indelfo['qr_gap_seq'][ich], indelfo['gl_gap_seq'][ich]
        qrcolor, glcolor = None, None
        if qrb in utils.gap_chars or glb in utils.gap_chars:
            qrcolor = 'light_blue'
            glcolor = 'light_blue'
        elif qrb in utils.ambiguous_bases:
            qrcolor = 'light_blue'
        elif glb in utils.ambiguous_bases:
            glcolor = 'light_blue'
        elif qrb != glb:
            qrcolor = 'red'
        qrprintstr.append(utils.color(qrcolor, qrb if qrb not in utils.gap_chars else '*'))  # change it to a start just cause that's what it originally was... at some point should switch to just leaving it whatever gap char it was
        glprintstr.append(utils.color(glcolor, glb if glb not in utils.gap_chars else '*'))
    qrprintstr = indelfo['reversed_seq'][ : indelfo['qrbounds'][0]] + ''.join(qrprintstr)
    glprintstr = ' ' * indelfo['qrbounds'][0] + ''.join(glprintstr)

    gene_str = ''
    gwidth = str(len('query'))
    if indelfo['v_gene'] is not None:
        gene_str = utils.color_gene(indelfo['v_gene'], width=int(gwidth), leftpad=True)
        gwidth = str(utils.len_excluding_colors(gene_str))
    dbg_str_list = [('  %' + gwidth + 's  %s  %s') % (gene_str, glprintstr, utils.color_gene(indelfo['j_gene']) if indelfo['j_gene'] is not None else ''),
                    ('  %' + gwidth + 's  %s') % ('query', qrprintstr)]
    for idl in indelfo['indels']:  # TODO for joint indels (both/all regions) positions are only right if you don't count stars in the query sequence (which might be what I want)
        dbg_str_list.append('%10s: %d base%s at %d (%s)' % (idl['type'], idl['len'], utils.plural(idl['len']), idl['pos'], idl['seqstr']))  # kind of dumb/confusing to subtract from the pos, but otherwise we have to figure out how to print out the usually-different-length non-matched qr and gl bits
    return '\n'.join(dbg_str_list)

# ----------------------------------------------------------------------------------------
def get_reversed_seq(qr_gap_seq, gl_gap_seq, v_5p_del_str, j_3p_del_str):
    reversed_match_seq = [(qrb if qrb not in utils.gap_chars else glb) for qrb, glb in zip(qr_gap_seq, gl_gap_seq) if glb not in utils.gap_chars]
    return v_5p_del_str + ''.join(reversed_match_seq) + j_3p_del_str

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

    # each position is the cigar code corresponding to that position in the alignment
    codestr = ''.join([length * code for code, length in cigars])

    # add each indel to <indelfo['indels']>, and build <tmp_indices> to keep track of what's going on at each position
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
    qr_gap_seq, gl_gap_seq = [], []
    iqr, igl = 0, 0
    for icode in range(len(codestr)):
        code = codestr[icode]
        qrb, glb = qrseq[iqr], glseq[igl]
        if code == 'M':
            qr_gap_seq += [qrb]
            gl_gap_seq += [glb]
        elif code == 'I':
            indelfo['indels'][tmp_indices[icode]]['seqstr'] += [qrb]  # and to the sequence of just this indel
            qr_gap_seq += [qrb]
            gl_gap_seq += ['.']
            igl -= 1
        elif code == 'D':
            indelfo['indels'][tmp_indices[icode]]['seqstr'] += [glb]  # and to the sequence of just this indel
            qr_gap_seq += ['.']
            gl_gap_seq += [glb]
            iqr -= 1
        else:
            raise Exception('unexpected cigar code %s' % code)
        iqr += 1
        igl += 1

    # convert character lists to strings (indels are rare enough that this probably isn't that much faster, but it just feels wrong not to)
    qr_gap_seq = ''.join(qr_gap_seq)
    gl_gap_seq = ''.join(gl_gap_seq)
    for ifo in indelfo['indels']:
        ifo['seqstr'] = ''.join(ifo['seqstr'])

    # at the start of this fcn we trimmed off the "non-matched" bits of the query and germline sequences, so now we have to account for them (it might be nicer to have it all done at once, but this is the way it is, for historical reasons) (where the definition of "non-matched" is a bit fuzzy depending on whether it's vsearch or ig-sw)
    for ifo in indelfo['indels']:
        ifo['pos'] += qrbounds[0]

    # NOTE gapped seqs do _not_ contain the v 5p and j 3p deletions or fv and jf insertions, because this makes it easier to combine indels from different regions later on
    for region in ['v', 'j']:
        indelfo[region + '_gene'] = gene if region == utils.get_region(gene) else None
    indelfo['qr_gap_seq'] = qr_gap_seq
    indelfo['gl_gap_seq'] = gl_gap_seq
    indelfo['qrbounds'] = qrbounds
    indelfo['reversed_seq'] = get_reversed_seq(qr_gap_seq, gl_gap_seq, full_qrseq[ : qrbounds[0]], full_qrseq[qrbounds[1] : ])

    if debug:
        print utils.pad_lines(get_dbg_str(indelfo, gene=gene), 0)

    return indelfo

# ----------------------------------------------------------------------------------------
def pad_indel_info(indelfo, leftstr, rightstr):
    # TODO update for any new keys
    indelfo['reversed_seq'] = leftstr + indelfo['reversed_seq'] + rightstr
    for indel in indelfo['indels']:
        indel['pos'] += len(leftstr)

# ----------------------------------------------------------------------------------------
def check_indelfo_consistency(glfo, line, debug=False):
    for iseq in range(len(line['unique_ids'])):
        check_single_sequence_indels(glfo, line, iseq, debug=debug)

# ----------------------------------------------------------------------------------------
def check_single_sequence_indels(glfo, line, iseq, debug=False):
    # TODO figure out what to do with this fcn

    # add a list of qr_gap_seqs and gl_gap_seqs, which are None if there's no indels
    # swich has_indels() to just checking a new bool in <line>

    indelfo = line['indelfos'][iseq]

    if not has_indels(indelfo):
        return

    # new_indelfo = get_indelfo_from_cigar(indelfo['cigarstr'], line['input_seqs'][iseq], indelfo['qrbounds'], glfo['seqs']['v'][line['v_gene']], indelfo['glbounds'], line['v_gene'])
    # if len(new_indelfo['indels']) != len(indelfo['indels']):
    #     print '%s different lengths %d %d' % (len(new_indelfo['indels']), len(indelfo['indels']))

    for iindel in range(len(indelfo['indels'])):
        ifo = indelfo['indels'][iindel]
        if debug:
            print '  len %d  pos %d  seqstr %s' % (ifo['len'], ifo['pos'], ifo['seqstr'])

# ----------------------------------------------------------------------------------------
        # new_ifo = new_indelfo['indels'][iindel]
        # if new_ifo['seqstr'] != ifo['seqstr']:
        #     print '\n%s inconsistent indel info for %s:' % (utils.color('red', 'error'), ':'.join(line['unique_ids']))
        #     utils.color_mutants(new_ifo['seqstr'], ifo['seqstr'], print_result=True, extra_str='    ', ref_label=ref_label + ' ')
        #     utils.print_reco_event(line)
        # else:
        #     if debug:
        #         print '  %s' % utils.color('green', 'ok')
# ----------------------------------------------------------------------------------------

        if ifo['type'] == 'insertion':
            deleted_str = line['input_seqs'][iseq][ifo['pos'] : ifo['pos'] + ifo['len']]
            ref_label = 'input seq'
        else:
            assert len(line['fv_insertion']) == indelfo['qrbounds'][0]
            gl_pos = ifo['pos'] - len(line['fv_insertion']) + indelfo['glbounds'][0]
            deleted_str = glfo['seqs']['v'][line['v_gene']][gl_pos : gl_pos + ifo['len']]
            ref_label = 'gl seq'
            # gl_pos = ifo['pos'] - len(line['fv_insertion'])
            # deleted_str = line['v_gl_seq'][gl_pos : gl_pos + ifo['len']]
            # deleted_str = line['indel_reversed_seqs'][iseq][ifo['pos'] : ifo['pos'] + ifo['len']]

        if deleted_str != ifo['seqstr']:
            print '%s inconsistent indel info for %s:' % (utils.color('red', 'error'), ':'.join(line['unique_ids']))
            utils.color_mutants(deleted_str, ifo['seqstr'], print_result=True, extra_str='    ', ref_label=ref_label + ' ')
            # utils.print_reco_event(line)

# ----------------------------------------------------------------------------------------
def combine_indels(vfo, jfo, full_qrseq):
    if vfo is None and jfo is None:
        assert False  # shouldn't be possible, since it'd require a d indel?
    elif jfo is None:
        print get_dbg_str(vfo)
        joint_indelfo = copy.deepcopy(vfo)
        joint_indelfo['v_gene'] = vfo['v_gene']
    elif vfo is None:
        print get_dbg_str(jfo)
        joint_indelfo = copy.deepcopy(jfo)
        joint_indelfo['j_gene'] = jfo['j_gene']
    else:
        print 'v'
        print get_dbg_str(vfo)
        print 'j'
        print get_dbg_str(jfo)

        ungapped_qr_central_str = full_qrseq[vfo['qrbounds'][1] : jfo['qrbounds'][0]]
        ungapped_gl_central_str = utils.ambiguous_bases[0] * (jfo['qrbounds'][0] - vfo['qrbounds'][1])  # the bit between the end of the v and the start of the j (we could instead use the d + insert stuff, but I'd rather keep it agnostic about that since all we're really trying to communicate here is where the shm indels are)
        combined_gapped_qr_seq = vfo['qr_gap_seq'] + ungapped_qr_central_str + jfo['qr_gap_seq']  # reminder: still don't contain v 5p or j 3p deletions or fv/jf insertions
        combined_gapped_gl_seq = vfo['gl_gap_seq'] + ungapped_gl_central_str + jfo['gl_gap_seq']
        assert len(combined_gapped_qr_seq) == len(combined_gapped_gl_seq)

        joint_indelfo = get_empty_indel()
        joint_indelfo['v_gene'] = vfo['v_gene']
        joint_indelfo['j_gene'] = jfo['j_gene']
        joint_indelfo['indels'] = copy.deepcopy(vfo['indels']) + copy.deepcopy(jfo['indels'])
        for iindel in range(len(vfo['indels']), len(joint_indelfo['indels'])):  # arg, kind of wish I didn't need this, but it's only really to keep the stupid 'pos' things right, and I wish I didn't even need them any more
            ifo = joint_indelfo['indels'][iindel]
            ifo['pos'] += utils.count_gap_chars(combined_gapped_qr_seq, ifo['pos'])  # NOTE this corresponds to the _new_ query position, i.e. in the reversed seq, _not_ the old one (this bullshit is why I shouldn't have used the 'pos' stuff in the first place)
        joint_indelfo['qr_gap_seq'] = combined_gapped_qr_seq
        joint_indelfo['gl_gap_seq'] = combined_gapped_gl_seq
        joint_indelfo['qrbounds'] = (vfo['qrbounds'][0], jfo['qrbounds'][1])
        joint_indelfo['reversed_seq'] = get_reversed_seq(combined_gapped_qr_seq, combined_gapped_gl_seq, full_qrseq[ : vfo['qrbounds'][0]], full_qrseq[jfo['qrbounds'][1] : ])
        assert 'N' not in joint_indelfo['reversed_seq']  # TODO remove this

    print 'combined'
    print get_dbg_str(joint_indelfo)
    return joint_indelfo
