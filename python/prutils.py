import copy
import sys

import utils

# ----------------------------------------------------------------------------------------
def get_uid_str(line, iseq, seed_uid):
    uid_width = max([len(uid) for uid in line['unique_ids']])
    fstr = '%' + str(uid_width) + 's'
    uidstr = fstr % line['unique_ids'][iseq]
    if seed_uid is not None and line['unique_ids'][iseq] == seed_uid:
        uidstr = utils.color('red', uidstr)
    return uidstr

# ----------------------------------------------------------------------------------------
def indel_shenanigans(outstrs, indels):  # NOTE similar to/overlaps with get_seq_with_indels_reinstated()
    # <outstrs> convention: [indels, d, vj, query]
    def is_qr(index):
        return index == 3
    def use_stars(ifo, index):
        if ifo['type'] == 'insertion': # for insertions, query sequence should *not* have stars
            return not is_qr(index)
        elif ifo['type'] == 'deletion':
            return is_qr(index)
        else:
            assert False
    def reinstate(seq, ifo, istr):
        indelstr = ifo['seqstr']
        if seq[ifo['pos']] not in utils.nukes + utils.ambiguous_bases:  # if this bit of the sequences is spaces, dots, or dashes, then we only want to insert spaces (note that this adds some arbitrariness on boundaries as to who gets the actual inserted string)
            indelstr = ' ' * len(ifo['seqstr'])
        elif use_stars(ifo, istr):
            indelstr = '*' * len(ifo['seqstr'])

        if ifo['type'] == 'deletion':
            return seq[ : ifo['pos']] + indelstr + seq[ifo['pos'] + ifo['len'] : ]
        else:
            return seq[ : ifo['pos']] + indelstr + seq[ifo['pos'] : ]

    for ifo in reversed(indels['indels']):
        outstrs = [reinstate(outstrs[istr], ifo, istr) for istr in range(len(outstrs))]

    return outstrs

# ----------------------------------------------------------------------------------------
def add_colors(outstrs, line):  # NOTE do *not* modify <line>
    # <outstrs> convention: [indels, d, vj, query]
    bluechars = utils.ambiguous_bases + ['*', '-']

    def ismuted(ch1, ch2):
        if ch1 in bluechars or ch2 in bluechars:
            return False
        if ch1 == ch2:
            return False
        return True

    # first color mutated bases and conserved codons in the query sequence
    codon_positions = [p for cpos in line['codon_positions'].values() for p in range(cpos, cpos + 3)]  # *all* the positions in both the codons
    qrseqlist = list(outstrs[-1])
    ipos = 0  # position in real (alphabetical) query sequence
    for inuke in range(len(qrseqlist)):
        if '*' in ''.join([outstrs[i][inuke] for i in range(3)]):  # if any of the germline lines have a star at this position, i.e. if we're in an shm insertion (if the query line has a star, it's an shm deletion, i.e. the star's position was actually there as a base in the hmm)
            continue
        glchars = [ostr[inuke] for ostr in outstrs[:3] if ostr[inuke] in utils.alphabet]
        if len(glchars) == 0:  # everybody's spaces, dashes, or dots (I think those are the only possibilities...)
            continue
        if len(glchars) > 1:
            raise Exception('more than one germline line has an alphabet character at %d: %s' % (inuke, glchars))
        if ismuted(qrseqlist[inuke], glchars[0]):
            qrseqlist[inuke] = utils.color('red', qrseqlist[inuke])
        if ipos in codon_positions:
            qrseqlist[inuke] = utils.color('reverse_video', qrseqlist[inuke])
        ipos += 1

    outstrs = [outstrs[i] for i in range(3)] + [''.join(qrseqlist)]

    # then color the blues in everybody
    for istr in range(len(outstrs)):
        if len(filter((bluechars).__contains__, outstrs[istr])) == 0:
            continue
        oslist = list(outstrs[istr])
        oslist = [utils.color('light_blue', ochar) if ochar in bluechars else ochar for ochar in oslist]
        outstrs[istr] = ''.join(oslist)

    return outstrs

# ----------------------------------------------------------------------------------------
def print_seq_in_reco_event(germlines, original_line, iseq, extra_str='', label='', one_line=False, seed_uid=None, check_line_integrity=False):
    """
    Print ascii summary of recombination event and mutation.
    If <one_line>, then skip the germline lines, and only print the final_seq line.
    """
    line = original_line
    if check_line_integrity:  # it's very important not to modify <line> -- this lets you verify that you aren't
        line = copy.deepcopy(original_line)  # copy that we can modify without changing <line>

    lengths = {r : line['lengths'][r] for r in utils.regions}  # copy so that we don't have to modify <line>
    glseqs = {r : line[r + '_gl_seq'] for r in utils.regions}  # copy so that we don't have to modify <line>

    # don't print a million dots if left-side v deletion is really big
    v_5p_del_str = '.'*line['v_5p_del']
    if line['v_5p_del'] > 50:
        v_5p_del_str = '...' + str(line['v_5p_del']) + '...'

    # if there isn't enough space for dots in the vj line, we add some dashes to everybody so things fit (very rare in heavy chain rearrangements, but pretty common in light chain)
    interior_length = len(line['vd_insertion']) + len(glseqs['d']) + len(line['dj_insertion'])  # length of the portion of the vj line that is normally taken up by dots (and spaces)
    if line['v_3p_del'] + line['j_5p_del'] > interior_length:  # not enough space
        v_3p_del_str = '.' + str(line['v_3p_del']) + '.'
        j_5p_del_str = '.' + str(line['j_5p_del']) + '.'
        gap_insertion_point = len(line['fv_insertion'] + glseqs['v'])
        gapstr = '-' * (len(v_3p_del_str + j_5p_del_str) - interior_length)
        extra_space_because_of_fixed_nospace = max(0, interior_length - len(v_3p_del_str + j_5p_del_str))
    else:
        v_3p_del_str = '.' * line['v_3p_del']
        j_5p_del_str = '.' * line['j_5p_del']
        gap_insertion_point = None
        gapstr = ''
        extra_space_because_of_fixed_nospace = 0

    eroded_seqs_dots = {
        'v' : glseqs['v'] + v_3p_del_str,
        'd' : '.'*line['d_5p_del'] + glseqs['d'] + '.'*line['d_3p_del'],
        'j' : j_5p_del_str + glseqs['j'] + '.'*line['j_3p_del'],
    }

    # build the three germline lines
    insert_line = ' ' * (len(line['fv_insertion']) + lengths['v'] + len(v_5p_del_str)) \
                  + line['vd_insertion'] + ' ' * lengths['d'] + line['dj_insertion'] \
                  + ' ' * (lengths['j'] + line['j_3p_del'] + len(line['jf_insertion']))
    germline_d_start = len(line['fv_insertion']) + lengths['v'] + len(line['vd_insertion']) - line['d_5p_del']
    germline_d_end = germline_d_start + len(germlines['d'][line['d_gene']])
    d_line = ' ' * (germline_d_start + len(v_5p_del_str)) \
             + eroded_seqs_dots['d'] \
             + ' ' * (len(glseqs['j']) + len(line['dj_insertion']) - line['d_3p_del'] + line['j_3p_del'] + len(line['jf_insertion']))
    germline_v_end = len(line['fv_insertion']) + len(glseqs['v']) + line['v_3p_del'] - 1  # position in the query sequence at which we find the last base of the v match. NOTE we subtract off the v_5p_del because we're *not* adding dots for that deletion (it's just too long)
    germline_j_start = germline_d_end + 1 - line['d_3p_del'] + len(line['dj_insertion']) - line['j_5p_del']
    vj_line = ' ' * len(line['fv_insertion']) + v_5p_del_str + eroded_seqs_dots['v'] + '.' * extra_space_because_of_fixed_nospace \
              + ' ' * (germline_j_start - germline_v_end - 2) + eroded_seqs_dots['j'] + ' ' * len(line['jf_insertion'])
    # and the query line
    qrseq_line = ' ' * len(v_5p_del_str) + line['seqs'][iseq] + ' ' * line['j_3p_del']

    if gap_insertion_point is not None:  # <gap_insertion_point> point is only right here as long as there's no colors in these lines... but there usually almost probably always aren't
        qrseq_line = qrseq_line[:gap_insertion_point] + gapstr + qrseq_line[gap_insertion_point:]
        insert_line = insert_line[:gap_insertion_point] + gapstr + insert_line[gap_insertion_point:]
        d_line = d_line[:gap_insertion_point] + gapstr + d_line[gap_insertion_point:]

    chain = utils.get_chain(line['v_gene'])
    if chain != 'h':
        assert lengths['d'] == 0 and len(line['vd_insertion']) == 0

    outstrs = [insert_line, d_line, vj_line, qrseq_line]
    if check_line_integrity and len(set([len(ostr) for ostr in outstrs])) > 1:
        raise Exception('outstrs not all the same length %s' % [len(ostr) for ostr in outstrs])
    outstrs = indel_shenanigans(outstrs, line['indelfos'][iseq])
    outstrs = add_colors(outstrs, line)
    if check_line_integrity and len(set([utils.len_excluding_colors(ostr) for ostr in outstrs])) > 1:
        raise Exception('outstrs not all the same length %s' % [utils.len_excluding_colors(ostr) for ostr in outstrs])
    suffixes = ['insert%s\n'       % ('s' if chain == 'h' else ''),
                '%s\n'             % (utils.color_gene(line['d_gene'])),
                '%s %s\n'          % (utils.color_gene(line['v_gene']), utils.color_gene(line['j_gene'])),
                '%s   %4.2f mut\n' % (get_uid_str(line, iseq, seed_uid), line['mut_freqs'][iseq])]
    outstrs = ['%s%s   %s' % (extra_str, ostr, suf) for ostr, suf in zip(outstrs, suffixes)]

    if label != '':
        offset = max(0, len(extra_str) - 2)  # skootch <label> this many positions leftward into <extra_str>
        outstrs[0] = outstrs[0][ : offset] + label + outstrs[0][utils.len_excluding_colors(label) + offset : ]  # NOTE this *replaces* the bases in <extra_str> with <label>, which is only fine if they're spaces

    if one_line:
        outstrs = outstrs[-1:]  # remove all except the query seq line
    elif chain != 'h':
        outstrs.pop(1)  # remove the d germline line

    print ''.join(outstrs),

    if check_line_integrity:
        if set(line.keys()) != set(original_line.keys()):
            raise Exception('ack 1')
        for k in line:
            if line[k] != original_line[k]:
                print 'key %s differs:\n  %s\n  %s ' % (k, line[k], original_line[k])
                raise Exception('')
