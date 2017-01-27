import sys

import utils

# ----------------------------------------------------------------------------------------
def process_position(original, final):  # optimizing this, and to a lesser extent get_query_line(), would speed up utils.print_reco_event() significantly
    # if original not in utils.expected_characters or final not in utils.expected_characters:
    #     raise Exception('one of %s %s not among expected characters' % (original, final))
    if original in utils.ambiguous_bases or final in utils.ambiguous_bases:
        return final

    if original != final:
        return utils.color('red', final)

    return final

# ----------------------------------------------------------------------------------------
def get_query_line(qseq, line, lengths, glseqs, indelfo=None):  # NOTE do not, on pain of death, modify <line>
    # build up the query sequence line, including colors for mutations and conserved codons
    j_right_extra = 0  # portion of query sequence to right of end of the j match
    n_inserted = 0
    qrseqlist = []
    if indelfo is not None:
        lastfo = indelfo['indels'][-1]  # if the "last" (arbitrary but necessary ordering) indel starts here
    for inuke in range(len(qseq)):
        # if we're at the position that the insertion started at (before we removed it)
        if indelfo is not None and lastfo['type'] == 'insertion':
            if inuke == lastfo['pos']:
                qrseqlist.append(lastfo['seqstr'])  # put the insertion back into the query sequence
                n_inserted += len(lastfo['seqstr'])
        if indelfo is not None and lastfo['type'] == 'deletion':
            if inuke - lastfo['pos'] >= 0 and inuke - lastfo['pos'] < lastfo['len']:  # if we're within the bases that we added to make up for the deletionlen
                qrseqlist.append('*')  # gets blue'd later on, so it can happen at the same time as the germline lines
                continue

        new_nuke = ''
        key = None
        ilocal = inuke
        if indelfo is not None:
            ilocal += n_inserted
        if ilocal < len(line['fv_insertion']):  # haven't got to start of v match yet, so just add on the query seq nuke
            pass
        else:
            ilocal -= len(line['fv_insertion'])
            if ilocal < lengths['v']:
                key = 'v'
            else:
                ilocal -= lengths['v']
                if ilocal < len(line['vd_insertion']):
                    key = 'vd_insertion'
                else:
                    ilocal -= len(line['vd_insertion'])
                    if ilocal < lengths['d']:
                        key = 'd'
                    else:
                        ilocal -= lengths['d']
                        if ilocal < len(line['dj_insertion']):
                            key = 'dj_insertion'
                        else:
                            ilocal -= len(line['dj_insertion'])
                            if ilocal < lengths['j']:
                                key = 'j'
                            else:
                                j_right_extra += 1

        if key is None:
            original = qseq[inuke]  # dummy value
        else:
            original = glseqs[key][ilocal] if key in glseqs else line[key][ilocal]
        new_nuke = process_position(original, qseq[inuke])

        for region, pos in line['codon_positions'].items():  # reverse video for the conserved codon positions
            if inuke >= pos and inuke < pos + 3:
                new_nuke = utils.color('reverse_video', new_nuke)  #'\033[7m' + new_nuke + '\033[0m'  # not sure why the hell I wasn't using color() here

        qrseqlist.append(new_nuke)

# ----------------------------------------------------------------------------------------
    SQ = ''.join(qrseqlist)
    for color_code in utils.Colors.values():
        SQ = SQ.replace(color_code, '')
    qrseqlist = list(SQ)
# ----------------------------------------------------------------------------------------
    return qrseqlist, j_right_extra

# ----------------------------------------------------------------------------------------
def handle_no_space(line, glseqs, qrseqlist):  # NOTE do not, on pain of death, modify <line>
    # if there isn't enough space for dots in the vj line, we add some blue dashes to everybody so things fit (very rare in heavy chain rearrangements, but pretty common in light chain)
    interior_length = len(line['vd_insertion']) + len(glseqs['d']) + len(line['dj_insertion'])  # length of the portion of the vj line that is normally taken up by dots (and spaces)
    if line['v_3p_del'] + line['j_5p_del'] > interior_length:  # not enough space
        v_3p_del_str = '.' + str(line['v_3p_del']) + '.'
        j_5p_del_str = '.' + str(line['j_5p_del']) + '.'
        extra_space_because_of_fixed_nospace = max(0, interior_length - len(v_3p_del_str + j_5p_del_str))

        gap_insertion_point = len(line['fv_insertion'] + glseqs['v'])
        gaps_to_add = len(v_3p_del_str + j_5p_del_str) - interior_length
        qrseqlist = qrseqlist[:gap_insertion_point] + gaps_to_add * [utils.color('blue', '-'), ] + qrseqlist[gap_insertion_point:]
    else:
        v_3p_del_str = '.' * line['v_3p_del']
        j_5p_del_str = '.' * line['j_5p_del']
        gap_insertion_point = None
        gaps_to_add = 0
        extra_space_because_of_fixed_nospace = 0

    return qrseqlist, gap_insertion_point, utils.color('blue', '-') * gaps_to_add, v_3p_del_str, j_5p_del_str, extra_space_because_of_fixed_nospace

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
    def use_stars(index, ifo):
        if ifo['type'] == 'insertion': # for insertions, query sequence should *not* have stars
            return not is_qr(index)
        elif ifo['type'] == 'deletion':
            return is_qr(index)
        else:
            assert False
    def reinstate(seq, ifo, stars=False):
        indelstr = ifo['seqstr']
        if seq[ifo['pos']] not in utils.nukes + utils.ambiguous_bases:  # if this bit of the sequences is spaces, dots, or dashes, then we only want to insert spaces (note that this adds some arbitrariness on boundaries as to who gets the actual inserted string)
            indelstr = ' ' * len(ifo['seqstr'])
        elif stars:
            indelstr = '*' * len(ifo['seqstr'])

        return seq[ : ifo['pos']] + indelstr + seq[ifo['pos'] : ]

    for ifo in reversed(indels['indels']):
        outstrs = [reinstate(outstrs[i], ifo, stars=use_stars(i, ifo)) for i in range(len(outstrs))]

    return outstrs

# ----------------------------------------------------------------------------------------
def color_query_seq(outstrs):
    # <outstrs> convention: [indels, d, vj, query]
    qseq = outstrs[-1]
    qseq = qseq[:35] + 'T' + qseq[36:]
    qrseqlist = list(qseq)
    if len(qrseqlist) != len(qseq):
        print qrseqlist
        print qseq
        assert False
    for inuke in range(len(qseq)):
        if '*' in ''.join([ostr[inuke] for ostr in outstrs]):  # if any of the four have a star at this position (i.e. if we're in an shm insertion or shm deletion)
            continue
        glchars = [ostr[inuke] for ostr in outstrs[:3] if ostr[inuke] in utils.alphabet]
        if len(glchars) > 1:
            raise Exception('more than one germline line has an alphabet character at %d: %s' % (inuke, glchars))
        if qseq[inuke] != glchars[0]:
            qrseqlist[inuke] = utils.color('red', qrseqlist[inuke])
    return outstrs[:3] + [''.join(qrseqlist)]
