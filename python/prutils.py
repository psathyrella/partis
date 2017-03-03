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
def indel_shenanigans(outstrs, colors, indels, iseq):  # NOTE similar to/overlaps with get_seq_with_indels_reinstated()
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
    def reinstate(ifo, istr):
        indelstr = ifo['seqstr']
        if outstrs[istr][ifo['pos']] not in utils.nukes + utils.ambiguous_bases:  # if this bit of the sequences is spaces, dots, or dashes, then we only want to insert spaces (note that this adds some arbitrariness on boundaries as to who gets the actual inserted string)
            indelstr = ' ' * len(ifo['seqstr'])
        elif use_stars(ifo, istr):
            indelstr = '*' * len(ifo['seqstr'])

        if ifo['type'] == 'deletion':
            outstrs[istr] = outstrs[istr][ : ifo['pos']] + indelstr + outstrs[istr][ifo['pos'] + ifo['len'] : ]
            colors[istr] = colors[istr][ : ifo['pos']] + [[] for _ in range(len(indelstr))] + colors[istr][ifo['pos'] + ifo['len'] : ]
        else:
            outstrs[istr] = outstrs[istr][ : ifo['pos']] + indelstr + outstrs[istr][ifo['pos'] : ]
            colors[istr] = colors[istr][ : ifo['pos']] + [[] for _ in range(len(indelstr))] + colors[istr][ifo['pos'] : ]
            if iseq > 0 and is_qr(istr):  # the germline lines are printed based on the *first* query sequence in a multi-seq alignment, and it'd be kinda hard to really account for all subsequent sequences' indels, so we compromise by blueing the insertions in the subsequent query sequences (deletions are already blue stars). Note that they still don't line up right.
                for inuke  in range(ifo['pos'], ifo['pos'] + ifo['len']):
                    colors[istr][inuke] += ['light_blue', 'reverse_video']

    for ifo in reversed(indels['indels']):
        for istr in range(len(outstrs)):
            reinstate(ifo, istr)

    return outstrs, colors

# ----------------------------------------------------------------------------------------
def add_colors(outstrs, colors, line):  # NOTE do *not* modify <line>
    # <outstrs> convention: [indels, d, vj, query]
    bluechars = utils.ambiguous_bases + ['*', '-']

    def ismuted(ch1, ch2):
        if ch1 in bluechars or ch2 in bluechars:
            return False
        if ch1 == ch2:
            return False
        return True

    # first add colors for mutated bases and conserved codons in the query sequence
    codon_positions = [p for cpos in line['codon_positions'].values() for p in range(cpos, cpos + 3)]  # *all* the positions in both the codons
    qrseq, qrcols = outstrs[-1], colors[-1]
    ipos = 0  # position in real (alphabetical) query sequence
    for inuke in range(len(qrseq)):
        glchars = ''.join([ostr[inuke] for ostr in outstrs[:3]])
        if '*' in glchars:  # if any of the germline lines have a star at this position, i.e. if we're in an shm insertion (if the query line has a star, it's an shm deletion, i.e. the star's position was actually there as a base in the hmm)
            continue
        alphagl = filter((utils.alphabet).__contains__, glchars)  # glchars that are also alphabet chars
        if len(alphagl) == 0:  # if none of the germline lines have an alphabet base, we're not within the coding region
            if qrseq[inuke] in utils.alphabet:  # fv insertions
                ipos += 1
            continue
        if len(alphagl) > 1:
            raise Exception('more than one germline line has an alphabet character at %d:\n  %s\n  %s\n  %s' % (inuke, glchars[0], glchars[1], glchars[2]))
        if ismuted(qrseq[inuke], alphagl[0]):
            qrcols[inuke].append('red')
        if ipos in codon_positions:
            qrcols[inuke].append('reverse_video')
        ipos += 1

    # then add the blues for everybody
    for istr in range(len(outstrs)):
        if len(filter((bluechars).__contains__, outstrs[istr])) == 0:  # skip anybody that doesn't have any blue chars
            continue
        blue_indices = [inuke for inuke in range(len(outstrs[istr])) if outstrs[istr][inuke] in bluechars]
        for inuke in blue_indices:
            colors[istr][inuke].append('light_blue')

    # and finally apply the colors from <colors>
    for istr in range(len(outstrs)):
        assert len(colors[istr]) == len(outstrs[istr])
        if colors.count([]) == len(colors):  # colorless line
            continue
        oslist = list(outstrs[istr])
        for inuke in range(len(oslist)):
            for col in colors[istr][inuke]:
                oslist[inuke] = utils.color(col, oslist[inuke])
        outstrs[istr] = ''.join(oslist)

    return outstrs

# ----------------------------------------------------------------------------------------
def print_seq_in_reco_event(original_line, iseq, extra_str='', label='', one_line=False, seed_uid=None, check_line_integrity=False):
    """
    Print ascii summary of recombination event and mutation.
    If <one_line>, then skip the germline lines, and only print the final_seq line.
    """
    line = original_line
    if check_line_integrity:  # it's very important not to modify <line> -- this lets you verify that you aren't
        line = copy.deepcopy(original_line)  # copy that we can modify without changing <line>

    delstrs = {d : '.' * line[d + '_del'] for d in utils.all_erosions}  # NOTE len(delstrs[<del>]) is not in general the same as len(line[<del>_del])
    if len(delstrs['v_5p']) > 50:  # don't print a million dots if left-side v deletion is really big
        delstrs['v_5p'] = '.%d.' % len(delstrs['v_5p'])

    # if there isn't enough space for dots in the vj line, we add some dashes to everybody so things fit (rare in heavy chain rearrangements, but pretty common in light chain)
    d_plus_inserts_length = len(line['vd_insertion'] + line['d_gl_seq'] + line['dj_insertion'])
    if line['v_3p_del'] + line['j_5p_del'] > d_plus_inserts_length:  # if dots for v and j interior deletions will be longer than <d_plus_inserts_length>
        delstrs['v_3p'] = '.%d.' % line['v_3p_del']
        delstrs['j_5p'] = '.%d.' % line['j_5p_del']
        gapstr = '-' * (len(delstrs['v_3p'] + delstrs['j_5p']) - d_plus_inserts_length)
        gap_insert_point = len(line['fv_insertion'] + delstrs['v_5p'] + line['v_gl_seq'])  # it doesn't really matter exactly where we put the blue dashes, as long as it's the same place in all four lines, but this is a good spot
        extra_space_because_of_fixed_nospace = max(0, d_plus_inserts_length - len(delstrs['v_3p'] + delstrs['j_5p']))  # if shortening the <delstrs> already over-compensated for the lack of space (i.e., if the number of dashes necessary is zero), then we need to add some dots to the vj line below
    else:
        gapstr = ''
        gap_insert_point = None
        extra_space_because_of_fixed_nospace = 0

    eroded_seqs_dots = {r : delstrs[r + '_5p'] + line[r + '_gl_seq'] + delstrs[r + '_3p'] for r in utils.regions}

    # build the three germline lines
    insert_line = ' ' * (len(line['fv_insertion']) + line['lengths']['v'] + len(delstrs['v_5p'])) \
                  + line['vd_insertion'] + ' ' * line['lengths']['d'] + line['dj_insertion'] \
                  + ' ' * (line['lengths']['j'] + line['j_3p_del'] + len(line['jf_insertion']))
    germline_d_start = len(line['fv_insertion']) + line['lengths']['v'] + len(line['vd_insertion']) - line['d_5p_del']
    germline_d_end = germline_d_start + line['d_5p_del'] + line['lengths']['d'] + line['d_3p_del']
    d_line = ' ' * (germline_d_start + len(delstrs['v_5p'])) \
             + eroded_seqs_dots['d'] \
             + ' ' * (len(line['j_gl_seq']) + len(line['dj_insertion']) - line['d_3p_del'] + line['j_3p_del'] + len(line['jf_insertion']))
    germline_v_end = len(line['fv_insertion']) + len(line['v_gl_seq']) + line['v_3p_del'] - 1  # position in the query sequence at which we find the last base of the v match. NOTE we subtract off the v_5p_del because we're *not* adding dots for that deletion (it's just too long)
    germline_j_start = germline_d_end + 1 - line['d_3p_del'] + len(line['dj_insertion']) - line['j_5p_del']
    vj_line = ' ' * len(line['fv_insertion']) + eroded_seqs_dots['v'] + '.' * extra_space_because_of_fixed_nospace \
              + ' ' * (germline_j_start - germline_v_end - 2) + eroded_seqs_dots['j'] + ' ' * len(line['jf_insertion'])
    # and the query line
    qrseq_line = ' ' * len(delstrs['v_5p']) + line['seqs'][iseq] + ' ' * line['j_3p_del']

    outstrs = [insert_line, d_line, vj_line, qrseq_line]

    if gap_insert_point is not None:
        for istr in [0, 1, 3]:  # everybody except the vj line, which already has the modified interior delstrs above
            outstrs[istr] = outstrs[istr][:gap_insert_point] + gapstr + outstrs[istr][gap_insert_point:]

    if len(set([len(ostr) for ostr in outstrs])) > 1:  # could put this in a bunch of different places, but things're probably most likely to get screwed up either when initally building the four lines, or dealing with the stupid gaps
        raise Exception('outstrs not all the same length %s' % [len(ostr) for ostr in outstrs])

    colors = [[[] for _ in range(len(ostr))] for ostr in outstrs]
    outstrs, colors = indel_shenanigans(outstrs, colors, line['indelfos'][iseq], iseq)
    outstrs = add_colors(outstrs, colors, line)

    suffixes = ['insert%s\n'       % ('s' if 'd' in utils.getregions(utils.get_locus(line['v_gene'])) else ''),
                '%s\n'             % (utils.color_gene(line['d_gene'])),
                '%s %s\n'          % (utils.color_gene(line['v_gene']), utils.color_gene(line['j_gene'])),
                '%s   %4.2f mut\n' % (get_uid_str(line, iseq, seed_uid), line['mut_freqs'][iseq])]
    outstrs = ['%s%s   %s' % (extra_str, ostr, suf) for ostr, suf in zip(outstrs, suffixes)]

    if label != '':
        offset = max(0, len(extra_str) - 2)  # skootch <label> this many positions leftward into <extra_str>
        outstrs[0] = outstrs[0][ : offset] + label + outstrs[0][utils.len_excluding_colors(label) + offset : ]  # NOTE this *replaces* the bases in <extra_str> with <label>, which is only fine if they're spaces

    if one_line:
        outstrs = outstrs[-1:]  # remove all except the query seq line
    elif 'd' not in utils.getregions(utils.get_locus(line['v_gene'])):
        outstrs.pop(1)  # remove the d germline line

    print ''.join(outstrs),

    if check_line_integrity:
        if set(line.keys()) != set(original_line.keys()):
            raise Exception('ack 1')
        for k in line:
            if line[k] != original_line[k]:
                print 'key %s differs:\n  %s\n  %s ' % (k, line[k], original_line[k])
                raise Exception('')
