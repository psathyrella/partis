from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import copy
import sys

from . import indelutils
from . import utils

# ----------------------------------------------------------------------------------------
def get_uid_str(line, iseq, queries_to_emphasize, duplicated_uids=None):
    uid_width = max([len(uid) for uid in line['unique_ids']])
    fstr = '%' + str(uid_width) + 's'
    uidstr = fstr % line['unique_ids'][iseq]
    if queries_to_emphasize is not None and line['unique_ids'][iseq] in queries_to_emphasize:
        uidstr = utils.color('red', uidstr)
    if duplicated_uids is not None and line['unique_ids'][iseq] in duplicated_uids:
        uidstr += ' ' + utils.color('red', 'duplicate: %d' % duplicated_uids[line['unique_ids'][iseq]])
    return uidstr

# ----------------------------------------------------------------------------------------
def check_outsr_lengths(line, outstrs, fix=False, debug=False):
    if len(set([len(ostr) for ostr in outstrs])) == 1:  # could put this in a bunch of different places, but things're probably most likely to get screwed up either when initally building the four lines, or dealing with the stupid gaps
        if debug:
            print('    outstr lengths ok')
        return

    if fix:
        max_len = max([len(ostr) for ostr in outstrs])
        if debug:
            print('      fixing outstr lengths: %s --> %d' % (' '.join(str(len(ostr)) for ostr in outstrs), max_len))
        for istr in range(len(outstrs)):
            if len(outstrs[istr]) < max_len:
                outstrs[istr] += ' ' * (max_len - len(outstrs[istr]))
        return

    print(':'.join(line['unique_ids']))
    for ostr in outstrs:
        print('%s%s%s' % (utils.color('red', 'x'), ostr, utils.color('red', 'x')))
    raise Exception('outstrs not all the same length %s' % [len(ostr) for ostr in outstrs])

# ----------------------------------------------------------------------------------------
def vj_ostr(ostrs):  # it would be nice to add a general fcn incorporating also the fcn below, rather than adding this, but whatever
    return ostrs[2]

# ----------------------------------------------------------------------------------------
def is_qr(index):
    return index == 3

# ----------------------------------------------------------------------------------------
def use_stars(itype, index):
    if itype == 'insertion': # for insertions, query sequence should *not* have stars
        return not is_qr(index)
    elif itype == 'deletion':
        return is_qr(index)
    else:
        assert False

# ----------------------------------------------------------------------------------------
def old_indel_shenanigans(line, iseq, outstrs, colors, debug=False):  # NOTE similar to/overlaps with get_seq_with_indels_reinstated()
    # <outstrs> convention: [inserts, d, vj, query]
    def reinstate(ifo, istr):
        indelstr = ifo['seqstr']

        if outstrs[istr][ifo['pos']] not in utils.nukes + utils.all_ambiguous_bases:  # if this bit of the sequences is spaces, dots, or dashes, then we only want to insert spaces (note that this adds some arbitrariness on boundaries as to who gets the actual inserted string)
            indelstr = ' ' * len(ifo['seqstr'])
        elif use_stars(ifo['type'], istr):
            indelstr = '*' * len(ifo['seqstr'])

        if ifo['type'] == 'deletion':
            outstrs[istr] = outstrs[istr][ : ifo['pos']] + indelstr + outstrs[istr][ifo['pos'] + ifo['len'] : ]
            colors[istr] = colors[istr][ : ifo['pos']] + [[] for _ in range(len(indelstr))] + colors[istr][ifo['pos'] + ifo['len'] : ]
        else:
            outstrs[istr] = outstrs[istr][ : ifo['pos']] + indelstr + outstrs[istr][ifo['pos'] : ]
            colors[istr] = colors[istr][ : ifo['pos']] + [[] for _ in range(len(indelstr))] + colors[istr][ifo['pos'] : ]
            if iseq > 0 and is_qr(istr):  # the germline lines are printed based on the *first* query sequence in a multi-seq alignment, and it'd be kinda hard to really account for all subsequent sequences' indels, so we compromise by blueing the insertions in the subsequent query sequences (deletions are already blue stars). Note that they still don't line up right.
                for inuke  in range(ifo['pos'], ifo['pos'] + ifo['len']):
                    colors[istr][inuke] += ['blue', 'reverse_video']

    ifo = line['indelfos'][iseq]
    for ifo in ifo['indels']:
        for istr in range(len(outstrs)):
            reinstate(ifo, istr)
        check_outsr_lengths(line, outstrs)

    return outstrs, colors

# ----------------------------------------------------------------------------------------
def indel_shenanigans(line, iseq, outstrs, colors, delstrs, debug=False):  # NOTE similar to/overlaps with get_seq_with_indels_reinstated()
    if debug:
        print('%s' % utils.color('blue', 'shenanigans'))

    def get_itype(qrchar, glchar):
        if qrchar not in utils.gap_chars and glchar not in utils.gap_chars:
            return None
        elif glchar not in utils.gap_chars:
            return 'deletion'
        elif qrchar not in utils.gap_chars:
            return 'insertion'
        else:
            return 'uh'  # i guess it's possible if there's overlapping indels, but I don't want to think about that a.t.m.
    def choose_char(ostrchar, istr, qrchar, glchar):  # <ostrchar> is the character that's there already
        if ostrchar not in utils.nukes + utils.all_ambiguous_bases:  # if this bit of the sequences is spaces, dots, or dashes, then we only want to insert spaces (note that this adds some arbitrariness on boundaries as to who gets the actual inserted string)
            return ' '
        elif use_stars(get_itype(qrchar, glchar), istr):
            return '*'
        elif qrchar not in utils.gap_chars and glchar not in utils.gap_chars:
            return ostrchar
        elif glchar not in utils.gap_chars:
            return glchar
        elif qrchar not in utils.gap_chars:
            return qrchar
        else:
            assert False

    ifo = line['indelfos'][iseq]

    # eh, fuck it, there's too many ways this can go wrong
    # print ' in     \'%s\'' % line['input_seqs'][iseq]
    # print ' qr     \'%s\'' % outstrs[-1]
    # print ' qr gap \'%s\'' % ifo['qr_gap_seq']
    # print ' gl gap \'%s\'' % ifo['gl_gap_seq']
    # assert len(ifo['qr_gap_seq']) == len(ifo['gl_gap_seq'])
    # # non_gap_len() is if there were dashes added because of no space, rstrip() is in case we fixed length problems in check_outsr_lengths()
    # if utils.non_gap_len(outstrs[-1].rstrip()) - line['v_5p_del'] - line['j_3p_del'] != len(ifo['qr_gap_seq']) - utils.gap_len(ifo['gl_gap_seq']):
    #     print '%d - %d - %d = %d != %d - %d = %d' % (utils.non_gap_len(outstrs[-1].rstrip()), line['v_5p_del'], line['j_3p_del'], utils.non_gap_len(outstrs[-1].rstrip()) - line['v_5p_del'] - line['j_3p_del'],
    #                                                  len(ifo['qr_gap_seq']), utils.gap_len(ifo['gl_gap_seq']), len(ifo['qr_gap_seq']) - utils.gap_len(ifo['gl_gap_seq']))
    #     assert False

    ipos, iqrgap, iglgap = 0, 0, 0  # line['v_5p_del']
    istop = len(outstrs[-1])
    new_outstrs, new_colors = [[] for _ in outstrs], [[] for _ in colors]
    while ipos < istop:
        outchars = [ostr[ipos] for ostr in outstrs]
        if ipos < len(delstrs['v_5p']) or ipos >= istop - line['j_3p_del'] or '-' in outchars or ' ' in outchars[-1]:  # '-' is to check for no-extra-space fix, and ' ' is in case we fixed length problems with check_outsr_lengths()
            for istr in range(len(outstrs)):
                new_outstrs[istr] += [outstrs[istr][ipos]]
                new_colors[istr] += [[]]
            ipos += 1
            continue

        if iqrgap >= len(ifo['qr_gap_seq']) or iglgap >= len(ifo['gl_gap_seq']):
            print(' ins     %3d  x%sx' % (ipos, outstrs[0]))
            print('  d           x%sx' % (outstrs[1]))
            print(' vj           x%sx' % (outstrs[2]))
            print(' qr           x%sx' % (outstrs[3]))
            print(' qr gap  %s  x%sx' % (utils.color('red' if iqrgap >= len(ifo['qr_gap_seq']) else None, '%3d' % iqrgap), ifo['qr_gap_seq']))
            print(' gl gap  %s  x%sx' % (utils.color('red' if iglgap >= len(ifo['gl_gap_seq']) else None, '%3d' % iglgap), ifo['gl_gap_seq']))
            print('%s index problem when adding indel info to print strings in prutils.indel_shenanigans() (probably due to overlapping indels)' % utils.color('red', 'error'))
            break
        # print ipos, iqrgap, iglgap, ifo['qr_gap_seq'][iqrgap], ifo['gl_gap_seq'][iglgap]
        if ifo['qr_gap_seq'][iqrgap] not in utils.gap_chars and ifo['gl_gap_seq'][iglgap] not in utils.gap_chars:
            for istr in range(len(outstrs)):
                new_outstrs[istr] += [outstrs[istr][ipos]]
                new_colors[istr] += [[]]
            ipos += 1
            iqrgap += 1
            iglgap += 1
        else:
            for istr in range(len(outstrs)):
                new_char = choose_char(outstrs[istr][ipos], istr, ifo['qr_gap_seq'][iqrgap], ifo['gl_gap_seq'][iglgap])
                # it would be nice to have a way to double check that we're at the right position here, but I can't figure one out (i think i'd have to totally rewrite this stuff so the indels got incorporated before constructing outstrs or something)
                # anyway, the below doesn't work if iseq > 0 or if there's more than one indel
                # if iseq == 0 and is_qr(istr) and new_char == '*':  # if this is a deletion, we can check that the bases in the gl gap seq are the same as in the germline (outstr) seq (this used to happen because i was using the v_5p len rather than the length in delstrs['v_5p'] above) (we can't check it if iseq > 0 since we don't expect the vj line to match up in that case)
                #     if vj_ostr(outstrs)[ipos] != vj_ostr(new_outstrs)[ipos]:  # i guess just call it a 'warning', since i think it might get triggered in d or j when things are ok
                #         print '    %s base at %d from gl gap seq %s different to base in gl outstr seq %s (indel is probably at wrong position) when printing ascii dbg' % (utils.color('red', 'warning'), ipos, vj_ostr(new_outstrs)[ipos], vj_ostr(outstrs)[ipos])
                new_outstrs[istr] += [new_char]
                new_colors[istr] += [[]]
                if iseq > 0 and is_qr(istr):  # the germline lines are printed based on the *first* query sequence in a multi-seq alignment, and it'd be kinda hard to really account for all subsequent sequences' indels, so we compromise by blueing the insertions in the subsequent query sequences (deletions are already blue stars). Note that they still don't line up right.
                    new_colors[istr][-1] += ['blue', 'reverse_video']
            if ifo['gl_gap_seq'][iglgap] not in utils.gap_chars:  # if this is a dot in the gl gap seq, then this base was an insertion in the query sequence, i.e. we removed it for the reversed seq, i.e. we don't want to increment <ipos>
                ipos += 1
            iqrgap += 1
            iglgap += 1

    outstrs = [''.join(ostr) for ostr in new_outstrs]
    colors = new_colors
    check_outsr_lengths(line, outstrs)

    if debug:
        for ostr in outstrs:
            print('    %s' % ostr)

    return outstrs, colors

# ----------------------------------------------------------------------------------------
def add_colors(outstrs, colors, line):  # NOTE do *not* modify <line>
    # <outstrs> convention: [indels, d, vj, query]
    bluechars = utils.all_ambiguous_bases + ['*', '-']

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
        alphagl = ''.join(filter((utils.alphabet).__contains__, glchars))  # glchars that are also alphabet chars
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
        if len(list(filter((bluechars).__contains__, outstrs[istr]))) == 0:  # skip anybody that doesn't have any blue chars
            continue
        blue_indices = [inuke for inuke in range(len(outstrs[istr])) if outstrs[istr][inuke] in bluechars]
        for inuke in blue_indices:
            colors[istr][inuke].append('blue')

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
def print_seq_in_reco_event(original_line, iseq, extra_str='', label='', one_line=False, queries_to_emphasize=None, duplicated_uids=None, check_line_integrity=False, uid_extra_str='', uid_extra_str_label=None):
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
    # UPDATE turning this off, now since I'm so frequently looking at multiple annotations on top of each other, it's more important to have them really line up vertically
    d_plus_inserts_length = len(line['vd_insertion'] + line['d_gl_seq'] + line['dj_insertion'])
    vj_delstr = ''
    if line['v_3p_del'] + line['j_5p_del'] > d_plus_inserts_length:  # if dots for v and j interior deletions will be longer than <d_plus_inserts_length>
        vj_del_strs = [('%s %d' % (d, line[d + '_del'])) for d in ['v_3p', 'j_5p'] if line[d + '_del'] > 0]
        vj_delstr = '  %s: %s' % (utils.color('blue', 'del'), ', '.join(vj_del_strs))
        delstrs['v_3p'] = '.' * d_plus_inserts_length
        delstrs['j_5p'] = ''
        gapstr = ''
        gap_insert_point = None
        extra_space_because_of_fixed_nospace = 0
        # old way, adding blue dashes:
        # delstrs['v_3p'] = '.%d.' % line['v_3p_del']
        # delstrs['j_5p'] = '.%d.' % line['j_5p_del']
        # gapstr = '-' * (len(delstrs['v_3p'] + delstrs['j_5p']) - d_plus_inserts_length)
        # gap_insert_point = len(line['fv_insertion'] + delstrs['v_5p'] + line['v_gl_seq'])  # it doesn't really matter exactly where we put the blue dashes, as long as it's the same place in all four lines, but this is a good spot
        # extra_space_because_of_fixed_nospace = max(0, d_plus_inserts_length - len(delstrs['v_3p'] + delstrs['j_5p']))  # if shortening the <delstrs> already over-compensated for the lack of space (i.e., if the number of dashes necessary is zero), then we need to add some dots to the vj line below
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
    check_outsr_lengths(line, outstrs, fix=True)  # I think the only way they can be different is if the d right side erosion is so long that it hangs over the right side of the j

    if gap_insert_point is not None:
        for istr in [0, 1, 3]:  # everybody except the vj line, which already has the modified interior delstrs above
            outstrs[istr] = outstrs[istr][:gap_insert_point] + gapstr + outstrs[istr][gap_insert_point:]

    check_outsr_lengths(line, outstrs, fix=True)

    colors = [[[] for _ in range(len(ostr))] for ostr in outstrs]
    if indelutils.has_indels(line['indelfos'][iseq]):
        # outstrs, colors = old_indel_shenanigans(line, iseq, outstrs, colors)
        outstrs, colors = indel_shenanigans(line, iseq, outstrs, colors, delstrs)
    outstrs = add_colors(outstrs, colors, line)

    mtpy = utils.get_multiplicity(line, None, iseq)
    mtpystr = ' ' if (mtpy == 1 or mtpy is None) else utils.color('blue', str(mtpy))
    uidstr = get_uid_str(line, iseq, queries_to_emphasize, duplicated_uids=duplicated_uids)
    vjlabelstr = '%s %s' % (utils.color_gene(line['v_gene']), utils.color_gene(line['j_gene']))
    if uid_extra_str != '':  # pad the uid str out to same length as vj str so the uid extra strs (e.g. extra_print_keys) are lined up with their labels
        uidstr += ' ' * (utils.len_excluding_colors(vjlabelstr) - utils.len_excluding_colors(uidstr))
    tmppad = utils.len_excluding_colors(uidstr) - utils.len_excluding_colors(vjlabelstr)
    suffixes = ['insert%s\n'       % ('s' if utils.has_d_gene(utils.get_locus(line['v_gene'])) else ''),
                '%s%s\n'             % (utils.color_gene(line['d_gene']), vj_delstr),
                '%s%s  %%shm  %s%s\n'        % (vjlabelstr, tmppad * ' ', utils.color('blue', 'N') if any((m is not None and m > 1) for m in utils.get_multiplicities(line)) else '', uid_extra_str_label if uid_extra_str_label is not None else ''),
                '%s  %4.1f  %s %s %s\n' % (uidstr, 100*line['mut_freqs'][iseq], mtpystr, uid_extra_str, utils.color('red', utils.is_functional_dbg_str(line, iseq)))]
    outstrs = ['%s%s   %s' % (extra_str, ostr, suf) for ostr, suf in zip(outstrs, suffixes)]

    if label != '':  # this doesn't really work if the edge of the removed string is the middle of a color code... but oh well, it doesn't really happen any more since I shortened the kbound label from waterer.py
        offset = max(0, len(extra_str) - 2)  # skootch <label> this many positions leftward into <extra_str>
        removed_str = outstrs[0][offset : offset + utils.len_excluding_colors(label)]
        outstrs[0] = outstrs[0][ : offset] + label + outstrs[0][utils.len_excluding_colors(label) + offset : ]  # NOTE this *replaces* the bases in <extra_str> with <label>, which is only fine if they're spaces
        if removed_str.strip() != '':
            print('%s%s (covered by label \'%s\')' % (' ' * offset, utils.color('red', removed_str), label))

    if one_line:
        outstrs = outstrs[-1:]  # remove all except the query seq line
    elif not utils.has_d_gene(utils.get_locus(line['v_gene'])) and len(vj_delstr) == 0:
        outstrs.pop(1)  # remove the d germline line

    print(''.join(outstrs), end='')

    if check_line_integrity:
        if set(line.keys()) != set(original_line.keys()):
            raise Exception('ack 1')
        for k in line:
            if line[k] != original_line[k]:
                print('key %s differs:\n  %s\n  %s ' % (k, line[k], original_line[k]))
                raise Exception('')
