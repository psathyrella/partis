#!/usr/bin/env python
# copied from https://github.com/briney/VRC01gH-GT3/blob/master/utils/vrc01.py
#  - then edited to work within partis, since original (e.g. abtools import) only worked on python 2.7
# filename: vrc01.py


#
# Copyright (c) 2016 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


import sys
import numpy as np
import time

from . import utils

# ----------------------------------------------------------------------------------------
# returns two lists (with entries for each input sequence): <shared_vals>: mutations shared with any vrc01-class sequence and <total_vals>: total number of mutations in input sequence
def vrc01_class_mutation_count(input_seqs, debug=True):
    start = time.time()
    shared_vals, total_vals = [], []
    vrc01_seqs = get_vrc01_class_sequences()
    vrc01_names = [s['name'] for s in vrc01_seqs]
    glvrc01 = get_vrc01_germline_sequence()
    glvrc01_name = glvrc01['name']
    msa_seqfos = utils.align_many_seqs([glvrc01] + vrc01_seqs + input_seqs, aa=True, extra_str='      ', debug=debug)
    for info in input_seqs:
        aln_in = utils.get_single_entry([s for s in msa_seqfos if s['name'] == info['name']])  # alignment of query sequence (info)
        aln_gl = utils.get_single_entry([s for s in msa_seqfos if s['name'] == glvrc01_name])  # alignment of germline sequence
        aln_vrc01s = [s for s in msa_seqfos if s['name'] in vrc01_names]  # alignment of vrc01-class (reference) seqs
        skip_chars = set(utils.gap_chars + utils.ambiguous_amino_acids)
        total_vals.append(sum([i != g for i, g in zip(aln_in['seq'], aln_gl['seq']) if i not in skip_chars and g not in skip_chars]))  # total number of differences between aln_in and aln_gl (if query base (<i>) is a gap, i guess that counts as a mutation?
        all_shared = {}  # all mutations shared with each vrc01-class seq
        if debug:
            print('  getting shared for %s:' % info['name'])
        for vrc01 in aln_vrc01s:
            _shared = []
            chk_tot = 0
            twidth = max(len(s['name']) for s in msa_seqfos)
            for q, g, v in zip(aln_in['seq'], aln_gl['seq'], vrc01['seq']):  # query, germline, and vrc01-class (ref) characters
                if any(c in skip_chars for c in [q, g, v]):  # if any of the three bases are either gaps or ambiguous, we can't know it's shared (or mutated)
                    _shared.append(False)
                    if all(c not in skip_chars for c in [q, g]) and q != g:  # if query and germline are non-skipped and different, it's a (non-shared) mutation [just for checking total above)
                        chk_tot += 1
                elif q == g:  # unmutated: query and germline are the same
                    _shared.append(False)
                elif q != g:
                    if q == v:  # same mutation as vrc01 (ref)
                        _shared.append(True)
                    else:
                        _shared.append(False)
                    chk_tot += 1
                else:
                    assert False
            if chk_tot != total_vals[-1]:  # this shouldn't happen, but i'm scared i forgot some case
                print('  %s different sum totals %d %d' % (utils.wrnstr(), total_vals[-1], chk_tot))
            if debug:
                print('      %s %s' % (utils.wfmt(vrc01['name'], twidth), ''.join(utils.color(None, 'x') if s else '-' for s in _shared)))
            all_shared[vrc01['name']] = _shared
        any_shared = 0
        for pos in zip(*all_shared.values()):  # loop over all positions in the aligned sequences
            if any(pos):  # if this input seq shares a mutation with any vrc01-class seq at this position...
                any_shared += 1  # ...then increment this
        shared_vals.append(any_shared)
        if debug:
            print('      shared: %d  total: %d' % (any_shared, total_vals[-1]))
    print('  vrc01 class mutation counts for %d seqs in %.1f sec' % (len(input_seqs), time.time() - start))
    return shared_vals, total_vals


# ----------------------------------------------------------------------------------------
def vrc01_class_mutation_positions(seqs):
# TODO
    assert False
    data = []
    input_seqs = [Sequence([s['seq_id'], s['vdj_aa']]) for s in seqs]
    input_names = [s.id for s in input_seqs]
    # get VRC01-class sequences
    hiv_seqs = get_vrc01_class_sequences()
    all_hiv_names = [s.id for s in hiv_seqs]
    # MSA
    seqs_for_alignment = input_seqs + hiv_seqs
    seqs_for_alignment.append(get_vrc01_germline_sequence(vgene_only=False))
    aln = muscle(seqs_for_alignment)
    aln_seqs = [seq for seq in aln if seq.id in input_names]
    aln_gl = [seq for seq in aln if seq.id == 'glVRC01'][0]
    aln_mins = [seq for seq in aln if seq.id in ['minVRC01', 'min12A21']]
    aln_hiv = [seq for seq in aln if seq.id in all_hiv_names]
    for seq in aln_seqs:
        seq_data = []
        for i, (s, g) in enumerate(zip(str(seq.seq), str(aln_gl.seq))):
            # if g == '-' and s == '-':
            if g == '-':
                continue
            min_residues = [seq[i] for seq in aln_mins]
            vrc01_residues = [seq[i] for seq in aln_hiv]
            if s == '-':
                seq_data.append(0)
            elif s == g:
                seq_data.append(0)
            elif s != g and s in min_residues:
                seq_data.append(2)
            elif s != g and s in vrc01_residues:
                seq_data.append(3)
            elif s != g and s not in vrc01_residues:
                seq_data.append(1)
            else:
                seq_data.append(0)
        data.append(np.asarray(seq_data))
    return np.asarray(data)


def get_vrc01_germline_sequence(vgene_only=True):
    if vgene_only:
        gl_vrc01 = ('glVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAR')
    else:
        gl_vrc01 = ('glVRC01', 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGKNSDYNWDFQHWGQGTLVTVSS')
    return {'name' : gl_vrc01[0], 'seq' : gl_vrc01[1]}


def get_vrc01_class_sequences(chain='heavy', vgene_only=True, only_include=None):
    if vgene_only:
        heavy = [('VRC01', 'QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTR'),
                 ('PGV04', 'QVQLVQSGSGVKKPGASVRVSCWTSEDIFERTELIHWVRQAPGQGLEWIGWVKTVTGAVNFGSPDFRQRVSLTRDRDLFTAHMDIRGLTQGDTATYFCAR'),
                 ('VRC-CH31', 'QVQLVQSGAAVRKPGASVTVSCKFAEDDDYSPYWVNPAPEHFIHFLRQAPGQQLEWLAWMNPTNGAVNYAWYLNGRVTATRDRSMTTAFLEVKSLRSDDTAVYYCAR'),
                 ('3BNC60', 'QVHLSQSGAAVTKPGASVRVSCEASGYKISDHFIHWWRQAPGQGLQWVGWINPKTGQPNNPRQFQGRVSLTRQASWDFDTYSFYMDLKAVRSDDTAIYFCAR'),
                 ('12A12', 'HLVQSGTQVKKPGASVRISCQASGYSFTDYVLHWWRQAPGQGLEWMGWIKPVYGARNYARRFQGRINFDRDIYREIAFMDLSGLRSDDTALYFCAR'),
                 ('PGV20', 'QVHLMQSGTEMKKPGASVRVTCQTSGYTFSDYFIHWLRQVPGRGFEWMGWMNPQWGQVNYARTFQGRVTMTRDVYREVAYLDLRSLTFADTAVYFCAR')]
        light = []
    else:
        heavy = [('VRC01', 'QVQLVQSGGQMKKPGESMRISCRASGYEFIDCTLNWIRLAPGKRPEWMGWLKPRGGAVNYARPLQGRVTMTRDVYSDTAFLELRSLTVDDTAVYFCTRGKNCDYNWDFEHWGRGTPVIVSS'),
                 ('PGV04', 'QVQLVQSGSGVKKPGASVRVSCWTSEDIFERTELIHWVRQAPGQGLEWIGWVKTVTGAVNFGSPDFRQRVSLTRDRDLFTAHMDIRGLTQGDTATYFCARQKFYTGGQGWYFDLWGRGTLIVVSS'),
                 ('VRC-CH31', 'QVQLVQSGAAVRKPGASVTVSCKFAEDDDYSPYWVNPAPEHFIHFLRQAPGQQLEWLAWMNPTNGAVNYAWYLNGRVTATRDRSMTTAFLEVKSLRSDDTAVYYCARAQKRGRSEWAYAHWGQGTPVVVSS'),
                 ('3BNC60', 'QVHLSQSGAAVTKPGASVRVSCEASGYKISDHFIHWWRQAPGQGLQWVGWINPKTGQPNNPRQFQGRVSLTRQASWDFDTYSFYMDLKAVRSDDTAIYFCARQRSDFWDFDVWGSGTQVTVSS'),
                 ('12A12', 'HLVQSGTQVKKPGASVRISCQASGYSFTDYVLHWWRQAPGQGLEWMGWIKPVYGARNYARRFQGRINFDRDIYREIAFMDLSGLRSDDTALYFCARDGSGDDTSWHLDPWGQGTLVIVSA'),
                 ('PGV20', 'QVHLMQSGTEMKKPGASVRVTCQTSGYTFSDYFIHWLRQVPGRGFEWMGWMNPQWGQVNYARTFQGRVTMTRDVYREVAYLDLRSLTFADTAVYFCARRMRSQDREWDFQHWGQGTRIIVSS')]
        light = []
    seqs = heavy if chain == 'heavy' else light
    if only_include is not None:
        if type(only_include) in [str, unicode]:
            only_include = [only_include, ]
        seqs = [s for s in seqs if s[0] in only_include]
    return [{'name' : n, 'seq' : s} for n, s in seqs]


def get_vrc01_class_mutations():
    vrc01_class = get_vrc01_class_sequences()
    glvrc01 = get_vrc01_germline_sequence()
    return list(set(_get_mutations(vrc01_class, glvrc01)))


def _get_mutations(seqs, standard):
    mutations = []
    msa_seqfos = utils.align_many_seqs([standard] + seqs, aa=True, debug=True)
    glfo, queryfos = msa_seqfos[0], msa_seqfos[1:]
    for qfo in queryfos:
        mutations.extend(_parse_mutations(qfo, glfo))
    return mutations


def _parse_mutations(qfo, glfo):
    muts = []
    tpos = 0
    for q, t in zip(qfo['seq'], glfo['seq']):
        if t == '-':
            continue
        tpos += 1
        if any([q == '-', q == t]):
            continue
        mut = '{}{}{}'.format(t, tpos, q)
        muts.append(mut)
    return muts
