from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import numpy as np
# from GCutils import hamming_distance
from random import randint
import sys
import collections

from . import utils
from . import treeutils

'''
copied from bcr-phylo/bin/COAR.py (yes it sucks to duplicate it, but there's too much pickle/ete dependent stuff in there, and i'm not even sure how to create the ete/pickle 'forest' files from the inference methods))
'''

# ----------------------------------------------------------------------------------------
def reconstruct_lineage(tree, node):
    lineage, uids = [], []  # uids is just for dbg atm
    while True:
        lineage.append(node.seq)
        uids.append(node.taxon.label)
        if node.parent_node is None:
            return lineage, uids
        node = node.parent_node

# ----------------------------------------------------------------------------------------
def find_node(tree, seq, uid):
    nodes = [n for n in tree.leaf_node_iter() if n.seq == seq]
    if len(nodes) > 1:
        nodes = [n for n in nodes if n.taxon.label == uid]
    dbgpair = ('seq', seq) if seq is not None else ('uid', uid)
    if len(nodes) == 0:
        raise Exception('couldn\'t find node with %s %s' % dbgpair)
    elif len(nodes) > 1:
        raise Exception('found multiple nodes with %s %s' % dbgpair)
    return nodes[0]

# ----------------------------------------------------------------------------------------
def align_lineages(node_t, tree_t, tree_i, gap_penalty_pct=0, known_root=True, allow_double_gap=False, test=False, debug=False):
    '''
    Standard implementation of a Needleman-Wunsch algorithm as described here:
    http://telliott99.blogspot.com/2009/08/alignment-needleman-wunsch.html
    https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
    And implemented here:
    https://github.com/alevchuk/pairwise-alignment-in-python/blob/master/alignment.py

    gap_penalty_pct is the gap penalty relative to the sequence length of the sequences on the tree.

    Set test to True to evalue the simple example from supplemental info in the paper.
    '''
    # ----------------------------------------------------------------------------------------
    def get_gap_penalties(len_t, len_i):
        # Gap penalty chosen not too large:
        gap_penalty = -1 * int((len(node_t.seq) / 100.0) * gap_penalty_pct)
        assert(gap_penalty <= 0)  # Penalties must be negative
        if gap_penalty == 0:  # If gap penalty is zero only gaps in the shortest sequence will be allowed
            assert(allow_double_gap is False)

        # Generate a score matrix:
        # Disallow gaps in the longest list:
        if allow_double_gap is False and len_t > len_i:
            # If true is longer than inferred allow gap only in inferred:
            gp_i = gap_penalty
            gp_j = -1 * float('inf')
        elif allow_double_gap is False and len_t < len_i:
            # If inferred is longer than true allow gap only in true:
            gp_i = -1 * float('inf')
            gp_j = gap_penalty
        elif allow_double_gap is False and len_t == len_i:
            # If lists are equally long no gaps are allowed:
            gp_i = -1 * float('inf')
            gp_j = -1 * float('inf')
        else:
            gp_i = gap_penalty
            gp_j = gap_penalty
        if debug:
            print('      gap penalties:  %.0f   i %.0f  j %.0f' % (gap_penalty, gp_i, gp_j))
        return gap_penalty, gp_i, gp_j
    # ----------------------------------------------------------------------------------------
    def check_test_values(alignment_score, max_penalty, align_t, align_i):
        correct_aln_score, correct_max_pen = -1, -3
        correct_align_t = ['TTT', 'ATT', 'AAT', 'AAA']
        correct_align_i = ['TTT', '-', 'TAT', 'AAA']
        if alignment_score != correct_aln_score:
            raise Exception('incorrect alignment score %.2f (should be %.2f)' % (alignment_score, correct_aln_score))
        if max_penalty != correct_max_pen:
            raise Exception('incorrect alignment score %.2f (should be %.2f)' % (max_penalty, correct_max_pen))
        if align_t != correct_align_t:
            raise Exception('incorrect true alignment %s (should be %s' % (align_t, correct_align_t))
        if align_i != correct_align_i:
            raise Exception('incorrect inferred alignment %s (should be %s' % (align_i, correct_align_i))
        print('  all test values ok')

    # ----------------------------------------------------------------------------------------
    gap_seq = '-'
    if test:  # example from supp info in paper
        lin_t = ['AAA', 'AAT', 'ATT', 'TTT']
        lin_i = ['AAA', 'TAT', 'TTT']
        uids_t = ['naive', 'a1', 'a2', 'leaf']
        uids_i = ['naive', 'a1', 'leaf']
    else:
        node_i = find_node(tree_i, node_t.seq, node_t.taxon.label)  # looks first by seq, then disambuguates with uid (if you only look by uid, if the inference swaps two nearby internal/leaf nodes with the same seq it'll crash)
        (lin_t, uids_t), (lin_i, uids_i) = [reconstruct_lineage(t, n) for t, n in [(tree_t, node_t), (tree_i, node_i)]]
    # One lineages must be longer than just the root and the terminal node
    if len(lin_t) <= 2 and len(lin_i) <= 2:
        return False
    if debug:
        print('      aligning lineage for true node %s: found inf node %s' % (uids_t[0], uids_i[0]))
        max_len = max(len(u) for u in uids_i + uids_t)
        print('              %s' % '  '.join(utils.wfmt(i, max_len, fmt='d', jfmt='-') for i in range(max(len(uids_t), len(uids_i)))))
        print('        true: %s' % '  '.join(utils.wfmt(u, max_len, jfmt='-') for u in uids_t))
        print('         inf: %s' % '  '.join(utils.wfmt(u, max_len, jfmt='-') for u in uids_i))

    len_t, len_i = [len(l) for l in [lin_t, lin_i]]
    gap_penalty, gp_i, gp_j = get_gap_penalties(len_t, len_i)

    sc_mat = np.zeros((len_t, len_i), dtype=np.float64)
    for i in range(len_t):
        for j in range(len_i):
            # Notice the score is defined by number of mismatches:
            #sc_mat[i, j] = len(lin_t[i]) - hamming_distance(lin_t[i], lin_i[j])
            sc_mat[i, j] = -1 * utils.hamming_distance(lin_t[i], lin_i[j])

    # Calculate the alignment scores:
    aln_sc = np.zeros((len_t+1, len_i+1), dtype=np.float64)
    for i in range(0, len_t+1):
        if known_root is True:
            aln_sc[i][0] = -1 * float('inf')
        else:
            aln_sc[i][0] = gp_i * i
    for j in range(0, len_i+1):
        if known_root is True:
            aln_sc[0][j] = -1 * float('inf')
        else:
            aln_sc[0][j] = gp_j * j
    aln_sc[0][0] = 0  # The top left is fixed to zero
    for i in range(1, len_t+1):
        for j in range(1, len_i+1):
            match = aln_sc[i-1][j-1] + sc_mat[i-1, j-1]
            gap_in_inferred = aln_sc[i-1][j] + gp_i
            gap_in_true = aln_sc[i][j-1] + gp_j
            aln_sc[i][j] = max(match, gap_in_inferred, gap_in_true)
    # Traceback to compute the alignment:
    align_t, align_i, asr_align, aln_ids = list(), list(), list(), {'true' : [], 'inf' : []}
    i, j = len_t, len_i
    alignment_score = aln_sc[i][j]
    while i > 0 and j > 0:
        sc_current = aln_sc[i][j]
        sc_diagonal = aln_sc[i-1][j-1]
        sc_up = aln_sc[i][j-1]
        sc_left = aln_sc[i-1][j]

        if sc_current == (sc_diagonal + sc_mat[i-1, j-1]):
            align_t.append(lin_t[i-1])
            align_i.append(lin_i[j-1])
            aln_ids['true'].append(uids_t[i-1])
            aln_ids['inf'].append(uids_i[j-1])
            i -= 1
            j -= 1
        elif sc_current == (sc_left + gp_i):
            align_t.append(lin_t[i-1])
            align_i.append(gap_seq)
            aln_ids['true'].append(uids_t[i-1])
            aln_ids['inf'].append('gap')
            i -= 1
        elif sc_current == (sc_up + gp_j):
            align_t.append(gap_seq)
            align_i.append(lin_i[j-1])
            aln_ids['true'].append('gap')
            aln_ids['inf'].append(uids_i[j-1])
            j -= 1

    # If space left fill it with gaps:
    while i > 0:
        asr_align.append(gp_i)
        align_t.append(lin_t[i-1])
        align_i.append(gap_seq)
        i -= 1
    while j > 0:
        asr_align.append(gp_j)
        align_t.append(gap_seq)
        align_i.append(lin_i[j-1])
        j -= 1

    max_penalty = 0
    for a, b in zip(align_t, align_i):
        if a == gap_seq or b == gap_seq:
            max_penalty += gap_penalty
        else:
            max_penalty += -len(a)
    # Notice that the root and the terminal node is excluded from this comparison.
    # by adding their length to the max_penalty:
    if known_root is True:
        max_penalty += 2 * len(lin_t[0])
    else:  # Or in the case of an unknown root, just add the terminal node
        max_penalty += len(lin_t[0])

    if debug:
        max_seq_len = max(len(s) for slist in [align_t, align_i] for s in slist)
        max_uid_len = max(len(u) for u in uids_i + uids_t)
        print('      aligned lineages:')
        print(('          index  hdist  %'+str(max_seq_len)+'s  %'+str(max_seq_len)+'s')  % ('true', 'inferred'))
        for uid_t, uid_i, seq_t, seq_i in zip(aln_ids['true'], aln_ids['inf'], align_t, align_i):
            def cfn(s): return utils.color('blue', s, width=max_seq_len, padside='right') if s == gap_seq else s
            str_t, str_i = cfn(seq_t), cfn(seq_i)
            hdstr = utils.color('blue', '-', width=3)
            if all(s != gap_seq for s in (seq_t, seq_i)):
                str_i, isnps = utils.color_mutants(seq_t, seq_i, return_isnps=True)
                hdstr = utils.color('red' if len(isnps) > 0 else None, '%3d' % len(isnps))
            def ustr(u): return utils.color('blue' if u=='gap' else None, utils.wfmt(u, max_uid_len))
            print('           %s  %s %s  %s   %s' % (hdstr, ustr(uid_t), ustr(uid_i), str_t, str_i))
        if alignment_score % 1 != 0 or max_penalty % 1 != 0:
            raise Exception('alignment_score score %s or max_penalty %s not integers, so need to fix dbg print in next line' % (alignment_score, max_penalty))
        print('      alignment score: %.0f   max penalty: %.0f' % (alignment_score, max_penalty))
        if test:
            check_test_values(alignment_score, max_penalty, align_t, align_i)

    return [align_t, align_i, alignment_score, max_penalty]

# ----------------------------------------------------------------------------------------
def COAR(true_tree, inferred_tree, known_root=True, allow_double_gap=False, debug=False):
    lineage_dists, n_skipped = list(), collections.OrderedDict([('inferred-internal', []), ('missing-leaf', []), ('aln-fail', [])])
    inf_leaf_nodes = [n.taxon.label for n in inferred_tree.leaf_node_iter()]  # just so we can skip true leaf nodes that were inferred to be internal
    for node_t in true_tree.leaf_node_iter():
        nlabel = node_t.taxon.label
        if debug:
            print('%s             %3d %s' % (nlabel, len(node_t.seq), node_t.seq))
        if nlabel not in inf_leaf_nodes:
            is_internal = any(n.taxon.label == nlabel for n in inferred_tree.preorder_node_iter())
            n_skipped['inferred-internal' if is_internal else 'missing-leaf'].append(nlabel)
            continue
        aln_res = align_lineages(node_t, true_tree, inferred_tree, known_root=known_root, allow_double_gap=allow_double_gap, debug=debug)
        if aln_res is False:  # Skip lineages less than three members long
            n_skipped['aln-fail'].append(nlabel)
            continue
        align_t, align_i, final_score, max_penalty = aln_res
        if max_penalty < 0:
            lineage_dists.append(final_score / float(max_penalty))
            if debug:
                print('    normalized dist: %.3f' % lineage_dists[-1])
        else:
            if debug:
                print('    max penalty not less than zero: %.3f' % max_penalty)

    if any(len(v) > 0 for v in n_skipped.values()):
        dstrs = {'missing-leaf' : '%d missing from inferred tree: %s', 'inferred-internal' : '%d were internal in inferred tree: %s', 'aln-fail' : '%d failed lineage alignment: %s'}
        print('    %s skipped %d / %d true leaf nodes in coar calculation (%s)' % (utils.wrnstr(), sum(len(v) for v in n_skipped.values()), len(list(true_tree.leaf_node_iter())), ', '.join((dstrs[k]%(len(n_skipped[k]), ' '.join(sorted(n_skipped[k])))) for k in n_skipped if len(n_skipped[k])>0)))

    if len(lineage_dists) == 0:  # max_penalty is 0 when all lineages have less than three members
        if debug:
            print('  all lineages shorter than 3, returning 0')
        return 0
    if debug:
        print('  mean over %d lineages: %.5f' % (len(lineage_dists), sum(lineage_dists) / float(len(lineage_dists))))
    return sum(lineage_dists) / float(len(lineage_dists))
