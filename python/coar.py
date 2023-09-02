import numpy as np
# from GCutils import hamming_distance
from random import randint
import sys

import utils
import treeutils

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
def find_node(tree, seq=None, uid=None):
# ----------------------------------------------------------------------------------------
# TODO TODO don't even need this anymore?
    assert [seq, uid].count(None) == 1  # specify either seq or uid
    if seq is not None:
        nodes = [n for n in tree.leaf_node_iter() if n.seq == seq] # and node.frequency > 0]
    elif uid is not None:
        nodes = [n for n in tree.leaf_node_iter() if n.taxon.label == seq] # and node.frequency > 0]
    else:
        assert False
    dbgpair = ('seq', seq) if seq is not None else ('uid', uid)
    if len(nodes) == 0:
        for l in tree.leaf_node_iter():
            print l.taxon.label, seq == l.seq, seq in l.seq, utils.color_mutants(seq, l.seq, align_if_necessary=True)
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
            print '      gap penalties:  %.0f   i %.0f  j %.0f' % (gap_penalty, gp_i, gp_j)
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
        print '  all test values ok'

    # ----------------------------------------------------------------------------------------
    gap_seq = '-'
    if test:  # example from supp info in paper
        lt = ['AAA', 'AAT', 'ATT', 'TTT']
        li = ['AAA', 'TAT', 'TTT']
        uids_t = ['naive', 'a1', 'a2', 'leaf']
        uids_i = ['naive', 'a1', 'leaf']
    else:
        nt = node_t
# TODO TODO why am i still using this?
        ni = find_node(tree_i, seq=node_t.seq) #, uid=node_t.taxon.label)
        (lt, uids_t), (li, uids_i) = [reconstruct_lineage(t, n) for t, n in [(tree_t, nt), (tree_i, ni)]]
    # One lineages must be longer than just the root and the terminal node
    if len(lt) <= 2 and len(li) <= 2:
        return False
    if debug:
        print '      aligning lineage for true node %s: found inf node %s' % (uids_t[0], uids_i[0])
        print '        true: %s' % ' '.join(uids_t)
        print '         inf: %s' % ' '.join(uids_i)

    len_t, len_i = [len(l) for l in [lt, li]]
    gap_penalty, gp_i, gp_j = get_gap_penalties(len_t, len_i)

    sc_mat = np.zeros((len_t, len_i), dtype=np.float64)
    for i in range(len_t):
        for j in range(len_i):
            # Notice the score is defined by number of mismatches:
            #sc_mat[i, j] = len(lt[i]) - hamming_distance(lt[i], li[j])
            sc_mat[i, j] = -1 * utils.hamming_distance(lt[i], li[j])

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
    align_t, align_i, asr_align = list(), list(), list()
    i, j = len_t, len_i
    alignment_score = aln_sc[i][j]
    while i > 0 and j > 0:
        sc_current = aln_sc[i][j]
        sc_diagonal = aln_sc[i-1][j-1]
        sc_up = aln_sc[i][j-1]
        sc_left = aln_sc[i-1][j]

        if sc_current == (sc_diagonal + sc_mat[i-1, j-1]):
            align_t.append(lt[i-1])
            align_i.append(li[j-1])
            i -= 1
            j -= 1
        elif sc_current == (sc_left + gp_i):
            align_t.append(lt[i-1])
            align_i.append(gap_seq)
            i -= 1
        elif sc_current == (sc_up + gp_j):
            align_t.append(gap_seq)
            align_i.append(li[j-1])
            j -= 1

    # If space left fill it with gaps:
    while i > 0:
        asr_align.append(gp_i)
        align_t.append(lt[i-1])
        align_i.append(gap_seq)
        i -= 1
    while j > 0:
        asr_align.append(gp_j)
        align_t.append(gap_seq)
        align_i.append(li[j-1])
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
        max_penalty += 2 * len(lt[0])
    else:  # Or in the case of an unknown root, just add the terminal node
        max_penalty += len(lt[0])

    if debug:
        max_len = max(len(s) for slist in [align_t, align_i] for s in slist)
        print '      aligned lineages:'
        print ('          hdist  %'+str(max_len)+'s  %'+str(max_len)+'s')  % ('true', 'inferred')
        for seq_t, seq_i in zip(align_t, align_i):
            def cfn(s): return utils.color('blue', s, width=max_len, padside='right') if s == gap_seq else s
            str_t, str_i = cfn(seq_t), cfn(seq_i)
            hdstr = utils.color('blue', '-', width=3)
            if all(s != gap_seq for s in (seq_t, seq_i)):
                str_i, isnps = utils.color_mutants(seq_t, seq_i, return_isnps=True)
                hdstr = utils.color('red' if len(isnps) > 0 else None, '%3d' % len(isnps))
            print '           %s  %s   %s' % (hdstr, str_t, str_i)
        if alignment_score % 1 != 0 or max_penalty % 1 != 0:
            raise Exception('alignment_score score %s or max_penalty %s not integers, so need to fix dbg print in next line' % (alignment_score, max_penalty))
        print '      alignment score: %.0f   max penalty: %.0f' % (alignment_score, max_penalty)
        if test:
            check_test_values(alignment_score, max_penalty, align_t, align_i)

    return [align_t, align_i, alignment_score, max_penalty]

# ----------------------------------------------------------------------------------------
def COAR(true_tree, inferred_tree, freq_weighting=False, known_root=True, allow_double_gap=False, debug=False):
    # assert False  # this doesn't work and hasn't been tested, need to go through it carefully before even attempting to run (basically just copied it from bcr-phylo)
# ----------------------------------------------------------------------------------------
# maybe makes sense to do this? (i think not, but saving it here just in case)
# ete_path = '/home/' + os.getenv('USER') + '/anaconda_ete/bin'
# bcr_phylo_path = os.getenv('PWD') + '/packages/bcr-phylo-benchmark'
# outbase = '/fh/fast/matsen_e/dralph/partis/paired-loci/bcr-phylo-antn/v0/seed-0/n-sim-seqs-per-generation-15/obs-times-150/simu/selection/simu/event-0/simu'
# cmd = '%s/bin/validation.py --outbase %s' % (bcr_phylo_path, outbase)
# utils.run_ete_script(cmd, ete_path)
# sys.exit()

# ----------------------------------------------------------------------------------------
    norm_lineage_dist = list()
    nlineages = 0
    for node_t in true_tree.leaf_node_iter():
# TODO TODO remove frequency stuff (?)
        if debug:
            print '%s             %3d %s' % (node_t.taxon.label, len(node_t.seq), node_t.seq)
        # if not node.frequency > 0:
        #     continue

        aln_res = align_lineages(node_t, true_tree, inferred_tree, known_root=known_root, allow_double_gap=allow_double_gap, debug=debug)
        if aln_res is False:  # Skip lineages less than three members long
            continue
        align_t, align_i, final_score, max_penalty = aln_res
        if freq_weighting is True:
            assert False  # don't have frequencies in this version (see old bcr-phylo implementation)
            total_max_penalty = max_penalty * node_t.frequency
            total_lineage_dist = final_score * node_t.frequency
        else:
# TODO TODO what's the point of total_max_penalty?
            total_max_penalty = max_penalty
            total_lineage_dist = final_score
        if total_max_penalty < 0:
            norm_lineage_dist.append(total_lineage_dist / total_max_penalty)
            if debug:
                print '    normalized dist: %.3f' % norm_lineage_dist[-1]

    if len(norm_lineage_dist) == 0:  # There can be total_max_penalty == 0 when all lineages have less than three members
        if debug:
            print '  all lineages shorter than 3, returning 0'
        return 0
    # Take the mean of the distances:
    mean_norm_lineage_dist = sum(norm_lineage_dist) / len(norm_lineage_dist)
    if debug:
        print '  mean over %d lineages: %.3f' % (len(norm_lineage_dist), mean_norm_lineage_dist)
    return mean_norm_lineage_dist
