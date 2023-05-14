from __future__ import print_function
import numpy as np
# from GCutils import hamming_distance
from random import randint

import utils

'''
copied from bcr-phylo/bin/COAR.py (yes it sucks to duplicate it, but there's too much pickle/ete dependent stuff in there, and i'm not even sure how to create the ete/pickle 'forest' files from the inference methods))
'''

# ----------------------------------------------------------------------------------------
def reconstruct_lineage(tree, node):
    lineage = list()
    while True:
        lineage.append(node.sequence)
        if node.up is None:
            return lineage
        node = node.up

# ----------------------------------------------------------------------------------------
def find_node_by_seq(tree, sequence):
    nodes = [node for node in tree.traverse() if node.sequence == sequence and node.frequency > 0]
    if len(nodes) > 1:
        nodes = [nodes[randint(0, len(nodes)-1)]]
    try:
        assert(len(nodes) == 1)
    except Exception as e:
        print('Nodes list:')
        print(nodes)
        print(sequence)
        print(tree)
        print([(node.frequency, node.name, node.sequence) for node in tree.traverse()])
        print(nodes[0])
        print(nodes[1])

        raise e
    return nodes[0]

# ----------------------------------------------------------------------------------------
def align_lineages(seq, tree_t, tree_i, gap_penalty_pct=0, known_root=True, allow_double_gap=False):
    '''
    Standard implementation of a Needleman-Wunsch algorithm as described here:
    http://telliott99.blogspot.com/2009/08/alignment-needleman-wunsch.html
    https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
    And implemented here:
    https://github.com/alevchuk/pairwise-alignment-in-python/blob/master/alignment.py

    gap_penalty_pct is the gap penalty relative to the sequence length of the sequences on the tree.
    '''
    nt = find_node_by_seq(tree_t, seq)
    lt = reconstruct_lineage(tree_t, nt)
    ni = find_node_by_seq(tree_i, seq)
    li = reconstruct_lineage(tree_i, ni)
    # One lineages must be longer than just the root and the terminal node
    if len(lt) <= 2 and len(li) <= 2:
        return False

    # Gap penalty chosen not too large:
    gap_penalty = -1 * int((len(seq) / 100.0) * gap_penalty_pct)
    assert(gap_penalty <= 0)  # Penalties must be negative
    if gap_penalty == 0:  # If gap penalty is zero only gaps in the shortes sequence will be allowed
        assert(allow_double_gap is False)

    # Generate a score matrix matrix:
    kt = len(lt)
    ki = len(li)
    # Disallow gaps in the longest list:
    if allow_double_gap is False and kt > ki:
        # If true is longer than inferred allow gap only in inferred:
        gap_penalty_i = gap_penalty
        gap_penalty_j = -1 * float('inf')
    elif allow_double_gap is False and kt < ki:
        # If inferred is longer than true allow gap only in true:
        gap_penalty_i = -1 * float('inf')
        gap_penalty_j = gap_penalty
    elif allow_double_gap is False and kt == ki:
        # If lists are equally long no gaps are allowed:
        gap_penalty_i = -1 * float('inf')
        gap_penalty_j = -1 * float('inf')
    else:
        gap_penalty_i = gap_penalty
        gap_penalty_j = gap_penalty

    sc_mat = np.zeros((kt, ki), dtype=np.float64)
    for i in range(kt):
        for j in range(ki):
            # Notice the score is defined by number of mismatches:
            #sc_mat[i, j] = len(lt[i]) - hamming_distance(lt[i], li[j])
            sc_mat[i, j] = -1 * hamming_distance(lt[i], li[j])

###    print(sc_mat)
    # Calculate the alignment scores:
    aln_sc = np.zeros((kt+1, ki+1), dtype=np.float64)
    for i in range(0, kt+1):
        if known_root is True:
            aln_sc[i][0] = -1 * float('inf')
        else:
            aln_sc[i][0] = gap_penalty_i * i
    for j in range(0, ki+1):
        if known_root is True:
            aln_sc[0][j] = -1 * float('inf')
        else:
            aln_sc[0][j] = gap_penalty_j * j
    aln_sc[0][0] = 0  # The top left is fixed to zero
###    print(aln_sc)
    for i in range(1, kt+1):
        for j in range(1, ki+1):
            match = aln_sc[i-1][j-1] + sc_mat[i-1, j-1]
            gap_in_inferred = aln_sc[i-1][j] + gap_penalty_i
            gap_in_true = aln_sc[i][j-1] + gap_penalty_j
            aln_sc[i][j] = max(match, gap_in_inferred, gap_in_true)
###    print(aln_sc)
    # Traceback to compute the alignment:
    align_t, align_i, asr_align = list(), list(), list()
    i, j = kt, ki
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
        elif sc_current == (sc_left + gap_penalty_i):
            align_t.append(lt[i-1])
            align_i.append('-')
            i -= 1
        elif sc_current == (sc_up + gap_penalty_j):
            align_t.append('-')
            align_i.append(li[j-1])
            j -= 1

    # If space left fill it with gaps:
    while i > 0:
        asr_align.append(gap_penalty_i)
        align_t.append(lt[i-1])
        align_i.append('-')
        i -= 1
    while j > 0:
        asr_align.append(gap_penalty_j)
        align_t.append('-')
        align_i.append(li[j-1])
        j -= 1

    max_penalty = 0
    for a, b in zip(align_t, align_i):
        if a == '-' or b == '-':
            max_penalty += gap_penalty
        else:
            max_penalty += -len(a)
    # Notice that the root and the terminal node is excluded from this comparison.
    # by adding their length to the max_penalty:
    if known_root is True:
        max_penalty += 2 * len(lt[0])
    else:  # Or in the case of an unknown root, just add the terminal node
        max_penalty += len(lt[0])

    return [align_t, align_i, alignment_score, max_penalty]

# ----------------------------------------------------------------------------------------
def COAR(true_tree, inferred_tree, freq_weigthing=False, known_root=True, allow_double_gap=False):
    norm_lineage_dist = list()
    nlineages = 0
#    for node in true_tree.tree.traverse():
    for node in true_tree.tree.iter_leaves():  # Iterate only through the leaves
        if not node.frequency > 0:
            continue

        aln_res = align_lineages(node.sequence, true_tree.tree, inferred_tree.tree, known_root=known_root, allow_double_gap=allow_double_gap)
        if aln_res is  False:  # Skip lineages less than three members long
            continue
        align_t, align_i, final_score, max_penalty = aln_res
        if freq_weigthing is True:
            total_max_penalty = max_penalty * node.frequency
            total_lineage_dist = final_score * node.frequency
            # Normalize with the max penalty:
            if total_max_penalty < 0:
                norm_lineage_dist.append(total_lineage_dist/total_max_penalty)
        else:
            total_max_penalty = max_penalty
            total_lineage_dist = final_score
            # Normalize with the max penalty:
            if total_max_penalty < 0:
                norm_lineage_dist.append(total_lineage_dist/total_max_penalty)

    if len(norm_lineage_dist) == 0:  # There can be total_max_penalty == 0 when all lineages have less than three members
        return 0
    # Take the mean of the distances:
    mean_norm_lineage_dist = sum(norm_lineage_dist) / len(norm_lineage_dist)
    return mean_norm_lineage_dist
