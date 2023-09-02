#!/usr/bin/env python
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import argparse
import colored_traceback.always
import json

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
import treeutils
import lbplotting
import coar

# ----------------------------------------------------------------------------------------
helpstr = """
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('true_tree_file', help='')
parser.add_argument('inferred_tree_file', help='')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()

_, tru_atn_list, _ = utils.read_output(args.true_tree_file)
_, inf_atn_list, _ = utils.read_output(args.inferred_tree_file)

# ----------------------------------------------------------------------------------------
def add_seqs_to_nodes(ttr, seqdict, tfn):
    for node in ttr.preorder_node_iter():
        node.seq = seqdict['naive' if node is ttr.seed_node else node.taxon.label]

# ----------------------------------------------------------------------------------------
def fix_seqs(atn_t, atn_i, tr_t, tr_i, seq_key='input_seqs', debug=False):  # inferred annotation has padded seqs, which means when the h and l seqs get smashed together (in paircluster.sumv() called from paircluster.make_fake_hl_pair_antns()) sometimes there's extra Ns that aren't in the true annotation
    # ----------------------------------------------------------------------------------------
    def combine_chain_seqs(uid, seq_i):
        new_seq_i = []
        for tch in 'hl':
            cseq_i = utils.per_seq_val(atn_i, '%s_seqs'%tch, uid)
            if cseq_i is None:
                if tch == 'h':
                    cseq_i = seq_i[ : cs_lens[tch]]
                else:
                    cseq_i = seq_i[cs_lens['h'] : ]
                cseq_i = cseq_i.strip('N')
            else:
                cseq_i = cseq_i.strip('N')
                if tch not in cs_lens:
                    cs_lens[tch] = len(cseq_i)  # keep track of the h/l seq lengths, so for inferred nodes where we don't know it, we can remove the same bases
                assert cs_lens[tch] == len(cseq_i)  # they should all be the same
            new_seq_i.append(cseq_i)
            print '  x', cseq_i, seq_i
        return ''.join(new_seq_i)
    # ----------------------------------------------------------------------------------------
    def check_seqs(uid, seq_i, seq_t, force=False):
        print uid
        if seq_t == seq_i:
            return False  # return whether we fixed it or not
        seq_i = combine_chain_seqs(uid, seq_i)
        if seq_t is None:
            assert force  # for seqs that are in inferred but not true, we already know we need to fix them (and how)
        else:
            if seq_t != seq_i:
                print uid
                utils.color_mutants(seq_t, seq_i, print_result=True, align_if_necessary=True)
                assert False
            seqs_t[uid] = seq_t
        seqs_i[uid] = seq_i
        if debug:
            print '    fixed seq for %s' % uid
        return True
    # ----------------------------------------------------------------------------------------
    leaf_ids_t = [l.taxon.label for l in tr_t.leaf_node_iter() if l.taxon.label in atn_t['unique_ids']]
    leaf_ids_i = [u for u in leaf_ids_t if u in atn_i['unique_ids']]  # inferred tree may swap internal/leaf nodes
    assert set(leaf_ids_t) == set(leaf_ids_i)  # maybe missing ones would be ok? but don't want to mess with it, and for now we assume below that they're the same
    seqs_t, seqs_i = [{u : utils.per_seq_val(atn, seq_key, u).strip('N') for u in atn['unique_ids']} for atn in (atn_t, atn_i)]
    seqs_t['naive'], seqs_i['naive'] = [a['naive_seq'].strip('N') for a in (atn_t, atn_i)]
    fixed, cs_lens = None, {}
    for uid in leaf_ids_i:
        tfx = check_seqs(uid, seqs_i[uid], seqs_t[uid])
        if fixed is None:
            fixed = tfx
        assert tfx == fixed  # if we fix one, we should fix all of them
    if fixed:
        for uid in [u for u in atn_i['unique_ids'] if u not in leaf_ids_i] + ['naive']:  # need to also fix any internal/inferred nodes
            check_seqs(uid, seqs_i[uid], seqs_t.get(uid), force=True)

    return seqs_t, seqs_i

# ----------------------------------------------------------------------------------------
def trnfn(u): return u + '_contig_igh+igk'
utils.translate_uids(tru_atn_list, trfcn=trnfn, expect_missing=True)
for atn_t in tru_atn_list:
    print '  starting true annotation with size %d' % len(atn_t['unique_ids'])
    atn_i = None
    for tatn in inf_atn_list:
        if len(set(atn_t['unique_ids']) & set(tatn['unique_ids'])) > 0:
            print '    found inferred annotation with %d / %d uids in common' % (len(set(atn_t['unique_ids']) & set(tatn['unique_ids'])), len(atn_t['unique_ids']))
            atn_i = tatn
            break
    if atn_i is None:
        raise Exception('couldn\'t find inferred annotation')
    dtree_t, dtree_i = [treeutils.get_dendro_tree(treestr=lbplotting.get_tree_in_line(l, is_true)) for is_true, l in [[True, atn_t], [False, atn_i]]]
    seqs_t, seqs_i = fix_seqs(atn_t, atn_i, dtree_t, dtree_i, debug=args.debug)
    if args.debug:
        for tstr, ttr in zip(['true', 'inf'], [dtree_t, dtree_i]):
            print '    %4s:' % tstr
            print utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=ttr, width=250))
    for ttr, seqdict, tfn in zip([dtree_t, dtree_i], [seqs_t, seqs_i], [args.true_tree_file, args.inferred_tree_file]):
        add_seqs_to_nodes(ttr, seqdict, tfn)
    coar.COAR(dtree_t, dtree_i, debug=False) #args.debug)
