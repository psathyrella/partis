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
args = parser.parse_args()

_, tru_atn_list, _ = utils.read_output(args.true_tree_file)
_, inf_atn_list, _ = utils.read_output(args.inferred_tree_file)

# ----------------------------------------------------------------------------------------
def add_seqs_to_nodes(ttr, tatn, tfn):
    for node in ttr.preorder_node_iter():
        if node is ttr.seed_node:
            node.seq = tatn['naive_seq'].strip('N')
        else:
            node.seq = utils.per_seq_val(tatn, 'input_seqs', node.taxon.label).strip('N')  # TODO not sure whether we want indel reversed or not
        if node.seq is None:
            raise Exception('no sequence for node %s from %s' % (node.taxon.label, tfn))

# ----------------------------------------------------------------------------------------
def trnfn(u): return u + '_contig_igh+igk'
utils.translate_uids(tru_atn_list, trfcn=trnfn, expect_missing=True)
for atn_t in tru_atn_list:
    atn_i = None
    for tatn in inf_atn_list:
        if len(set(atn_t['unique_ids']) & set(tatn['unique_ids'])) > 0:
            print '  found inferred annotation with %d / %d uids in common' % (len(set(atn_t['unique_ids']) & set(tatn['unique_ids'])), len(atn_t['unique_ids']))
            atn_i = tatn
            break
    if atn_i is None:
        raise Exception('couldn\'t find inferred annotation')
    dtree_t, dtree_i = [treeutils.get_dendro_tree(treestr=lbplotting.get_tree_in_line(l, is_true)) for is_true, l in [[True, atn_t], [False, atn_i]]]
    # print utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=dtree, width=250))
    for ttr, tatn, tfn in zip([dtree_t, dtree_i], [atn_t, atn_i], [args.true_tree_file, args.inferred_tree_file]):
        add_seqs_to_nodes(ttr, tatn, tfn)
    coar.COAR(dtree_t, dtree_i)
