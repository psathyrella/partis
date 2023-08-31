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
def add_seqs_to_nodes(ttr, tatn, tfn):
    for node in ttr.preorder_node_iter():
        if node is ttr.seed_node:
# TODO figure out how to remove all the .strip('N') calls
            node.seq = tatn['naive_seq'].strip('N')
        else:
            node.seq = utils.per_seq_val(tatn, 'input_seqs', node.taxon.label).strip('N')  # TODO not sure whether we want indel reversed or not
        if node.seq is None:
            raise Exception('no sequence for node %s from %s' % (node.taxon.label, tfn))

# ----------------------------------------------------------------------------------------
def fix_seqs(atn_t, atn_i):  # inferred annotation has padded seqs, which means when the h and l seqs get smashed together (in paircluster.sumv() called from paircluster.make_fake_hl_pair_antns()) sometimes there's extra Ns that aren't in the true annotation
    # ----------------------------------------------------------------------------------------
    def rm_nnn(tatn, iseq, ireplace):
        print '    removing NNN from seq for %s' % tatn['unique_ids'][iseq]
        tseq = tatn['input_seqs'][iseq].strip('N')
        if tseq[ireplace : ireplace + 3] != 'NNN':
            print 'wtf', tseq[ireplace : ireplace + 3]
        tatn['input_seqs'][iseq] = tseq[ : ireplace] + tseq[ireplace + 3 : ]
    # ----------------------------------------------------------------------------------------
    ireplace = None  # inferred seq index at which to replace NNN
    for inf_index, inf_id in enumerate(atn_i['unique_ids']):
        tru_indices = [i for i, u in enumerate(atn_t['unique_ids']) if u == inf_id]
        if len(tru_indices) == 0:  # probably an inferred ancestral sequence (i.e. with different name in true and inferred annotations)
            if ireplace is None:
                raise Exception('couldn\'t find true index')  # maybe could just return/continue here?
            else:  # already figured out what we need to change, so chagne it
                rm_nnn(atn_i, inf_index, ireplace)
                continue
        tru_index = utils.get_single_entry(tru_indices)
        seq_t, seq_i = [s.strip('N') for s in (atn_t['input_seqs'][tru_index], atn_i['input_seqs'][inf_index])]
        if seq_t == seq_i:
            return  # don't need to fix anything
        else:
            if 'NNN' in seq_i and seq_i.count('NNN') == 1 and seq_i.replace('NNN', '') == seq_t:
                if ireplace is None:
                    ireplace = seq_i.find('NNN')
                    print '    found ireplace %d: %s %s' % (ireplace, seq_i[ireplace : ireplace + 3], seq_i)
                rm_nnn(atn_i, inf_index, ireplace)
            else:
                print 'NNN' in seq_i, seq_i.count('NNN'), seq_i.replace('NNN', '') == seq_t
                print inf_id, atn_t['unique_ids'][tru_index]
                utils.color_mutants(seq_t, seq_i, print_result=True, align_if_necessary=True)
                continue
                raise Exception('seqs don\'t match\n  %s\n  %s' % (seq_t, seq_i))
    if ireplace is not None:  # have to also do naive seq
        tseq = atn_i['naive_seq']
        atn_i['naive_seq'] = tseq[ : ireplace] + tseq[ireplace + 3 : ]

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
            fix_seqs(atn_t, atn_i)
            break
    if atn_i is None:
        raise Exception('couldn\'t find inferred annotation')
    dtree_t, dtree_i = [treeutils.get_dendro_tree(treestr=lbplotting.get_tree_in_line(l, is_true)) for is_true, l in [[True, atn_t], [False, atn_i]]]
# TODO make a lookup table between t/i nodes/seqs
    if args.debug:
        for tstr, ttr in zip(['true', 'inf'], [dtree_t, dtree_i]):
            print '    %4s:' % tstr
            print utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=ttr, width=250))
    for ttr, tatn, tfn in zip([dtree_t, dtree_i], [atn_t, atn_i], [args.true_tree_file, args.inferred_tree_file]):
        add_seqs_to_nodes(ttr, tatn, tfn)
    coar.COAR(dtree_t, dtree_i, debug=args.debug)
