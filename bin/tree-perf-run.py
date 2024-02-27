#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import sys
import csv
from io import open
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import argparse
import colored_traceback.always
import json
import dendropy

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir) # + '/python')

import python.utils as utils
import python.glutils as glutils
import python.treeutils as treeutils
import python.lbplotting as lbplotting
import python.coar as coar

# ----------------------------------------------------------------------------------------
helpstr = """
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('--true-tree-file', required=True, help='partis yaml file with true annotations from which to extract true trees')
parser.add_argument('--inferred-tree-file', required=True, help='partis yaml file with inferred annotations and inferred trees')
parser.add_argument('--outdir')
parser.add_argument('--n-procs', type=int, help='NOTE not used, just putting here for consistency with other scripts')
parser.add_argument('--overwrite', action='store_true', help='NOTE just for compatibility, not used atm')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()

_, tru_atn_list, _ = utils.read_output(args.true_tree_file)
_, inf_atn_list, _ = utils.read_output(args.inferred_tree_file)

naive_name = 'naive'

# ----------------------------------------------------------------------------------------
def add_seqs_to_nodes(ttr, seqdict, tfn):
    for node in ttr.preorder_node_iter():
        node.seq = seqdict[naive_name if node is ttr.seed_node else node.taxon.label]

# ----------------------------------------------------------------------------------------
def fix_seqs(atn_t, atn_i, tr_t, tr_i, seq_key='input_seqs', debug=False):  # inferred annotation has padded seqs, which means when the h and l seqs get smashed together (in paircluster.sumv() called from paircluster.make_fake_hl_pair_antns()) sometimes there's extra Ns that aren't in the true annotation
    # ----------------------------------------------------------------------------------------
    def combine_chain_seqs(uid, seq_i):  # basic idea is that we need to remove any N padding from waterer.py, then repad but just for translation
        new_seq_i = []
        for tch in 'hl':
            cseq_i = utils.per_seq_val(atn_i, '%s_seqs'%tch, uid)
            if cseq_i is None:  # inferred ancestral seqs won't have h/l seqs set (maybe also naive?)
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
            cseq_i = utils.pad_seq_for_translation(atn_i, cseq_i)
            new_seq_i.append(cseq_i)
        return ''.join(new_seq_i).strip('N')
    # ----------------------------------------------------------------------------------------
    def check_seqs(uid, seq_i, seq_t, force=False, dont_fix=False):
        if debug:
            print('      fixing %s:' % uid, end='')
        if seq_t == seq_i:
            if debug:
                print('   nothing to fix (seqs already the same)')
            return False  # return whether we fixed it or not
        seq_i = combine_chain_seqs(uid, seq_i)
        if seq_t is None:
            assert force  # for seqs that are in inferred but not true, we already know we need to fix them (and how)
        else:
            if seq_t != seq_i and not dont_fix:
                print('%s tried to fix %s but seqs still different:' % (utils.wrnstr(), uid))
                utils.color_mutants(seq_t, seq_i, print_result=True, align_if_necessary=True, ref_label='true ', seq_label='inf ')
                assert False
            seqs_t[uid] = seq_t
        seqs_i[uid] = seq_i
        if debug:
            print('    successfully fixed')
        return True
    # ----------------------------------------------------------------------------------------
    leaf_ids_t = [l.taxon.label for l in tr_t.leaf_node_iter() if l.taxon.label in atn_t['unique_ids']]
    leaf_ids_i = [u for u in leaf_ids_t if u in atn_i['unique_ids']]  # inferred tree may swap internal/leaf nodes
    if set(leaf_ids_i) != set(leaf_ids_t):
        only_true, only_inf = set(leaf_ids_t) - set(leaf_ids_i), set(leaf_ids_i) - set(leaf_ids_t)
        print('    %s inferred leaf ids not the same as true leaf ids when trying to fix seqs (this is probably ok, since the coar calculation will probably skip them).\n      %d extra true: %s\n      %d extra inf: %s' % (utils.wrnstr(), len(only_true), ' '.join(only_true), len(only_inf), ' '.join(only_inf)))
    common_leaf_ids = set(leaf_ids_t) & set(leaf_ids_i)  # maybe missing ones would be ok? but don't want to mess with it, and for now we assume below that they're the same
    seqs_t, seqs_i = [{u : utils.per_seq_val(atn, seq_key, u).strip('N') for u in atn['unique_ids']} for atn in (atn_t, atn_i)]
    seqs_t[naive_name], seqs_i[naive_name] = [a['naive_seq'].strip('N') for a in (atn_t, atn_i)]
    fixed, cs_lens = None, {}
    for uid in common_leaf_ids:
        tfx = check_seqs(uid, seqs_i[uid], seqs_t[uid])
        if fixed is None:
            fixed = tfx
        assert tfx == fixed  # if we fix one, we should fix all of them
    if fixed:
        for uid in [u for u in atn_i['unique_ids'] if u not in leaf_ids_i] + [naive_name]:  # need to also fix any internal/inferred nodes
            check_seqs(uid, seqs_i[uid], seqs_t.get(uid), force=True, dont_fix=uid==naive_name)

    return seqs_t, seqs_i

# ----------------------------------------------------------------------------------------
def trnfn(u): return u + '_contig_igh+igk'
utils.translate_uids(tru_atn_list, trfcn=trnfn, expect_missing=True)
jvals = {'coar' : [], 'rf' : []}
for atn_t in tru_atn_list:
    print('  starting true annotation with size %d' % len(atn_t['unique_ids']))
    atn_i = None
    for tatn in inf_atn_list:
        common_ids = set(atn_t['unique_ids']) & set(tatn['unique_ids'])
        if len(common_ids) > 0:
            estr = '' if not args.debug else ' (missing %d: %s)' % (len(atn_t['unique_ids']) - len(common_ids), ' '.join(sorted(set(atn_t['unique_ids']) - common_ids)))
            print('    found inferred annotation with %d / %d uids in common%s' % (len(common_ids), len(atn_t['unique_ids']), estr))
            atn_i = tatn
            break
    if atn_i is None:
        raise Exception('couldn\'t find inferred annotation (looked in %d inferred annotations)' % len(inf_atn_list))
    dtree_t, dtree_i = [treeutils.get_dendro_tree(treestr=lbplotting.get_tree_in_line(l, is_true)) for is_true, l in [[True, atn_t], [False, atn_i]]]
    if args.debug:
        for tstr, ttr in zip(['true', 'inf'], [dtree_t, dtree_i]):
            print('    %4s:' % tstr)
            print(utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=ttr, width=250)))
    seqs_t, seqs_i = fix_seqs(atn_t, atn_i, dtree_t, dtree_i, debug=args.debug)
    for ttr, seqdict, tfn in zip([dtree_t, dtree_i], [seqs_t, seqs_i], [args.true_tree_file, args.inferred_tree_file]):
        add_seqs_to_nodes(ttr, seqdict, tfn)
    cval = coar.COAR(dtree_t, dtree_i, debug=args.debug)
    jvals['coar'].append(cval)
    dtree_t, dtree_i = treeutils.sync_taxon_namespaces(dtree_t, dtree_i, only_leaves=True)
    jvals['rf'].append(dendropy.calculate.treecompare.robinson_foulds_distance(dtree_t, dtree_i))

if args.outdir is None:
    print('  %s no --outdir specified, so not writing anything' % utils.wrnstr())
    sys.exit(0)

ofn = '%s/tree-perf-vals.yaml' % args.outdir
print('  writing tree perf values to %s' % ofn)
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
utils.jsdump(ofn, jvals)