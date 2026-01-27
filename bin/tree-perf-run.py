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

import partis.utils as utils
import partis.glutils as glutils
import partis.treeutils as treeutils
import partis.lbplotting as lbplotting
import partis.coar as coar

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
parser.add_argument('--metrics', default='coar:rf:mrca')
parser.add_argument('--n-procs', type=int, help='NOTE not used, just putting here for consistency with other scripts')
parser.add_argument('--overwrite', action='store_true', help='NOTE just for compatibility, not used atm')
parser.add_argument('--itree', type=int, help='only run on tree/annotation with this index')
parser.add_argument('--debug', type=int, default=0)
args = parser.parse_args()
args.metrics = utils.get_arg_list(args.metrics, choices=['coar', 'rf', 'mrca'])

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
                # n_lstrip, n_rstrip = len(cseq_i) - len(cseq_i.lstrip('N')), len(cseq_i) - len(cseq_i.rstrip('N'))  # if this starts causing problems again, it might be worth doing something like this to keep track of n bases removed from each side, and making sure it's the same for all seqs
                if tch not in cs_lens:
                    cs_lens[tch] = len(cseq_i)  # keep track of the h/l seq lengths, so for inferred nodes where we don't know it, we can remove the same bases
                assert cs_lens[tch] == len(cseq_i)  # they should all be the same
            cseq_i = utils.pad_seq_for_translation(atn_i, cseq_i)
            new_seq_i.append(cseq_i)
        return ''.join(new_seq_i).strip('N')
    # ----------------------------------------------------------------------------------------
    def check_seqs(uid, seq_i, seq_t, fix_counts, force=False, dont_fix=False):  # check/fix that any nodes that are in both trees have the same sequence
        fix_counts['total'] += 1
        if seq_t == seq_i:
            return False  # return whether we fixed it or not
        # utils.color_mutants(seq_t, seq_i, align_if_necessary=True, print_result=True)
        seq_i = combine_chain_seqs(uid, seq_i)
        if seq_t is None:
            assert force  # for seqs that are in inferred but not true, we already know we need to fix them (and how)
        else:
            if seq_t != seq_i and not dont_fix:
                print('%s tried to fix %s but seqs still different:' % (utils.wrnstr(), uid))
                utils.color_mutants(seq_t, seq_i, print_result=True, align_if_necessary=True, ref_label='true ', seq_label='inf ')
                assert False  # NOTE if you stop crashing here, you probably need to increment something in fix_counts
            seqs_t[uid] = seq_t
        seqs_i[uid] = seq_i
        fix_counts['fixed'].append(uid)
        return True
    # ----------------------------------------------------------------------------------------
    def check_all_lengths(seqs_t, seqs_i):  # check/fix that all seqs in both trees have the same length
        lens_t, lens_i = [list(set(len(s) for s in slist.values())) for slist in [seqs_t, seqs_i]]
        true_len = utils.get_single_entry(lens_t)
        if len(lens_i) == 1 and lens_i[0] == true_len:
            return
        tseq = list(seqs_t.values())[0]
        for uid in [u for u, s in seqs_i.items() if len(s) != true_len]:
            utils.color_mutants(tseq, seqs_i[uid], align_if_necessary=True, print_result=True, extra_str='        ', ref_label='arb. true ', seq_label=uid+' ')
            _, new_seq = utils.align_seqs(tseq, seqs_i[uid])  # i added this to fix a case that i ended up fixing a different (much better) way, but it might be useful in future, so leaving here
            seqs_i[uid] = new_seq.replace('-', utils.ambig_base)  # UGH
            utils.color_mutants(tseq, seqs_i[uid], align_if_necessary=True, print_result=True, extra_str='        ', ref_label='arb. true ', seq_label=uid+' ')
        raise Exception('different sequence lengths (probably from inferred internal nodes), see previous lines')
    # ----------------------------------------------------------------------------------------
    leaf_ids_t = [l.taxon.label for l in tr_t.leaf_node_iter() if l.taxon.label in atn_t['unique_ids']]
    leaf_ids_i = [u for u in leaf_ids_t if u in atn_i['unique_ids']]  # inferred tree may swap internal/leaf nodes
    if set(leaf_ids_i) != set(leaf_ids_t):
        only_true, only_inf = set(leaf_ids_t) - set(leaf_ids_i), set(leaf_ids_i) - set(leaf_ids_t)
        print('    %s inferred leaf ids not the same as true leaf ids when trying to fix seqs (this is probably ok, since the coar calculation will probably skip them).\n      %d extra true: %s\n      %d extra inf: %s' % (utils.wrnstr(), len(only_true), ' '.join(only_true), len(only_inf), ' '.join(only_inf)))
    common_leaf_ids = set(leaf_ids_t) & set(leaf_ids_i)  # maybe missing ones would be ok? but don't want to mess with it, and for now we assume below that they're the same
    seqs_t, seqs_i = [{u : utils.per_seq_val(atn, seq_key, u).strip('N') for u in atn['unique_ids']} for atn in (atn_t, atn_i)]
    seqs_t[naive_name], seqs_i[naive_name] = [a['naive_seq'].strip('N') for a in (atn_t, atn_i)]
    fixed, cs_lens, fix_counts = None, {}, {'fixed' : [], 'total' : 0}
    for uid in common_leaf_ids:
        tfx = check_seqs(uid, seqs_i[uid], seqs_t[uid], fix_counts)
        if fixed is None:
            fixed = tfx
        assert tfx == fixed  # if we fix one, we should fix all of them
    if fixed:
        for uid in [u for u in atn_i['unique_ids'] if u not in leaf_ids_i] + [naive_name]:  # need to also fix any internal/inferred nodes
            check_seqs(uid, seqs_i[uid], seqs_t.get(uid), fix_counts, force=True, dont_fix=uid==naive_name)
    print('    no nodes needed fixing (all seqs already the same for common true/inferred nodes)' if len(fix_counts['fixed'])==0 else '    fixed %d / %d nodes' % (len(fix_counts['fixed']), fix_counts['total']))
    check_all_lengths(seqs_t, seqs_i)
    if debug and len(fix_counts['fixed']) > 0:
        print('      fixed seqs: %s' % ' '.join(sorted(fix_counts['fixed'])))

    return seqs_t, seqs_i

# ----------------------------------------------------------------------------------------
def get_n_parsimony_trees(n_clusters):
# other way to get this number:
#              with open('gctree_base.inference.parsimony_forest.p', 'rb') as fh:
#                  forest = pickle.load(fh)
#              n_parsimony_trees = forest._forest.count_histories()
    n_ptree_list = []
    for iclust in range(n_clusters):
        logfn = '%s/%s/iclust-%d/log' % (os.path.dirname(args.inferred_tree_file), os.path.basename(args.inferred_tree_file).replace('-annotations.yaml', ''), iclust)
        out, err = utils.simplerun('grep "number of trees with integer branch lengths:" %s ' % logfn, shell=True, return_out_err=True, debug=False)
        n_ptree_list.append(int(out.split()[-1]))
    return n_ptree_list

# ----------------------------------------------------------------------------------------
# don't need this now that i'm using --simultaneous-true-clonal-seqs (yes, ick)
def trnfn(u): return u + '_contig_igh+igk'
utils.translate_uids(tru_atn_list, trfcn=trnfn, expect_missing=True)

# ----------------------------------------------------------------------------------------
jvals = {'coar' : [], 'rf' : [], 'mrca' : []}
for itree, atn_t in enumerate(tru_atn_list):
    if args.itree is not None and itree != args.itree:
        continue
    print('  %d: starting true annotation with size %d' % (itree, len(atn_t['unique_ids'])))
    atn_i = None
    for tatn in inf_atn_list:
        common_ids = set(atn_t['unique_ids']) & set(tatn['unique_ids'])
        if len(common_ids) > 0:
            estr = '' if not args.debug else ' (missing %d: %s)' % (len(atn_t['unique_ids']) - len(common_ids), ' '.join(sorted(set(atn_t['unique_ids']) - common_ids)))
            print('    found inferred annotation with %d / %d uids in common%s' % (len(common_ids), len(atn_t['unique_ids']), estr))
            atn_i = tatn
            break
    if atn_i is None:
        raise Exception('couldn\'t find inferred annotation for true annotation (looked in %d inferred annotations, maybe try uncommenting translation above): %s' % (len(inf_atn_list), ' '.join(atn_t['unique_ids'])))
    dtree_t, dtree_i = [treeutils.get_dendro_tree(treestr=lbplotting.get_tree_in_line(l, is_true)) for is_true, l in [[True, atn_t], [False, atn_i]]]
    if args.debug:
        for tstr, ttr in zip(['true', 'inf'], [dtree_t, dtree_i]):
            print('    %4s:' % tstr)
            print(utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=ttr, width=250)))  # , label_fcn=lambda l: l.replace('_contig_igh+igk', '')
    seqs_t, seqs_i = fix_seqs(atn_t, atn_i, dtree_t, dtree_i, debug=args.debug)
    for ttr, seqdict, tfn in zip([dtree_t, dtree_i], [seqs_t, seqs_i], [args.true_tree_file, args.inferred_tree_file]):
        add_seqs_to_nodes(ttr, seqdict, tfn)
    if 'coar' in args.metrics:
        jvals['coar'].append(coar.COAR(dtree_t, dtree_i, known_root=False, debug=args.debug))
    if 'mrca' in args.metrics:
        jvals['mrca'].append(treeutils.mrca_dist(dtree_t, dtree_i, debug=args.debug))
    if 'rf' in args.metrics:
        dts_t, dts_i = treeutils.sync_taxon_namespaces(dtree_t, dtree_i, only_leaves=True) #, debug=True)
        # this is weighted (i.e. depends on edge length), could also use unweighted (fcn symmetric_difference()) [from /loc/dralph/.local/lib/python3.6/site-packages/dendropy/calculate/treecompare.py]
        jvals['rf'].append(dendropy.calculate.treecompare.weighted_robinson_foulds_distance(dts_t, dts_i))
        # print(treeutils.get_ete_rf(dtree_t, dtree_i)

# if os.path.basename(args.inferred_tree_file).split('-')[0] == 'gctree':
#     jvals['n-pars-trees'] = get_n_parsimony_trees(len(tru_atn_list))

if args.outdir is None:
    print('  %s no --outdir specified, so not writing anything' % utils.wrnstr())
    sys.exit(0)

ofn = '%s/tree-perf-vals.yaml' % args.outdir
print('  writing tree perf values to %s' % ofn)
if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)
utils.jsdump(ofn, jvals)
