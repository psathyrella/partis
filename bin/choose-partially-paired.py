#!/usr/bin/env python
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import argparse
import colored_traceback.always
import operator

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
from clusterpath import ClusterPath
import seqfileopener

helpstr = """
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('igh_fname')
parser.add_argument('igk_fname')
parser.add_argument('igl_fname')
parser.add_argument('--input-metafname')
parser.add_argument('--n-to-choose', type=int, default=3)
parser.add_argument('--choose-paired', action='store_true')
args = parser.parse_args()

cpaths, antn_lists = {}, {}
for ltmp, fn in zip(['igh', 'igk', 'igl'], [args.igh_fname, args.igk_fname, args.igl_fname]):
    _, antn_lists[ltmp], cpaths[ltmp] = utils.read_output(fn)
    for tline in antn_lists[ltmp]:
        tline['paired-uids'] = [[] for _ in tline['unique_ids']]
    if args.input_metafname is not None:
        seqfileopener.read_single_input_metafo(args.input_metafname, antn_lists[ltmp], debug=True)

antn_dicts = {l : utils.get_annotation_dict(alist) for l, alist in antn_lists.items()}
for hclust in sorted(cpaths['igh'].best(), key=len, reverse=True):
    hline = antn_dicts['igh'][':'.join(hclust)]
    h_paired_ids, all_pids = zip(*[(u, pids[0]) for u, pids in zip(hline['unique_ids'], hline['paired-uids']) if len(pids)==1])
    if len(all_pids) == 0:
        continue
    lclusts = []
    for ltmp in ['igk', 'igl']:
        for lc in cpaths[ltmp].best():
            if len(set(lc) & set(all_pids)) > 0:
                lclusts.append((ltmp, lc))
    if len(lclusts) != 1:
        print '    %s couldn\'t find unique light cluster for paired ids %s (from heavy ids %s)' % (utils.color('yellow', 'warning'), ' '.join(all_pids), ' '.join(h_paired_ids))
        continue
    l_locus, lclust = lclusts[0]
    lline = antn_dicts[l_locus][':'.join(lclust)]

    # add aa-cdist (it's probably usually already there, but it's easy to add, and should always end up the same)
    import treeutils
    chosen_ids = []
    for tline in [hline, lline]:
        tline['tree-info'] = {'lb' : {}}
        treeutils.add_cdists_to_lbfo(tline, tline['tree-info']['lb'], 'cons-dist-aa')
        cids, cdists = zip(*sorted(tline['tree-info']['lb']['cons-dist-aa'].items(), key=operator.itemgetter(1), reverse=True))
        chosen_ids += cids[:args.n_to_choose]
        print '    %s: chose %d with aa-cdists %s' % (tline['loci'][0], len(cids), ' '.join(str(d) for d in cdists))

    for tline in [hline, lline]:
        utils.print_reco_event(tline, extra_print_keys=['paired-uids', 'cons-dist-aa'], queries_to_emphasize=chosen_ids)
    sys.exit()
