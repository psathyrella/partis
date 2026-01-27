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
import operator

import partis.utils as utils
import partis.glutils as glutils
from partis.clusterpath import ClusterPath
import partis.seqfileopener as seqfileopener
import partis.indelutils as indelutils
import partis.treeutils as treeutils

# ----------------------------------------------------------------------------------------
def addseq(ltmp, tline, uid, iclust):
    if any(uid==s['name'] for s in chosen_seqs[ltmp]):  # don't add it twice
        return
    if indelutils.has_indels_line(tline, tline['unique_ids'].index(uid)):
        indel_warning_strs.append('  %s shm indels in chosen seq %s, which means you need to decide by hand whether you want to choose the input or indel-reversed seq (indel-reversed is written to output file' % (utils.color('yellow', 'warning'), uid))
    chosen_seqs[ltmp].append({
        'name' : uid,
        'seq' : utils.per_seq_val(tline, 'seqs', uid),
        'locus' : ltmp,
        'igh_iclust' : iclust,
        'aa-cdist' : treeutils.smvals(tline, 'cons-dist-aa', uid=uid)
    })

# ----------------------------------------------------------------------------------------
def translate_paired_ids(ltmp, pids):
    return ['%s-%s-%s' % (args.sample_prefix, ltmp, u) for u in pids]

# ----------------------------------------------------------------------------------------
helpstr = """
I think this was an old script to kind of do a hackey version of approximate bulk pairing (but not really sure, would need to read through it more carefully, atm i'm adding this late and i forget).
see usage: datascripts/meta/qa013-synth/run.sh
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('igh_fname')
parser.add_argument('igk_fname')
parser.add_argument('igl_fname')
parser.add_argument('--input-metafnames')
parser.add_argument('--sample-prefix', default='QA013-10x-pre', help='str that needs to be prepended to paired uid to match the uid in \'paired-uids\' (necessary because when we merge samples in datascripts/preprocess.py we don\'t know which sample each paired id is from)')
parser.add_argument('--outfname')
parser.add_argument('--n-largest-clusters', type=int, default=3)
parser.add_argument('--n-to-choose', type=int, default=2)
parser.add_argument('--choose-paired', action='store_true')
# parser.add_argument('--paired-sample-prefix')
# parser.add_argument('--n-max-queries', type=int, default=-1)  # just for testing DAMMIT can't do this, you have to read all the other chains to find the right cluster
args = parser.parse_args()
args.input_metafnames = utils.get_arg_list(args.input_metafnames)

cpaths, antn_lists = {}, {}
for ltmp, fn in zip(['igh', 'igk', 'igl'], [args.igh_fname, args.igk_fname, args.igl_fname]):
    _, antn_lists[ltmp], cpaths[ltmp] = utils.read_output(fn) #, n_max_queries=args.n_max_queries)
    for tline in antn_lists[ltmp]:
        tline['paired-uids'] = [[] for _ in tline['unique_ids']]
    if args.input_metafnames is not None:
        seqfileopener.read_input_metafo(args.input_metafnames, antn_lists[ltmp])

chosen_seqs = {l : [] for l in utils.sub_loci('ig')}
indel_warning_strs = []
lp_antn_pairs = []  # somewhat similar to paircluster.find_cluster_pairs()
antn_dicts = {l : utils.get_annotation_dict(alist) for l, alist in antn_lists.items()}
sorted_hclusters = sorted(cpaths['igh'].best(), key=len, reverse=True)[:args.n_largest_clusters]
print('  choosing seqs from %d largest igh clusters with sizes %s' % (args.n_largest_clusters, ' '.join(str(len(c)) for c in sorted_hclusters)))
print('             igh    N igh   light  l clust     N chosen')
print('    iclust   size  paired   locus  size    aa-cdist   paired')
for iclust, hclust in enumerate(sorted_hclusters):
    print('    %3d    %5d' % (iclust, len(hclust)), end=' ')
    hline = antn_dicts['igh'][':'.join(hclust)]
    tid_lists = [(u, pids[0]) for u, pids in zip(hline['unique_ids'], hline['paired-uids']) if len(pids)==1]  # NOTE this doesn't check that the pairing info is reciprocal (i.e. that the paired-uids in the light chain correpsond to the h seqs)
    if len(tid_lists) == 0:
        non_zero_pids = [pids for pids in hline['paired-uids'] if len(pids) > 0]
        print('         no uniquely-paired seqs (paired-uids lengths: %s, +%d unpaired)' % (' '.join(str(n) for n in sorted((len(pids) for pids in non_zero_pids), reverse=True)), len(hline['unique_ids']) - len(non_zero_pids)))
        continue
    h_paired_ids, l_paired_ids = list(zip(*tid_lists))
    print(' %3d' % len(h_paired_ids), end=' ')

    lclusts = []
    for ltmp in ['igk', 'igl']:
        l_tmp_ids = translate_paired_ids(ltmp, l_paired_ids)
        for lc in cpaths[ltmp].best():
            if len(set(lc) & set(l_tmp_ids)) > 0:
                lclusts.append((ltmp, lc))
    if len(lclusts) != 1:
        print('         couldn\'t find unique light cluster (found %d) for %d paired ids (from %d heavy ids)' % (len(lclusts), len(l_paired_ids), len(h_paired_ids)))
        continue
    l_locus, lclust = lclusts[0]
    lline = antn_dicts[l_locus][':'.join(lclust)]
    lp_antn_pairs.append((hline, lline))
    print('      %3s  %4d' % (l_locus, len(lclust)), end=' ')

    # add aa-cdist (it's probably usually already there, but it's easy to add, and should always end up the same)
    import partis.treeutils as treeutils
    tmpids = {}
    for ltmp, tline in zip(('igh', l_locus), (hline, lline)):
        tline['tree-info'] = {'lb' : {}}
        treeutils.add_cdists_to_lbfo(tline, tline['tree-info']['lb'], 'cons-dist-aa')
        tmpids[ltmp], _ = list(zip(*sorted(list(tline['tree-info']['lb']['cons-dist-aa'].items()), key=operator.itemgetter(1), reverse=True)))
        tmpids[ltmp] = tmpids[ltmp][:args.n_to_choose]
        for uid in tmpids[ltmp]:
            addseq(ltmp, tline, uid, iclust)
    print('     %2d %2d' % (len(tmpids['igh']), len(tmpids[l_locus])), end=' ')

    if args.choose_paired:
        for ltmp, tline, cids in zip(('igh', l_locus), (hline, lline), (h_paired_ids, translate_paired_ids(l_locus, l_paired_ids))):
            for uid in cids:
                addseq(ltmp, tline, uid, iclust)
        print('     %2d %2d' % (len(h_paired_ids), len(l_paired_ids)), end=' ')
    print('')

if len(indel_warning_strs) > 0:
    print('\n'.join(indel_warning_strs))

all_chosen_ids = [s['name'] for sfos in chosen_seqs.values() for s in sfos]
for hline, lline in lp_antn_pairs:
    print('%s   sizes %d %d  chose %d %d' % (utils.color('green', '-->'), len(hline['unique_ids']), len(lline['unique_ids']), len([u for u in all_chosen_ids if u in hline['unique_ids']]), len([u for u in all_chosen_ids if u in hline['unique_ids']])))
    for tline in [hline, lline]:
        utils.print_reco_event(tline, extra_print_keys=['cons-dist-aa', 'paired-uids'], queries_to_emphasize=all_chosen_ids, extra_str='        ')

if args.outfname is not None:
    print('  writing %d chosen seqs to %s' % (len(chosen_seqs), args.outfname))
    utils.mkdir(args.outfname, isfile=True)
    with open(args.outfname, utils.csv_wmode()) as ofile:
        writer = csv.DictWriter(ofile, sorted(chosen_seqs['igh'][0].keys()))  # NOTE dammit this is way too similar to treeutils.combine_selection_metrics(), i need to maybe split the csv writing code out of there?
        writer.writeheader()
        for ltmp, seqfos in chosen_seqs.items():
            for sfo in seqfos:
                writer.writerow(sfo)
