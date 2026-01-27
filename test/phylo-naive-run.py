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

import partis.utils as utils
import partis.glutils as glutils
from partis.clusterpath import ClusterPath
import partis.paircluster as paircluster
from partis.performanceplotter import PerformancePlotter

# ----------------------------------------------------------------------------------------
helpstr = """
Run 'infer-trees' with --phylo-naive-inference-fuzz set to get phyl-inferred naive sequence (i.e. mask uncertain bits of partis naive with Ns, then pass with observed seqs to phylo method.
Then read resulting annotation, which has desired sequence as an observed seq called 'phylo-naive', and compare to simulation truth (if set).
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('--input-partition-dir', required=True, help='dir with partis partition output')
parser.add_argument('--outdir', required=True, help='if set, read simulation from this file and compare to results')
parser.add_argument('--simdir', help='if set, read simulation from this file and compare to results')
parser.add_argument('--phylo-naive-inference-fuzz', type=int, default=1)
parser.add_argument('--tree-inference-method', required=True)
parser.add_argument('--n-procs')  # not used atm
parser.add_argument('--debug', action='store_true')

args = parser.parse_args()

for ltmp in utils.sub_loci('ig'):
    tfn = paircluster.paired_fn(args.simdir, ltmp, suffix='.yaml')
    if not os.path.exists(tfn):
        print('  %s: no simulation file, skipping' % ltmp)
        continue
    pfn = '%s/partition-%s.yaml' % (args.outdir, ltmp)
    phylo_inf_fname = '%s/%s-annotations.yaml' % (utils.getprefix(pfn), args.tree_inference_method)
    if os.path.exists(phylo_inf_fname):
        print('  %s output exists, skipping: %s' % (ltmp, phylo_inf_fname))
        continue
    cmd = './bin/partis partition --simultaneous-true-clonal-seqs --is-simu --infname %s --parameter-dir %s/parameters/%s --outfname %s --locus %s' % (tfn, args.input_partition_dir, ltmp, pfn, ltmp)
    utils.simplerun(cmd)
    cmd = './bin/partis infer-trees --debug 1 --outfname %s --tree-inference-method %s --phylo-naive-inference-fuzz %s' % (pfn, args.tree_inference_method, args.phylo_naive_inference_fuzz)
    utils.simplerun(cmd)

    glfo, inf_antns, cpath = utils.read_output(phylo_inf_fname)

    if args.simdir is not None:
        true_glfo, true_antns, true_cpath = utils.read_output(tfn)
        perfplotter = PerformancePlotter(args.tree_inference_method)
        def get_reco_info(antn_list):
            return {u : utils.synthesize_single_seq_line(l, i, dont_deep_copy=True) for l in antn_list for i, u in enumerate(l['unique_ids'])}  # dont_deep_copy is fine, we're not modifying anything
        reco_info = get_reco_info(true_antns)
        inf_rinfo = get_reco_info(inf_antns)  # single-seq annotations from inferred annotations
        print('  missing %d / %d uids from inf' % (len(set(reco_info) - set(inf_rinfo)), len(reco_info)))
        printed_true_naives = set()
        for inf_line in inf_antns:
            phylo_naive_seq = utils.per_seq_val(inf_line, 'seqs', 'phylo-naive')
            if phylo_naive_seq is None:
                raise Exception('couldn\'t find phylo naive seq')
            tline = reco_info[inf_line['unique_ids'][0]]  # first one usually always be there since inferred ancestral seqs are appended
            if args.debug and phylo_naive_seq != inf_line['naive_seq'] and tline['naive_seq'] not in printed_true_naives:
                utils.print_reco_event(tline, extra_str='         ')
                # utils.print_reco_event(inf_line, extra_str='         ')
                utils.color_mutants(tline['naive_seq'], inf_line['naive_seq'], seq_label='partis   ', ref_label='true    ', print_result=True, align_if_necessary=True, print_n_snps=True)
                utils.color_mutants(tline['naive_seq'], phylo_naive_seq, seq_label='phylo    ', print_result=True, only_print_seq=True, align_if_necessary=True, print_n_snps=True)
                printed_true_naives.add(tline['naive_seq'])
            inf_line['naive_seq'] = phylo_naive_seq  # replace naive seq with phylo-inferred seq NOTE we now have an inconsistent annotation since this doesn't update e.g. deletions (how could we? phylo method doesn't know about those)
            for iseq, uid in enumerate(inf_line['unique_ids']):
                if uid not in reco_info:
                    continue
                perfplotter.evaluate(reco_info[uid], utils.synthesize_single_seq_line(inf_line, iseq), simglfo=true_glfo)
        perfplotter.plot('%s/%s'%(args.outdir, ltmp), only_csv=True)
