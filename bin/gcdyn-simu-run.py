#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import glob
import sys
import csv
from io import open
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import copy
import argparse
import colored_traceback.always
import json

import partis.utils as utils
import partis.paircluster as paircluster
import partis.glutils as glutils
from partis.clusterpath import ClusterPath
import partis.treeutils as treeutils
import partis.indelutils as indelutils
from partis.event import RecombinationEvent

# ----------------------------------------------------------------------------------------
def get_replay_naive_antn(glfos, ltmp, add_empty_mature_keys=False, debug=False):
    antn, naive_seq = {}, []
    for region in utils.getregions(ltmp):
        gene, seq = utils.get_single_entry(list(glfos[ltmp]['seqs'][region].items()))
        antn['%s_gene'%region] = gene
        if ltmp == 'igh' and region == 'j' and args.rm_last_ighj_base:
            if debug:
                print('    removing last ighj base: %s --> %s' % (seq, seq[:-1]))
            seq = seq[:-1]
        naive_seq.append(seq)
    for dstr in utils.all_erosions:
        antn['%s_del'%dstr] = 0
    if not utils.has_d_gene(ltmp):
        antn['d_gene'] = glutils.dummy_d_genes[ltmp]
        antn['d_5p_del'] = 1
    for bstr in utils.all_boundaries:
        antn['%s_insertion'%bstr] = ''
    antn['naive_seq'] = ''.join(naive_seq)
    antn['invalid'] = False
    if add_empty_mature_keys:  # keys that need to be replaced in mature annotation
        antn['unique_ids'] = []
        antn['seqs'] = []
        antn['input_seqs'] = antn['seqs']
        antn['indelfos'] = []
        antn['paired-uids'] = []
    else:
        antn['unique_ids'] = ['replay_naive']
        antn['seqs'] = [''.join(naive_seq)]
        antn['input_seqs'] = antn['seqs']
        antn['indelfos'] = [indelutils.get_empty_indel()]
        antn['paired-uids'] = [[]]
    if debug:
        utils.add_implicit_info(glfos[ltmp], antn)
        utils.print_reco_event(antn)
    return antn

# ----------------------------------------------------------------------------------------
def get_uid(name, ltmp):
    return '%s-%s' % (name, ltmp)

# ----------------------------------------------------------------------------------------
def process_tree(glfos, treelist, tree_sfos, leaf_meta, itree, lp_infos, lpair):
    antns = {}
    for ltmp in lpair:
        antns[ltmp] = get_replay_naive_antn(glfos, ltmp, add_empty_mature_keys=True, debug=args.debug)
        for tkey in args.meta_columns:
            antns[ltmp][utils.input_metafile_keys[tkey]] = []

    hloc = utils.heavy_locus(args.ig_or_tr)
    lloc = utils.get_single_entry([l for l in glfos if l!=hloc])
    h_seq_len = len(antns[hloc]['naive_seq'])
    for sfo in tree_sfos[itree]:
        joint_naive_seq = ''.join(antns[l]['naive_seq'] for l in lpair)
        if len(sfo['seq']) == len(joint_naive_seq) - 1:
            raise Exception('seq read from gcdyn file has len %d, one less than replay naive seq %d, you probably need to set --rm-last-ighj-base' % (len(sfo['seq']), len(joint_naive_seq)))
        assert len(sfo['seq']) == len(joint_naive_seq)
        if 'naive' in sfo['name']:
            assert sfo['seq'] == joint_naive_seq
            continue
        dumb_offset = 0
        if args.rm_last_ighj_base:
            dumb_offset = 1
        antns[hloc]['seqs'].append(sfo['seq'][:h_seq_len + dumb_offset])
        antns[lloc]['seqs'].append(sfo['seq'][h_seq_len:])
    
        for ltmp in lpair:
            assert antns[ltmp]['input_seqs'][-1] == antns[ltmp]['seqs'][-1]
            antns[ltmp]['unique_ids'].append(get_uid(sfo['name'], ltmp))
            antns[ltmp]['indelfos'].append(indelutils.get_empty_indel())
            other_locus = utils.get_single_entry([l for l in lpair if l!=ltmp])
            antns[ltmp]['paired-uids'].append([get_uid(sfo['name'], other_locus)])
            for tkey in args.meta_columns:
                antns[ltmp][utils.input_metafile_keys[tkey]].append(leaf_meta[sfo['name']][tkey])

    for ltmp in lpair:
        dtree = treeutils.get_dendro_tree(treestr=treelist[itree])
        treeutils.translate_labels(dtree, [(s['name'], get_uid(s['name'], ltmp)) for s in tree_sfos[itree] if 'naive' not in s['name']], expect_missing=True)
        antns[ltmp]['tree'] = dtree.as_string(schema='newick').strip()
        tmp_event = RecombinationEvent(glfos[ltmp])  # I don't want to move the function out of event.py right now
        tmp_event.set_reco_id(antns[ltmp], irandom=itree)  # not sure that setting <irandom> here actually does anything
        utils.add_implicit_info(glfos[ltmp], antns[ltmp])  # easiest way to add codon_positions, which we want to write to file

    for ltmp in lpair:
        lp_infos[lpair]['antn_lists'][ltmp].append(antns[ltmp])

    if args.debug:
        for ltmp in sorted(glfos):
            utils.print_reco_event(antns[ltmp], extra_str='        ')
    
# ----------------------------------------------------------------------------------------
def mfname():  # have to look for old meta file name for backwards compatibility
    mfn = '%s/meta.csv'%gcd_dir
    if os.path.exists(mfn):  # if new name exists, return that
        return mfn
    return '%s/leaf-meta.csv'%gcd_dir  # otherwise return old name

# ----------------------------------------------------------------------------------------
def run_gcdyn():
    if os.path.exists('%s/encoded-trees.npy'%gcd_dir):
        print('    gcdyn output exists in %s' % gcd_dir)
        return
    cmds = ['#!/bin/bash']
    cmds += utils.mamba_cmds('gcdyn')
    gcmd = 'gcd-simulate --sample-internal-nodes --label-leaf-internal-nodes --outdir %s --xshift-values 2.5 --xscale-values 5 --yscale-values 1000000' % gcd_dir
    if args.n_sub_procs is not None:
        gcmd += ' --n-sub-procs %d' % args.n_sub_procs
    if args.seed is not None:
        gcmd += ' --seed %d' % args.seed
    if args.n_sim_events is not None:
        gcmd += ' --n-trials %d' % args.n_sim_events
    if args.obs_times is not None:
        gcmd += ' --time-to-sampling-values %d' % args.obs_times
    cmds += [gcmd]
    cmdfname = '%s/run.sh' % gcd_dir
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname=cmdfname)

# ----------------------------------------------------------------------------------------
def process_output():
    glfos = {}
    for ltmp in utils.loci:
        if os.path.exists('%s/%s' % (args.replay_germline_dir, ltmp)):
            glfos[ltmp] = glutils.read_glfo(args.replay_germline_dir, ltmp)

    seqfos = utils.read_fastx('%s/seqs.fasta'%gcd_dir)
    treelist = treeutils.get_treestrs_from_file('%s/trees.nwk'%gcd_dir)
    lmetalines = utils.csvlines(mfname())
    leaf_meta = {l['name'] : {'affinity' : float(l['affinity'])} for l in lmetalines}
    print('  read %d trees, %d seqs  (plus leaf metafo) from %s' % (len(treelist), len(seqfos), gcd_dir))
    tree_sfos, all_uids = {}, set()  # collect up the seqfos for each tree
    for sfo in seqfos:
        if sfo['name'].count('-') == 1:  # old-style, and naive seq
            itree, sname = sfo['name'].split('-')
        elif sfo['name'].count('-') == 2:  # new-style
            itree, listr, sname = sfo['name'].split('-')  # listr is either 'leaf' or 'mrca'
        else:
            assert False
        itree = int(itree)
        if itree not in tree_sfos:
            tree_sfos[itree] = []
        tree_sfos[itree].append(sfo)
        if sfo['name'] in all_uids:
            raise Exception('found uid %s twice' % sfo['name'])
        all_uids.add(sfo['name'])
    print('      %d tree seqfos with lengths: %s' % (len(tree_sfos), ' '.join(str(len(slist)) for slist in sorted(tree_sfos.values(), key=len, reverse=True))))
    if sorted(tree_sfos) != list(range(len(treelist))):
        raise Exception('tree indices from sequence names didn\'t match number of trees')

    lpairs = [tuple(sorted(glfos))]  # could use utils.locus_pairs() if i update to more than one lpair
    lp_infos = {lp : {'antn_lists' : {l : [] for l in glfos}, 'glfos' : {l : g for l, g in glfos.items()}, 'cpaths' : {}} for lp in lpairs}
    for itree in range(len(treelist)):
        process_tree(glfos, treelist, tree_sfos, leaf_meta, itree, lp_infos, utils.get_single_entry(lpairs))

    print('  writing annotations to %s' % args.outdir)
    headers = utils.simulation_headers + args.meta_columns
    def ofn_fn(locus, lpair=None, joint=None):
        return paircluster.paired_fn(args.outdir, locus, lpair=lpair, suffix='.yaml')
    paircluster.write_lpair_output_files(lpairs, lp_infos, ofn_fn, headers=headers)
    glfos, antn_lists, _ = paircluster.concat_heavy_chain(lpairs, lp_infos)  # per-locus glfos with concat'd heavy chain
    paircluster.write_concatd_output_files(glfos, antn_lists, ofn_fn, headers)
    outfos, metafos = paircluster.get_combined_outmetafos(antn_lists, extra_meta_headers=[utils.input_metafile_keys[k] for k in args.meta_columns])
    paircluster.write_combined_fasta_and_meta(args.outdir+'/all-seqs.fa', args.outdir+'/meta.yaml', outfos, metafos)

# ----------------------------------------------------------------------------------------
helpstr = """
run gcdyn simulation, then process results into partis-format paired output dir, for example:
    ./bin/gcdyn-simu-run.py --gcd-dir <gcd> --n-sub-procs 10 --seed 0 --n-trials 10 --rm-last-ighj-base
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('--outdir', required=True)
parser.add_argument('--actions', default='run:process')
parser.add_argument('--input-simu-dir', help='if set, only run \'process\' action, reading gcdyn simu from this dir and writing partis files to --outdir')
parser.add_argument('--replay-germline-dir', default='datascripts/meta/taraki-gctree-2021-10/germlines', help='dir with gcreplay germline sequences')
parser.add_argument('--rm-last-ighj-base', action='store_true', help='sometimes the ighj gene has an extra G at the end, sometimes not, this says to remove it from the seqs read from --replay-germline-dir')
parser.add_argument('--n-sub-procs', type=int)
parser.add_argument('--seed', type=int)
parser.add_argument('--n-sim-events', type=int)
parser.add_argument('--obs-times', type=int)
parser.add_argument('--meta-columns', default='affinity')
parser.add_argument('--ig-or-tr', default='ig')
parser.add_argument('--debug', action='store_true')
args = parser.parse_args()
args.meta_columns = utils.get_arg_list(args.meta_columns, choices=utils.input_metafile_keys.keys())
if args.input_simu_dir is None:
    gcd_dir = '%s/gcdyn' % args.outdir
    args.actions = utils.get_arg_list(args.actions)
else:
    gcd_dir = args.input_simu_dir
    args.actions = ['process']
if 'run' in args.actions:
    run_gcdyn()
if 'process' in args.actions:
    process_output()
