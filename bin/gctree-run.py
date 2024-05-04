#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import numpy
import csv
import yaml
import time
import colored_traceback.always
import argparse
import subprocess
import sys
import os
import dendropy
import json
from io import open

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir) #'./python')
import python.utils as utils
import python.glutils as glutils
import python.treeutils as treeutils

# ----------------------------------------------------------------------------------------
def get_inf_int_name(gname):  # <gname> is just an integer, which won't be unique and will break things
    return '%s-%s' % (args.inf_int_label, gname)

# ----------------------------------------------------------------------------------------
def gctofn(ft):
    ftstrs = {
        'tree' : 'gctree.out.inference.1.nk',
        'seqs' : 'gctree.out.inference.1.fasta',
        'dnapars' : 'outfile',
    }
    return '%s/%s' % (args.outdir, ftstrs[ft])

# ----------------------------------------------------------------------------------------
def fofn(ft):
    assert ft in ['tree', 'seqs']
    return '%s/%s%s' % (args.outdir, ft if ft=='tree' else 'inferred-%s'%ft, '.nwk' if ft=='tree' else '.fa')

# ----------------------------------------------------------------------------------------
def idfn():
    return 'idmap.txt'

# ----------------------------------------------------------------------------------------
def install():
    cmds = ['#!/bin/bash']
    cmds += utils.mamba_cmds(args.env_label, only_prep=True)
    cmds += ['micromamba create -n %s python=3.9' % args.env_label]  # 3.10 currently has problems with ete
    cmds += ['micromamba activate %s' % args.env_label]
    cmds += ['micromamba install -c bioconda phylip']
    cmds += ['micromamba install -c conda-forge gctree click']
    # micromamba remove -n gctree --all  # to nuke it and start over
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname='/tmp/tmprun.sh', debug=True)

# ----------------------------------------------------------------------------------------
def update():
    cmds = ['#!/bin/bash']
    cmds += utils.mamba_cmds(args.env_label)
    cmds += ['micromamba update phylip gctree click']
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname='/tmp/tmprun.sh', debug=True)

# ----------------------------------------------------------------------------------------
def add_mfo(tcmd, mfn):
    kdict = {'frame' : 'frame', 'h_frame' : 'frame', 'l_frame' : 'frame2', 'l_offset' : 'chain_split'}  # translates from metafo dict to gctree command line args
    with open(args.metafname) as mfile:
        metafo = json.load(mfile)
    for tk, tc in kdict.items():
        if tk in metafo:
            tcmd += ' --%s %d' % (tc, metafo[tk])
    return tcmd

# ----------------------------------------------------------------------------------------
def run_gctree():
    # ----------------------------------------------------------------------------------------
    def get_gctree_cmd():
        tcmd = '%s/bin/xvfb-run -a gctree infer outfile abundances.csv --root %s --verbose --idlabel' % (utils.get_partis_dir(), args.root_label)  # --idlabel writes the output fasta file
        if not args.base_model:
            tcmd += ' --mutability %s/HS5F_Mutability.csv --substitution %s/HS5F_Substitution.csv' % (args.data_dir, args.data_dir)
        if args.ranking_coeffs is not None:
            tcmd += ' --ranking_coeffs %s' % (' '.join(c for c in args.ranking_coeffs))
        if args.branching_process_ranking_coeff is not None:
            tcmd += ' --branching_process_ranking_coeff %d' % args.branching_process_ranking_coeff
        if os.path.exists(args.metafname):
            tcmd = add_mfo(tcmd, args.metafname)
        return tcmd
    # ----------------------------------------------------------------------------------------
    def get_cmds():
        cmds = ['#!/bin/bash']
        cmds += utils.mamba_cmds(args.env_label)
        if args.run_help:
            cmds += ['gctree infer -h']
            return cmds
        if not os.path.exists(args.infname):
            raise Exception('--infname %s doesn\'t exist' % args.infname)
        cmds += ['cd %s' % args.outdir]
        if args.input_forest_dir is None:
            ofn = '%s/outfile' % args.outdir  # dnapars output file (this is what takes the longest to make
            if os.path.exists(ofn) and os.stat(ofn).st_size > 0:
                print('    dnapars output already exists, not rerunning: %s' % ofn)
            else:
                if os.path.exists(ofn) and os.stat(ofn).st_size == 0:
                    print('    removing zero length dnapars output %s' % ofn)
                    utils.prep_dir(args.outdir, wildlings=['outfile', 'outtree'], allow_other_files=True)  # phylip barfs like a mfer if its outputs exist (probably you'll get a KeyError 'naive')
                cmds += ['deduplicate %s --root %s --abundance_file abundances.csv --idmapfile %s > deduplicated.phylip' % (args.infname, args.root_label, idfn())]
                cmds += ['mkconfig deduplicated.phylip dnapars > dnapars.cfg']
                cmds += ['dnapars < dnapars.cfg > dnapars.log']  # NOTE if things fail, look in dnaparse.log (but it's super verbose so we can't print it to std out by default)
        else:
            print('    --input-forest-dir: copying abundance, idmap, and forest files from %s' % args.input_forest_dir)
            cmds += ['cp %s/{abundances.csv,%s,outfile} %s/' % (args.input_forest_dir, idfn(), args.outdir)]
        if not args.only_write_forest:
            cmds.append(get_gctree_cmd())
        return cmds
    # ----------------------------------------------------------------------------------------
    if not args.run_help and utils.output_exists(args, gctofn('dnapars' if args.only_write_forest else 'tree')):
        return

    cmds = get_cmds()  # also preps dir + other stuff

    utils.simplerun('\n'.join(cmds) + '\n', cmdfname=args.outdir + '/run.sh', print_time='gctree', debug=True, dryrun=args.dry_run)
    if args.run_help:
        sys.exit()

# ----------------------------------------------------------------------------------------
def parse_output():
    if utils.output_exists(args, fofn('seqs')):
        return

    # read translations (this only includes input sequences, not inferred intermediates)
    idm_trns = {}
    with open('%s/idmap.txt' % args.outdir) as idfile:
        reader = csv.DictReader(idfile, fieldnames=('name', 'orig_names'))
        for line in reader:
            if line['orig_names'] == '':
                continue
            idm_trns[line['name']] = line['orig_names'].split(':')

    # read fasta (mostly for inferred intermediate seqs)
    seqfos = utils.read_fastx(gctofn('seqs'), look_for_tuples=True)
    if any(s['name']=='' for s in seqfos):
        n_removed = len([s for s in seqfos if s['name']==''])
        seqfos = [s for s in seqfos if s['name']!='']
        print('  %s removed %d seqs with zero-length names \'\' (I\'m *not* sure this is the right thing to do, but it just kicked this error when I was doing the python 3 conversion)' % (utils.wrnstr(), n_removed))
    nfos = [s for s in seqfos if s['name']==args.root_label]
    if len(nfos) != 1:
        print('  %s expected 1 naive seq with label \'%s\' but found %d: %s  (in %s)' % (utils.wrnstr(), args.root_label, len(nfos), ' '.join(n['name'] for n in nfos), gctofn('seqs')))
    seqfos = [s for s in seqfos if s['name'] != args.root_label]  # don't want naive seq in final fasta
    seq_len = numpy.mean([len(s['seq']) for s in seqfos])
    if not args.expand_all_nodes:  # also remove input seqs (well, gctree's new names for input seqs), unless we're expanding all nodes, in which case we need the gctree-named-nodes as fake new internal nodes
        seqfos = [s for s in seqfos if s['name'] not in idm_trns]
    if len(seqfos) == 0:
        print('  %s no inferred sequences (all seqs read from gctree output were input seqs' % utils.wrnstr())
    inf_int_trns = []
    for sfo in seqfos:
        inf_int_trns.append((sfo['name'], get_inf_int_name(sfo['name'])))
        sfo['name'] = get_inf_int_name(sfo['name'])

    # read tree
    dtree = treeutils.get_dendro_tree(treefname=gctofn('tree'), debug=args.debug)
    dtree.scale_edges(1. / seq_len)
    dtree.seed_node.taxon.label = args.root_label
    for gname, onames in idm_trns.items():
        node = dtree.find_node_with_taxon_label(gname)
        if node is None:
            raise Exception('couldn\'t find node with name \'%s\' in tree from gctree in %s' % (gname, gctofn('tree')))
        if args.debug and len(onames) > 1:
            print('    abundance > 1 for %s: %d (%s)' % (gname, len(onames), ' '.join(onames)))
        for onm in onames:
            if node.taxon.label == gname and not args.expand_all_nodes:
                node.taxon.label = onm
                if args.debug and len(onames) > 1:
                    print('        setting node to %s' % onm)
                continue
            new_taxon = dendropy.Taxon(onm)  # add duplicates as children with zero-length edges
            dtree.taxon_namespace.add_taxon(new_taxon)
            new_node = dendropy.Node(taxon=new_taxon)
            node.add_child(new_node)
            new_node.edge_length = 0
            if args.debug and len(onames) > 1:
                print('        adding child node %s' % onm)
    treeutils.translate_labels(dtree, inf_int_trns, expect_missing=True, debug=args.debug)

    if args.debug:
        print('    final tree:')
        print(treeutils.get_ascii_tree(dendro_tree=dtree, extra_str='      ', width=350))
    with open(fofn('tree'), 'w') as ofile:
        ofile.write('%s\n' % treeutils.as_str(dtree))
    utils.write_fasta(fofn('seqs'), seqfos)

# ----------------------------------------------------------------------------------------
def convert_pickle_tree():
    assert False  # doesn't work yet
    cmd = '%s/bin/read-bcr-phylo-trees.py --pickle-tree-file %s --newick-tree-file %s/tree.nwk' % (utils.get_partis_dir(), args.infname, args.outdir)
    # assert False  # this also needs updating
    # utils.run_ete_script(cmd, None, conda_path=args.condapath, conda_env='ete3', pyversion='3')

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='run:parse')
parser.add_argument('--infname')
parser.add_argument('--metafname', help='if you need --frame (v region doesn\'t start at first position) or --chain_split and --frame2 (heavy/light chain smooshed together), pass the info in json format with this arg (see code above for format).')
parser.add_argument('--outdir')
parser.add_argument('--only-write-forest', action='store_true', help='only run preparatory steps for gctree, i.e. up through dnapars, to write parsimony forest')
parser.add_argument('--input-forest-dir', help='If set, skips preparatory steps (see --only-write-forest), and looks for \'abundance.csv\' and parsimony forest file (\'outfile\') in the specified dir')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--base-model', action='store_true', help='By default, we pass gctree info for the s5f mutation model; if this is set, we don\'t, and it instead use the base model.')
parser.add_argument('--ranking-coeffs', nargs='+', help='see gctree help')
parser.add_argument('--branching-process-ranking-coeff', type=int, help='see gctree help')
parser.add_argument('--env-label', default='gctree')
parser.add_argument('--root-label', default='naive')
parser.add_argument('--data-dir', default='%s/data/s5f'%utils.get_partis_dir())
parser.add_argument('--inf-int-label', default='inf', help='base name for inferred intermediate seqs (numerical name is appended with -')
parser.add_argument('--expand-all-nodes', action='store_true', help='Gctree collapses duplicate observed seqs into nodes with new names and abundance N > 1. By default, we expand these such that the node is named for one of the observed seqs, and add N-1 (zero-length) children. If this arg is set, however, we leave the node and add N (zero-length) children.')
parser.add_argument('--run-help', action='store_true', help='run gctree help')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--dry-run', action='store_true')

args = parser.parse_args()
if args.only_write_forest and args.input_forest_dir:
    raise Exception('doesn\'t make sense to specify both')
args.actions = utils.get_arg_list(args.actions, choices=['install', 'update', 'run', 'parse', 'convert-pickle-tree'])
args.infname = utils.fpath(args.infname)
args.outdir = utils.fpath(args.outdir)

if 'install' in args.actions:
    install()
if 'update' in args.actions:
    update()
if 'run' in args.actions:
    run_gctree()
if 'parse' in args.actions:
    parse_output()
if 'convert-pickle-tree' in args.actions:
    convert_pickle_tree()
