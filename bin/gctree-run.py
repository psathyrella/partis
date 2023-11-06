#!/usr/bin/env python
from __future__ import absolute_import, division
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

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/test', '')
sys.path.insert(1, partis_dir) #'./python')
import python.utils as utils
import python.glutils as glutils
import python.treeutils as treeutils

# ----------------------------------------------------------------------------------------
def get_inf_int_name(gname):  # <gname> is just an integer, which won't be unique and will break things
    return '%s-%s' % (args.inf_int_label, gname)

# ----------------------------------------------------------------------------------------
def getpathcmd():
    cmds = ['#!/bin/bash']
    cmds += ['. %s/etc/profile.d/conda.sh' % args.condapath]  # NOTE have to update conda (using the old version in the next two lines) in order to get this to work
    cmds += ['export PATH=%s/bin:$PATH' % args.condapath]
    # cmds += ['export PYTHONNOUSERSITE=True']  # otherwise it finds the pip-installed packages in .local and breaks (see https://github.com/conda/conda/issues/448)
    return cmds

# ----------------------------------------------------------------------------------------
def gctofn(ft):
    assert ft in ['tree', 'seqs']
    return '%s/gctree.out.inference.1%s' % (args.outdir, '.nk' if ft=='tree' else '.fasta')

# ----------------------------------------------------------------------------------------
def fofn(ft):
    assert ft in ['tree', 'seqs']
    return '%s/%s%s' % (args.outdir, ft if ft=='tree' else 'inferred-%s'%ft, '.nwk' if ft=='tree' else '.fa')

# ----------------------------------------------------------------------------------------
def idfn():
    return 'idmap.txt'

# ----------------------------------------------------------------------------------------
def install():
    cmds = getpathcmd()

    args.env_label = 'gctree'
    install_dir = partis_dir + '/packages'
    # cmds += ['conda update -n base -c defaults conda']  # this updates conda
    cmds += ['conda create -n %s python=3.9' % args.env_label]
    cmds += ['conda activate %s' % args.env_label]
    cmds += ['conda install -c conda-forge gctree']
    cmds += ['conda install -c bioconda phylip']
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname='/tmp/tmprun.sh', debug=True)

# ----------------------------------------------------------------------------------------
def run_gctree(infname):
    if not args.run_help and utils.output_exists(args, gctofn('tree')):
        return

    cmds = getpathcmd()
    cmds += ['conda activate %s' % args.env_label]
    if args.run_help:
        cmds += ['gctree infer -h']
    else:
        if not os.path.exists(args.infname):
            raise Exception('--infname %s doesn\'t exist' % args.infname)
        utils.prep_dir(args.outdir, wildlings=['outfile', 'outtree'], allow_other_files=True)  # phylip barfs like a mfer if its outputs exist (probably you'll get a KeyError 'naive')
        # cmds += ['gctree &>/dev/null && dnapars &>/dev/null || echo "command \'gctree\' or \'dnapars\' not in path, maybe need to run: gctree-run.py --actions install"']  # ick doesn't quite work, whatever
        cmds += ['cd %s' % args.outdir]
        cmds += ['deduplicate %s --root %s --abundance_file abundances.csv --idmapfile %s > deduplicated.phylip' % (args.infname, args.root_label, idfn())]
        cmds += ['mkconfig deduplicated.phylip dnapars > dnapars.cfg']
        cmds += ['dnapars < dnapars.cfg > dnapars.log']  # NOTE if things fail, look in dnaparse.log (but it's super verbose so we can't print it to std out by default)
        cmds += ['%s/bin/xvfb-run -a gctree infer outfile abundances.csv --root %s --frame 1 --verbose --idlabel' % (utils.get_partis_dir(), args.root_label)]  # --idlabel writes the output fasta file
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
    nfos = [s for s in seqfos if s['name']==args.root_label]
    if len(nfos) != 1:
        print '  %s expected 1 naive seq with label \'%s\' but found %d: %s  (in %s)' % (utils.wrnstr(), args.root_label, len(nfos), ' '.join(n['name'] for n in nfos), gctofn('seqs'))
    seqfos = [s for s in seqfos if s['name'] != args.root_label]  # don't want naive seq in final fasta
    seqfos = [s for s in seqfos if s['name'] not in idm_trns]  # also remove input seqs
    inf_int_trns = []
    for sfo in seqfos:
        inf_int_trns.append((sfo['name'], get_inf_int_name(sfo['name'])))
        sfo['name'] = get_inf_int_name(sfo['name'])

    # read tree
    dtree = treeutils.get_dendro_tree(treefname=gctofn('tree'), debug=args.debug)
    dtree.scale_edges(1. / numpy.mean([len(s['seq']) for s in seqfos]))
    dtree.seed_node.taxon.label = args.root_label
    for gname, onames in idm_trns.items():
        node = dtree.find_node_with_taxon_label(gname)
        if node is None:
            raise Exception('couldn\'t find node with name \'%s\' in tree from gctree in %s' % (gname, gctofn('tree')))
        for onm in onames:
            if node.taxon.label == gname:
                node.taxon.label = onm
                continue
            new_taxon = dendropy.Taxon(onm)  # add duplicates as children with zero-length edges
            dtree.taxon_namespace.add_taxon(new_taxon)
            new_node = dendropy.Node(taxon=new_taxon)
            node.add_child(new_node)
            new_node.edge_length = 0
    treeutils.translate_labels(dtree, inf_int_trns, debug=args.debug)

    if args.debug:
        print '    final tree:'
        print treeutils.get_ascii_tree(dendro_tree=dtree, extra_str='      ', width=350)
    with open(fofn('tree'), 'w') as ofile:
        ofile.write('%s\n' % treeutils.as_str(dtree))
    utils.write_fasta(fofn('seqs'), seqfos)

# ----------------------------------------------------------------------------------------
def convert_pickle_tree():
    assert False  # doesn't work yet
    cmd = '%s/bin/read-bcr-phylo-trees.py --pickle-tree-file %s --newick-tree-file %s/tree.nwk' % (utils.get_partis_dir(), args.infname, args.outdir)
    utils.run_ete_script(cmd, None, conda_path=args.condapath, conda_env='ete3', pyversion='3')

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='run:parse')
parser.add_argument('--infname')
parser.add_argument('--outdir')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--condapath', default=os.getenv('HOME') + '/miniconda3')
parser.add_argument('--env-label', default='gctree')
parser.add_argument('--root-label', default='naive')
parser.add_argument('--inf-int-label', default='inf', help='base name for inferred intermediate seqs (numerical name is appended with -')
parser.add_argument('--run-help', action='store_true', help='run gctree help')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--dry-run', action='store_true')
# parser.add_argument('--ete-path', default='/home/%s/anaconda_ete/bin' % os.getenv('USER') if os.getenv('USER') is not None else None)

args = parser.parse_args()
args.actions = utils.get_arg_list(args.actions, choices=['install', 'run', 'parse', 'convert-pickle-tree'])
args.infname = utils.fpath(args.infname)
args.outdir = utils.fpath(args.outdir)

if 'install' in args.actions:
    install()
if 'run' in args.actions:
    run_gctree(args.infname)
if 'parse' in args.actions:
    parse_output()
if 'convert-pickle-tree' in args.actions:
    convert_pickle_tree()
