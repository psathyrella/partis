#!/usr/bin/env python
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
sys.path.insert(1, './python')
import utils
import glutils
import treeutils

# ----------------------------------------------------------------------------------------
def getpathcmd():
    cmds = ['#!/bin/bash']
    cmds += ['. %s/etc/profile.d/conda.sh' % args.condapath]  # NOTE have to update conda (using the old version in the next two lines) in order to get this to work
    cmds += ['export PATH=%s/bin:$PATH' % args.condapath]
    # cmds += ['export PYTHONNOUSERSITE=True']  # otherwise it finds the pip-installed packages in .local and breaks (see https://github.com/conda/conda/issues/448)
    return cmds

# ----------------------------------------------------------------------------------------
def gctofn():
    return '%s/gctree.out.inference.1.nk' % args.outdir

# ----------------------------------------------------------------------------------------
def fofn():
    return '%s/%s' % (args.outdir, args.outfile_basename)

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
    if utils.output_exists(args, gctofn()):
        return

    cmds = getpathcmd()
    cmds += ['conda activate %s' % args.env_label]
    if args.run_help:
        cmds += ['gctree infer -h']
    else:
        if not os.path.exists(args.infname):
            raise Exception('--infname %s doesn\'t exist' % args.infname)
        cmds += ['cd %s' % args.outdir]
        cmds += ['rm -v outfile']  # dnaparse barfs if it's output exists
        cmds += ['deduplicate %s --root %s --abundance_file abundances.csv --idmapfile %s > deduplicated.phylip' % (args.infname, args.root_label, idfn())]
        cmds += ['mkconfig deduplicated.phylip dnapars > dnapars.cfg']
        cmds += ['dnapars < dnapars.cfg > dnapars.log']  # NOTE if things fail, look in dnaparse.log (but it's super verbose so we can't print it to std out by default)
        cmds += ['%s/bin/xvfb-run -a gctree infer outfile abundances.csv --root %s --frame 1 --verbose' % (utils.get_partis_dir(), args.root_label)]
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname=args.outdir + '/run.sh', print_time='gctree', debug=True, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def parse_output(debug=False):
    if utils.output_exists(args, fofn()):
        return
    translations = {}
    with open('%s/idmap.txt' % args.outdir) as idfile:
        reader = csv.DictReader(idfile, fieldnames=('name', 'orig_names'))
        for line in reader:
            translations[line['name']] = line['orig_names'].split(':')
    dtree = treeutils.get_dendro_tree(treefname=gctofn(), debug=True)
    dtree.seed_node.taxon.label = args.root_label
    for gname, onames in translations.items():
        node = dtree.find_node_with_taxon_label(gname)
        if node is None:
            raise Exception('couldn\'t find node with name \'%s\' in tree from gctree in %s' % (gname, gctofn()))
        for onm in onames:
            if node.taxon.label == gname:
                node.taxon.label = onm
                continue
            new_taxon = dendropy.Taxon(onm)  # add duplicates as children with zero-length edges
            dtree.taxon_namespace.add_taxon(new_taxon)
            new_node = dendropy.Node(taxon=new_taxon)
            node.add_child(new_node)
            new_node.edge_length = 0
    if debug:
        print treeutils.get_ascii_tree(dendro_tree=dtree, extra_str='      ', width=350)
    with open(fofn(), 'w') as ofile:
        ofile.write('%s\n' % treeutils.as_str(dtree))

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--install', action='store_true')
parser.add_argument('--infname')
parser.add_argument('--outdir')
parser.add_argument('--outfile-basename', default='tree.nwk')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--condapath', default=os.getenv('HOME') + '/miniconda3')
parser.add_argument('--env-label', default='gctree')
parser.add_argument('--root-label', default='naive')
parser.add_argument('--run-help', action='store_true', help='run gctree help')
parser.add_argument('--dry-run', action='store_true')
args = parser.parse_args()

if args.install:
    install()
    sys.exit()

# ----------------------------------------------------------------------------------------
if args.infname[0] != '/':
    args.infname = '%s/%s' % (os.getcwd(), args.infname)
run_gctree(args.infname)
parse_output()
