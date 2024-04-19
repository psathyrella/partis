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

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/projects', '')
sys.path.insert(1, partis_dir) #'./python')
import python.utils as utils
import python.glutils as glutils
import python.treeutils as treeutils
from python.clusterpath import ClusterPath

# ----------------------------------------------------------------------------------------
def install():
    cmds = ['#!/bin/bash']
    cmds += utils.mamba_cmds(args.env_label, only_prep=True)
    cmds += ['micromamba create -n %s' % args.env_label]
    cmds += ['micromamba activate %s' % args.env_label]  #  python= 3.6 and 3.9 failed, so i let it choose, it chose 3.5 which seems to work
    cmds += ['micromamba install -c bioconda changeo']
    # micromamba remove --all -n args.env_label  # to nuke it and start over
    cmds += ['cd packages']
    cmds += ['git clone https://bitbucket.org/kleinstein/igphyml']
    cmds += ['cd igphyml']
    cmds += ['./make_phyml_omp']
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname='/tmp/tmprun.sh', debug=True)


# ----------------------------------------------------------------------------------------
def igp_ofn(ft, modelstr='hlp'):
    assert ft in ['tree', 'seqs']
    return '%s/input/0.fasta_igphyml_asr_%s.%s' % (args.outdir, modelstr, 'nex' if ft=='tree' else 'fasta')

# ----------------------------------------------------------------------------------------
def run_igphyml():
    afn = '%s/input.tsv' % args.outdir
    rfn = '%s/input_lineages.tsv' % args.outdir
    if utils.all_outputs_exist(args, [afn, rfn, igp_ofn('tree', modelstr='gy'), igp_ofn('seqs', modelstr='gy'), igp_ofn('tree'), igp_ofn('seqs')], outlabel='igphyml'):
        return

    if utils.output_exists(args, afn):  # don't put this in the run file below since partis probably won't work in the igphyml mamba env
        print('    airr tsv exists, not rerunning: %s' % afn)
    else:
        cmd = '%s/bin/parse-output.py %s %s --airr-output' % (partis_dir, args.infname, afn)
        utils.simplerun(cmd)

    cmds = ['#!/bin/bash']
    cmds += utils.mamba_cmds(args.env_label)
    cmds += ['BuildTrees.py -d %s --log %s/build-trees.log --collapse' % (afn, args.outdir)]
    # doesn't correct for shm biases:
    cmds += ['%s/packages/igphyml/src/igphyml --repfile %s --ASR --ASRc 0.9 -m GY --outrep %s/input_lineages_gy.tsv --run_id gy --threads %d' % (partis_dir, rfn, args.outdir, args.n_procs)]
    # corrects for shm biases:
    cmds += ['%s/packages/igphyml/src/igphyml --repfile %s/input_lineages_gy.tsv --ASR --ASRc 0.9 -m HLP --optimize lr --run_id hlp --threads %d' % (partis_dir, args.outdir, args.n_procs)]
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname=args.outdir + '/run.sh', print_time='igphyml', debug=True, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def fofn(ft):
    assert ft in ['tree', 'seqs']
    return '%s/%s%s' % (args.outdir, ft if ft=='tree' else 'inferred-%s'%ft, '.nwk' if ft=='tree' else '.fa')

# ----------------------------------------------------------------------------------------
def parse_output():
    if utils.all_outputs_exist(args, [fofn('tree'), fofn('seqs')], outlabel='final parsed'):
        return
    utils.makelink(args.outdir, igp_ofn('seqs'), fofn('seqs'), debug=True)
    dtree = treeutils.get_dendro_tree(treefname=igp_ofn('tree'), schema='nexus', debug=True)
    with open(fofn('tree'), 'w') as tfile:
        tfile.write(dtree.as_string(schema='newick'))

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='run:parse')
parser.add_argument('--infname', help='partis yaml format input file')
parser.add_argument('--outdir')
parser.add_argument('--dry-run', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--env-label', default='igphyml')
parser.add_argument('--n-procs', type=int, default=1)

args = parser.parse_args()
args.actions = utils.get_arg_list(args.actions, choices=['run', 'parse', 'install'])
args.infname = utils.fpath(args.infname)
args.outdir = utils.fpath(args.outdir)

if 'install' in args.actions:
    install()
if 'run' in args.actions:
    run_igphyml()
if 'parse' in args.actions:
    parse_output()
