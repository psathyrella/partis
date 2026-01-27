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

import partis.utils as utils
import partis.glutils as glutils
import partis.treeutils as treeutils
from partis.clusterpath import ClusterPath

# ----------------------------------------------------------------------------------------
def install():
    cmds = ['#!/bin/bash']
    cmds += utils.mamba_cmds(args.env_label, only_prep=True)
    cmds += ['micromamba create -n %s python=3.9' % args.env_label]  # 3.10 currently has problems with ete
    cmds += ['micromamba activate %s' % args.env_label]
    # cmds += ['micromamba install -c bioconda phylip']
    cmds += ['pip install hilary']
    # micromamba remove -n gctree --all  # to nuke it and start over
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname='/tmp/tmprun.sh', debug=True)

# ----------------------------------------------------------------------------------------
def update():
    cmds = ['#!/bin/bash']
    cmds += utils.mamba_cmds(args.env_label)
    cmds += ['pip install --upgrade hilary']
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname='/tmp/tmprun.sh', debug=True)

# ----------------------------------------------------------------------------------------
def hofn():  # same as the current default
    return '%s/hilary_results/inferred_full_method_input.tsv' % args.outdir  # "input" has to match afn below (hilary doesn't seem to havey any options to control output file)

# ----------------------------------------------------------------------------------------
def fofn():
    return '%s/partition.yaml' % args.outdir

# ----------------------------------------------------------------------------------------
def run_hilary():
    if utils.output_exists(args, hofn(), outlabel='hilary'):
        return

    afn = '%s/input.tsv' % args.outdir
    if utils.output_exists(args, afn):
        print('    airr tsv exists, not rerunning: %s' % afn)
    else:
        cmd = '%s/bin/parse-output.py %s %s --airr-output' % (partis_dir, args.infname, afn)
        utils.simplerun(cmd)

    cmds = ['#!/bin/bash']
    cmds += utils.mamba_cmds(args.env_label)
    cmds += ['infer-lineages full-method %s --silent --result-folder %s' % (afn, os.path.dirname(hofn()))]  # --silent is supposed to kill the progress bar, but doesn't actually seem to)  --verbose
    # --kappa-file   # path to igk file, which tells it to do paired clustering
    if args.overwrite:
        cmds[-1] += ' --override'  # don't really need to set this, but it's probably nice
    utils.simplerun('\n'.join(cmds) + '\n', cmdfname=args.outdir + '/run.sh', print_time='hilary', debug=True, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def parse_output():
    if utils.output_exists(args, fofn(), outlabel='final parsed'):
        return
    glfo, _, _ = utils.read_output(args.infname, skip_annotations=True)
    glfo, annotation_list, cpath = utils.read_airr_output(hofn(), locus=args.locus, glfo=glfo, clone_id_field='family')  # not sure why they don't use clone_id
    utils.write_annotations(fofn(), glfo, annotation_list, utils.annotation_headers, partition_lines=cpath.get_partition_lines())

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='run:parse')
parser.add_argument('--infname', help='partis yaml format input file')
parser.add_argument('--outdir')
parser.add_argument('--locus', default='igh')
parser.add_argument('--dry-run', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--env-label', default='hilary')

args = parser.parse_args()
args.actions = utils.get_arg_list(args.actions, choices=['run', 'parse', 'install', 'update'])
args.infname = utils.fpath(args.infname)
args.outdir = utils.fpath(args.outdir)

if 'install' in args.actions:
    install()
if 'update' in args.actions:
    update()
if 'run' in args.actions:
    run_hilary()
if 'parse' in args.actions:
    parse_output()
