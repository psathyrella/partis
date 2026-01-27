#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
import sys
import argparse
import random
import colored_traceback.always
import os

import partis.utils as utils
import partis.mds as mds

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('infname')
parser.add_argument('--n-clusters', type=int)  # if not set, it just runs mds (i.e. without k-means clustering)
parser.add_argument('--n-components', type=int, default=2)
parser.add_argument('--plotdir')
parser.add_argument('--plotname')
parser.add_argument('--title')
parser.add_argument('--leg-title')
parser.add_argument('--queries-to-include')
parser.add_argument('--workdir', default='/tmp/dralph/mds/' + str(random.randint(0, 999999)))
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--aligned', action='store_true')
args = parser.parse_args()
args.queries_to_include = utils.get_arg_list(args.queries_to_include, key_val_pairs=True)
if args.title is not None:
    args.title = args.title.replace('@', ' ')  # this is kind of hackey

seqfos = utils.read_fastx(args.infname)
color_scale_vals = {}
for sfo in seqfos:
    if len(sfo['infostrs']) == 2:
        color_scale_vals[sfo['name']] = int(sfo['infostrs'][1])
if len(color_scale_vals) == 0:
    color_scale_vals = None

# mds.run_sklearn_mds(args.n_components, args.n_clusters, seqfos, args.seed, plotdir=args.plotdir)
mds.run_bios2mds(args.n_components, args.n_clusters, seqfos, args.workdir, args.seed,
                 aligned=args.aligned, plotdir=args.plotdir, plotname=args.plotname,
                 queries_to_include=args.queries_to_include, color_scale_vals=color_scale_vals, title=args.title, leg_title=args.leg_title)
