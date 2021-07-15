#!/usr/bin/env python
import argparse
import os
import sys
import colored_traceback.always
import numpy
import math

sys.path.insert(1, './python')
import utils
import paircluster

# ----------------------------------------------------------------------------------------
in_param_dir = '_output/paired-simulation/parameters'  # TODO

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='simulate:cache-parameters:partition')
parser.add_argument('--n-sim-events', type=int, default=10)
parser.add_argument('--n-leaf-list', default='1') #'2:3:4:10') #1 5; do10)
parser.add_argument('--cells-per-drop-list', default='1') #0.8 2 3 5 10; do #1.1; do
parser.add_argument('--allowed-cdr3-lengths', default='30:33:36:42:45:48')
parser.add_argument('--mutation-multiplier', type=float, default=1)
parser.add_argument('--n-procs', type=int, default=10)
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--version', default='v0')
parser.add_argument('--label', default='test')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--no-plots', action='store_true')
parser.add_argument('--extra-args')
args = parser.parse_args()

args.actions = utils.get_arg_list(args.actions, choices=['simulate', 'cache-parameters', 'partition', 'merge-paired-partitions', 'get-selection-metrics'])
args.n_leaf_list = utils.get_arg_list(args.n_leaf_list, intify=True)
args.cells_per_drop_list = utils.get_arg_list(args.cells_per_drop_list, floatify=True)

# ----------------------------------------------------------------------------------------
def run_simu(ncells, nleaf):
    if utils.all_outputs_exist(args, paircluster.paired_dir_fnames('%s/simu'%outdir, suffix='.yaml'), debug=False):
        print '    simulation output exists %s' % ('%s/simu'%outdir)
        return
    cmd = './bin/partis simulate --paired-loci --simulate-from-scratch --random-seed %d --paired-outdir %s/simu --mean-cells-per-droplet %f' % (args.seed, outdir, ncells)  #  --parameter-dir %s in_param_dir
    cmd += ' --n-sim-events %d --n-leaves %d --n-procs %d --no-per-base-mutation --allowed-cdr3-lengths %s' % (args.n_sim_events, nleaf, args.n_procs, args.allowed_cdr3_lengths)
    cmd += ' --mutation-multiplier %.2f --constant-number-of-leaves' % args.mutation_multiplier
    if args.extra_args is not None:
        cmd += ' %s' % args.extra_args
    utils.simplerun(cmd, logfname='%s/simu.log'%outdir, dryrun=args.dry)

# ----------------------------------------------------------------------------------------
def outpath(action):
    if action == 'cache-parameters':
        outstr = 'parameters'
    elif action in ['partition', 'merge-paired-partitions']:  # i guess it makes sense to use the same file for both?
        outstr = 'partition-igh.yaml'  # partition-igh.yaml is about the last file to be written, so it's probably ok to use for this
    elif action == 'get-selection-metrics':
        outstr = 'igh+igk/partition-igh-selection-metrics.yaml'
    else:
        assert False
    return '%s/inferred/%s' % (outdir, outstr)

# ----------------------------------------------------------------------------------------
for ncells in args.cells_per_drop_list:
    for nleaf in args.n_leaf_list:
        ncstr = ('%.0f'%ncells) if int(ncells)==ncells else '%.1f'%ncells
        outlabel = '%s%s%s-%s' % (args.label, '-'+str(nleaf) if len(args.n_leaf_list)>1 else '', '-'+ncstr if len(args.cells_per_drop_list)>1 else '', args.version)
        outdir = '%s/partis/paired-loci/%s' % (os.getenv('fs'), outlabel)
        out_param_dir = '%s/params' % outdir

        for action in args.actions:
            if action == 'simulate':
                run_simu(ncells, nleaf)
                continue
            if utils.output_exists(args, outpath(action)):
                continue
            cmd = './bin/partis %s --paired-loci --paired-indir %s/simu --input-metafname %s/simu/meta.yaml --paired-outdir %s/inferred' % (action, outdir, outdir, outdir)
            cmd += ' --random-seed %d' % args.seed
            cmd += ' --n-procs %d' % args.n_procs
            if action != 'get-selection-metrics':  # it just breaks here because i don't want to set --simultaneous-true-clonal-seqs (but maybe i should?)
                cmd += ' --is-simu'
            if action != 'cache-parameters' and not args.no_plots:
        	cmd += ' --plotdir paired-outdir'
                if action in ['partition', 'merge-paired-partitions']:
                    cmd += ' --no-partition-plots' #--no-mds-plots' #
            if args.extra_args is not None:
                cmd += ' %s' % args.extra_args
            utils.simplerun(cmd, logfname='%s/%s.log'%(outdir, action), dryrun=args.dry)
        break

# bd=_output/cells-per-drop
# subd=inferred/plots
# ./bin/compare-plotdirs.py --outdir ~/Dropbox/tmp-plots/cells-per-drop \
#      --normalize --translegend=-0.2:-0.2 \
#      --plotdirs $bd-1.0/$subd:$bd-1.2/$subd:$bd-1.7/$subd:$bd-2.0/$subd:$bd-3.0/$subd \
#      --names 1.0:1.2:1.7:2.0:3.0
