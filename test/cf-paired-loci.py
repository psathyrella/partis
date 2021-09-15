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
# in_param_dir = '_output/paired-simulation/parameters'  # TODO

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='simulate:cache-parameters:partition')  # can also be merge-paired-partitions and get-selection-metrics
parser.add_argument('--base-outdir', default='%s/partis/paired-loci'%os.getenv('fs'))
parser.add_argument('--n-sim-events', type=int, default=10)
parser.add_argument('--n-leaves-list', default='1') #'2:3:4:10') #1 5; do10)
parser.add_argument('--mean-cells-per-droplet-list', default='None') #, default='1') #0.8 2 3 5 10; do #1.1; do
# TODO add fraction-of-reads-to-remove
parser.add_argument('--allowed-cdr3-lengths-list', default='30,45:30,33,36,42,45,48')
parser.add_argument('--mutation-multiplier', type=float, default=1)
parser.add_argument('--n-procs', type=int, default=10)
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--version', default='v0')
parser.add_argument('--label', default='test')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--no-plots', action='store_true')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--extra-args')
parser.add_argument('--zip-vars', help='colon-separated list of variables for which to pair up values sequentially, rather than doing all combinations')
args = parser.parse_args()
args.scan_vars = {'simu' : ['n-leaves', 'mean-cells-per-droplet', 'allowed-cdr3-lengths']}
for act in ['cache-parameters', 'partition']:
    args.scan_vars[act] = args.scan_vars['simu']
args.str_list_vars = ['allowed-cdr3-lengths']
args.svartypes = {'int' : ['n-leaves', 'allowed-cdr3-lengths'], 'float' : []}  # 'mean-cells-per-droplet' # i think can't float() this since we want to allow None as a value

args.actions = utils.get_arg_list(args.actions, choices=['simulate', 'cache-parameters', 'partition', 'merge-paired-partitions', 'get-selection-metrics'])

utils.get_scanvar_arg_lists(args)

# ----------------------------------------------------------------------------------------
def outpath(action):
    if action == 'cache-parameters':
        return 'parameters'
    elif action in ['partition', 'merge-paired-partitions']:  # i guess it makes sense to use the same file for both?
        return 'partition-igh.yaml'  # partition-igh.yaml is about the last file to be written, so it's probably ok to use for this
    elif action == 'get-selection-metrics':
        return 'igh+igk/partition-igh-selection-metrics.yaml'
    else:
        assert False

# ----------------------------------------------------------------------------------------
def run_simu():
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars['simu'])
    # cmdfos = []
    print '  simu: running %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    if args.debug:
        print '   %s' % ' '.join(varnames)
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print '   %s' % ' '.join(vstrs)
        outdir = '%s/simu' % utils.svoutdir(args, varnames, vstrs, 'simu')
        if utils.all_outputs_exist(args, paircluster.paired_dir_fnames(outdir, suffix='.yaml'), debug=False):
            print '    simulation output exists %s' % outdir
            n_already_there += 1
            continue
        cmd = './bin/partis simulate --paired-loci --simulate-from-scratch --random-seed %d --paired-outdir %s %s' % (args.seed, outdir, ' '.join(base_args))  #  --parameter-dir %s in_param_dir
        cmd += ' --n-sim-events %d --n-procs %d --no-per-base-mutation --mutation-multiplier %.2f --constant-number-of-leaves' % (args.n_sim_events, args.n_procs, args.mutation_multiplier)
        if args.extra_args is not None:
            cmd += ' %s' % args.extra_args
        for vname, vstr in zip(varnames, vstrs):
            vstr_for_cmd = vstr
            cmd += ' --%s %s' % (vname, vstr_for_cmd)
        utils.simplerun(cmd, logfname='%s.log'%outdir, dryrun=args.dry)

# ----------------------------------------------------------------------------------------
def run_partis(action):
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars[action])
    # cmdfos = []
    print '  %s: running %d combinations of: %s' % (action, len(valstrs), ' '.join(varnames))
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print '   %s' % ' '.join(vstrs)

        outdir = '%s/inferred' % utils.svoutdir(args, varnames, vstrs, action)
        if utils.output_exists(args, '%s/%s' % (outdir, outpath(action))):
            n_already_there += 1
            continue

        cmd = './bin/partis %s --paired-loci --paired-indir %s/simu --paired-outdir %s' % (action, utils.svoutdir(args, varnames, vstrs, 'simu'), outdir)
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
        utils.simplerun(cmd, logfname='%s-%s.log'%(outdir, action), dryrun=args.dry)

for action in args.actions:
    if action == 'simulate':
        run_simu()
    elif action in ['cache-parameters', 'partition']:
        run_partis(action)

# bd=_output/cells-per-drop
# subd=inferred/plots
# ./bin/compare-plotdirs.py --outdir ~/Dropbox/tmp-plots/cells-per-drop \
#      --normalize --translegend=-0.2:-0.2 \
#      --plotdirs $bd-1.0/$subd:$bd-1.2/$subd:$bd-1.7/$subd:$bd-2.0/$subd:$bd-3.0/$subd \
#      --names 1.0:1.2:1.7:2.0:3.0
