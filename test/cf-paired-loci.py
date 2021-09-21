#!/usr/bin/env python
import argparse
import os
import sys
import colored_traceback.always
import numpy
import math
import collections

sys.path.insert(1, './python')
import utils
import paircluster
import scanplot
import plotting

# ----------------------------------------------------------------------------------------
# in_param_dir = '_output/paired-simulation/parameters'  # TODO

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='simulate:cache-parameters:partition:plot')  # can also be merge-paired-partitions and get-selection-metrics
parser.add_argument('--base-outdir', default='%s/partis/paired-loci'%os.getenv('fs'))
parser.add_argument('--n-sim-events-per-proc', type=int, default=10)
parser.add_argument('--n-leaves-list', default='1') #'2:3:4:10') #1 5; do10)
parser.add_argument('--n-replicates', default=1, type=int)
parser.add_argument('--iseeds', help='if set, only run these replicate indices (i.e. these corresponds to the increment *above* the random seed)')
parser.add_argument('--mean-cells-per-droplet-list', default='None') #, default='1') #0.8 2 3 5 10; do #1.1; do
# TODO add fraction-of-reads-to-remove
parser.add_argument('--allowed-cdr3-lengths-list', default='30,45:30,33,36,42,45,48')
parser.add_argument('--mutation-multiplier', type=float, default=1)
parser.add_argument('--n-procs', type=int, default=10)
parser.add_argument('--random-seed', default=0, type=int, help='note that if --n-replicates is greater than 1, this is only the random seed of the first replicate')
parser.add_argument('--version', default='v0')
parser.add_argument('--label', default='test')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--no-plots', action='store_true')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--extra-args')
parser.add_argument('--plot-metrics', default='partis', help='NOTE these are methods, but in tree metric script + scanplot they\'re metrics, so we have to call them metrics here')
parser.add_argument('--perf-metrics', default='precision:sensitivity:f1')
parser.add_argument('--zip-vars', help='colon-separated list of variables for which to pair up values sequentially, rather than doing all combinations')
parser.add_argument('--final-plot-xvar', default='mean-cells-per-droplet', help='variable to put on the x axis of the final comparison plots')
parser.add_argument('--pvks-to-plot', help='only plot these line/legend values when combining plots')
parser.add_argument('--plot-metric-extra-strs', help='extra strs for each metric in --plot-metrics (i.e. corresponding to what --extra-plotstr was set to during get-tree-metrics for that metric)')
parser.add_argument('--dont-plot-extra-strs', action='store_true', help='while we still use the strings in --plot-metric-extra-strs to find the right dir to get the plot info from, we don\'t actually put the str in the plot (i.e. final plot versions where we don\'t want to see which dtr version it is)')
parser.add_argument('--combo-extra-str', help='extra label for combine-plots action i.e. write to combined-%s/ subdir instead of combined/')
parser.add_argument('--legend-var')
args = parser.parse_args()
args.scan_vars = {'simu' : ['seed', 'n-leaves', 'mean-cells-per-droplet', 'allowed-cdr3-lengths']}
for act in ['cache-parameters', 'partition']:
    args.scan_vars[act] = args.scan_vars['simu']
args.str_list_vars = ['allowed-cdr3-lengths']
args.svartypes = {'int' : ['n-leaves', 'allowed-cdr3-lengths'], 'float' : []}  # 'mean-cells-per-droplet' # i think can't float() this since we want to allow None as a value

args.actions = utils.get_arg_list(args.actions, choices=['simulate', 'cache-parameters', 'partition', 'merge-paired-partitions', 'get-selection-metrics', 'plot', 'combine-plots'])
args.plot_metrics = utils.get_arg_list(args.plot_metrics)
args.plot_metric_extra_strs = utils.get_arg_list(args.plot_metric_extra_strs)
if args.plot_metric_extra_strs is None:
    args.plot_metric_extra_strs = ['' for _ in args.plot_metrics]
if len(args.plot_metrics) != len(args.plot_metric_extra_strs):
    raise Exception('--plot-metrics %d not same length as --plot-metric-extra-strs %d' % (len(args.plot_metrics), len(args.plot_metric_extra_strs)))
args.pvks_to_plot = utils.get_arg_list(args.pvks_to_plot)
args.perf_metrics = utils.get_arg_list(args.perf_metrics)

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
# TODO
        outdir = '%s/simu' % utils.svoutdir(args, varnames, vstrs, 'simu')
        if utils.all_outputs_exist(args, paircluster.paired_dir_fnames(outdir, suffix='.yaml'), debug=False):
            print '    simulation output exists %s' % outdir
            n_already_there += 1
            continue
        cmd = './bin/partis simulate --paired-loci --simulate-from-scratch --paired-outdir %s %s' % (outdir, ' '.join(base_args))  #  --parameter-dir %s in_param_dir
        cmd += ' --n-procs %d --no-per-base-mutation --mutation-multiplier %.2f --constant-number-of-leaves' % (args.n_procs, args.mutation_multiplier)
        if args.n_sim_events_per_proc is not None:
            cmd += ' --n-sim-events %d' % args.n_sim_events_per_proc
        if args.extra_args is not None:
            cmd += ' %s' % args.extra_args
        for vname, vstr in zip(varnames, vstrs):
            vstr_for_cmd = vstr
            cmd += ' --%s %s' % (vname, vstr_for_cmd)
        utils.simplerun(cmd, logfname='%s.log'%outdir, dryrun=args.dry)

# ----------------------------------------------------------------------------------------
def odir(args, varnames, vstrs, action):
    return '%s/%s' % (utils.svoutdir(args, varnames, vstrs, action), action if action=='simu' else 'inferred')

# ----------------------------------------------------------------------------------------
def ofname(args, varnames, vstrs, action):
    return '%s/%s' % (odir(args, varnames, vstrs, action), outpath(action))

# ----------------------------------------------------------------------------------------
def run_partis(action):
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars[action])
    # cmdfos = []
    print '  %s: running %d combinations of: %s' % (action, len(valstrs), ' '.join(varnames))
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print '   %s' % ' '.join(vstrs)

        if utils.output_exists(args, ofname(args, varnames, vstrs, action)):
            n_already_there += 1
            continue

        cmd = './bin/partis %s --paired-loci --paired-indir %s --paired-outdir %s' % (action, odir(args, varnames, vstrs, 'simu'), odir(args, varnames, vstrs, action))
        cmd += ' --n-procs %d' % args.n_procs
        if action != 'get-selection-metrics':  # it just breaks here because i don't want to set --simultaneous-true-clonal-seqs (but maybe i should?)
            cmd += ' --is-simu'
        if action != 'cache-parameters' and not args.no_plots:
            cmd += ' --plotdir paired-outdir'
            if action in ['partition', 'merge-paired-partitions']:
                cmd += ' --no-partition-plots' #--no-mds-plots' #
        if args.extra_args is not None:
            cmd += ' %s' % args.extra_args
        utils.simplerun(cmd, logfname='%s-%s.log'%(odir(args, varnames, vstrs, action), action), dryrun=args.dry)

# ----------------------------------------------------------------------------------------
import random
random.seed(args.random_seed)
numpy.random.seed(args.random_seed)

for action in args.actions:
    if action == 'simulate':
        run_simu()
    elif action in ['cache-parameters', 'partition']:
        run_partis(action)
    elif action in ['plot', 'combine-plots'] and not args.dry:
        _, varnames, val_lists, valstrs = utils.get_var_info(args, args.scan_vars['partition'])
        def fnfcn(varnames, vstrs, tmet, x_axis_label): return ofname(args, varnames, vstrs, 'partition')
        if action == 'plot':
            print 'plotting %d combinations of %d variable%s (%s) with %d families per combination to %s' % (len(valstrs), len(varnames), utils.plural(len(varnames)), ', '.join(varnames), 1 if args.n_sim_events_per_proc is None else args.n_sim_events_per_proc, scanplot.get_comparison_plotdir(args, None, None))
            for method in args.plot_metrics:  # NOTE in cf-tree-metrics.py these are metrics, but here they're more like different methods
                utils.prep_dir(scanplot.get_comparison_plotdir(args, method, None), subdirs=args.perf_metrics, wildlings=['*.html', '*.svg', '*.yaml'])
                for pmetr in args.perf_metrics:
                    print '    ', pmetr
                    scanplot.make_plots(args, args.scan_vars['partition'], action, method, pmetr, plotting.legends.get(pmetr, pmetr), args.final_plot_xvar, fnfcn, debug=args.debug)
            for method in args.plot_metrics:
                plotting.make_html(scanplot.get_comparison_plotdir(args, method, None), n_columns=3)
        elif action == 'combine-plots':
            utils.prep_dir(scanplot.get_comparison_plotdir(args, 'combined', None), wildlings=['*.html', '*.svg'])
            for pmetr in args.perf_metrics:
                print '    ', pmetr
                scanplot.make_plots(args, args.scan_vars['partition'], action, None, pmetr, plotting.legends.get(pmetr, pmetr), args.final_plot_xvar, fnfcn, debug=args.debug)
            plotting.make_html(scanplot.get_comparison_plotdir(args, 'combined', None), n_columns=2)
        else:
            assert False

# bd=_output/cells-per-drop
# subd=inferred/plots
# ./bin/compare-plotdirs.py --outdir ~/Dropbox/tmp-plots/cells-per-drop \
#      --normalize --translegend=-0.2:-0.2 \
#      --plotdirs $bd-1.0/$subd:$bd-1.2/$subd:$bd-1.7/$subd:$bd-2.0/$subd:$bd-3.0/$subd \
#      --names 1.0:1.2:1.7:2.0:3.0
