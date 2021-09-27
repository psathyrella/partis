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
partition_types = ['single', 'joint']
all_perf_metrics = ['precision', 'sensitivity', 'f1']

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='simu:cache-parameters:partition:plot')  # can also be merge-paired-partitions and get-selection-metrics
parser.add_argument('--base-outdir', default='%s/partis/paired-loci'%os.getenv('fs'))
parser.add_argument('--n-sim-events-per-proc', type=int, default=10)
parser.add_argument('--n-leaves-list', default='1') #'2:3:4:10') #1 5; do10)
parser.add_argument('--n-replicates', default=1, type=int)
parser.add_argument('--iseeds', help='if set, only run these replicate indices (i.e. these corresponds to the increment *above* the random seed)')
parser.add_argument('--mean-cells-per-droplet-list') #, default='None')
parser.add_argument('--fraction-of-reads-to-remove-list')
parser.add_argument('--allowed-cdr3-lengths-list') #, default='30,45:30,33,36,42,45,48')
parser.add_argument('--n-genes-per-region-list')
parser.add_argument('--n-sim-alleles-per-gene-list')
parser.add_argument('--scratch-mute-freq-list') #, type=float, default=1)
parser.add_argument('--mutation-multiplier-list') #, type=float, default=1)
parser.add_argument('--n-max-procs', type=int, help='Max number of *child* procs (see --n-sub-procs). Default (None) results in no limit.')
parser.add_argument('--n-sub-procs', type=int, default=1, help='Max number of *grandchild* procs (see --n-max-procs)')
parser.add_argument('--random-seed', default=0, type=int, help='note that if --n-replicates is greater than 1, this is only the random seed of the first replicate')
parser.add_argument('--version', default='v0')
parser.add_argument('--label', default='test')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--no-plots', action='store_true')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--simu-extra-args')
parser.add_argument('--inference-extra-args')
parser.add_argument('--plot-metrics', default='partis', help='NOTE these are methods, but in tree metric script + scanplot they\'re metrics, so we have to call them metrics here')
parser.add_argument('--perf-metrics', default=':'.join(all_perf_metrics))
parser.add_argument('--zip-vars', help='colon-separated list of variables for which to pair up values sequentially, rather than doing all combinations')
parser.add_argument('--final-plot-xvar', help='variable to put on the x axis of the final comparison plots')
parser.add_argument('--pvks-to-plot', help='only plot these line/legend values when combining plots')
parser.add_argument('--plot-metric-extra-strs', help='extra strs for each metric in --plot-metrics (i.e. corresponding to what --extra-plotstr was set to during get-tree-metrics for that metric)')
parser.add_argument('--dont-plot-extra-strs', action='store_true', help='while we still use the strings in --plot-metric-extra-strs to find the right dir to get the plot info from, we don\'t actually put the str in the plot (i.e. final plot versions where we don\'t want to see which dtr version it is)')
parser.add_argument('--combo-extra-str', help='extra label for combine-plots action i.e. write to combined-%s/ subdir instead of combined/')
parser.add_argument('--legend-var')
parser.add_argument('--workdir')  # default set below
args = parser.parse_args()
args.scan_vars = {'simu' : ['seed', 'n-leaves', 'scratch-mute-freq', 'mutation-multiplier', 'mean-cells-per-droplet', 'fraction-of-reads-to-remove', 'allowed-cdr3-lengths', 'n-genes-per-region', 'n-sim-alleles-per-gene']}
for act in ['cache-parameters', 'partition']:
    args.scan_vars[act] = args.scan_vars['simu']
args.str_list_vars = ['allowed-cdr3-lengths', 'n-genes-per-region', 'n-sim-alleles-per-gene']
args.svartypes = {'int' : ['n-leaves', 'allowed-cdr3-lengths', 'n-genes-per-region'], 'float' : ['n-sim-alleles-per-gene']}  # 'mean-cells-per-droplet' # i think can't float() this since we want to allow None as a value

args.actions = utils.get_arg_list(args.actions, choices=['simu', 'cache-parameters', 'partition', 'merge-paired-partitions', 'get-selection-metrics', 'plot', 'combine-plots'])
args.plot_metrics = utils.get_arg_list(args.plot_metrics)
args.plot_metric_extra_strs = utils.get_arg_list(args.plot_metric_extra_strs)
if args.plot_metric_extra_strs is None:
    args.plot_metric_extra_strs = ['' for _ in args.plot_metrics]
if len(args.plot_metrics) != len(args.plot_metric_extra_strs):
    raise Exception('--plot-metrics %d not same length as --plot-metric-extra-strs %d' % (len(args.plot_metrics), len(args.plot_metric_extra_strs)))
args.pvks_to_plot = utils.get_arg_list(args.pvks_to_plot)
args.perf_metrics = utils.get_arg_list(args.perf_metrics, choices=all_perf_metrics)

utils.get_scanvar_arg_lists(args)
if args.final_plot_xvar is None:  # set default value based on scan vars
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars['simu'])
    args.final_plot_xvar = [v for v in varnames if v != 'seed'][0]

# ----------------------------------------------------------------------------------------
def odir(args, varnames, vstrs, action):
    return '%s/%s' % (utils.svoutdir(args, varnames, vstrs, action), action if action=='simu' else 'inferred')

# ----------------------------------------------------------------------------------------
def ofname(args, varnames, vstrs, action, locus=None, single_chain=False, single_file=False):
    outdir = odir(args, varnames, vstrs, action)
    if action == 'cache-parameters':
        return '%s/parameters' % outdir
    assert action in ['simu', 'partition', 'merge-paired-partitions']
    if single_file:
        assert locus is None
        locus = 'igh'
    assert locus is not None
    return paircluster.paired_fn(outdir, locus, suffix='.yaml', actstr=None if action=='simu' else 'partition', single_chain=single_chain)

# ----------------------------------------------------------------------------------------
def run_simu(action):  # TODO this (and run_partis()) should really be combined with their equivalent fcns in cf-tree-metrics.py NOTE actually all four could be combined into one fcn i think?
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars['simu'])
    cmdfos = []
    print '  simu: running %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    if args.debug:
        print '   %s' % ' '.join(varnames)
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print '   %s' % ' '.join(vstrs)
        ofn = ofname(args, varnames, vstrs, action, single_file=True)
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            continue
        cmd = './bin/partis simulate --paired-loci --simulate-from-scratch --paired-outdir %s %s' % (odir(args, varnames, vstrs, action), ' '.join(base_args))  #  --parameter-dir %s in_param_dir
        cmd += ' --no-per-base-mutation'
        if args.n_sub_procs > 1:
            cmd += ' --n-procs %d' % args.n_sub_procs
        if args.n_sim_events_per_proc is not None:
            cmd += ' --n-sim-events %d' % args.n_sim_events_per_proc
        if args.simu_extra_args is not None:
            cmd += ' %s' % args.simu_extra_args
        for vname, vstr in zip(varnames, vstrs):
            vstr_for_cmd = vstr
            cmd += ' --%s %s' % (vname, vstr_for_cmd)
        # utils.simplerun(cmd, logfname='%s.log'%odir(args, varnames, vstrs, action), dryrun=args.dry)
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : odir(args, varnames, vstrs, action),
            'workdir' : '%s/bcr-phylo-work/%d' % (args.workdir, icombo),
        }]
    if n_already_there > 0:
        print '      %d / %d skipped (outputs exist, e.g. %s)' % (n_already_there, len(valstrs), ofn)
    if len(cmdfos) > 0:
        if args.dry:
            print '    %s' % '\n    '.join(cfo['cmd_str'] for cfo in cmdfos)
        else:
            print '      starting %d jobs' % len(cmdfos)
            utils.run_cmds(cmdfos, debug='write:simu.log', n_max_procs=args.n_max_procs, allow_failure=True)

# ----------------------------------------------------------------------------------------
def run_partis(action):
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars[action])
    cmdfos = []
    print '  %s: running %d combinations of: %s' % (action, len(valstrs), ' '.join(varnames))
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print '   %s' % ' '.join(vstrs)

        ofn = ofname(args, varnames, vstrs, action, single_file=True)
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            continue

        cmd = './bin/partis %s --paired-loci --paired-indir %s --paired-outdir %s' % (action, odir(args, varnames, vstrs, 'simu'), odir(args, varnames, vstrs, action))
        if args.n_sub_procs > 1:
            cmd += ' --n-procs %d' % args.n_sub_procs
        if action != 'get-selection-metrics':  # it just breaks here because i don't want to set --simultaneous-true-clonal-seqs (but maybe i should?)
            cmd += ' --is-simu'
        if action != 'cache-parameters' and not args.no_plots:
            cmd += ' --plotdir paired-outdir'
            if action in ['partition', 'merge-paired-partitions']:
                cmd += ' --no-partition-plots' #--no-mds-plots' #
        if args.inference_extra_args is not None:
            cmd += ' %s' % args.inference_extra_args
        # utils.simplerun(cmd, logfname='%s-%s.log'%(odir(args, varnames, vstrs, action), action), dryrun=args.dry)
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : odir(args, varnames, vstrs, action),
            'workdir' : '%s/partis-work/%d' % (args.workdir, icombo),
        }]

    if n_already_there > 0:
        print '      %d / %d skipped (outputs exist, e.g. %s)' % (n_already_there, len(valstrs), ofn)
    if len(cmdfos) > 0:
        print '      %s %d jobs' % ('--dry: would start' if args.dry else 'starting', len(cmdfos))
        if args.dry:
            print '  first command: %s' % cmdfos[0]['cmd_str']
        else:
            utils.run_cmds(cmdfos, debug='write:%s.log'%action, n_max_procs=args.n_max_procs, allow_failure=True)

# ----------------------------------------------------------------------------------------
import random
random.seed(args.random_seed)
numpy.random.seed(args.random_seed)
if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work')))

for action in args.actions:
    if action == 'simu':
        run_simu(action)
    elif action in ['cache-parameters', 'partition']:
        run_partis(action)
    elif action in ['plot', 'combine-plots'] and not args.dry:
        _, varnames, val_lists, valstrs = utils.get_var_info(args, args.scan_vars['partition'])
        if action == 'plot':
            print 'plotting %d combinations of %d variable%s (%s) with %d families per combination to %s' % (len(valstrs), len(varnames), utils.plural(len(varnames)), ', '.join(varnames), 1 if args.n_sim_events_per_proc is None else args.n_sim_events_per_proc, scanplot.get_comparison_plotdir(args, None))
            fnames = {meth : {pmetr : [[] for _ in partition_types] for pmetr in args.perf_metrics} for meth in args.plot_metrics}
            for method in args.plot_metrics:  # NOTE in cf-tree-metrics.py these are [selection] metrics, but here they're [clustering] methods
                cfpdir = scanplot.get_comparison_plotdir(args, method)
                utils.prep_dir(cfpdir, subdirs=all_perf_metrics, wildlings=['*.html', '*.svg', '*.yaml'])
                for ipt, ptntype in enumerate(partition_types):
                    for ltmp in utils.sub_loci('ig'):
                        def fnfcn(varnames, vstrs, tmet, x_axis_label): return ofname(args, varnames, vstrs, 'partition', locus=ltmp, single_chain=ptntype=='single')  # NOTE tmet (e.g. 'precision') and x_axis_label (e.g. 'precision') aren't used for this fcn atm
                        for pmetr in args.perf_metrics:
                            print '  %12s  %6s partition: %3s %s' % (method, ptntype.replace('single', 'single chain'), ltmp, pmetr)
                            scanplot.make_plots(args, args.scan_vars['partition'], action, method, pmetr, plotting.legends.get(pmetr, pmetr), args.final_plot_xvar, fnfcn=fnfcn, locus=ltmp, ptntype=ptntype, fnames=fnames[method][pmetr][ipt], debug=args.debug)
            for method in args.plot_metrics:
                for pmetr in args.perf_metrics:
                    pmcdir = cfpdir + '/' + pmetr
                    fnames[method][pmetr] = [[f.replace(pmcdir, '') for f in flist] for flist in fnames[method][pmetr]]
                    plotting.make_html(pmcdir, n_columns=3, fnames=fnames[method][pmetr])
        elif action == 'combine-plots':
            cfpdir = scanplot.get_comparison_plotdir(args, 'combined')
            utils.prep_dir(cfpdir, wildlings=['*.html', '*.svg'])
            for pmetr in args.perf_metrics:
                print '    ', pmetr
                for ptntype in partition_types:
                    for ltmp in utils.sub_loci('ig'):
                        scanplot.make_plots(args, args.scan_vars['partition'], action, None, pmetr, plotting.legends.get(pmetr, pmetr), args.final_plot_xvar, locus=ltmp, ptntype=ptntype, debug=args.debug)
            plotting.make_html(cfpdir, n_columns=3)
        else:
            assert False

# bd=_output/cells-per-drop
# subd=inferred/plots
# ./bin/compare-plotdirs.py --outdir ~/Dropbox/tmp-plots/cells-per-drop \
#      --normalize --translegend=-0.2:-0.2 \
#      --plotdirs $bd-1.0/$subd:$bd-1.2/$subd:$bd-1.7/$subd:$bd-2.0/$subd:$bd-3.0/$subd \
#      --names 1.0:1.2:1.7:2.0:3.0
