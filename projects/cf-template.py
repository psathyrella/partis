#!/usr/bin/env python2
import argparse
import os
import sys
import colored_traceback.always
import numpy
import math
import collections
import multiprocessing
import yaml

sys.path.insert(1, './python')
import utils
import paircluster
import scanplot
import plotting
import clusterpath

# ----------------------------------------------------------------------------------------
# This is the template/minimal implementation (not necessarily functional) for scan scripts, final versions of which include test/{cf-{germline-inference,paired-loci,tree-metrics}.py and projects/cf-gcdyn.py
# Each of those scripts run most of the steps for a paper, and they (roughly speaking) run a bunch of different methods/metrics and compare the results while scanning over a bunch of different variables.
# No real documentation, although help msgs are there.

# ----------------------------------------------------------------------------------------
script_base = os.path.basename(__file__).replace('cf-', '').replace('.py', '')
all_perf_metrics = ['f1'] #'precision', 'sensitivity', 'f1', 'time-reqd', 'naive-hdist', 'cln-frac']  # pcfrac-*: pair info cleaning correct fraction, cln-frac: collision fraction
after_actions = ['cache-parameters', 'partition']  # actions that come after simulation (e.g. partition)
plot_actions = []  # these are any actions that don't require running any new action, and instead are just plotting stuff that was run by a previous action (e.g. single-chain-partis in cf-paired-loci) (note, does not include 'plot' or 'combine-plots')

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='simu:cache-parameters:partition:plot')
parser.add_argument('--base-outdir', default='%s/partis/%s'%(os.getenv('fs'), script_base))
# parser.add_argument('--n-sim-events-list', default='10', help='N sim events in each repertoire/"proc"/partis simulate run')
parser.add_argument('--dataset-in-list', help='list of samples from --data-in-cfg')
parser.add_argument('--data-in-cfg', help='instead of running simulation to get input for subsequent actions, use input files specified in this yaml cfg (the actual samples to use are specified with --dataset-in-list)')
parser.add_argument('--n-replicates', default=1, type=int)
parser.add_argument('--iseeds', help='if set, only run these replicate indices (i.e. these corresponds to the increment *above* the random seed)')
parser.add_argument('--n-max-procs', type=int, help='Max number of *child* procs (see --n-sub-procs). Default (None) results in no limit.')
parser.add_argument('--n-sub-procs', type=int, default=1, help='Max number of *grandchild* procs (see --n-max-procs)')
parser.add_argument('--random-seed', default=0, type=int, help='note that if --n-replicates is greater than 1, this is only the random seed of the first replicate')
parser.add_argument('--single-light-locus')
# scan fwk stuff (mostly):
parser.add_argument('--version', default='v0')
parser.add_argument('--label', default='test')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--paired-loci', action='store_true', help='same as partis arg')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--test', action='store_true', help='don\'t parallelize \'plot\' action')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--simu-extra-args')
parser.add_argument('--inference-extra-args')
parser.add_argument('--plot-metrics', default='partition', help='these can be either methods (paired-loci) or metrics (tree metrics), but they\'re called metrics in scanplot so that\'s what they have to be called everywhere')
parser.add_argument('--perf-metrics', default=':'.join(all_perf_metrics), help='performance metrics (i.e. usually y axis) that we want to plot vs the scan vars. Only necessary if --plot-metrics are actually methods (as in cf-paired-loci.py).')
parser.add_argument('--zip-vars', help='colon-separated list of variables for which to pair up values sequentially, rather than doing all combinations')
parser.add_argument('--final-plot-xvar', help='variable to put on the x axis of the final comparison plots')
parser.add_argument('--legend-var', help='non-default "component" variable (e.g. obs-frac) to use to label different lines in the legend')
parser.add_argument('--x-legend-var', help='derived variable with which to label the x axis (e.g. mfreq [shm %] when --final-plot-x-var is scratch-mute-freq)')
parser.add_argument('--pvks-to-plot', help='only plot these line/legend values when combining plots')
parser.add_argument('--use-val-cfgs', action='store_true', help='use plotting.val_cfgs dict (we can\'t always use it)')
parser.add_argument('--plot-metric-extra-strs', help='extra strs for each metric in --plot-metrics (i.e. corresponding to what --extra-plotstr was set to during get-tree-metrics for that metric)')
parser.add_argument('--dont-plot-extra-strs', action='store_true', help='while we still use the strings in --plot-metric-extra-strs to find the right dir to get the plot info from, we don\'t actually put the str in the plot (i.e. final plot versions where we don\'t want to see which dtr version it is)')
parser.add_argument('--combo-extra-str', help='extra label for combine-plots action i.e. write to combined-%s/ subdir instead of combined/')
parser.add_argument('--make-hist-plots', action='store_true')
parser.add_argument('--empty-bin-range', help='remove empty bins only outside this range (will kick a warning if this ignores any non-empty bins)')
parser.add_argument('--workdir')  # default set below
args = parser.parse_args()
args.scan_vars = {
    'simu' : ['seed', 'dataset-in'],
}
for act in after_actions + plot_actions:
    if act not in args.scan_vars:
        args.scan_vars[act] = []
    args.scan_vars[act] = args.scan_vars['simu'] + args.scan_vars[act]
args.str_list_vars = []  #  scan vars that are colon-separated lists (e.g. allowed-cdr3-lengths)
args.recurse_replace_vars = []  # scan vars that require weird more complex parsing (e.g. allowed-cdr3-lengths, see cf-paired-loci.py)
args.bool_args = []  # need to keep track of bool args separately (see utils.add_to_scan_cmd())

args.actions = utils.get_arg_list(args.actions, choices=['simu', 'plot', 'combine-plots', ] + after_actions + plot_actions)
args.plot_metrics = utils.get_arg_list(args.plot_metrics)
args.zip_vars = utils.get_arg_list(args.zip_vars)
args.plot_metric_extra_strs = utils.get_arg_list(args.plot_metric_extra_strs)
if args.plot_metric_extra_strs is None:
    args.plot_metric_extra_strs = ['' for _ in args.plot_metrics]
if len(args.plot_metrics) != len(args.plot_metric_extra_strs):
    raise Exception('--plot-metrics %d not same length as --plot-metric-extra-strs %d' % (len(args.plot_metrics), len(args.plot_metric_extra_strs)))
args.pvks_to_plot = utils.get_arg_list(args.pvks_to_plot)
args.perf_metrics = utils.get_arg_list(args.perf_metrics, choices=all_perf_metrics)
args.iseeds = utils.get_arg_list(args.iseeds, intify=True)
args.empty_bin_range = utils.get_arg_list(args.empty_bin_range, floatify=True)

if args.dataset_in_list is not None and args.data_in_cfg is None:
    raise Exception('have to set --data-in-cfg if --dataset-in-list is set')
if args.data_in_cfg is not None:
    assert args.n_replicates == 1
    with open(args.data_in_cfg) as dfile:
        args.data_in_cfg = yaml.load(dfile, Loader=yaml.CLoader)
    if 'simu' in args.actions and not args.dry:
        print '  %s turning on --dry, since the links still get made, and not sure what happens if it isn\'t set' % utils.wrnstr()
        args.dry = True

utils.get_scanvar_arg_lists(args)
if args.final_plot_xvar is None:  # set default value based on scan vars
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars['simu'])
    svars = [v for v in varnames if v != 'seed']
    args.final_plot_xvar = svars[0] if len(svars) > 0 else 'seed'  # if we're not scanning over any vars, i'm not sure what we should use
if len(args.dataset_in_list) < 2:
    raise Exception('at least for now you need to specify at least two samples')

# ----------------------------------------------------------------------------------------
def odir(args, varnames, vstrs, action):
    actstr = action
    if action in ['cache-parameters', 'partition', 'single-chain-partis']:
        actstr = 'partis'
    return '%s/%s' % (utils.svoutdir(args, varnames, vstrs, action), actstr)

# ----------------------------------------------------------------------------------------
def get_fnfcn(method, locus, pmetr):  # have to adjust signature before passing to fcn in scavars
    def tmpfcn(varnames, vstrs): return ofname(args, varnames, vstrs, method, locus=locus, single_file=True)
    return tmpfcn

# ----------------------------------------------------------------------------------------
def ofname(args, varnames, vstrs, action, single_file=False, locus=None):
    outdir = odir(args, varnames, vstrs, action)
    if action == 'cache-parameters':
        if not args.paired_loci or single_file:
            outdir += '/parameters'
        if locus is None and not single_file:
            return outdir
    if args.paired_loci and single_file and locus is None:
        locus = 'igk'
    if args.paired_loci:
        ofn = outdir
        if single_file:
            assert locus is not None
            if action == 'cache-parameters':
                ofn = '%s/%s' % (outdir, locus)
            else:
                ofn = paircluster.paired_fn(outdir, locus, suffix='.yaml', actstr=None if action=='simu' else 'partition')
    else:
        if action == 'cache-parameters':
            ofn = outdir
        else:
            ofn = '%s/%s.yaml' % (outdir, action)
    if action == 'cache-parameters' and single_file:
        # ofn += '/hmm/germline-sets/%s/%sv.fasta' % (locus, locus)
        ofn += '/hmm/all-mean-mute-freqs.csv'
    return ofn

# ----------------------------------------------------------------------------------------
def get_cmd(action, base_args, varnames, vlists, vstrs):
    cmd = './bin/partis %s%s' % (action.replace('simu', 'simulate'), ' --paired-loci' if args.paired_loci else '')
    if args.paired_loci or action != 'cache-parameters':
        cmd += ' --%s %s' % ('paired-outdir' if args.paired_loci else 'outfname',  ofname(args, varnames, vstrs, action))
    if action == 'simu':
        cmd += ' --simulate-from-scratch'
        if args.simu_extra_args is not None:
            cmd += ' %s' % args.simu_extra_args
    else:
        cmd += ' --%s %s' % ('paired-indir' if args.paired_loci else 'infname', ofname(args, varnames, vstrs, 'simu'))
        if args.data_in_cfg is None:
            cmd += ' --is-simu'
        if not args.paired_loci:
            cmd += ' --parameter-dir %s' % ofname(args, varnames, vstrs, 'cache-parameters')
        if action != 'cache-parameters':
            cmd += ' --refuse-to-cache-parameters'
        if args.inference_extra_args is not None:
            cmd += ' %s' % args.inference_extra_args
    return cmd

# ----------------------------------------------------------------------------------------
def plot_loci():
    if args.paired_loci:
        if args.single_light_locus is None:
            return utils.sub_loci('ig')
        else:
            return [utils.heavy_locus('ig'), args.single_light_locus]
    else:
        return [None]

# ----------------------------------------------------------------------------------------
# would be nice to combine this also with fcns in cf-tree-metrics.py (and put it in scanplot)
def run_scan(action):
    base_args, varnames, val_lists, valstrs = utils.get_var_info(args, args.scan_vars[action])
    cmdfos = []
    print '  %s: running %d combinations of: %s' % (utils.color('blue_bkg', action), len(valstrs), ' '.join(varnames))
    if args.debug:
        print '   %s' % ' '.join(varnames)
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print '   %s' % ' '.join(vstrs)

        single_file, locus = (True, None) if args.data_in_cfg is None else (True, args.data_in_cfg[vstrs[0]]['locus'])  # locus may need fixing for data_in_cfg
        ofn = ofname(args, varnames, vstrs, action, single_file=single_file, locus=locus)
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            continue

        if action == 'simu' and args.data_in_cfg is not None:
            if varnames != ['dataset-in']:  # at least for now
                raise Exception('at least for now, if --data-in-cfg/--dataset-in-list are set, \'dataset-in\' must be a scan var (and be the only one), but got %s' % (varnames))
            sample = vstrs[0]
            if sample not in args.data_in_cfg:
                raise Exception('sample \'%s\' not found in meta info with choices: %s' % (sample, args.data_in_cfg.keys()))
            inpath = args.data_in_cfg[sample]['infname']
            if inpath[0] != '/':
                inpath = os.path.realpath(inpath)
            if args.paired_loci and not os.path.isdir(inpath) and os.path.exists(inpath):  # if the data isn't paired (but we're running paired loci), we need to make the paired loci structure UPDATE maybe don't need this any more?
                utils.mkdir(inpath, isfile=True)
                link_name = ofname(args, varnames, vstrs, action, locus=args.data_in_cfg[sample]['locus'])
                assert args.data_in_cfg[sample]['locus'] == 'igh'  # would need to update things if not
                utils.write_empty_annotations(ofname(args, varnames, vstrs, action, locus='igk'), 'igk')  # make a fake light file (ick)
            else:  # NOTE haven't run this, it may need testing
                link_name = '%s%s' % (odir(args, varnames, vstrs, 'simu'), '' if args.paired_loci else '/simu.yaml')
            if not os.path.exists(ofn):
                utils.makelink(os.path.dirname(link_name), inpath, link_name, debug=True)  # , dryrun=True

        cmd = get_cmd(action, base_args, varnames, val_lists, vstrs)
        # utils.simplerun(cmd, logfname='%s-%s.log'%(odir(args, varnames, vstrs, action), action), dryrun=args.dry)
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : odir(args, varnames, vstrs, action),
            'workdir' : '%s/partis-work/%d' % (args.workdir, icombo),
        }]

    utils.run_scan_cmds(args, cmdfos, '%s.log'%action, len(valstrs), n_already_there, ofn)

# ----------------------------------------------------------------------------------------
import random
random.seed(args.random_seed)
numpy.random.seed(args.random_seed)
if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work')))

for action in args.actions:
    if action in ['simu', ] + after_actions:
        run_scan(action)
    elif action in ['plot', 'combine-plots']:
        if args.dry:
            print '  --dry: not plotting'
            continue
        _, varnames, val_lists, valstrs = utils.get_var_info(args, args.scan_vars[after_actions[0]])
        if action == 'plot':
            print 'plotting %d combinations of %d variable%s (%s) to %s' % (len(valstrs), len(varnames), utils.plural(len(varnames)), ', '.join(varnames), scanplot.get_comparison_plotdir(args, None))
            fnames = {meth : {pmetr : [] for pmetr in args.perf_metrics} for meth in args.plot_metrics}
            procs = []
            for method in args.plot_metrics:
                for ltmp in plot_loci():
                    for pmetr in args.perf_metrics:
                        utils.prep_dir(scanplot.get_comparison_plotdir(args, method) + '/' + pmetr, wildlings=['*.html', '*.svg', '*.yaml'])  # , subdirs=args.perf_metrics
                    for pmetr in args.perf_metrics:
                        print '  %12s %s%s' % (method, pmetr, '' if ltmp is None else ' %s'%ltmp)
                        arglist, kwargs = (args, args.scan_vars[after_actions[0]], action, method, pmetr, args.final_plot_xvar), {'fnfcn' : get_fnfcn(method, ltmp, pmetr), 'fnames' : fnames[method][pmetr], 'script_base' : script_base, 'debug' : args.debug}
                        if args.test:
                            scanplot.make_plots(*arglist, **kwargs)
                        else:
                            procs.append(multiprocessing.Process(target=scanplot.make_plots, args=arglist, kwargs=kwargs))
            if not args.test:
                utils.run_proc_functions(procs)
            for method in args.plot_metrics:
                for pmetr in args.perf_metrics:
                    pmcdir = scanplot.get_comparison_plotdir(args, method) + '/' + pmetr
                    fnames[method][pmetr] = [[f.replace(pmcdir+'/', '') for f in fnames[method][pmetr]]]
                    plotting.make_html(pmcdir, n_columns=3, fnames=fnames[method][pmetr])  # this doesn't work unless --test is set since multiprocessing uses copies of <fnames>, but whatever, just run combine-plots
        elif action == 'combine-plots':
            cfpdir = scanplot.get_comparison_plotdir(args, 'combined')
            utils.prep_dir(cfpdir, wildlings=['*.html', '*.svg'])
            fnames = [[] for _ in args.perf_metrics]
            for ipm, pmetr in enumerate(args.perf_metrics):
                print '    ', pmetr
                for ltmp in plot_loci():
                    scanplot.make_plots(args, args.scan_vars[after_actions[0]], action, None, pmetr, args.final_plot_xvar, fnames=fnames[ipm], debug=args.debug)
            plotting.make_html(cfpdir, fnames=fnames)
        else:
            raise Exception('unsupported action %s' % action)
    else:
        raise Exception('unsupported action %s' % action)
