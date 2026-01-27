#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import os
import sys
import colored_traceback.always
import numpy
import math
import collections
import multiprocessing
import yaml

import partis.utils as utils
import partis.paircluster as paircluster
import partis.scanplot as scanplot
import partis.plotting as plotting
import partis.clusterpath as clusterpath

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
# parser.add_argument('--n-sim-events-list', default='10', help='N sim events in each repertoire/"proc"/partis simulate run')
utils.add_scanvar_args(parser, script_base, all_perf_metrics)
args = parser.parse_args()
args.scan_vars = {
    'simu' : ['seed', 'dataset-in'],
}
args.str_list_vars = []  #  scan vars that are colon-separated lists (e.g. allowed-cdr3-lengths)
args.recurse_replace_vars = []  # scan vars that require weird more complex parsing (e.g. allowed-cdr3-lengths, see cf-paired-loci.py)
args.bool_args = []  # need to keep track of bool args separately (see add_to_scan_cmd())
utils.process_scanvar_args(args, after_actions, plot_actions, all_perf_metrics)

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
    print('  %s: running %d combinations of: %s' % (utils.color('blue_bkg', action), len(valstrs), ' '.join(varnames)))
    if args.debug:
        print('   %s' % ' '.join(varnames))
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print('   %s' % ' '.join(vstrs))

        single_file, locus = (True, None) if args.data_in_cfg is None else (True, args.data_in_cfg[vstrs[0]]['locus'])  # locus may need fixing for data_in_cfg
        ofn = ofname(args, varnames, vstrs, action, single_file=single_file, locus=locus)
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            continue

        if action == 'simu' and args.data_in_cfg is not None:
            utils.make_scanvar_data_links(args, varnames, vstrs, action, ofn, ofname, odir)

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
        plt_scn_vars = args.scan_vars['partition']  # NOTE may need to switch this from 'partition'
        if args.dry:
            print('  --dry: not plotting')
            continue
        _, varnames, val_lists, valstrs = utils.get_var_info(args, plt_scn_vars)
        if action == 'plot':
            print('plotting %d combinations of %d variable%s (%s) to %s' % (len(valstrs), len(varnames), utils.plural(len(varnames)), ', '.join(varnames), scanplot.get_comparison_plotdir(args, None)))
            fnames = {meth : {pmetr : [] for pmetr in args.perf_metrics} for meth in args.plot_metrics}
            procs = []
            for method in args.plot_metrics:
                for ltmp in plot_loci():
                    for pmetr in args.perf_metrics:
                        utils.prep_dir(scanplot.get_comparison_plotdir(args, method) + '/' + pmetr, wildlings=['*.html', '*.svg', '*.yaml'])  # , subdirs=args.perf_metrics
                    for pmetr in args.perf_metrics:
                        print('  %12s %s%s' % (method, pmetr, '' if ltmp is None else ' %s'%ltmp))
                        arglist, kwargs = (args, plt_scn_vars, action, method, pmetr, args.final_plot_xvar), {'fnfcn' : get_fnfcn(method, ltmp, pmetr), 'fnames' : fnames[method][pmetr], 'script_base' : script_base, 'debug' : args.debug}
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
                print('    ', pmetr)
                for ltmp in plot_loci():
                    scanplot.make_plots(args, plt_scn_vars, action, None, pmetr, args.final_plot_xvar, fnames=fnames[ipm], debug=args.debug)
            plotting.make_html(cfpdir, fnames=fnames)
        else:
            raise Exception('unsupported action %s' % action)
    else:
        raise Exception('unsupported action %s' % action)
