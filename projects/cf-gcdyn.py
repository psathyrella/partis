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
import glob
import copy

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
dl_metrics = ['%s-%s-%s' % (p, s, m) for s in ['train', 'test'] for m in ['bias', 'variance', 'mae'] for p in ['xscale', 'xshift']] + ['xscale-train-vs-test-mae']
all_perf_metrics = ['max-abundances', 'distr-abundances', 'distr-hdists', 'all-dl', 'all-test-dl'] + dl_metrics #'precision', 'sensitivity', 'f1', 'time-reqd', 'naive-hdist', 'cln-frac']  # pcfrac-*: pair info cleaning correct fraction, cln-frac: collision fraction
# 'check-dl': make simulation using the dl-inferred parameters (so you can then run 'replay-plot-ckdl')
after_actions = ['merge-simu', 'replay-plot', 'dl-infer', 'data', 'check-dl', 'replay-plot-ckdl', 'dl-infer-merged', 'write-partis-simu-files', 'partis']  # actions that come after simulation (e.g. partition)
plot_actions = ['group-expts']  # these are any actions that don't require running any new action, and instead are just plotting stuff that was run by a previous action (e.g. single-chain-partis in cf-paired-loci) (note, does not include 'plot' or 'combine-plots')
merge_actions = ['merge-simu', 'dl-infer-merged']  # actions that act on all scanned values at once (i.e. only run once, regardless of how many scan vars/values)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='don\'t use')  # old defaults: simu:merge-simu:process:plot (partis is just to make tree plots, which aren't really slow but otoh we don't need for every seed/variable combo)
parser.add_argument('--birth-response-list')
parser.add_argument('--death-response-list')
parser.add_argument('--death-values-list')
parser.add_argument('--death-range-list')
parser.add_argument('--xscale-values-list')
parser.add_argument('--xshift-values-list')
parser.add_argument('--yscale-values-list')
parser.add_argument('--yshift-values-list')
parser.add_argument('--xscale-range-list')
parser.add_argument('--xshift-range-list')
parser.add_argument('--yscale-range-list')
parser.add_argument('--yshift-range-list')
parser.add_argument('--x-ceil-start-range-list')
parser.add_argument('--x-ceil-start-values-list')
parser.add_argument('--initial-birth-rate-range-list')
parser.add_argument('--carry-cap-values-list')
parser.add_argument('--carry-cap-range-list')
parser.add_argument('--capacity-method-list')
parser.add_argument('--nonsense-phenotype-value-list')
parser.add_argument('--init-population-values-list')
parser.add_argument('--time-to-sampling-range-list')
parser.add_argument('--n-seqs-range-list')
parser.add_argument('--n-trials-list')
parser.add_argument('--dl-bundle-size-list', help='size of bundles during dl inference (must be equal to or less than simulation bundle size)')
parser.add_argument('--model-type-list')
parser.add_argument('--epochs-list')
parser.add_argument('--batch-size-list')
parser.add_argument('--dropout-rate-list')
parser.add_argument('--learning-rate-list')
parser.add_argument('--ema-momentum-list')
parser.add_argument('--prebundle-layer-cfg-list')
parser.add_argument('--dont-scale-params-list')
parser.add_argument('--params-to-predict-list')
parser.add_argument('--n-trees-per-expt-list', help='Number of per-tree predictions to group together and average over during the \'group-expts\' action (see also --dl-bundle-size-list)')
parser.add_argument('--simu-bundle-size-list', help='Number of trees to simulate with each chosen set of parameter values, in each simulation subprocess (see also --n-trees-per-expt-list')
parser.add_argument('--data-samples-list', help='List of data samples to run on. Don\'t need to set this, it can use glob on --data-dir')
utils.add_scanvar_args(parser, script_base, all_perf_metrics, default_plot_metric='replay-plot')
parser.add_argument('--check-dl-n-trials', type=int, default=120)
parser.add_argument('--dl-extra-args')
parser.add_argument('--gcddir', default='%s/work/partis/projects/gcdyn'%os.getenv('HOME'))
# parser.add_argument('--gcreplay-data-dir', default='/fh/fast/matsen_e/%s/gcdyn/gcreplay-observed'%os.getenv('USER'))
parser.add_argument('--gcreplay-germline-dir', default='datascripts/meta/taraki-gctree-2021-10/germlines')  # only used for partis partition step, which I'm not using atm
parser.add_argument('--dl-model-dir')
parser.add_argument('--dl-model-label-str')
parser.add_argument('--dl-model-vars')
parser.add_argument('--tree-inference-method', choices=['iqtree', 'gctree'], help='if set, we both run tree inference with this method on simulation, and use those inferred trees during training')
parser.add_argument('--data-dir')
parser.add_argument('--check-dl-input', action='store_true', help='for actions that would normally use input from \'simu\' action (e.g. dl-infer), instead use the simulation from \'check-dl\' action')
parser.add_argument('--iqtree-data-dir')
parser.add_argument('--iqtree-version')
parser.add_argument('--gcreplay-dir')
args = parser.parse_args()
args.scan_vars = {
    'simu' : ['seed', 'birth-response', 'death-response', 'death-values', 'death-range', 'xscale-values', 'xshift-values', 'yscale-values', 'yshift-values', 'xscale-range', 'xshift-range', 'yscale-range', 'yshift-range', 'x-ceil-start-range', 'x-ceil-start-values', 'initial-birth-rate-range', 'carry-cap-values', 'carry-cap-range', 'capacity-method', 'nonsense-phenotype-value', 'init-population-values', 'time-to-sampling-range', 'n-seqs-range', 'n-trials', 'simu-bundle-size'],
    'dl-infer' : ['dl-bundle-size', 'model-type', 'epochs', 'batch-size', 'dropout-rate', 'learning-rate', 'ema-momentum', 'prebundle-layer-cfg', 'dont-scale-params', 'params-to-predict'],
    'data' : ['data-samples', 'model-type', 'carry-cap-values', 'init-population-values', 'death-values'],
}
args.scan_vars['group-expts'] = copy.deepcopy(args.scan_vars['dl-infer'])
args.scan_vars['check-dl'] = copy.deepcopy(args.scan_vars['data']) + [v for v in args.scan_vars['dl-infer'] if v not in args.scan_vars['data']]
check_dl_args = [v for v in args.scan_vars['simu'] if not any(p in v for p in ['xscale', 'xshift', 'yscale', 'x-ceil'])]  # UGH
args.scan_vars['replay-plot-ckdl'] = copy.deepcopy(args.scan_vars['check-dl'])
args.str_list_vars = ['death-values', 'death-range', 'xscale-values', 'xshift-values', 'yscale-values', 'yshift-values', 'xscale-range', 'xshift-range', 'yscale-range', 'yshift-range', 'x-ceil-start-range', 'x-ceil-start-values', 'initial-birth-rate-range', 'time-to-sampling-range', 'carry-cap-values', 'carry-cap-range', 'capacity-method', 'init-population-values', 'n-seqs-range', 'params-to-predict']  #  scan vars that are colon-separated lists (e.g. allowed-cdr3-lengths)
args.recurse_replace_vars = []  # scan vars that require weird more complex parsing (e.g. allowed-cdr3-lengths, see cf-paired-loci.py)
args.bool_args = ['dont-scale-params']  # need to keep track of bool args separately (see utils.add_to_scan_cmd())
if 'data' in args.actions:
    assert args.data_dir is not None
    if args.data_samples_list is None:
        args.data_samples_list = ':'.join(os.path.basename(d) for d in glob.glob('%s/*' % args.data_dir) if os.path.isdir(d))
    print('  running on %d data sample%s from %s' % (len(args.data_samples_list.split(':')), utils.plural(len(args.data_samples_list.split(':'))), args.data_dir))
data_actions = ['data', 'check-dl', 'replay-plot-ckdl']
if args.check_dl_input:  # this is ugly, but if we're running dl-infer on check-dl input (after data inference), we need to do this to make sure the order matches the other actions
    data_actions.append('dl-infer')
    args.scan_vars['dl-infer'] = copy.deepcopy(args.scan_vars['check-dl']) + [v for v in args.scan_vars['dl-infer'] if v not in args.scan_vars['check-dl']]
utils.process_scanvar_args(args, after_actions, plot_actions, all_perf_metrics, data_actions=data_actions)
args.dl_model_vars = utils.get_arg_list(args.dl_model_vars)
if args.inference_extra_args is not None:
    raise Exception('not used atm')
if 'all-dl' in args.perf_metrics:
    args.perf_metrics += dl_metrics
    args.perf_metrics.remove('all-dl')
if 'all-test-dl' in args.perf_metrics:
    args.perf_metrics += [m for m in dl_metrics if 'test' in m]
    args.perf_metrics.remove('all-test-dl')
if 'group-expts' in args.plot_metrics and any('train' in m for m in args.perf_metrics):
    acts_to_remove = [m for m in args.perf_metrics if 'train' in m]
    print('    can\'t plot any training metrics for \'group-expts\' so removing: %s' % ' '.join(acts_to_remove))
    args.perf_metrics = [m for m in args.perf_metrics if 'train' not in m]
if args.params_to_predict_list is not None:
    before_metrics = copy.copy(args.perf_metrics)
    args.perf_metrics = [m for m in args.perf_metrics if any(p in m for plist in args.params_to_predict_list for p in plist)]
    print('    restricting perf metrics to correspond to --params-to-predict (removed %s)' % (' '.join(m for m in before_metrics if m not in args.perf_metrics)))

# ----------------------------------------------------------------------------------------
def odir(args, varnames, vstrs, action):
    return utils.svoutdir(args, varnames, vstrs, action)

# ----------------------------------------------------------------------------------------
def ofname(args, varnames, vstrs, action, ftype='npy', current_action=None):  # <action> is the action that create[ds] the file we're looking for, whereas <current_action> is the action we're now running
    if action == 'merge-simu':
        varnames = []
        vstrs = []
    if action in ['simu', 'merge-simu', 'check-dl']:
        ftstrs = {  # copied from gcdyn/scripts/multi_simulation.py
            'fasta' : 'seqs',
            'nwk' : 'trees',
            'npy' : 'encoded-trees',
            'pkl' : 'responses',
        }
        sfx = '%s.%s' % (ftstrs.get(ftype, 'simu'), ftype)
    elif 'replay-plot' in action:
        sfx = 'diff-vals.yaml'
    elif action in ['dl-infer', 'dl-infer-merged', 'group-expts']:
        if args.dl_model_dir is None and args.dl_model_label_str is None:
            sfx = 'test.csv'
        else:
            sfx = 'infer.csv'
    elif action == 'write-partis-simu-files':
        sfx = 'igh.yaml'
    elif action == 'partis':
        sfx = 'partition.yaml'
    elif action == 'data':
        sfx = 'infer.csv'
    else:
        assert False
    if args.tree_inference_method is not None and 'dl-infer' in [action, current_action]:  # NOTE not sure which actions really want this added, but don't want to test all of them atm
        sfx = '%s/%s' % (args.tree_inference_method, sfx)
    return '%s/%s/%s' % (odir(args, varnames, vstrs, action), action, sfx)

# ----------------------------------------------------------------------------------------
def add_mamba_cmds(cmd):
    if os.getenv('CONDA_PREFIX') is None:  # if we're already in a conda/mamba env, shouldn't need to add the commands
        cmd = utils.mamba_cmds('gcdyn') + [cmd]
        cmd = ' && '.join(cmd)
    return cmd

# ----------------------------------------------------------------------------------------
def get_cmd(action, base_args, varnames, vstrs, all_simdirs=None):
    # ----------------------------------------------------------------------------------------
    def add_scan_args(cmd, skip_fcn=None):  # using nargs='+' syntax for these rather than partis-style colons
        for barg in base_args:  # NOTE don't modify base_args here
            bname, bstr = barg.split(' ') if len(barg.split(' '))==2 else (barg, '')  # first clause: arg with value, second: bool-type arg without value
            bname = bname.replace('--', '')
            if skip_fcn is not None and skip_fcn(bname):
                continue
            if bname in args.str_list_vars:
                bstr = bstr.replace(':', ' ')
            cmd += ' --%s %s' % (bname, bstr)
        for vname, vstr in zip(varnames, vstrs):
            if skip_fcn is not None and skip_fcn(vname):
                continue
            if len(varnames)==1 and vstr is None:  # used force_val when getting var info
                continue
            if vname in args.str_list_vars:
                vstr = vstr.replace(':', ' ')
            cmd = utils.add_to_scan_cmd(args, vname, vstr, cmd)
        return cmd
    # ----------------------------------------------------------------------------------------
    def add_model_dir(cmd):
        if args.dl_model_dir is not None:
            cmd += ' --model-dir %s' % args.dl_model_dir
        if args.dl_model_label_str is not None:
            tvnames, tvstrs = varnames, vstrs
            if args.dl_model_vars is not None:
                tvnames, tvstrs = zip(*[(v, s) for v, s in zip(varnames, vstrs) if v in args.dl_model_vars])
            tmpdir = os.path.dirname(ofname(args, tvnames, tvstrs, 'dl-infer'))  # location at which we'd write dl training results if we were training here
            assert '%s/%s'%(args.label, args.version) in tmpdir
            cmd += ' --model-dir %s' % tmpdir.replace('%s/%s'%(args.label, args.version), args.dl_model_label_str)
        return cmd
    # ----------------------------------------------------------------------------------------
    odr = os.path.dirname(ofname(args, varnames, vstrs, action))  # NOTE that this isn't the same as odir() (ugh)
    if action in ['simu', 'check-dl', 'merge-simu']:
        cmd = 'gcd-simulate' if action in ['simu', 'check-dl'] else 'python %s/scripts/%s.py' % (args.gcddir, 'combine-simu-files.py')
        if action in ['simu', 'check-dl']:
            cmd += ' --outdir %s --debug 1' % odr  #  --debug 1
            # --make-plots
            if args.tree_inference_method is not None:
                cmd += ' --tree-inference-method %s' % args.tree_inference_method
            if args.test:
                cmd += ' --test'
            cmd = add_scan_args(cmd, skip_fcn=lambda v: v not in args.scan_vars[action] or action=='check-dl' and v not in check_dl_args)
        else:
            cmd += ' %s --outdir %s' % (' '.join(all_simdirs), os.path.dirname(ofname(args, [], [], action)))
        if args.n_sub_procs > 1:
            cmd += ' --n-sub-procs %d' % args.n_sub_procs
        if action == 'check-dl':
            cmd += ' --dl-prediction-file %s/plots/median-curve-%s.csv' % (os.path.dirname(ofname(args, varnames, vstrs, 'dl-infer' if args.data_samples_list is None else 'data')), 'test' if args.data_samples_list is None else 'infer')
            trials_per_proc = 10 #4 # 10
            if args.n_sub_procs > 1 and '--n-trials' in cmd:
                trials_per_proc = int(int(utils.get_val_from_arglist(cmd.split(), '--n-trials')) / args.n_sub_procs)  # keep the same trials/proc as for regular simulation action
            cmd = utils.replace_in_argstr(cmd, '--n-sub-procs', str(int(args.check_dl_n_trials / trials_per_proc)))
            cmd = utils.replace_in_argstr(cmd, '--n-trials', str(args.check_dl_n_trials))
        if action in ['simu', 'check-dl'] and args.simu_extra_args is not None:
            cmd += ' %s' % args.simu_extra_args
        cmd = add_mamba_cmds(cmd)
    elif 'replay-plot' in action:
        #  --min-seqs-per-gc 70 --max-seqs-per-gc 70 --n-max-simu-trees 61  # don't want these turned on as long as e.g. N sampled seqs is varying a lot in simulation
        cmd = './projects/replay-plot.py --simu-like-dir %s --outdir %s --plot-labels iqt-data:simu:simu-iqtree --normalize --short-legends --n-max-simu-trees 120 --legend-plots abundances' % (os.path.dirname(ofname(args, varnames, vstrs, 'check-dl' if action=='replay-plot-ckdl' else 'simu')), odr)  #  --n-max-simu-trees 85  # :bst-data
        if args.iqtree_version is not None:
            cmd += ' --iqtree-version %s' % args.iqtree_version
        if args.gcreplay_dir is not None:
            cmd += ' --gcreplay-dir %s' % args.gcreplay_dir
        if args.iqtree_data_dir is not None:
            cmd += ' --iqtree-data-dir %s' % args.iqtree_data_dir
    elif action in ['dl-infer', 'dl-infer-merged', 'group-expts']:
        # NOTE use 'dl-infer' with --dl-model-dir or --dl-model-label-str if the *input* simulation is the same as the *output* label (but different label for model)
        #   but use 'data' below (with --data-dir and either --dl-model-dir or --dl-model-label-str) if they're different, i.e. you have three dirs: one for model, one for input, and one for output
        if 'dl-infer' in action:  # could be 'dl-infer' or 'dl-infer-merged'
            ofn = ofname(args, varnames, vstrs, 'merge-simu' if action=='dl-infer-merged' else ('check-dl' if args.check_dl_input else 'simu'), ftype='npy', current_action=action)
            cmd = 'gcd-dl %s --is-simu --indir %s --outdir %s' % ('train' if args.dl_model_dir is None and args.dl_model_label_str is None else 'infer', os.path.dirname(ofn), odr)
            if args.dl_extra_args is not None:
                cmd += ' %s' % args.dl_extra_args
            # if args.check_dl_input:  # turn this on if you want the many-curves plot even if simulation GCs don't all have the same parameters
            #     cmd += ' --force-many-plot'
            cmd = add_model_dir(cmd)
        else:
            cmd = 'python %s/scripts/group-gcdyn-expts.py --indir %s --outdir %s' % (args.gcddir, os.path.dirname(ofname(args, varnames, vstrs, 'dl-infer')), odr)
            cmd = add_mamba_cmds(cmd)
        cmd = add_scan_args(cmd, skip_fcn=lambda v: v in args.scan_vars['simu'] or v not in args.scan_vars[action] or action=='group-expts' and v in args.scan_vars['dl-infer'] or args.check_dl_input and v=='data-samples')  # ick
    elif action == 'write-partis-simu-files':  # didn't really use/test this (found a way around it), but want to keep it around
        cmd = './bin/gcdyn-simu-run.py --actions process --rm-last-ighj-base --input-simu-dir %s --outdir %s --n-sub-procs %d' % (os.path.dirname(ofname(args, varnames, vstrs, 'simu')), odr, args.n_sub_procs)
    elif action == 'partis':
        assert False  # needs updating
        cmd = './bin/partis partition --species mouse --infname %s --input-partition-fname %s --treefname %s' % (ofname(args, varnames, vstrs, 'simu', ftype='fasta'), ofname(args, varnames, vstrs, 'write-true-partition'), ofname(args, varnames, vstrs, 'simu', ftype='nwk'))
        cmd += ' --parameter-dir %s/parameters --outfname %s' % (odr, ofname(args, varnames, vstrs, action))
        cmd += ' --initial-germline-dir %s --no-insertions-or-deletions --min-selection-metric-cluster-size 3' % args.gcreplay_germline_dir
        cmd += ' --plotdir %s/plots --partition-plot-cfg trees' % odr
        # if args.inference_extra_args is not None:
        #     cmd += ' %s' % args.inference_extra_args
    elif action == 'data':
        cmd = 'gcd-dl infer --indir %s/%s --outdir %s --discard-extra-trees' % (args.data_dir, utils.vlval(args, vstrs, varnames, 'data-samples'), odr)
        cmd = add_model_dir(cmd)
        cmd = add_scan_args(cmd, skip_fcn=lambda v: v in ['data-samples'])
        if args.dl_extra_args is not None:
            cmd += ' %s' % args.dl_extra_args
        if args.dl_bundle_size_list is not None:  # kind of weird to use it for this as well as a scan var, but whatever
            assert False  # I *think* this is now handled by adding the add_scan_args() call just above, but it needs to be checked
            assert len(args.dl_bundle_size_list) == 1
            cmd += ' --dl-bundle-size %d' % int(args.dl_bundle_size_list[0])
    else:
        assert False
    return cmd

# ----------------------------------------------------------------------------------------
# would be nice to combine this also with fcns in cf-tree-metrics.py (and put it in scanplot)
def run_scan(action):
    # ----------------------------------------------------------------------------------------
    def init_cmd(local_varnames, vstrs, ofn, icombo):
        cmd = get_cmd(action, base_args, local_varnames, vstrs, all_simdirs=all_simdirs)
        cmdfos.append({
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : odir(args, local_varnames, vstrs, action) if args.tree_inference_method is None else os.path.dirname(ofn),
            'workdir' : '%s/partis-work/%d' % (args.workdir, icombo),
        })
    # ----------------------------------------------------------------------------------------
    base_args, varnames, val_lists, valstrs = utils.get_var_info(args, args.scan_vars[action], action=action)
    cmdfos = []
    print('  %s: running %d combination%s of: %s' % (utils.color('blue_bkg', action), len(valstrs), utils.plural(len(valstrs)), ' '.join(varnames)))
    if any(len(s)==0 for s in valstrs):  # shouldn't need this any more, but don't want to remove it yet
        raise Exception('zero length val str in valstrs (probably don\'t have any variables with more than one value): %s (varnames: %s)' % (valstrs, varnames))
    if args.debug:
        print('   %s' % ' '.join(varnames))
    n_already_there, n_missing_input, ifn = 0, 0, None
    all_simdirs = []
    for icombo, vstrs in enumerate(valstrs):
        # if vstrs[varnames.index('n-trials')] != '1000': # or vstrs[varnames.index('model-type')] != 'sigmoid':
        #     continue
        # if vstrs[varnames.index('model-type')] != 'sigmoid': # or vstrs[varnames.index('yscale-values')] not in ['5', '10', '20']:
        #     continue
        # if vstrs[varnames.index('carry-cap-values')] != '500': # or vstrs[varnames.index('yscale-values')] not in ['5', '10', '20']:
        #     continue
        # if vstrs[varnames.index('seed')] != '0': # or vstrs[varnames.index('simu-bundle-size')] != '1':
        if 'data-samples' in varnames and vstrs[varnames.index('data-samples')] != 'combo-trees': # or vstrs[varnames.index('simu-bundle-size')] != '1':
            continue
        if args.debug:
            print('   %s' % ' '.join(vstrs))

        ofn = ofname(args, varnames, vstrs, action) if action not in merge_actions else ofname(args, [], [], action)  #, single_file=True)
        if utils.output_exists(args, ofn, debug=False):
            n_already_there += 1
            continue

        if action in ['replay-plot', 'replay-plot-ckdl'] + merge_actions:
            ifn = ofname(args, varnames, vstrs, 'check-dl' if action=='replay-plot-ckdl' else 'simu')
            if not os.path.exists(ifn):  # do *not* use this, it'll delete it if --overwrite is set: utils.output_exists(args, ifn, debug=False):
                n_missing_input += 1
                continue
            if action in merge_actions:
                all_simdirs.append(os.path.dirname(ifn))
                continue
        if (action in ['check-dl', 'replay-plot-ckdl'] or action=='dl-infer' and args.check_dl_input) and utils.vlval(args, vstrs, varnames, 'model-type') != 'sigmoid':  # 'replay-plot-ckdl' doesn't actually get here cause it's missing input, so is skipped above
            # print('    skipping %s for non-sigmoid model' % action)
            continue

        init_cmd(varnames, vstrs, ofn, icombo)

    if action in merge_actions and len(all_simdirs) > 0:
        if action == 'merge-simu' and n_missing_input > 0:
            print('    %s writing merged simulation file despite missing some of its input files' % utils.wrnstr())
        init_cmd([], [], ofn, 0)

    utils.run_scan_cmds(args, cmdfos, '%s.log'%action, len(valstrs), n_already_there, ofn, n_missing_input=n_missing_input, single_ifn=ifn, shell=any('micromamba activate' in c['cmd_str'] for c in cmdfos))

# ----------------------------------------------------------------------------------------
def get_fnfcn(method, pmetr):  # have to adjust signature before passing to fcn in scavars
    def tmpfcn(varnames, vstrs): return ofname(args, varnames, vstrs, method)
    return tmpfcn

# ----------------------------------------------------------------------------------------
import random
random.seed(args.random_seed)
numpy.random.seed(args.random_seed)
if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work')))

for action in args.actions:
    if action in ['simu', 'data', ] + after_actions + plot_actions:
        run_scan(action)
    elif action in ['plot', 'combine-plots']:
        plt_scn_vars = args.scan_vars['group-expts']
        if args.dry:
            print('  --dry: not plotting')
            continue
        _, varnames, val_lists, valstrs = utils.get_var_info(args, plt_scn_vars, action='group-expts')
        if action == 'plot':
            print('plotting %d combinations of %d variable%s (%s) to %s' % (len(valstrs), len(varnames), utils.plural(len(varnames)), ', '.join(varnames), scanplot.get_comparison_plotdir(args, None)))
            fnames = {meth : {pmetr : [] for pmetr in args.perf_metrics} for meth in args.plot_metrics}
            procs = []
            for method in args.plot_metrics:
                for pmetr in args.perf_metrics:
                    utils.prep_dir(scanplot.get_comparison_plotdir(args, method) + '/' + pmetr, wildlings=['*.html', '*.svg', '*.yaml'])  # , subdirs=args.perf_metrics
                for pmetr in args.perf_metrics:
                    print('  %12s %s' % (method, pmetr))
                    arglist, kwargs = (args, plt_scn_vars, action, method, pmetr, args.final_plot_xvar), {'fnfcn' : get_fnfcn(method, pmetr), 'fnames' : fnames[method][pmetr], 'script_base' : script_base, 'debug' : args.debug}
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
                scanplot.make_plots(args, plt_scn_vars, action, None, pmetr, args.final_plot_xvar, fnames=fnames[ipm], make_legend=True, debug=args.debug)
            plotting.make_html(cfpdir, fnames=fnames)
        else:
            raise Exception('unsupported action %s' % action)
    else:
        raise Exception('unsupported action %s' % action)
