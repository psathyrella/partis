#!/usr/bin/env python
import argparse
import os
import sys
import colored_traceback.always
import numpy
import math
import collections
import multiprocessing

sys.path.insert(1, './python')
import utils
import paircluster
import scanplot
import plotting
import clusterpath

# ----------------------------------------------------------------------------------------
partition_types = ['single', 'joint']
all_perf_metrics = ['precision', 'sensitivity', 'f1', 'pcfrac-corr', 'pcfrac-mis', 'pcfrac-un', 'time-reqd', 'cln-frac']  # pcfrac-*: pair info cleaning correct fraction, cln-frac: collision fraction
synth_actions = ['synth-%s'%a for a in ['distance-0.03', 'reassign-0.10', 'singletons-0.40', 'singletons-0.20']]
ptn_actions = ['partition', 'partition-lthresh', 'vsearch-partition', 'vjcdr3-0.9', 'scoper', 'mobille'] + synth_actions  # using the likelihood (rather than hamming-fraction) threshold makes basically zero difference
def is_single_chain(action):
    return 'synth-' in action or 'vjcdr3-' in action or action in ['scoper', 'mobille']

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='simu:cache-parameters:partition:plot')  # can also be merge-paired-partitions and get-selection-metrics
parser.add_argument('--merge-paired-partitions', action='store_true', help='for partis partition actions, don\'t re-partition, just merge paired partitions')
parser.add_argument('--base-outdir', default='%s/partis/paired-loci'%os.getenv('fs'))
parser.add_argument('--n-sim-events-list', default='10', help='N sim events in each repertoire/"proc"/partis simulate run')
parser.add_argument('--n-leaves-list', default='1') #'2:3:4:10') #1 5; do10)
parser.add_argument('--constant-number-of-leaves-list')
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
parser.add_argument('--single-light-locus')
parser.add_argument('--prep', action='store_true', help='only for mobille run script atm')
# scan fwk stuff (mostly):
parser.add_argument('--version', default='v0')
parser.add_argument('--label', default='test')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--make-plots', action='store_true')
parser.add_argument('--test', action='store_true', help='don\'t parallelize \'plot\' action')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--simu-extra-args')
parser.add_argument('--inference-extra-args')
parser.add_argument('--plot-metrics', default='partition', help='NOTE these are methods, but in tree metric script + scanplot they\'re metrics, so we have to call them metrics here')
parser.add_argument('--perf-metrics', default='precision:sensitivity:f1') #':'.join(all_perf_metrics))
parser.add_argument('--zip-vars', help='colon-separated list of variables for which to pair up values sequentially, rather than doing all combinations')
parser.add_argument('--final-plot-xvar', help='variable to put on the x axis of the final comparison plots')
parser.add_argument('--legend-var', help='non-default "component" variable (e.g. obs-frac) to use to label different lines in the legend')
parser.add_argument('--x-legend-var', help='derived variable with which to label the x axis (e.g. mfreq [shm %] when --final-plot-x-var is scratch-mute-freq)')
parser.add_argument('--pvks-to-plot', help='only plot these line/legend values when combining plots')
parser.add_argument('--plot-metric-extra-strs', help='extra strs for each metric in --plot-metrics (i.e. corresponding to what --extra-plotstr was set to during get-tree-metrics for that metric)')
parser.add_argument('--dont-plot-extra-strs', action='store_true', help='while we still use the strings in --plot-metric-extra-strs to find the right dir to get the plot info from, we don\'t actually put the str in the plot (i.e. final plot versions where we don\'t want to see which dtr version it is)')
parser.add_argument('--combo-extra-str', help='extra label for combine-plots action i.e. write to combined-%s/ subdir instead of combined/')
parser.add_argument('--workdir')  # default set below
args = parser.parse_args()
args.scan_vars = {'simu' : ['seed', 'n-leaves', 'constant-number-of-leaves', 'scratch-mute-freq', 'mutation-multiplier', 'mean-cells-per-droplet', 'fraction-of-reads-to-remove', 'allowed-cdr3-lengths', 'n-genes-per-region', 'n-sim-alleles-per-gene', 'n-sim-events']}
for act in ['cache-parameters'] + ptn_actions:
    args.scan_vars[act] = args.scan_vars['simu']
args.str_list_vars = ['allowed-cdr3-lengths', 'n-genes-per-region', 'n-sim-alleles-per-gene']
args.bool_args = ['constant-number-of-leaves']  # NOTE different purpose to svartypes below (this isn't to convert all the values to the proper type, it's just to handle flag-type args
# NOTE ignoring svartypes atm, which may actually work?
# args.svartypes = {'int' : ['n-leaves', 'allowed-cdr3-lengths', 'n-sim-events'], 'float' : ['scratch-mute-freq', 'mutation-multiplier']}  # 'mean-cells-per-droplet' # i think can't float() this since we want to allow None as a value
# and these i think we can't since we want to allow blanks, 'n-genes-per-region' 'n-sim-alleles-per-gene'
# args.float_var_digits = 2  # ick

args.actions = utils.get_arg_list(args.actions, choices=['simu', 'cache-parameters', 'merge-paired-partitions', 'get-selection-metrics', 'plot', 'combine-plots'] + ptn_actions)
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

utils.get_scanvar_arg_lists(args)
if args.final_plot_xvar is None:  # set default value based on scan vars
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars['simu'])
    svars = [v for v in varnames if v != 'seed']
    args.final_plot_xvar = svars[0] if len(svars) > 0 else 'seed'  # if we're not scanning over any vars, i'm not sure what we should use

# ----------------------------------------------------------------------------------------
def odir(args, varnames, vstrs, action):
    return '%s/%s' % (utils.svoutdir(args, varnames, vstrs, action), 'partis' if action in ['cache-parameters', 'partition'] else action)

# ----------------------------------------------------------------------------------------
def ofname(args, varnames, vstrs, action, locus=None, single_chain=False, single_file=False, logfile=False, pcleancsv=False):
    outdir = odir(args, varnames, vstrs, action)
    if action == 'cache-parameters' and not logfile:
        outdir += '/parameters'
        if locus is None and not single_file:
            return outdir
    if single_file:
        assert locus is None
        locus = 'igk'
    assert locus is not None
    if logfile:
        ofn = '%s/%s%s.log' % (outdir, 'work/%s/'%locus if action=='mobille' else '',  action)
    elif pcleancsv:
        ofn = '%s/true-pair-clean-performance.csv' % outdir
    elif action == 'cache-parameters':
        ofn = '%s/%s' % (outdir, locus)
        if single_file:
            # ofn += '/hmm/germline-sets/%s/%sv.fasta' % (locus, locus)
            ofn += '/hmm/all-mean-mute-freqs.csv'
    else:
        ofn = paircluster.paired_fn(outdir, locus, suffix='.yaml', actstr=None if action=='simu' else 'partition', single_chain=single_chain or is_single_chain(action))
    return ofn

# ----------------------------------------------------------------------------------------
# distance-based partitions get made by running partis, but then we make the other types here
def make_synthetic_partition(action, varnames, vstrs):
    for ltmp in plot_loci():
        _, _, true_cpath = utils.read_output(paircluster.paired_fn(odir(args, varnames, vstrs, 'simu'), ltmp, suffix='.yaml'))
        _, mistype, misfrac = action.split('-')
        new_partition = utils.generate_incorrect_partition(true_cpath.best(), float(misfrac), mistype)
        new_cpath = clusterpath.ClusterPath(partition=new_partition)
        new_cpath.calculate_missing_values(true_partition=true_cpath.best())
        ofn = ofname(args, varnames, vstrs, action, locus=ltmp, single_chain=True)
        utils.write_annotations(ofn, None, [], utils.annotation_headers, partition_lines=new_cpath.get_partition_lines())
        print '    %s: wrote synthetic partition to %s' % (ltmp, ofn)

# ----------------------------------------------------------------------------------------
def get_cmd(action, base_args, varnames, vstrs, synth_frac=None):
    if action == 'scoper':
        cmd = './test/scoper-run.py --indir %s --outdir %s --simdir %s' % (ofname(args, varnames, vstrs, 'cache-parameters'), odir(args, varnames, vstrs, action), odir(args, varnames, vstrs, 'simu'))
        return cmd
    if action == 'mobille':
        cmd = './test/mobille-run.py --simdir %s --outdir %s --id-str %s --base-imgt-outdir %s' % (odir(args, varnames, vstrs, 'simu'), odir(args, varnames, vstrs, action), '_'.join('%s-%s'%(n, s) for n, s in zip(varnames, vstrs)), '%s/%s/imgt-output' % (args.base_outdir, args.label))
        if args.prep:
            cmd += ' --prep'
            # then do this by hand, and submit to imgt/high vquest by hand, then download results and put them in the proper dir (run mobille run script to get dir)
            # tar czf /path/somewhere/to/rsync/imgt-input.tgz /fh/local/dralph/partis/paired-loci/vs-shm/v2/seed-*/scratch-mute-freq-*/mobille/work/*/imgt-input/*.fa
        return cmd
    actstr = action
    if 'synth-distance-' in action or action == 'vsearch-partition' or action == 'partition-lthresh':
        actstr = 'partition'
    if 'vjcdr3-' in action:
        actstr = 'annotate'
    if args.merge_paired_partitions:
        actstr = 'merge-paired-partitions'
    cmd = './bin/partis %s --paired-loci --paired-outdir %s' % (actstr.replace('simu', 'simulate'), odir(args, varnames, vstrs, action))
    if args.n_sub_procs > 1:
        cmd += ' --n-procs %d' % args.n_sub_procs
    if action == 'simu':
        cmd += ' --simulate-from-scratch --no-per-base-mutation %s' % ' '.join(base_args)
        if args.single_light_locus is not None:
            cmd += ' --single-light-locus %s' % args.single_light_locus
        if args.simu_extra_args is not None:
            cmd += ' %s' % args.simu_extra_args
        for vname, vstr in zip(varnames, vstrs):
            if vname in args.bool_args:  # in cf-tree-metrics this was handled for e.g. context dependence in bcr-phylo-run, but anyway this should probably be moved to scanvar stuff in utils
                if vstr == '0':
                    continue
                elif vstr == '1':
                    vstr = ''
                else:
                    assert False
            cmd += ' --%s %s' % (vname, vstr)
    else:
        cmd += ' --paired-indir %s' % odir(args, varnames, vstrs, 'simu')
        if action == 'vsearch-partition':
            cmd += ' --naive-vsearch'
        if 'synth-distance-' in action:
            synth_hfrac = float(action.replace('synth-distance-', ''))
            cmd += ' --synthetic-distance-based-partition --naive-hamming-bounds %.2f:%.2f' % (synth_hfrac, synth_hfrac)
        if 'vjcdr3-' in action:
            cmd += ' --annotation-clustering --annotation-clustering-threshold %.2f' % float(action.split('-')[1])
        if action == 'partition-lthresh':
            cmd += ' --paired-naive-hfrac-threshold-type likelihood'
        if action != 'get-selection-metrics':  # it just breaks here because i don't want to set --simultaneous-true-clonal-seqs (but maybe i should?)
            cmd += ' --is-simu'
        if action != 'cache-parameters':
            cmd += ' --refuse-to-cache-parameters'
        if 'synth-distance-' in action or 'vjcdr3-' in action or action == 'vsearch-partition' or action == 'partition-lthresh':
            cmd += ' --parameter-dir %s' % ofname(args, varnames, vstrs, 'cache-parameters')
        if action in ptn_actions and 'vjcdr3-' not in action and not args.make_plots:
            cmd += ' --dont-calculate-annotations'
        if args.make_plots:
            cmd += ' --plotdir paired-outdir'
            if action in ptn_actions:
                cmd += ' --no-partition-plots' #--no-mds-plots' #
        if args.inference_extra_args is not None:
            cmd += ' %s' % args.inference_extra_args

    return cmd

# ----------------------------------------------------------------------------------------
# TODO combine this also with fcns in cf-tree-metrics.py (and put it in scanplot)
def run_scan(action):
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars[action])
    cmdfos = []
    print '  %s: running %d combinations of: %s' % (utils.color('blue_bkg', action), len(valstrs), ' '.join(varnames))
    if args.debug:
        print '   %s' % ' '.join(varnames)
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print '   %s' % ' '.join(vstrs)

        ofn = ofname(args, varnames, vstrs, action, single_file=True)
        if args.merge_paired_partitions:
            assert action in ['partition', 'vsearch-partition']
        else:
            if utils.output_exists(args, ofn, debug=False):
                n_already_there += 1
                continue

        if 'synth-reassign-' in action or 'synth-singletons-' in action:
            make_synthetic_partition(action, varnames, vstrs)
            continue

        cmd = get_cmd(action, base_args, varnames, vstrs)
        # utils.simplerun(cmd, logfname='%s-%s.log'%(odir(args, varnames, vstrs, action), action), dryrun=args.dry)
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : ofn,
            'logdir' : odir(args, varnames, vstrs, action),
            'workdir' : '%s/partis-work/%d' % (args.workdir, icombo),
        }]

    utils.run_scan_cmds(args, cmdfos, '%s.log'%(action if not args.merge_paired_partitions else 'merge-paired-partitions'), len(valstrs), n_already_there, ofn)

# ----------------------------------------------------------------------------------------
def plot_loci():
    if args.single_light_locus is None:
        return utils.sub_loci('ig')
    else:
        return [utils.heavy_locus('ig'), args.single_light_locus]

# ----------------------------------------------------------------------------------------
def get_fnfcn(method, locus, ptntype, pmetr):
    def tmpfcn(varnames, vstrs): return ofname(args, varnames, vstrs, method, locus=locus, single_chain=ptntype=='single', logfile=pmetr=='time-reqd', pcleancsv='pcfrac-' in pmetr)
    return tmpfcn

# ----------------------------------------------------------------------------------------
def get_pdirfcn(locus):
    def tmpfcn(varnames, vstrs): return ofname(args, varnames, vstrs, 'cache-parameters', locus=locus)
    return tmpfcn

# ----------------------------------------------------------------------------------------
def get_ptntypes(method):
    ptypes = partition_types
    if is_single_chain(method):
        ptypes = [t for t in ptypes if t != 'joint']  # this is probably needlessly general
    return ptypes

# ----------------------------------------------------------------------------------------
import random
random.seed(args.random_seed)
numpy.random.seed(args.random_seed)
if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work')))

for action in args.actions:
    if action in ['simu', 'cache-parameters'] + ptn_actions:
        run_scan(action)
    elif action in ['plot', 'combine-plots']:
        if args.dry:
            print '  --dry: not plotting'
            continue
        _, varnames, val_lists, valstrs = utils.get_var_info(args, args.scan_vars['partition'])
        if action == 'plot':
            print 'plotting %d combinations of %d variable%s (%s) to %s' % (len(valstrs), len(varnames), utils.plural(len(varnames)), ', '.join(varnames), scanplot.get_comparison_plotdir(args, None))
            fnames = {meth : {pmetr : [[] for _ in partition_types] for pmetr in args.perf_metrics} for meth in args.plot_metrics}
            procs = []
            for method in args.plot_metrics:  # NOTE in cf-tree-metrics.py these are [selection] metrics, but here they're [clustering] methods
                for pmetr in args.perf_metrics:
                    utils.prep_dir(scanplot.get_comparison_plotdir(args, method) + '/' + pmetr, wildlings=['*.html', '*.svg', '*.yaml'])  # , subdirs=args.perf_metrics
                for ipt, ptntype in enumerate(get_ptntypes(method)):
                    for ltmp in plot_loci():
                        for pmetr in args.perf_metrics:
                            if 'pcfrac-' in pmetr and (ptntype != 'joint' or ltmp != 'igh'):
                                continue
                            print '  %12s  %6s partition: %3s %s' % (method, ptntype.replace('single', 'single chain'), ltmp, pmetr)
                            arglist, kwargs = (args, args.scan_vars['partition'], action, method, pmetr, plotting.legends.get(pmetr, pmetr), args.final_plot_xvar), {'fnfcn' : get_fnfcn(method, ltmp, ptntype, pmetr), 'locus' : ltmp, 'ptntype' : ptntype, 'fnames' : fnames[method][pmetr][ipt], 'pdirfcn' : get_pdirfcn(ltmp), 'debug' : args.debug}
                            if args.test:
                                scanplot.make_plots(*arglist, **kwargs)
                            else:
                                procs.append(multiprocessing.Process(target=scanplot.make_plots, args=arglist, kwargs=kwargs))
            if not args.test:
                utils.run_proc_functions(procs)
            for method in args.plot_metrics:
                for pmetr in args.perf_metrics:
                    pmcdir = scanplot.get_comparison_plotdir(args, method) + '/' + pmetr
                    fnames[method][pmetr] = [[f.replace(pmcdir, '') for f in flist] for flist in fnames[method][pmetr]]
                    plotting.make_html(pmcdir, n_columns=3, fnames=fnames[method][pmetr])  # this doesn't work unless --test is set, and i'm not sure why, but whatever, just run combine-plots
        elif action == 'combine-plots':
            cfpdir = scanplot.get_comparison_plotdir(args, 'combined')
            utils.prep_dir(cfpdir, wildlings=['*.html', '*.svg'])
            fnames, iplot = [[] for _ in args.perf_metrics], 0
            for ipm, pmetr in enumerate(args.perf_metrics):
                print '    ', pmetr
                for ptntype in partition_types:
                    for ltmp in plot_loci():
                        if 'pcfrac-' in pmetr and (ptntype != 'joint' or ltmp != 'igh'):
                            continue
                        scanplot.make_plots(args, args.scan_vars['partition'], action, None, pmetr, plotting.legends.get(pmetr, pmetr), args.final_plot_xvar, locus=ltmp, ptntype=ptntype, fnames=fnames[ipm], make_legend=iplot==0, debug=args.debug)
                        iplot += 1
            fnames += [[os.path.dirname(fnames[0][0]) + '/legend.svg']]
            plotting.make_html(cfpdir, n_columns=3 if len(plot_loci())==3 else 4, fnames=fnames)
        else:
            raise Exception('unsupported action %s' % action)
    else:
        raise Exception('unsupported action %s' % action)

# bd=_output/cells-per-drop
# subd=inferred/plots
# ./bin/compare-plotdirs.py --outdir ~/Dropbox/tmp-plots/cells-per-drop \
#      --normalize --translegend=-0.2:-0.2 \
#      --plotdirs $bd-1.0/$subd:$bd-1.2/$subd:$bd-1.7/$subd:$bd-2.0/$subd:$bd-3.0/$subd \
#      --names 1.0:1.2:1.7:2.0:3.0
