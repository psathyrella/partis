#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import operator
import os
import sys
import yaml
import colored_traceback.always
import numpy
import subprocess
import multiprocessing
import copy
from io import open

# ----------------------------------------------------------------------------------------
def ura_vals(xvar):  # list of whether we're using relative affinity values
    if xvar == 'affinity' and args.include_relative_affy_plots:
        return [False, True]  # i.e. [don'e use relative affy, do use relative affy]
    else:
        return [False]
# ----------------------------------------------------------------------------------------
def get_n_generations(ntl, tau):  # NOTE duplicates code in treeutils.calculate_max_lbi()
    return max(1, int(args.seq_len * tau * ntl))

# ----------------------------------------------------------------------------------------
def get_outfname(outdir):
    return '%s/vals.yaml' % outdir

# ----------------------------------------------------------------------------------------
def make_lb_bound_plots(args, outdir, metric, btype, parsed_info, print_results=False):
    def make_plot(log, parsed_info):
        fig, ax = plotting.mpl_init()
        for lbt in sorted(parsed_info[metric][btype], reverse=True):
            n_gen_list, lb_list = zip(*sorted(list(parsed_info[metric][btype][lbt].items()), key=operator.itemgetter(0)))
            if lbt == 1. / args.seq_len:  # add a horizontal line corresponding to the asymptote for tau = 1/seq_len
                ax.plot(ax.get_xlim(), (lb_list[-1], lb_list[-1]), linewidth=3, alpha=0.7, color='darkred', linestyle='--') #, label='1/seq len')
            ax.plot(n_gen_list, lb_list, label='%.4f' % lbt, alpha=0.7, linewidth=4)
            if print_results and log == '':
                print('      %7.4f    %6.4f' % (lbt, lb_list[-1]))
        plotname = 'tau-vs-n-gen-vs-%s-%s' % (btype, metric)
        if log == '':
            ybounds = None
            leg_loc = (0.1, 0.5)
        else:
            plotname += '-log'
            if metric == 'lbi':
                ybounds = (2*min(parsed_info[metric][btype]), 3*ax.get_ylim()[1])
            else:
                ybounds = None
            leg_loc = (0.04, 0.57)
        plotting.mpl_finish(ax, outdir, plotname, log=log, xbounds=(min(n_gen_list), max(n_gen_list)), ybounds=ybounds,
                            xlabel='N generations', ylabel='%s %s' % (btype.capitalize(), metric.upper()), leg_title='tau', leg_prop={'size' : 12}, leg_loc=leg_loc)

    if print_results:
        print('%s:     tau    %s %s' % (btype, btype, metric))

    for log in ['', 'y']:
        make_plot(log, parsed_info)

# ----------------------------------------------------------------------------------------
def calc_lb_bounds(args, n_max_gen_to_plot=4, lbt_bounds=(0.001, 0.005), print_results=False):
    print_results = True
    btypes = ['min', 'max']

    outdir = '%s/lb-tau-normalization/%s' % (args.base_outdir, args.label)

    parsed_info = {m : {b : {} for b in btypes} for m in args.only_metrics}
    for lbt in args.lb_tau_list:
        if args.make_plots and (lbt < lbt_bounds[0] or lbt > lbt_bounds[1]):
            print('    skipping tau %.4f outside of bounds [%.4f, %.4f] for bound plotting' % (lbt, lbt_bounds[0], lbt_bounds[1]))
            continue

        gen_list = args.n_generations_list
        if gen_list is None:
            gen_list = [get_n_generations(ntl, lbt) for ntl in args.n_tau_lengths_list]
        if args.lb_tau_list.index(lbt) == 0 or args.n_tau_lengths_list is not None:  # if --n-tau-lengths-list is set, they could be different for each tau
            print(' seq len: %d   N gen list: %s' % (args.seq_len, ' '.join(str(n) for n in gen_list)))
        print('   %s %.4f' % (utils.color('green', 'lb tau'), lbt))
        for n_gen in gen_list:
            if args.debug:
                print('     %s %d  %s %.4f' % (utils.color('purple', 'N gen'), n_gen, utils.color('purple', 'lb tau'), lbt))

            this_outdir = '%s/n_gen-%d-lbt-%.4f' % (outdir, n_gen, lbt)  # if for some reason I want to specify --n-tau-lengths-list instead of --n-generations-list, this makes the output path structure still correspond to n generations, but that's ok since that's what the trees do

            if args.make_plots:  # just let it crash if you forgot to actually run it first
                with open(get_outfname(this_outdir)) as outfile:
                    info = yaml.load(outfile, Loader=yaml.Loader)
                for metric in args.only_metrics:
                    for btype in btypes:
                        if lbt not in parsed_info[metric][btype]:
                            parsed_info[metric][btype][lbt] = {}
                        parsed_info[metric][btype][lbt][n_gen] = info[metric][btype][metric]  # it's weird to have metric as the key twice, but it actually makes sense since we're subverting the normal calculation infrastructure to only run one metric or the other at a time (i.e. the righthand key is pulling out the metric we want from the lb info that, in principle, usually has both; while the lefthand key is identifying a run during which we were only caring about that metric)
                continue
            elif utils.output_exists(args, get_outfname(this_outdir)):
                continue

            print('         running n gen %d' % n_gen)

            if not os.path.exists(this_outdir):
                os.makedirs(this_outdir)

            lbvals = treeutils.calculate_lb_bounds(args.seq_len, lbt, n_generations=n_gen, n_offspring=args.max_lb_n_offspring, only_metrics=args.only_metrics, btypes=btypes, debug=args.debug)

            with open(get_outfname(this_outdir), 'w') as outfile:
                yamlfo = {m : {b : {k : v for k, v in lbvals[m][b].items() if k != 'vals'} for b in btypes} for m in args.only_metrics}  # writing these to yaml is really slow, and they're only used for plotting below
                yaml.dump(yamlfo, outfile)

            if n_gen > n_max_gen_to_plot:
                continue

            # this is really slow on large trees
            plotdir = this_outdir + '/plots'
            utils.prep_dir(plotdir, wildlings='*.svg')
            for metric in args.only_metrics:
                for btype in btypes:
                    if lbvals[metric][btype]['vals'] is None:
                        continue
                    cmdfos = [lbplotting.get_lb_tree_cmd(lbvals[metric][btype]['vals']['tree'], '%s/%s-%s-tree.svg' % (plotdir, metric, btype), metric, 'affinities', args.workdir, metafo=lbvals[metric][btype]['vals'], tree_style='circular')]
                    utils.run_cmds(cmdfos, clean_on_success=True, shell=True, debug='print')

    if args.make_plots:
        print('  writing plots to %s' % outdir)
        for metric in args.only_metrics:
            for btype in btypes:
                if 'lbr' in metric and btype == 'min':  # it's just zero, and confuses the log plots
                    continue
                if len(parsed_info[metric][btype]) == 0:
                    print('nothing to do (<parsed_info> empty)')
                    continue
                make_lb_bound_plots(args, outdir, metric, btype, parsed_info, print_results=print_results)

# ----------------------------------------------------------------------------------------
def get_bcr_phylo_outdir(varnames, vstr):
    return utils.svoutdir(args, varnames, vstr, 'simu') + '/bcr-phylo'

# ----------------------------------------------------------------------------------------
def get_simfname(varnames, vstr):
    return '%s/selection/simu/mutated-simu.yaml' % get_bcr_phylo_outdir(varnames, vstr)

# ----------------------------------------------------------------------------------------
def get_parameter_dir(varnames, vstr):
    return '%s/selection/partis/params' % get_bcr_phylo_outdir(varnames, vstr)

# ----------------------------------------------------------------------------------------
def get_tree_metric_outdir(varnames, vstr, metric_method=None):  # metric_method is only set if it's neither lbi nor lbr
    return utils.svoutdir(args, varnames, vstr, 'get-tree-metrics') + '/' + ('partis' if metric_method is None else metric_method)

# ----------------------------------------------------------------------------------------
def get_partition_fname(varnames, vstr, action, metric_method=None):  # if action is 'bcr-phylo', we want the original partition output file, but if it's 'get-tree-metrics', we want the copied one, that had tree metrics added to it (and which is in the e.g. tau subdir) UPDATE no longer modifying output files by default, so no longer doing the copying thing
    if action == 'bcr-phylo' or metric_method is not None:  # for non-lb metrics (i.e. if metric_method is set) we won't modify the partition file, so can just read the one in the bcr-phylo dir
        outdir = '%s/selection/partis' % get_bcr_phylo_outdir(varnames, vstr)
    else:
        outdir = get_tree_metric_outdir(varnames, vstr)
    return '%s/partition.yaml' % outdir

# ----------------------------------------------------------------------------------------
def get_tree_metric_plotdir(varnames, vstr, metric_method=None, extra_str=''):
    return  '%s/%splots' % (get_tree_metric_outdir(varnames, vstr, metric_method=metric_method), extra_str+'-' if extra_str != '' else '')

# ----------------------------------------------------------------------------------------
def get_dtr_model_dir(varnames, vstr, extra_str=''):
    plotdir = get_tree_metric_plotdir(varnames, vstr, metric_method='dtr', extra_str=extra_str)
    if plotdir.split('/')[-1] == 'plots':
        assert extra_str == ''  # i think?
        delim = '/'
    elif plotdir.split('-')[-1] == 'plots':  # i.e. args.extra_plotstr was set when we trained
        assert extra_str != ''  # i think?
        delim = '-'
    else:
        assert False
    plstr = delim + 'plots'
    assert plotdir.count(plstr) == 1
    return plotdir.replace(plstr, delim + 'dtr-models')

# ----------------------------------------------------------------------------------------
def rel_affy_str(use_relative_affy, metric):
    return '-relative' if use_relative_affy and metric == 'lbi' else ''

# ----------------------------------------------------------------------------------------
def get_tm_fname(varnames, vstr, metric, x_axis_label, use_relative_affy=False, cg=None, tv=None, extra_str=''):  # note that there are separate svg files for each iclust, but info for all clusters are written to the same yaml file (but split apart with separate info for each cluster)
    if metric == 'dtr':
        assert cg is not None and tv is not None
    # if metric in ['lbi', 'lbr']:  # NOTE using <metric> and <metric_method> for slightly different but overlapping things: former is the actual metric name, whereas setting the latter says we want a non-lb metric (necessary because by default we want to calculate lbi and lbr, but also be able treat lbi and lbr separately when plotting)
    #     plotdir = get_tree_metric_plotdir(varnames, vstr, extra_str=extra_str)
    #     old_path = '%s/true-tree-metrics/%s-vs-%s-true-tree-ptiles%s.yaml' % (plotdir, metric, x_axis_label, rel_affy_str(use_relative_affy, metric))  # just for backwards compatibility, could probably remove at some point (note: not updating this when I'm adding non-lb metrics like shm)
    #     if os.path.exists(old_path):
    #         return old_path
    # else:
    plotdir = get_tree_metric_plotdir(varnames, vstr, metric_method=metric, extra_str=extra_str)
    return treeutils.tmfname(plotdir, metric, x_axis_label, cg=cg, tv=tv, use_relative_affy=use_relative_affy)

# ----------------------------------------------------------------------------------------
def get_all_tm_fnames(varnames, vstr, metric_method=None, extra_str=''):
    if metric_method is None:
        return [get_tm_fname(varnames, vstr, mtmp, xatmp, use_relative_affy=use_relative_affy, extra_str=extra_str)
                for mtmp, cfglist in lbplotting.lb_metric_axis_cfg(args.metric_method)
                for xatmp, _ in cfglist
                for use_relative_affy in ura_vals(xatmp)]  # arg wow that's kind of complicated and ugly
    elif metric_method == 'dtr':
        if args.train_dtr:  # training
            return [treeutils.dtrfname(get_dtr_model_dir(varnames, vstr, extra_str=extra_str), cg, tvar) for cg in treeutils.cgroups for tvar in treeutils.dtr_targets[cg]]
        else:  # testing
            return [get_tm_fname(varnames, vstr, metric_method, lbplotting.getptvar(tv), cg=cg, tv=tv, use_relative_affy=use_relative_affy, extra_str=extra_str) for cg in treeutils.cgroups for tv in treeutils.dtr_targets[cg] for use_relative_affy in ura_vals(tv)]
    else:
        return [get_tm_fname(varnames, vstr, metric_method, 'n-ancestor' if metric_method in treeutils.daffy_metrics else 'affinity', extra_str=extra_str)]  # this hard coding sucks, and it has to match some stuff in treeutils.calculate_non_lb_tree_metrics()

# ----------------------------------------------------------------------------------------
def run_bcr_phylo(args):  # also caches parameters
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars['simu'])
    cmdfos = []
    print('  bcr-phylo: running %d combinations of: %s' % (len(valstrs), ' '.join(varnames)))
    if args.debug:
        print('   %s' % ' '.join(varnames))
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print('   %s' % ' '.join(vstrs))
        outdir = get_bcr_phylo_outdir(varnames, vstrs)
        if utils.output_exists(args, get_partition_fname(varnames, vstrs, 'bcr-phylo'), offset=8, debug=args.debug):
            n_already_there += 1
            continue
        cmd = './bin/bcr-phylo-run.py --actions %s --dont-get-tree-metrics --base-outdir %s %s' % (args.bcr_phylo_actions, outdir, ' '.join(base_args))
        for vname, vstr in zip(varnames, vstrs):
            vstr_for_cmd = vstr
            if vname == 'parameter-variances':
                vstr_for_cmd = vstr_for_cmd.replace('_c_', ':')  # necessary so we can have multiple different parameters with variances for each bcr-phylo-run.py cmd
                cmd += ' --%s %s' % (vname, vstr_for_cmd)
            else:
                cmd = utils.add_to_scan_cmd(args, vname, vstr, cmd)
            if 'context' in vname:
                cmd += ' --restrict-available-genes'
        if args.no_scan_parameter_variances is not None:
            cmd += ' --parameter-variances %s' % args.no_scan_parameter_variances  # we don't parse through this at all here, which means it's the same for all combos of variables (which I think makes sense -- we probably don't even really want to vary most variables if this is set)
        if args.n_sim_events_per_proc is not None:
            cmd += ' --n-sim-events %d' % args.n_sim_events_per_proc
        if args.n_max_queries is not None:
            cmd += ' --n-max-queries %d' % args.n_max_queries
        if args.dont_observe_common_ancestors:
            cmd += ' --dont-observe-common-ancestors'
        if args.overwrite:
            cmd += ' --overwrite'
        if args.only_csv_plots:
            cmd += ' --only-csv-plots'
        if args.n_sub_procs > 1:
            cmd += ' --n-procs %d' % args.n_sub_procs
        if args.sub_slurm:
            cmd += ' --slurm'
        if args.simu_extra_args is not None:
            cmd += ' %s' % args.simu_extra_args
        # cmd += ' --debug 2'
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : get_partition_fname(varnames, vstrs, 'bcr-phylo'),
            'logdir' : outdir,
            'workdir' : '%s/bcr-phylo-work/%d' % (args.workdir, icombo),
        }]
    if n_already_there > 0:
        print('      %d / %d skipped (outputs exist, e.g. %s)' % (n_already_there, len(valstrs), get_partition_fname(varnames, vstrs, 'bcr-phylo')))
    if len(cmdfos) > 0:
        if args.dry:
            print('    %s' % '\n    '.join(cfo['cmd_str'] for cfo in cmdfos))
        else:
            print('      starting %d jobs' % len(cmdfos))
            utils.run_cmds(cmdfos, debug='write:bcr-phylo.log', batch_system='slurm' if args.slurm else None, n_max_procs=args.n_max_procs, allow_failure=True)

# ----------------------------------------------------------------------------------------
def get_tree_metrics(args):
    _, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars['get-tree-metrics'])  # can't use base_args a.t.m. since it has the simulation/bcr-phylo args in it
    cmdfos = []
    print('  get-tree-metrics (%s): running %d combinations of: %s' % (args.metric_method, len(valstrs), ' '.join(varnames)))
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print('   %s' % ' '.join(vstrs))

        if utils.all_outputs_exist(args, get_all_tm_fnames(varnames, vstrs, metric_method=args.metric_method, extra_str=args.extra_plotstr), outlabel='get-tree-metrics', offset=8, debug=args.debug):
            n_already_there += 1
            continue

        if not args.dry:
            tmpoutdir = get_tree_metric_outdir(varnames, vstrs, metric_method=args.metric_method)
            if not os.path.isdir(tmpoutdir):
                os.makedirs(tmpoutdir)

        tmpdir = get_tree_metric_plotdir(varnames, vstrs, metric_method=args.metric_method, extra_str=args.extra_plotstr)  # if we're running lbi/lbr as a --metric-method to avoid running partis get-selection-metrics

        # old way, but don't want to delete it atm:
        # tmpdir = get_tree_metric_plotdir(varnames, vstrs, metric_method=None if args.metric_method in ['lbi', 'lbr'] else args.metric_method, extra_str=args.extra_plotstr)  # if we're running lbi/lbr as a --metric-method to avoid running partis get-selection-metrics
        # if args.metric_method is None:
        #     cmd = './bin/partis get-selection-metrics --is-simu --infname %s --plotdir %s --outfname %s --selection-metric-fname %s' % (get_simfname(varnames, vstrs), tmpdir,
        #                                                                                                                                 get_partition_fname(varnames, vstrs, 'bcr-phylo'), utils.insert_before_suffix('-selection-metrics', get_partition_fname(varnames, vstrs, 'get-tree-metrics')))  # we don't actually use the --selection-metric-fname for anything, but if we don't set it then all the different get-tree-metric jobs write their output files to the same selection metric file in the bcr-phylo dir
        #     cmd += ' --seed %s' % args.random_seed  # NOTE second/commented version this is actually wrong: vstrs[varnames.index('seed')]  # there isn't actually a reason for different seeds here (we want the different seeds when running bcr-phylo), but oh well, maybe it's a little clearer this way
        #     cmd += ' --selection-metrics-to-calculate lbi:lbr'  # TODO it would be better to just always use smetric-run.py, which you can do now, but i don't want to break backwards compatibility
        #     if args.no_tree_plots:
        #         assert False  # needs to be updated to modify plot cfg arg
        #         cmd += ' --ete-path None'
        #     # if args.n_sub_procs > 1:  # TODO get-tree-metrics doesn't paralellize anything atm
        #     #     cmd += ' --n-procs %d' % args.n_sub_procs
        # else:

        cmd = './bin/smetric-run.py --metric-method %s --infname %s --base-plotdir %s' % (args.metric_method, get_simfname(varnames, vstrs), tmpdir)
        if args.metric_method == 'dtr':
            if args.train_dtr and args.overwrite:  # make sure no training files exist, since we don\'t want treeutils.train_dtr_model() to overwrite existing ones (since training can be really slow)
                assert set([os.path.exists(f) for f in get_all_tm_fnames(varnames, vstrs, metric_method=args.metric_method, extra_str=args.extra_plotstr)]) == set([False])
            cmd += ' --action %s' % ('train' if args.train_dtr else 'test')
            cmd += ' --dtr-path %s' % (args.dtr_path if args.dtr_path is not None else get_dtr_model_dir(varnames, vstrs, extra_str=args.extra_plotstr))  # if --dtr-path is set, we're reading the model from there; otherwise we write a new model to the normal/auto location for these parameters (i.e. the point of --dtr-path is to point at the location from a different set of parameters)
            if args.dtr_cfg is not None:
                cmd += ' --dtr-cfg %s' % args.dtr_cfg
        if not args.only_csv_plots and not args.no_tree_plots:
            cmd += ' --make-tree-plots'
        if args.lb_tau_list is not None:
            cmd += ' --lb-tau %s' % utils.vlval(args, vstrs, varnames, 'lb-tau')
            if len(args.lb_tau_list) > 1:
                cmd += ' --dont-normalize-lbi'
        if args.only_csv_plots:
            cmd += ' --only-csv-plots'
        if args.n_max_queries is not None:
            cmd += ' --n-max-queries %d' % args.n_max_queries
        cmd += ' --min-selection-metric-cluster-size 5'  # if n per gen is small, sometimes the clusters are a bit smaller than 10, but we don't really want to skip any clusters here (especially because it confuses the plotting above)
        if args.include_relative_affy_plots:
            cmd += ' --include-relative-affy-plots'

        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : get_all_tm_fnames(varnames, vstrs, metric_method=args.metric_method, extra_str=args.extra_plotstr)[0],
            'workdir' : tmpdir,
        }]

    # TODO use utils.run_scan_cmds()
    if n_already_there > 0:
        print('      %d / %d skipped (outputs exist, e.g. %s)' % (n_already_there, len(valstrs), get_all_tm_fnames(varnames, vstrs, metric_method=args.metric_method, extra_str=args.extra_plotstr)[0]))
    if len(cmdfos) > 0:
        print('      %s %d jobs' % ('--dry: would start' if args.dry else 'starting', len(cmdfos)))
        if args.dry:
            print('  first command: %s' % cmdfos[0]['cmd_str'])
        else:
            logstr = 'get-tree-metrics'
            if args.metric_method == 'dtr':
                logstr += '-train' if args.train_dtr else '-test'
            utils.run_cmds(cmdfos, debug='write:%s.log'%logstr, batch_system='slurm' if args.slurm else None, n_max_procs=args.n_max_procs, allow_failure=True)

# ----------------------------------------------------------------------------------------
all_actions = ['get-lb-bounds', 'bcr-phylo', 'get-tree-metrics', 'plot', 'combine-plots']
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--actions', default=':'.join(a for a in all_actions if a not in ['get-lb-bounds', 'combine-plots']))
parser.add_argument('--bcr-phylo-actions', default='simu:cache-parameters:partition')
parser.add_argument('--test', action='store_true')
parser.add_argument('--carry-cap-list', default='1000')
parser.add_argument('--n-sim-seqs-per-gen-list', default='30:50:75:100:150:200', help='colon-separated list of comma-separated lists of the number of sequences for bcr-phylo to sample at the times specified by --obs-times-list')
parser.add_argument('--n-sim-events-per-proc', type=int, help='number of rearrangement events to simulate in each process (default is set in bin/bcr-phylo-run.py)')
parser.add_argument('--obs-times-list', default='125,150', help='colon-separated list of comma-separated lists of bcr-phylo observation times')
parser.add_argument('--lb-tau-list') #, default='0.0005:0.001:0.002:0.003:0.004:0.008:0.012')
parser.add_argument('--target-distance-list')
parser.add_argument('--metric-for-target-distance-list') #, default='aa')  # it would be nice to not set defaults here, since it clutters up the bcr-phylo simulator.py calls, but this insulates us against the defaults in bcr-phylo simulator.py changing at some point
parser.add_argument('--leaf-sampling-scheme-list') #, default='uniform-random')
parser.add_argument('--target-count-list') #, default='1')
parser.add_argument('--n-target-clusters-list')  # NOTE do *not* set a default here, since in bcr-phylo simulator.py the default is None
parser.add_argument('--min-target-distance-list')
parser.add_argument('--context-depend-list')
parser.add_argument('--multifurcating-tree-list')
parser.add_argument('--paratope-positions-list')
parser.add_argument('--affinity-measurement-error-list')
parser.add_argument('--metric-method', choices=['shm', 'fay-wu-h', 'cons-dist-nuc', 'cons-dist-aa', 'delta-lbi', 'lbi', 'aa-lbi', 'lbr', 'lbf', 'aa-lbr', 'aa-lbf', 'dtr', 'cons-lbi'], help='method/metric to compare to/correlate with affinity (for use with get-tree-metrics action). If not set, run partis to get lb metrics. UPDATE can now run this way also for plain/nuc lbi and lbr.')
parser.add_argument('--plot-metrics', default='lbi:lbr')  # don't add dtr until it can really run with default options (i.e. model files are really settled)
parser.add_argument('--plot-metric-extra-strs', help='extra strs for each metric in --plot-metrics (i.e. corresponding to what --extra-plotstr was set to during get-tree-metrics for that metric)')
parser.add_argument('--dont-plot-extra-strs', action='store_true', help='while we still use the strings in --plot-metric-extra-strs to find the right dir to get the plot info from, we don\'t actually put the str in the plot (i.e. final plot versions where we don\'t want to see which dtr version it is)')
parser.add_argument('--combo-extra-str', help='extra label for combine-plots action i.e. write to combined-%s/ subdir instead of combined/')
parser.add_argument('--daffy-distr-hist-limit', default='frac:0.03', help='lower bound type:value to use for new distribution performance hists for delta affinity metrics (aa-lbr/lbr). frac:<f> takes the top 100*<f> percent, N:<n> takes the top <n> seqs, val:<v> takes all seqs with metric value above <v>. Default set below.')
parser.add_argument('--abs-affy-distr-hist-limit', default='frac:0.10', help='same as --daffy-distr-hist-limit, but for absolute affinity')
parser.add_argument('--pvks-to-plot', help='only plot these line/legend values when combining plots')
parser.add_argument('--use-val-cfgs', action='store_true', help='use plotting.val_cfgs dict (we can\'t always use it)')
parser.add_argument('--train-dtr', action='store_true')
parser.add_argument('--dtr-path', help='Path from which to read decision tree regression training data. If not set (and --metric-method is dtr), we use a default (see --train-dtr).')
parser.add_argument('--dtr-cfg', help='yaml file with dtr training parameters (read by treeutils.calculate_non_lb_tree_metrics()). If not set, default parameters are taken from treeutils.py')
parser.add_argument('--selection-strength-list') #, default='1.0')
parser.add_argument('--no-scan-parameter-variances', help='Configures parameter variance among families (see bcr-phylo-run.py help for details). Use this version if you only want *one* combination, i.e. if you\'re not scanning across variable combinations: all the different variances go into one bcr-phylo-run.py run (this could be subsumed into the next arg, but for backwards compatibility/cmd line readability it\'s nice to keep it).')
parser.add_argument('--parameter-variances-list', help='Configures parameter variance among families (see bcr-phylo-run.py help for details). Use this version for scanning several combinations. Colons \':\' separate different bcr-phylo-run.py runs, while \'_c_\' separate parameter-variances for multiple variables within each bcr-phylo-run.py run.')
parser.add_argument('--dont-observe-common-ancestors', action='store_true')
parser.add_argument('--simu-extra-args')
parser.add_argument('--zip-vars', help='colon-separated list of variables for which to pair up values sequentially, rather than doing all combinations')
parser.add_argument('--seq-len', default=400, type=int)
parser.add_argument('--n-replicates', default=1, type=int)
parser.add_argument('--iseeds', help='if set, only run these replicate indices (i.e. these corresponds to the increment *above* the random seed)')
parser.add_argument('--n-max-procs', type=int, help='Max number of *child* procs (see --n-sub-procs). Default (None) results in no limit.')
parser.add_argument('--n-sub-procs', type=int, default=1, help='Max number of *grandchild* procs (see --n-max-procs)')
parser.add_argument('--n-max-queries', type=int, help='stop after reading this many queries from whatever is input file for this step (NOTE doesn\'t necessarily work for every action)')
parser.add_argument('--random-seed', default=0, type=int, help='note that if --n-replicates is greater than 1, this is only the random seed of the first replicate')
parser.add_argument('--base-outdir', default='%s/partis/tree-metrics' % os.getenv('fs', default=os.getenv('HOME')))
parser.add_argument('--label', default='test')
parser.add_argument('--extra-plotstr', default='', help='if set, put plots resulting from \'get-tree-metrics\' into a separate subdir using this string, rather than just plots/ (e.g. for plotting with many different dtr versions)')
parser.add_argument('--include-relative-affy-plots', action='store_true')
parser.add_argument('--only-csv-plots', action='store_true', help='only write csv/yaml versions of plots (for future parsing), and not the actual svg files (which is slow)')
parser.add_argument('--no-tree-plots', action='store_true', help='don\'t make any of the tree plots, which are slow')
parser.add_argument('--overwrite', action='store_true')  # not really propagated to everything I think
parser.add_argument('--debug', action='store_true')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--slurm', action='store_true', help='run child procs with slurm')
parser.add_argument('--sub-slurm', action='store_true', help='run grandchild procs with slurm')
parser.add_argument('--workdir')  # default set below
parser.add_argument('--final-plot-xvar', help='variable to put on the x axis of the final comparison plots')  # , default='lb-tau'
parser.add_argument('--legend-var', help='non-default "component" variable (e.g. obs-frac) to use to label different lines in the legend')
parser.add_argument('--x-legend-var', help='derived variable with which to label the x axis (e.g. mfreq [shm percent] when --final-plot-x-var is scratch-mute-freq)')
parser.add_argument('--partis-dir', default=os.getcwd(), help='path to main partis install dir')
parser.add_argument('--make-hist-plots', action='store_true')
# specific to get-lb-bounds:
parser.add_argument('--n-tau-lengths-list', help='set either this or --n-generations-list')
parser.add_argument('--n-generations-list', default='4:5:6:7:8:9:10:12', help='set either this or --n-tau-lengths-list')  # going to 20 uses a ton of memory, not really worth waiting for
parser.add_argument('--max-lb-n-offspring', default=2, type=int, help='multifurcation number for max lb calculation')
parser.add_argument('--only-metrics', default='lbi:lbr', help='which (of lbi, lbr) metrics to do lb bound calculation')
parser.add_argument('--make-plots', action='store_true', help='note: only for get-lb-bounds')
args = parser.parse_args()

args.scan_vars = {'simu' : ['carry-cap', 'n-sim-seqs-per-gen', 'obs-times', 'seed', 'target-distance', 'metric-for-target-distance', 'selection-strength', 'leaf-sampling-scheme', 'target-count', 'n-target-clusters', 'min-target-distance', 'context-depend', 'paratope-positions', 'affinity-measurement-error', 'parameter-variances', 'multifurcating-tree'],}
args.scan_vars['get-tree-metrics'] = args.scan_vars['simu'] + ['lb-tau']
args.str_list_vars = ['n-sim-seqs-per-gen', 'obs-times', 'context-depend', 'n-sim-seqs-per-gen', 'obs-times', 'context-depend']
args.bool_args = ['multifurcating-tree']

sys.path.insert(1, args.partis_dir + '/python')
try:
    import utils
    import treeutils
    import plotting
    import lbplotting
    import scanplot
    from hist import Hist
except ImportError as e:
    print(e)
    raise Exception('couldn\'t import from main partis dir \'%s\' (set with --partis-dir)' % args.partis_dir)

args.actions = utils.get_arg_list(args.actions, choices=all_actions)
args.zip_vars = utils.get_arg_list(args.zip_vars)
# TODO use utils.get_scanvar_arg_lists()?
args.carry_cap_list = utils.get_arg_list(args.carry_cap_list, intify=True, forbid_duplicates=args.zip_vars is None or 'carry-cap' not in args.zip_vars)  # if we're zipping the var, we have to allow duplicates, but then check for them again after we've done combos in get_var_info()
args.n_sim_seqs_per_gen_list = utils.get_arg_list(args.n_sim_seqs_per_gen_list, list_of_lists=True, intify=True, forbid_duplicates=args.zip_vars is None or 'n-sim-seqs-per-gen' not in args.zip_vars)
args.obs_times_list = utils.get_arg_list(args.obs_times_list, list_of_lists=True, intify=True, forbid_duplicates=args.zip_vars is None or 'obs-times' not in args.zip_vars)
if args.lb_tau_list == 'auto':
    assert 'get-lb-bounds' in args.actions
    args.lb_tau_list = '%f' % utils.round_to_n_digits(1. / args.seq_len, 2)  # (1. / args.seq_len) #
    print('  setting --lb-tau-list to 1 / %d = %s' % (args.seq_len, args.lb_tau_list))
args.lb_tau_list = utils.get_arg_list(args.lb_tau_list, floatify=True, forbid_duplicates=True)
args.target_distance_list = utils.get_arg_list(args.target_distance_list, intify=True)
args.metric_for_target_distance_list = utils.get_arg_list(args.metric_for_target_distance_list, forbid_duplicates=True, choices=['aa', 'nuc', 'aa-sim-ascii', 'aa-sim-blosum'])
args.leaf_sampling_scheme_list = utils.get_arg_list(args.leaf_sampling_scheme_list, forbid_duplicates=True, choices=['uniform-random', 'affinity-biased', 'high-affinity'])  # WARNING 'high-affinity' gets called 'perfect' in the legends and 'affinity-biased' gets called 'high affinity'
args.target_count_list = utils.get_arg_list(args.target_count_list, forbid_duplicates=True)
args.n_target_clusters_list = utils.get_arg_list(args.n_target_clusters_list, forbid_duplicates=True)
args.min_target_distance_list = utils.get_arg_list(args.min_target_distance_list, forbid_duplicates=True)
args.context_depend_list = utils.get_arg_list(args.context_depend_list, forbid_duplicates=True)
args.multifurcating_tree_list = utils.get_arg_list(args.multifurcating_tree_list, forbid_duplicates=True)
args.paratope_positions_list = utils.get_arg_list(args.paratope_positions_list, forbid_duplicates=True, choices=['all', 'cdrs'])
args.affinity_measurement_error_list = utils.get_arg_list(args.affinity_measurement_error_list, forbid_duplicates=True, floatify=True)
args.parameter_variances_list = utils.get_arg_list(args.parameter_variances_list, forbid_duplicates=True)
args.plot_metrics = utils.get_arg_list(args.plot_metrics)
args.plot_metric_extra_strs = utils.get_arg_list(args.plot_metric_extra_strs)
if args.plot_metric_extra_strs is None:
    args.plot_metric_extra_strs = ['' for _ in args.plot_metrics]
if len(args.plot_metrics) != len(args.plot_metric_extra_strs):
    raise Exception('--plot-metrics %d not same length as --plot-metric-extra-strs %d' % (len(args.plot_metrics), len(args.plot_metric_extra_strs)))
args.pvks_to_plot = utils.get_arg_list(args.pvks_to_plot)

args.daffy_distr_hist_limit = utils.get_arg_list(args.daffy_distr_hist_limit) #, key_val_pairs=True)
assert args.daffy_distr_hist_limit[0] in ['N', 'frac', 'val']
args.daffy_distr_hist_limit[1] = float(args.daffy_distr_hist_limit[1])

args.abs_affy_distr_hist_limit = utils.get_arg_list(args.abs_affy_distr_hist_limit) #, key_val_pairs=True)
assert args.abs_affy_distr_hist_limit[0] in ['N', 'frac', 'val']
args.abs_affy_distr_hist_limit[1] = float(args.abs_affy_distr_hist_limit[1])

args.selection_strength_list = utils.get_arg_list(args.selection_strength_list, floatify=True, forbid_duplicates=True)
args.n_tau_lengths_list = utils.get_arg_list(args.n_tau_lengths_list, floatify=True)
args.n_generations_list = utils.get_arg_list(args.n_generations_list, intify=True)
args.only_metrics = utils.get_arg_list(args.only_metrics)
args.iseeds = utils.get_arg_list(args.iseeds, intify=True)
if [args.n_tau_lengths_list, args.n_generations_list].count(None) != 1:
    raise Exception('have to set exactly one of --n-tau-lengths, --n-generations')
if args.final_plot_xvar is None:  # set default value based on scan vars
    base_args, varnames, _, valstrs = utils.get_var_info(args, args.scan_vars['simu'])
    svars = [v for v in varnames if v != 'seed']
    args.final_plot_xvar = svars[0] if len(svars) > 0 else 'seed'  # if we're not scanning over any vars, i'm not sure what we should use

import random
random.seed(args.random_seed)  # somehow this is necessary to get the same results, even though I'm not using the module anywhere directly
numpy.random.seed(args.random_seed)

if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work')))

# ----------------------------------------------------------------------------------------
for action in args.actions:
    if action == 'get-lb-bounds':
        calc_lb_bounds(args)
    elif action == 'bcr-phylo':
        run_bcr_phylo(args)
    elif action == 'get-tree-metrics':
        assert args.metric_method is None  # don't use it the old way (requires running multiple times on the command line
        for metric in args.plot_metrics:
            args.metric_method = metric  # NOTE this is hackey, but it's probably better not to try to rewrie the functions to take the metric as an argument, since e.g. get_all_tm_fnames() needs it -- basically lots of code assumes that setting --metric-method means something, so leave it that way)
            get_tree_metrics(args)
        args.metric_method = None
    elif action in ['plot', 'combine-plots'] and not args.dry:
        assert args.extra_plotstr == ''  # only use --extra-plotstr for get-tree-metrics, for this use --plot-metric-extra-strs (because we in general have multiple --plot-metrics when we're here)
        assert args.metric_method is None  # when plotting, you should only be using --plot-metrics
        _, varnames, val_lists, valstrs = utils.get_var_info(args, args.scan_vars['get-tree-metrics'])
        print('plotting %d combinations of %d variable%s (%s) with %d families per combination to %s' % (len(valstrs), len(varnames), utils.plural(len(varnames)), ', '.join(varnames), 1 if args.n_sim_events_per_proc is None else args.n_sim_events_per_proc, scanplot.get_comparison_plotdir(args, None)))
        pchoice = 'per-seq'
        if action == 'plot':
            procs = []
            for metric, estr in zip(args.plot_metrics, args.plot_metric_extra_strs):
                utils.prep_dir(scanplot.get_comparison_plotdir(args, metric, extra_str=estr), subdirs=[pchoice], wildlings=['*.html', '*.svg', '*.yaml'])
                cfg_list = lbplotting.single_lbma_cfg_vars(metric)
                cfg_list = lbplotting.add_lbma_cfg_permutations(cfg_list, include_relative_affy_plots=args.include_relative_affy_plots, make_distr_hists=True)
                for ptvar, ptlabel, use_relative_affy, distr_hists in cfg_list:
                    print('  %s  %-s %-13s%-s' % (utils.color('blue', metric), utils.color('purple', estr, width=20, padside='right') if estr != '' else 20*' ', ptvar, utils.color('green', '(relative)') if use_relative_affy else ''))
                    for cgroup in treeutils.cgroups:
                        print('    %-12s  %15s  %s' % (pchoice, cgroup, ptvar))
                        arglist, kwargs = (args, args.scan_vars['get-tree-metrics'], action, metric, ptvar, args.final_plot_xvar), {'ptilelabel' : ptlabel, 'fnfcn' : get_tm_fname, 'per_x' : pchoice, 'choice_grouping' : cgroup, 'use_relative_affy' : use_relative_affy, 'distr_hists' : distr_hists, 'metric_extra_str' : estr, 'debug' : args.debug}
                        if args.test:
                            scanplot.make_plots(*arglist, **kwargs)
                        else:
                            procs.append(multiprocessing.Process(target=scanplot.make_plots, args=arglist, kwargs=kwargs))
            if not args.test:
                utils.run_proc_functions(procs)
            for metric, estr in zip(args.plot_metrics, args.plot_metric_extra_strs):
                plotting.make_html(scanplot.get_comparison_plotdir(args, metric, per_x=pchoice, extra_str=estr), n_columns=2)
        elif action == 'combine-plots':
            utils.prep_dir(scanplot.get_comparison_plotdir(args, 'combined'), subdirs=[pchoice], wildlings=['*.html', '*.svg'])
            cfg_list = set([ppair for mtmp in args.plot_metrics for ppair in lbplotting.single_lbma_cfg_vars(mtmp)])  # I don't think we care about the order
            cfg_list = lbplotting.add_lbma_cfg_permutations(cfg_list, include_relative_affy_plots=args.include_relative_affy_plots, make_distr_hists=True)
            iplot = 0
            for ptvar, ptlabel, use_relative_affy, distr_hists in cfg_list:
                print(ptvar)
                for cgroup in treeutils.cgroups:
                    print('  ', cgroup)
                    scanplot.make_plots(args, args.scan_vars['get-tree-metrics'], action, None, ptvar, args.final_plot_xvar, ptilelabel=ptlabel, per_x=pchoice, choice_grouping=cgroup, use_relative_affy=use_relative_affy, distr_hists=distr_hists, make_legend=iplot==0)
                    iplot += 1
            plotting.make_html(scanplot.get_comparison_plotdir(args, 'combined', per_x=pchoice), n_columns=2)
        else:
            assert False
