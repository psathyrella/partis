#!/usr/bin/env python
import argparse
import operator
import os
import sys
import yaml
import json
import colored_traceback.always
import collections
import numpy
import math
import subprocess
import multiprocessing

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
            n_gen_list, lb_list = zip(*sorted(parsed_info[metric][btype][lbt].items(), key=operator.itemgetter(0)))
            if lbt == 1. / args.seq_len:  # add a horizontal line corresponding to the asymptote for tau = 1/seq_len
                ax.plot(ax.get_xlim(), (lb_list[-1], lb_list[-1]), linewidth=3, alpha=0.7, color='darkred', linestyle='--') #, label='1/seq len')
            ax.plot(n_gen_list, lb_list, label='%.4f' % lbt, alpha=0.7, linewidth=4)
            if print_results and log == '':
                print '      %7.4f    %6.4f' % (lbt, lb_list[-1])
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
        print '%s:     tau    %s %s' % (btype, btype, metric)

    for log in ['', 'y']:
        make_plot(log, parsed_info)

# ----------------------------------------------------------------------------------------
def calc_lb_bounds(args, n_max_gen_to_plot=4, lbt_bounds=(0.001, 0.005), print_results=False):
    print_results = True
    btypes = ['min', 'max']

    outdir = '%s/lb-tau-normalization/%s' % (args.base_outdir, args.label)

    parsed_info = {m : {b : {} for b in btypes} for m in args.only_metrics}
    for lbt in args.lb_tau_list:
        if lbt < lbt_bounds[0] or lbt > lbt_bounds[1]:
            print '    skipping tau outside of bounds for bound plotting'
            continue

        gen_list = args.n_generations_list
        if gen_list is None:
            gen_list = [get_n_generations(ntl, lbt) for ntl in args.n_tau_lengths_list]
        if args.lb_tau_list.index(lbt) == 0 or args.n_tau_lengths_list is not None:  # if --n-tau-lengths-list is set, they could be different for each tau
            print '      n gen: %s' % ' '.join(str(n) for n in gen_list)
        print '   lb tau %.4f' % lbt
        for n_gen in gen_list:

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

            print '         running %d' % n_gen

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
                    cmdfos = [lbplotting.get_lb_tree_cmd(lbvals[metric][btype]['vals']['tree'], '%s/%s-%s-tree.svg' % (plotdir, metric, btype), metric, 'affinities', args.ete_path, args.workdir, metafo=lbvals[metric][btype]['vals'], tree_style='circular')]
                    utils.run_cmds(cmdfos, clean_on_success=True, shell=True, debug='print')

    if args.make_plots:
        for metric in args.only_metrics:
            for btype in btypes:
                if metric == 'lbr' and btype == 'min':  # it's just zero, and confuses the log plots
                    continue
                if len(parsed_info[metric][btype]) == 0:
                    print 'nothing to do (<parsed_info> empty)'
                    continue
                make_lb_bound_plots(args, outdir, metric, btype, parsed_info, print_results=print_results)

# ----------------------------------------------------------------------------------------
def get_outdir(varnames, vstr, svtype):
    assert len(varnames) == len(vstr)
    outdir = [args.base_outdir, args.label]
    for vn, vstr in zip(varnames, vstr):
        if vn not in args.scan_vars[svtype]:  # e.g. lb tau, which is only for lb calculation
            continue
        outdir.append('%s-%s' % (vn, vstr))
    return '/'.join(outdir)

# ----------------------------------------------------------------------------------------
def get_bcr_phylo_outdir(varnames, vstr):
    return get_outdir(varnames, vstr, 'simu') + '/bcr-phylo'

# ----------------------------------------------------------------------------------------
def get_simfname(varnames, vstr):
    return '%s/selection/simu/mutated-simu.yaml' % get_bcr_phylo_outdir(varnames, vstr)

# ----------------------------------------------------------------------------------------
def get_parameter_dir(varnames, vstr):
    return '%s/selection/partis/params' % get_bcr_phylo_outdir(varnames, vstr)

# ----------------------------------------------------------------------------------------
def get_tree_metric_outdir(varnames, vstr):
    return get_outdir(varnames, vstr, 'get-tree-metrics') + '/partis'

# ----------------------------------------------------------------------------------------
def get_partition_fname(varnames, vstr, action):  # if action is 'bcr-phylo', we want the original partition output file, but if it's 'get-tree-metrics', we want the copied one, that had tree metrics added to it (and which is in the e.g. tau subdir)
    outdir = '%s/selection/partis' % (get_bcr_phylo_outdir(varnames, vstr)) if action == 'bcr-phylo' else get_tree_metric_outdir(varnames, vstr)
    return '%s/partition.yaml' % outdir

# ----------------------------------------------------------------------------------------
def get_tree_metric_plotdir(varnames, vstr):
    return get_tree_metric_outdir(varnames, vstr) + '/plots'

# ----------------------------------------------------------------------------------------
def get_tree_metric_fname(varnames, vstr, metric='lbi', x_axis_label='affinity', use_relative_affy=True):  # note that there are separate svg files for each iclust, but info for all clusters are written to the same yaml file (which is all we need here)
    # note that we set all the defaults just so when we're checking for output before running, we can easily get one of the file names
    old_path = '%s/true-tree-metrics/%s-vs-%s-true-tree-ptiles%s.yaml' % (get_tree_metric_plotdir(varnames, vstr), metric, x_axis_label, '-relative' if use_relative_affy and metric == 'lbi' else '')
    if os.path.exists(old_path):
        print 'exists', old_path
        return old_path
    vs_str = '%s-vs%s-%s' % (metric, '-relative' if use_relative_affy and metric == 'lbi' else '', x_axis_label)
    return '%s/true-tree-metrics/%s/%s-ptiles/%s-true-tree-ptiles-all-clusters.yaml' % (get_tree_metric_plotdir(varnames, vstr), metric, vs_str, vs_str)

# ----------------------------------------------------------------------------------------
def get_comparison_plotdir():
    return '%s/%s/plots' % (args.base_outdir, args.label)

# ----------------------------------------------------------------------------------------
def getsargval(sv):  # ick this name sucks
    def dkey(sv):
        return sv.replace('-', '_') + '_list'
    if sv == 'seed':
        riter = range(args.n_replicates) if args.iseed is None else [args.iseed]
        return [args.random_seed + i for i in riter]
    else:
        return args.__dict__[dkey(sv)]

# ----------------------------------------------------------------------------------------
def get_vlval(vlists, varnames, vname):  # ok this name also sucks, but they're doing complicated things while also needing really short names...
    if vname in varnames:
        return vlists[varnames.index(vname)]
    else:
        assert len(getsargval(vname))  # um, I think?
        return getsargval(vname)[0]

# ----------------------------------------------------------------------------------------
def get_var_info(args, scan_vars):
    def handle_var(svar, val_lists, valstrs):
        convert_fcn = str if svar in ['carry-cap', 'seed', 'lb-tau'] else lambda vlist: ':'.join(str(v) for v in vlist)
        sargv = getsargval(svar)
        if len(sargv) > 1 or (svar == 'seed' and args.iseed is not None):  # if --iseed is set, then we know there must be more than one replicate, but/and we also know the fcn will only be returning one of 'em
            varnames.append(svar)
            val_lists = [vlist + [sv] for vlist in val_lists for sv in sargv]
            valstrs = [vlist + [convert_fcn(sv)] for vlist in valstrs for sv in sargv]
        else:
            base_args.append('--%s %s' % (svar, convert_fcn(sargv[0])))
        return val_lists, valstrs

    base_args = []
    varnames = []
    val_lists, valstrs = [[]], [[]]
    for svar in scan_vars:
        val_lists, valstrs = handle_var(svar, val_lists, valstrs)

    return base_args, varnames, val_lists, valstrs

# ----------------------------------------------------------------------------------------
def make_plots(args, metric, ptilestr, ptilelabel, xvar='lb-tau', min_ptile_to_plot=75., debug=False):
    vlabels = {
        'obs_frac' : 'fraction sampled',
        'n-sim-seqs-per-gen' : 'N/gen',
        'obs-times' : 't obs',
    }
    # ----------------------------------------------------------------------------------------
    def get_obs_frac(vlists, varnames):
        obs_times = get_vlval(vlists, varnames, 'obs-times')
        n_per_gen_vals = get_vlval(vlists, varnames, 'n-sim-seqs-per-gen')
        if len(obs_times) == len(n_per_gen_vals):  # note that this duplicates logic in bcr-phylo simulator.py
            n_sampled = sum(n_per_gen_vals)
        elif len(n_per_gen_vals) == 1:
            n_sampled = len(obs_times) * n_per_gen_vals[0]
        else:
            assert False
        n_total = get_vlval(vlists, varnames, 'carry-cap')  # note that this is of course the number alive at a given time, and very different from the total number that ever lived
        obs_frac = n_sampled / float(n_total)
        dbgstr = '    %-12s %-12s   %-5d     %8s / %-4d = %.3f' % (' '.join(str(o) for o in obs_times), ' '.join(str(n) for n in n_per_gen_vals), n_total,
                                                                   ('(%s)' % ' + '.join(str(n) for n in n_per_gen_vals)) if len(obs_times) == len(n_per_gen_vals) else ('%d * %d' % (len(obs_times), n_per_gen_vals[0])),
                                                                   n_total, n_sampled / float(n_total))
        return obs_frac, dbgstr
    # ----------------------------------------------------------------------------------------
    def pvkeystr(vlists, varnames, obs_frac):
        def valstr(vname):
            vval = obs_frac if vname == 'obs_frac' else get_vlval(vlists, varnames, vname)
            if vname == 'obs_frac':
                return '%.4f' % obs_frac
            else:
                def strfcn(x):
                    return str(x)  # TODO
                if isinstance(vval, list):
                    return ', '.join(strfcn(v) for v in vval)
                else:
                    return strfcn(vval)
        pvnames = sorted(set(varnames) - set(['seed', xvar]))
        if pvnames == ['n-sim-seqs-per-gen']:  # if this is the only thing that's different between different runs (except for the x variable and seed/replicate) then we want to use obs_frac
            pvnames = ['obs_frac']
        pvkey = ', '.join(valstr(vn) for vn in pvnames)  # key identifying each line of a different color
        pvlabel = ', '.join(vlabels.get(vn, vn) for vn in pvnames)
        return pvkey, pvlabel
    # ----------------------------------------------------------------------------------------
    def get_diff_vals(yamlfo, iclust=None):
        ytmpfo = yamlfo
        if iclust is not None:
            if 'iclust-%d' % iclust not in yamlfo:
                print '    %s requested per-cluster ptile vals, but they\'re not in the yaml file (probably just an old file)' % utils.color('yellow', 'warning')
            ytmpfo = yamlfo['iclust-%d' % iclust]
        return [abs(pafp - afp) for lbp, afp, pafp in zip(ytmpfo['lb_ptiles'], ytmpfo[yval_key], ytmpfo['perfect_vals']) if lbp > min_ptile_to_plot]
    # ----------------------------------------------------------------------------------------
    def add_plot_vals(yamlfo, vlists, varnames, obs_frac, iclust=None):
        diff_vals = get_diff_vals(yamlfo, iclust=iclust)
        if len(diff_vals) == 0:
            missing_vstrs['empty'].append(vstrs)  # empty may be from empty list in yaml file, or may be from none of them being above <min_ptile_to_plot>
            return
        diff_to_perfect = numpy.mean(diff_vals)
        tau = get_vlval(vlists, varnames, xvar)
        pvkey, pvlabel = pvkeystr(vlists, varnames, obs_frac)  # key identifying each line in the plot, each with a different color, (it's kind of ugly to get the label here but not use it til we plot, but oh well)
        if args.n_replicates > 1:  # need to average over the replicates/clusters (NOTE I'm not really sure this will work if there's only one replicate but more than one event per proc)
            if args.n_sim_events_per_proc is not None:
                if pvkey not in plotvals:
                    plotvals[pvkey] = {('%d-%d' % (i, j)) : [] for i in getsargval('seed') for j in range(args.n_sim_events_per_proc)}
                ikey = '%d-%d' % (vlists[varnames.index('seed')], iclust)
                plotvals[pvkey][ikey].append((tau, diff_to_perfect))  # TODO arg, is this key wrong?
            else:
                if pvkey not in plotvals:
                    plotvals[pvkey] = {i : [] for i in getsargval('seed')}
                plotvals[pvkey][vlists[varnames.index('seed')]].append((tau, diff_to_perfect))  # TODO arg, is this key wrong?
        else:
            if pvkey not in plotvals:
                plotvals[pvkey] = []
            plotvals[pvkey].append((tau, diff_to_perfect))
    # ----------------------------------------------------------------------------------------
    def get_varname_str():
        return ''.join('%10s' % vlabels.get(v, v) for v in varnames)
    def get_varval_str(vstrs):
        return ''.join('%10s' % v for v in vstrs)

    # ----------------------------------------------------------------------------------------
    debug = True
    _, varnames, val_lists, valstrs = get_var_info(args, args.scan_vars['get-tree-metrics'])
    plotvals, errvals = collections.OrderedDict(), collections.OrderedDict()
    print '  plotting %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    if debug:
        print '%s   | obs times    N/gen        carry cap       fraction sampled' % get_varname_str()
    missing_vstrs = {'missing' : [], 'empty' : []}
    pvlabel = '?'  # arg ick ugh
    for vlists, vstrs in zip(val_lists, valstrs):
        obs_frac, dbgstr = get_obs_frac(vlists, varnames)
        if debug:
            print '%s   | %s' % (get_varval_str(vstrs), dbgstr)
        yfname = get_tree_metric_fname(varnames, vstrs, metric, ptilestr)  # why is this called vstrs rather than vstr?
        try:
            with open(yfname) as yfile:
                yamlfo = json.load(yfile)  # too slow with yaml
        except IOError:  # os.path.exists() is too slow with this many files
            missing_vstrs['missing'].append(vstrs)
            continue
        # the perfect line is higher for lbi, but lower for lbr, hence the abs(). Occasional values can go past/better than perfect, so maybe it would make sense to reverse sign for lbi/lbr rather than taking abs(), but I think this is better
        yval_key = 'mean_%s_ptiles' % ('affy' if ptilestr == 'affinity' else ptilestr)  # arg, would've been nice if that was different
        if args.n_sim_events_per_proc is not None:
            iclusts_in_file = sorted([int(k.split('-')[1]) for k in yamlfo if 'iclust-' in k])  # if there's info for each cluster, it's in sub-dicts under 'iclust-N' (older files won't have it)
            missing_iclusts = [i for i in range(args.n_sim_events_per_proc) if i not in iclusts_in_file]
            if len(missing_iclusts) > 0:
                print '  %s missing %d iclusts (i = %s) from file' % (utils.color('red', 'error'), len(missing_iclusts), ' '.join(str(i) for i in missing_iclusts))
            for iclust in iclusts_in_file:
                add_plot_vals(yamlfo, vlists, varnames, obs_frac, iclust=iclust)
        else:
            add_plot_vals(yamlfo, vlists, varnames, obs_frac)

    # print info about missing and empty results
    for mkey, vstrs_list in missing_vstrs.items():
        if len(vstrs_list) == 0:
            continue
        print '  %s:' % mkey
        print '     %s' % get_varname_str()
        for vstrs in vstrs_list:
            print '      %s' % get_varval_str(vstrs)

    # average over the replicates/clusters
    if (args.n_replicates > 1 or args.n_sim_events_per_proc is not None) and len(plotvals) > 0:
        if debug:
            print '  averaging over %d replicates' % args.n_replicates,
            if args.n_sim_events_per_proc is not None:
                print 'times %d clusters per proc:' % args.n_sim_events_per_proc,
            print ''
            tmplen = str(max(len(pvkey) for pvkey in plotvals))
            print ('    %'+tmplen+'s   N used  N expected') % 'pvkey'
        for pvkey, ofvals in plotvals.items():
            mean_vals, err_vals = [], []
            ofvals = {i : vals for i, vals in ofvals.items() if len(vals) > 0}  # remove zero-length ones (which should correspond to 'missing')
            assert len(set([len(ofvals[i]) for i in ofvals])) == 1
            n_used = []  # just for dbg
            for ipair in range(len(ofvals.values()[0])):  # NOTE if this first one is ever empty this will probably break
                tau = [ofvals[i][ipair][0] for i in ofvals]  # ick, now I wish I hadn't done it as a 2-tuple
                assert len(set(tau)) == 1  # all of 'em better have the same tau
                tau = tau[0]
                ltmp = [ofvals[i][ipair][1] for i in ofvals]
                mean_vals.append((tau, numpy.mean(ltmp)))
                err_vals.append((tau, numpy.std(ltmp, ddof=1) / math.sqrt(len(ltmp))))
                n_used.append(len([ofvals[i][ipair][1] for i in ofvals]))
            plotvals[pvkey] = mean_vals
            errvals[pvkey] = err_vals
            if debug:
                n_expected = args.n_replicates
                if args.n_sim_events_per_proc is not None:
                    n_expected *= args.n_sim_events_per_proc
                print ('    %'+tmplen+'s    %s   %4d%s') % (pvkey, ('%4d' % n_used[0]) if len(set(n_used)) == 1 else utils.color('red', ' '.join(str(n) for n in set(n_used))), n_expected, '' if n_used[0] == n_expected else utils.color('red', ' <--'))

    fig, ax = plotting.mpl_init()
    lb_taus, xticklabels = None, None
    for pvkey in plotvals:
        lb_taus, diffs_to_perfect = zip(*plotvals[pvkey])
        xticklabels = [str(t) for t in lb_taus]
        markersize = None if len(lb_taus) > 1 else 15
        if pvkey in errvals:
            _, yerrs = zip(*errvals[pvkey])  # first item would just be the same as <lb_taus>
            ax.errorbar(lb_taus, diffs_to_perfect, yerr=yerrs, alpha=0.7, markersize=markersize, linewidth=2, marker='.')  #, title='position ' + str(position))
        else:
            ax.plot(lb_taus, diffs_to_perfect, label=pvkey, alpha=0.7, linewidth=4)
    ax.plot([1./args.seq_len, 1./args.seq_len], ax.get_ylim(), linewidth=3, alpha=0.7, color='darkred', linestyle='--') #, label='1/seq len')
    plotting.mpl_finish(ax, get_comparison_plotdir(),
                        '%s-%s-ptiles-obs-frac-vs-lb-tau' % (ptilestr, metric),
                        xlabel=xvar.replace('-', ' '),
                        ylabel='mean %s to perfect\nfor %s ptiles in [%.0f, 100]' % ('percentile' if ptilelabel == 'affinity' else ptilelabel, metric.upper(), min_ptile_to_plot),
                        title=metric.upper(), leg_title=pvlabel, leg_prop={'size' : 12}, leg_loc=(0.04, 0.67),
                        xticks=lb_taus, xticklabels=xticklabels, xticklabelsize=16,
    )

# ----------------------------------------------------------------------------------------
def run_bcr_phylo(args):  # also caches parameters
    base_args, varnames, _, valstrs = get_var_info(args, args.scan_vars['simu'])
    cmdfos = []
    print '  bcr-phylo: running %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        outdir = get_bcr_phylo_outdir(varnames, vstrs)
        if utils.output_exists(args, get_partition_fname(varnames, vstrs, 'bcr-phylo'), offset=8, debug=args.debug):
            n_already_there += 1
            continue
        cmd = './bin/bcr-phylo-run.py --actions simu:cache-parameters:partition --dont-get-tree-metrics --base-outdir %s %s' % (outdir, ' '.join(base_args))
        for vname, vstr in zip(varnames, vstrs):
            cmd += ' --%s %s' % (vname, vstr)
        if args.n_sim_events_per_proc is not None:
            cmd += ' --n-sim-events %d' % args.n_sim_events_per_proc
        if args.overwrite:
            cmd += ' --overwrite'
        if args.only_csv_plots:
            cmd += ' --only-csv-plots'
        # cmd += ' --debug 1'
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : get_partition_fname(varnames, vstrs, 'bcr-phylo'),
            'logdir' : outdir,
            'workdir' : '%s/bcr-phylo-work/%d' % (args.workdir, icombo),
        }]
    if n_already_there > 0:
        print '      %d / %d skipped (outputs exist, e.g. %s)' % (n_already_there, len(valstrs), get_partition_fname(varnames, vstrs, 'bcr-phylo'))
    if len(cmdfos) > 0:
        print '      starting %d jobs' % len(cmdfos)
        utils.run_cmds(cmdfos, debug='write:bcr-phylo.log', batch_system='slurm' if args.slurm else None, n_max_procs=args.n_max_procs, proc_limit_str='bin/bcr-phylo-run')

# ----------------------------------------------------------------------------------------
def get_tree_metrics(args):
    _, varnames, _, valstrs = get_var_info(args, args.scan_vars['get-tree-metrics'])  # can't use base_args a.t.m. since it has the simulation/bcr-phylo args in it
    cmdfos = []
    print '  get-tree-metrics: running %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print '   %s' % ' '.join(vstrs)

        if utils.output_exists(args, get_tree_metric_fname(varnames, vstrs), outlabel='get-tree-metrics', offset=8, debug=args.debug):
            n_already_there += 1
            continue

        if not os.path.isdir(get_tree_metric_outdir(varnames, vstrs)):
            os.makedirs(get_tree_metric_outdir(varnames, vstrs))
        subprocess.check_call(['cp', get_partition_fname(varnames, vstrs, 'bcr-phylo'), get_partition_fname(varnames, vstrs, 'get-tree-metrics')])
        cmd = './bin/partis get-tree-metrics --is-simu --infname %s --plotdir %s --outfname %s' % (get_simfname(varnames, vstrs), get_tree_metric_plotdir(varnames, vstrs), get_partition_fname(varnames, vstrs, 'get-tree-metrics'))
        cmd += ' --lb-tau %s --lbr-tau-factor 1' % get_vlval(vstrs, varnames, 'lb-tau')
        cmd += ' --dont-normalize-lbi'
        cmd += ' --seed %s' % args.random_seed  # NOTE second/commented version this is actually wrong: vstrs[varnames.index('seed')]  # there isn't actually a reason for different seeds here (we want the different seeds when running bcr-phylo), but oh well, maybe it's a little clearer this way
        if args.only_csv_plots:
            cmd += ' --only-csv-plots'
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : get_tree_metric_fname(varnames, vstrs),
            'workdir' : get_tree_metric_outdir(varnames, vstrs),
        }]
    if n_already_there > 0:
        print '      %d / %d skipped (outputs exist, e.g. %s)' % (n_already_there, len(valstrs), get_tree_metric_fname(varnames, vstrs))
    if len(cmdfos) > 0:
        print '      starting %d jobs' % len(cmdfos)
        utils.run_cmds(cmdfos, debug='write:get-tree-metrics.log', batch_system='slurm' if args.slurm else None, n_max_procs=args.n_max_procs, proc_limit_str='bin/partis')

# ----------------------------------------------------------------------------------------
all_actions = ['get-lb-bounds', 'bcr-phylo', 'get-tree-metrics', 'plot']
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--actions', default=':'.join(a for a in all_actions if a != 'get-lb-bounds'))
parser.add_argument('--carry-cap-list', default='1000')
parser.add_argument('--n-sim-seqs-per-gen-list', default='30:50:75:100:150:200', help='colon-separated list of comma-separated lists of the number of sequences for bcr-phylo to sample at the times specified by --obs-times-list')
parser.add_argument('--n-sim-events-per-proc', type=int, help='number of rearrangement events to simulate in each process (default is set in bin/bcr-phylo-run.py)')
parser.add_argument('--obs-times-list', default='125,150', help='colon-separated list of comma-separated lists of bcr-phylo observation times')
parser.add_argument('--lb-tau-list', default='0.0005:0.001:0.002:0.0025:0.003:0.004:0.005:0.008:0.012')
parser.add_argument('--n-tau-lengths-list', help='set either this or --n-generations-list')
parser.add_argument('--n-generations-list', default='4:5:6:7:8:9:10:12', help='set either this or --n-tau-lengths-list')  # going to 20 uses a ton of memory, not really worth waiting for
parser.add_argument('--max-lb-n-offspring', default=2, type=int, help='multifurcation number for max lb calculation')
parser.add_argument('--seq-len', default=400, type=int)
parser.add_argument('--n-replicates', default=1, type=int)
parser.add_argument('--iseed', type=int, help='if set, only run this replicate index (i.e. this corresponds to the increment *above* the random seed)')
parser.add_argument('--n-max-procs', type=int)  # NOTE that with slurm this thinks there's twice as many jobs as there are
parser.add_argument('--only-metrics', default='lbi:lbr', help='which (of lbi, lbr) metrics to do lb bound calculation')
parser.add_argument('--random-seed', default=0, type=int, help='note that if --n-replicates is greater than 1, this is only the random seed of the first replicate')
parser.add_argument('--base-outdir', default='%s/partis/tree-metrics' % os.getenv('fs', default=os.getenv('HOME')))
parser.add_argument('--label', default='test')
parser.add_argument('--make-plots', action='store_true')
parser.add_argument('--only-csv-plots', action='store_true')
parser.add_argument('--overwrite', action='store_true')  # not really propagated to everything I think
parser.add_argument('--debug', action='store_true')
parser.add_argument('--slurm', action='store_true')
parser.add_argument('--workdir')  # default set below
parser.add_argument('--partis-dir', default=os.getcwd(), help='path to main partis install dir')
parser.add_argument('--ete-path', default=('/home/%s/anaconda_ete/bin' % os.getenv('USER')) if os.getenv('USER') is not None else None)
args = parser.parse_args()

args.scan_vars = {
    'simu' : ['carry-cap', 'n-sim-seqs-per-gen', 'obs-times', 'seed'],
    'get-tree-metrics' : ['carry-cap', 'n-sim-seqs-per-gen', 'obs-times', 'seed', 'lb-tau'],
}

sys.path.insert(1, args.partis_dir + '/python')
try:
    import utils
    import treeutils
    import plotting
    import lbplotting
except ImportError as e:
    print e
    raise Exception('couldn\'t import from main partis dir \'%s\' (set with --partis-dir)' % args.partis_dir)

args.actions = utils.get_arg_list(args.actions, choices=all_actions)
args.carry_cap_list = utils.get_arg_list(args.carry_cap_list, intify=True)
args.n_sim_seqs_per_gen_list = utils.get_arg_list(args.n_sim_seqs_per_gen_list, list_of_lists=True, intify=True)
args.obs_times_list = utils.get_arg_list(args.obs_times_list, list_of_lists=True, intify=True)
args.lb_tau_list = utils.get_arg_list(args.lb_tau_list, floatify=True)
args.n_tau_lengths_list = utils.get_arg_list(args.n_tau_lengths_list, floatify=True)
args.n_generations_list = utils.get_arg_list(args.n_generations_list, intify=True)
args.only_metrics = utils.get_arg_list(args.only_metrics)
if [args.n_tau_lengths_list, args.n_generations_list].count(None) != 1:
    raise Exception('have to set exactly one of --n-tau-lengths, --n-generations')
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
        get_tree_metrics(args)
    elif action == 'plot':
        utils.prep_dir(get_comparison_plotdir(), wildlings='*.svg')
        procs = [multiprocessing.Process(target=make_plots, args=(args, metric, ptilestr, ptilelabel))  # time is almost entirely due to file open + json.load
                 for metric, ptilestr, ptilelabel in lbplotting.lb_metric_axis_stuff]
        utils.run_proc_functions(procs)
        # for metric, ptilestr, ptilelabel in lbplotting.lb_metric_axis_stuff:
        #     make_plots(args, metric, ptilestr, ptilelabel)
        #     break
