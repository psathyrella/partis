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
                    cmdfos = [plotting.get_lb_tree_cmd(lbvals[metric][btype]['vals']['tree'], '%s/%s-%s-tree.svg' % (plotdir, metric, btype), metric, 'affinities', args.ete_path, args.workdir, metafo=lbvals[metric][btype]['vals'], tree_style='circular')]
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
        if vn not in args.scan_vars[svtype]:  # e.g. lb tau, which is only for partitioning
            continue
        outdir.append('%s-%s' % (vn, vstr))
    return '/'.join(outdir)

# ----------------------------------------------------------------------------------------
def get_bcr_phylo_outdir(varnames, vstr):
    return get_outdir(varnames, vstr, 'simu') + '/bcr-phylo'

# ----------------------------------------------------------------------------------------
def get_bcr_phylo_outfname(varnames, vstr):
    return '%s/selection/simu/mutated-simu.yaml' % get_bcr_phylo_outdir(varnames, vstr)

# ----------------------------------------------------------------------------------------
def get_parameter_dir(varnames, vstr):
    return '%s/selection/partis/params' % get_bcr_phylo_outdir(varnames, vstr)

# ----------------------------------------------------------------------------------------
def get_partition_outdir(varnames, vstr):
    return get_outdir(varnames, vstr, 'partition') + '/partis'

# ----------------------------------------------------------------------------------------
def get_partition_fname(varnames, vstr):
    return '%s/partitions.yaml' % get_partition_outdir(varnames, vstr)

# ----------------------------------------------------------------------------------------
def get_partition_plotdir(varnames, vstr):
    return get_partition_outdir(varnames, vstr) + '/plots'

# ----------------------------------------------------------------------------------------
def get_plotdir():
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
def make_plots(args, metric, use_relative_affy=True, min_ptile_to_plot=75.):  # have to go lower than 85. for small sample sizes
    _, varnames, val_lists, valstrs = get_var_info(args, args.scan_vars['partition'])
    plotvals = collections.OrderedDict()
    print '  plotting %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    missing_vstrs = []
    for vlists, vstrs in zip(val_lists, valstrs):
        n_per_gen = vlists[varnames.index('n-sim-seqs-per-gen')]
        assert len(n_per_gen) == 1
        n_per_gen = n_per_gen[0]
        assert len(args.obs_times_list) == 1
        n_obs = n_per_gen * len(args.obs_times_list[0])
        assert len(args.carry_cap_list) == 1
        n_tot = args.carry_cap_list[0]
        obs_frac = n_obs / float(n_tot)

        yfname = '%s/true-tree-metrics/%s-vs-%s-true-tree-ptiles%s.yaml' % (get_partition_plotdir(varnames, vstrs), metric,
                                                                            'affinity' if metric == 'lbi' else 'n-ancestors',
                                                                            '-relative' if use_relative_affy and metric == 'lbi' else '')
        try:
            with open(yfname) as yfile:
                info = json.load(yfile)  # too slow with yaml
        except IOError:  # os.path.exists() is too slow with this many files
            missing_vstrs.append(vstrs)
            continue
        # the perfect line is higher for lbi, but lower for lbr, hence the abs(). Occasional values can go past/better than perfect, so maybe it would make sense to reverse sign for lbi/lbr rather than taking abs(), but I think this is better
        yval_key = 'mean_%s_ptiles' % ('affy' if metric == 'lbi' else 'n_ancestor')
        diff_to_perfect = numpy.mean([abs(pafp - afp) for lbp, afp, pafp in zip(info['lb_ptiles'], info[yval_key], info['perfect_vals']) if lbp > min_ptile_to_plot])
        if math.isnan(diff_to_perfect):
            raise Exception('  empty value list above min ptile to plot of %f -- probably just reduce this value' % min_ptile_to_plot)
            # continue  # can't just continue, it crashes further down, which I think means things have to be the same length and you can't just skip

        tau = vlists[varnames.index('lb-tau')]
        if args.n_replicates > 1:  # need to average over the replicates
            if obs_frac not in plotvals:
                plotvals[obs_frac] = {i : [] for i in getsargval('seed')}
            plotvals[obs_frac][vlists[varnames.index('seed')]].append((tau, diff_to_perfect))
        else:
            if obs_frac not in plotvals:
                plotvals[obs_frac] = []
            plotvals[obs_frac].append((tau, diff_to_perfect))

    if len(missing_vstrs) > 0:
        print '  missing:  %s' % '   '.join(varnames)
        for vstrs in missing_vstrs:
            print '      %s' % '  '.join('%7s' % v for v in vstrs)
        sys.exit()

    if args.n_replicates > 1:  # need to average over the replicates
        for obs_frac, ofvals in plotvals.items():
            mean_vals = []
            for ipair in range(len(ofvals[0])):  # note that 0 is a dict key, not an index
                tau = [ofvals[i][ipair][0] for i in ofvals]  # ick, now I wish I hadn't done it as a 2-tuple
                assert len(set(tau)) == 1  # all of 'em better have the same tau
                tau = tau[0]
                mean_vals.append((tau, numpy.mean([ofvals[i][ipair][1] for i in ofvals])))
            plotvals[obs_frac] = mean_vals

    fig, ax = plotting.mpl_init()
    for obs_frac in plotvals:
        lb_taus, diffs_to_perfect = zip(*plotvals[obs_frac])
        ax.plot(lb_taus, diffs_to_perfect, label='%.4f' % obs_frac, alpha=0.7, linewidth=4)
    ax.plot([1./args.seq_len, 1./args.seq_len], ax.get_ylim(), linewidth=3, alpha=0.7, color='darkred', linestyle='--') #, label='1/seq len')
    plotting.mpl_finish(ax, get_plotdir(), '%s-affy-ptiles-obs-frac-vs-lb-tau' % metric, xlabel='tau', ylabel='mean %s to perfect\nfor %s ptiles in [%.0f, 100]' % ('percentile' if metric == 'lbi' else 'N ancestors', metric.upper(), min_ptile_to_plot),
                        title=metric.upper(), leg_title='fraction sampled', leg_prop={'size' : 12}, leg_loc=(0.04, 0.67),
                        xticks=lb_taus, xticklabels=[str(t) for t in lb_taus], xticklabelsize=16,
    )

# ----------------------------------------------------------------------------------------
def run_bcr_phylo(args):  # also caches parameters
    base_args, varnames, _, valstrs = get_var_info(args, args.scan_vars['simu'])
    cmdfos = []
    print '  running %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        print '   %s' % ' '.join(vstrs)
        outdir = get_bcr_phylo_outdir(varnames, vstrs)
        if utils.output_exists(args, get_parameter_dir(varnames, vstrs) + '/hmm/hmms', offset=8):
            n_already_there += 1
            continue
        cmd = './bin/bcr-phylo-run.py --actions simu:cache-parameters --base-outdir %s %s' % (outdir, ' '.join(base_args))
        for vname, vstr in zip(varnames, vstrs):
            cmd += ' --%s %s' % (vname, vstr)
        if args.overwrite:
            cmd += ' --overwrite'
        # cmd += ' --debug 1'
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : get_bcr_phylo_outfname(varnames, vstrs),
            'logdir' : outdir,
            'workdir' : '%s/bcr-phylo-work/%d' % (args.workdir, icombo),
        }]
        print '     %s %s' % (utils.color('red', 'run'), cmd)
    if n_already_there > 0:
        print '      %d / %d skipped (outputs exist, e.g. %s)' % (n_already_there, len(valstrs), get_parameter_dir(varnames, vstrs))
    if len(cmdfos) > 0:
        print '      starting %d jobs' % len(cmdfos)
    utils.run_cmds(cmdfos, debug='write:bcr-phylo.log', batch_system='slurm' if args.slurm else None, n_max_procs=args.n_max_procs, proc_limit_str='bin/bcr-phylo-run')

# ----------------------------------------------------------------------------------------
def partition(args):
    _, varnames, _, valstrs = get_var_info(args, args.scan_vars['partition'])  # can't use base_args a.t.m. since it has the simulation/bcr-phylo args in it
    cmdfos = []
    print '  running %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    n_already_there = 0
    for icombo, vstrs in enumerate(valstrs):
        if args.debug:
            print '   %s' % ' '.join(vstrs)

        partition_fname = get_partition_fname(varnames, vstrs)
        if args.action == 'partition' and utils.output_exists(args, partition_fname, outlabel='partition', offset=8, debug=args.debug):
            n_already_there += 1
            continue

        cmd = './bin/partis %s --is-simu --infname %s --plotdir %s --outfname %s' % (args.action, get_bcr_phylo_outfname(varnames, vstrs), get_partition_plotdir(varnames, vstrs), partition_fname)
        if args.action == 'partition':
            cmd += ' --parameter-dir %s --n-final-clusters 1 --get-tree-metrics' % (get_parameter_dir(varnames, vstrs))
        cmd += ' --lb-tau %s' % vstrs[varnames.index('lb-tau')]
        cmd += ' --seed %s' % args.random_seed  # NOTE second/commented version this is actually wrong: vstrs[varnames.index('seed')]  # there isn't actually a reason for different seeds here (we want the different seeds when running bcr-phylo), but oh well, maybe it's a little clearer this way
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : partition_fname,
            'workdir' : get_partition_outdir(varnames, vstrs),
        }]
    if n_already_there > 0:
        print '      %d / %d skipped (outputs exist, e.g. %s)' % (n_already_there, len(valstrs), partition_fname)
    if len(cmdfos) > 0:
        print '      starting %d jobs' % len(cmdfos)
        utils.run_cmds(cmdfos, debug='write:partition.log', batch_system='slurm' if args.slurm else None, n_max_procs=args.n_max_procs, proc_limit_str='bin/partis')

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('action', choices=['get-lb-bounds', 'run-bcr-phylo', 'partition', 'get-tree-metrics', 'plot'])
parser.add_argument('--carry-cap-list', default='1000')
parser.add_argument('--n-sim-seqs-per-gen-list', default='30:50:75:100:150:200', help='colon-separated list of comma-separated lists of the number of sequences for bcr-phylo to sample at the times specified by --obs-times-list')
parser.add_argument('--obs-times-list', default='125,150', help='colon-separated list of comma-separated lists of bcr-phylo observation times')
parser.add_argument('--lb-tau-list', default='0.0005:0.001:0.002:0.0025:0.003:0.004:0.005:0.008:0.012')
parser.add_argument('--n-tau-lengths-list', help='set either this or --n-generations-list')
parser.add_argument('--n-generations-list', default='4:5:6:7:8:9:10:12', help='set either this or --n-tau-lengths-list')  # going to 20 uses a ton of memory, not really worth waiting for
parser.add_argument('--max-lb-n-offspring', default=2, type=int, help='multifurcation number for max lb calculation')
parser.add_argument('--seq-len', default=400, type=int)
parser.add_argument('--n-replicates', default=1, type=int)
parser.add_argument('--iseed', type=int, help='if set, only run this replicate index')
parser.add_argument('--n-max-procs', type=int)
parser.add_argument('--only-metrics', default='lbi:lbr', help='which (of lbi, lbr) metrics to do lb bound calculation')
parser.add_argument('--random-seed', default=0, type=int, help='note that if --n-replicates is greater than 1, this is only the random seed of the first replicate')
parser.add_argument('--base-outdir', default='%s/partis/tree-metrics' % os.getenv('fs', default=os.getenv('HOME')))
parser.add_argument('--label', default='test')
parser.add_argument('--make-plots', action='store_true')
parser.add_argument('--overwrite', action='store_true')  # not really propagated to everything I think
parser.add_argument('--debug', action='store_true')
parser.add_argument('--slurm', action='store_true')
parser.add_argument('--workdir')  # default set below
parser.add_argument('--partis-dir', default=os.getcwd(), help='path to main partis install dir')
parser.add_argument('--ete-path', default=('/home/%s/anaconda_ete/bin' % os.getenv('USER')) if os.getenv('USER') is not None else None)
args = parser.parse_args()

args.scan_vars = {
    'simu' : ['carry-cap', 'n-sim-seqs-per-gen', 'obs-times', 'seed'],
    'partition' : ['carry-cap', 'n-sim-seqs-per-gen', 'obs-times', 'seed', 'lb-tau'],
}

sys.path.insert(1, args.partis_dir + '/python')
try:
    import utils
    import treeutils
    import plotting
except ImportError as e:
    print e
    raise Exception('couldn\'t import from main partis dir \'%s\' (set with --partis-dir)' % args.partis_dir)

args.carry_cap_list = utils.get_arg_list(args.carry_cap_list, intify=True)
args.n_sim_seqs_per_gen_list = utils.get_arg_list(args.n_sim_seqs_per_gen_list, list_of_lists=True, intify=True)
args.obs_times_list = utils.get_arg_list(args.obs_times_list, list_of_lists=True, intify=True)
args.lb_tau_list = utils.get_arg_list(args.lb_tau_list, floatify=True)
args.n_tau_lengths_list = utils.get_arg_list(args.n_tau_lengths_list, floatify=True)
args.n_generations_list = utils.get_arg_list(args.n_generations_list, intify=True)
args.only_metrics = utils.get_arg_list(args.only_metrics)
if [args.n_tau_lengths_list, args.n_generations_list].count(None) != 1:
    raise Exception('have to set exactly one of --n-tau-lengths, --n-generations')

if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work')))

# ----------------------------------------------------------------------------------------
if args.action == 'get-lb-bounds':
    calc_lb_bounds(args)
elif args.action == 'run-bcr-phylo':
    run_bcr_phylo(args)
elif args.action in ['partition', 'get-tree-metrics']:
    partition(args)
elif args.action == 'plot':
    procs = [multiprocessing.Process(target=make_plots, args=(args, metric))  # time is almost entirely due to file open + json.load
             for metric in treeutils.lb_metrics]
    utils.run_proc_functions(procs)
    # for metric in treeutils.lb_metrics:
    #     make_plots(args, metric)
