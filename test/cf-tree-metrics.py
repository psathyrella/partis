#!/usr/bin/env python
import argparse
import operator
import os
import sys
import yaml
import colored_traceback.always
import collections
import numpy
import math

# ----------------------------------------------------------------------------------------
def get_n_generations(ntl, tau):  # NOTE duplicates code in treeutils.get_max_lbi()
    return max(1, int(args.seq_len * tau * ntl))

# ----------------------------------------------------------------------------------------
def get_outfname(outdir):
    return '%s/vals.yaml' % outdir

# ----------------------------------------------------------------------------------------
def calc_max_lbi(args):
    if args.overwrite:
        raise Exception('not implemented')

    outdir = '%s/lb-tau-optimization/%s' % (args.base_outdir, args.label)

    parsed_info = {}
    for lbt in args.lb_tau_list:

        gen_list = args.n_generations_list
        if gen_list is None:
            gen_list = [get_n_generations(ntl, lbt) for ntl in args.n_tau_lengths_list]
        print '   lb tau %.4f' % lbt
        print '      n gen: %s' % ' '.join(str(n) for n in gen_list)
        for n_gen in gen_list:

            # if ntl is not None:
            #     this_outdir = '%s/ XXX %s/n-tau-%.2f-lbt-%.4f' % (args.base_outdir, args.label, ntl, lbt)
            # elif n_gen is not None:
            this_outdir = '%s/n_gen-%d-lbt-%.4f' % (outdir, n_gen, lbt)

            if os.path.exists(get_outfname(this_outdir)):
                if args.make_plots:
                    with open(get_outfname(this_outdir)) as outfile:
                        info = yaml.load(outfile, Loader=yaml.Loader)
                    if lbt not in parsed_info:
                        parsed_info[lbt] = {}
                    parsed_info[lbt][n_gen] = info['max']['lbi']
                else:
                    print '         output exists, skipping: %s' % get_outfname(this_outdir)
                continue

            if not os.path.exists(this_outdir):
                os.makedirs(this_outdir)

            # lbvals = treeutils.get_min_lbi(args.seq_len, args. XXX lb_tau)
            max_name, max_lbi, lbvals = treeutils.get_max_lbi(args.seq_len, lbt, n_generations=n_gen)
            # TODO maybe should write tree + lb values to file here?

            with open(get_outfname(this_outdir), 'w') as outfile:
                yaml.dump({'max' : {'name' : max_name, 'lbi' : max_lbi}}, outfile)

            plotdir = this_outdir + '/plots'
            utils.prep_dir(plotdir, wildlings='*.svg')
            cmdfos = [plotting.get_lb_tree_cmd(lbvals['tree'], '%s/tree.svg' % plotdir, 'lbi', 'affinities', args.ete_path, args.workdir, metafo=lbvals, tree_style='circular')]
            utils.run_cmds(cmdfos, clean_on_success=True, shell=True, debug='print')

    if args.make_plots:
        fig, ax = plotting.mpl_init()
        for lbt in sorted(parsed_info, reverse=True):
            n_gen_list, max_lbi_list = zip(*sorted(parsed_info[lbt].items(), key=operator.itemgetter(0)))
            ax.plot(n_gen_list, max_lbi_list, label='%.4f' % lbt, alpha=0.7, linewidth=4)
        plotting.mpl_finish(ax, outdir, 'tau-vs-n-gen-vs-max-lbi', xlabel='N generations', ylabel='Max LBI', leg_title='tau', leg_prop={'size' : 12}, leg_loc=(0.04, 0.67))

        # there's got to be a way to get a log plot without redoing everything, but I'm not sure what it is
        fig, ax = plotting.mpl_init()
        for lbt in sorted(parsed_info, reverse=True):
            n_gen_list, max_lbi_list = zip(*sorted(parsed_info[lbt].items(), key=operator.itemgetter(0)))
            ax.plot(n_gen_list, max_lbi_list, label='%.4f' % lbt, alpha=0.7, linewidth=4)
        plotting.mpl_finish(ax, outdir, 'tau-vs-n-gen-vs-max-lbi-log', log='y', xlabel='N generations', ylabel='Max LBI', leg_title='tau', leg_prop={'size' : 12}, leg_loc=(0.04, 0.67))

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
        return [args.random_seed + i for i in range(args.n_replicates)]
    else:
        return args.__dict__[dkey(sv)]

# ----------------------------------------------------------------------------------------
def get_var_info(args, scan_vars):
    def handle_var(svar, val_lists, valstrs):
        convert_fcn = str if svar in ['carry-cap', 'seed', 'lb-tau'] else lambda vlist: ':'.join(str(v) for v in vlist)
        if len(getsargval(svar)) > 1:
            varnames.append(svar)
            val_lists = [vlist + [sv] for vlist in val_lists for sv in getsargval(svar)]
            valstrs = [vlist + [convert_fcn(sv)] for vlist in valstrs for sv in getsargval(svar)]
        else:
            base_args.append('--%s %s' % (svar, convert_fcn(getsargval(svar)[0])))
        return val_lists, valstrs

    base_args = []
    varnames = []
    val_lists, valstrs = [[]], [[]]
    for svar in scan_vars:
        val_lists, valstrs = handle_var(svar, val_lists, valstrs)

    return base_args, varnames, val_lists, valstrs

# ----------------------------------------------------------------------------------------
def make_plots(args, use_relative_affy=True, min_ptile_to_plot=85.):  # have to go lower than 85. for small sample sizes
    _, varnames, val_lists, valstrs = get_var_info(args, args.scan_vars['partition'])
    plotvals = collections.OrderedDict()
    print '  plotting %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    for vlists, vstrs in zip(val_lists, valstrs):
        n_per_gen = vlists[varnames.index('n-sim-seqs-per-gen')]
        assert len(n_per_gen) == 1
        n_per_gen = n_per_gen[0]
        assert len(args.obs_times_list) == 1
        n_obs = n_per_gen * len(args.obs_times_list[0])
        assert len(args.carry_cap_list) == 1
        n_tot = args.carry_cap_list[0]
        obs_frac = n_obs / float(n_tot)

        yfname = '%s/true-tree-metrics/lbi-vs-affinity-true-tree-ptiles%s.yaml' % (get_partition_plotdir(varnames, vstrs), '-relative' if use_relative_affy else '')
        with open(yfname) as yfile:
            info = yaml.load(yfile, Loader=yaml.Loader)
        diff_to_perfect = numpy.mean([pafp - afp for lbp, afp, pafp in zip(info['lb_ptiles'], info['mean_affy_ptiles'], info['perfect_vals']) if lbp > min_ptile_to_plot])
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
    plotting.mpl_finish(ax, get_plotdir(), 'affy-ptiles-obs-frac-vs-lb-tau', xlabel='tau', ylabel='mean difference to perfect for lbi ptiles in [%.0f, 100]' % min_ptile_to_plot,
                        leg_title='fraction sampled', leg_prop={'size' : 12}, leg_loc=(0.04, 0.67),
                        xticks=lb_taus, xticklabels=[str(t) for t in lb_taus], xticklabelsize=16,
    )

# ----------------------------------------------------------------------------------------
def run_bcr_phylo(args):  # also caches parameters
    base_args, varnames, _, valstrs = get_var_info(args, args.scan_vars['simu'])
    cmdfos = []
    print '  running %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    for icombo, vstrs in enumerate(valstrs):
        print '   %s' % ' '.join(vstrs)
        outdir = get_bcr_phylo_outdir(varnames, vstrs)
        if utils.output_exists(args, get_parameter_dir(varnames, vstrs) + '/hmm/hmms', offset=8):
            continue
        cmd = './bin/bcr-phylo-run.py --actions simu:cache-parameters --base-outdir %s %s' % (outdir, ' '.join(base_args))
        for vname, vstr in zip(varnames, vstrs):
            cmd += ' --%s %s' % (vname, vstr)
        if args.n_procs > 1:
            cmd += ' --n-procs %d' % args.n_procs  # i think this only gets used for partitioning
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
    utils.run_cmds(cmdfos, debug='write', batch_system='slurm' if args.slurm else None)

# ----------------------------------------------------------------------------------------
def partition(args):
    _, varnames, _, valstrs = get_var_info(args, args.scan_vars['partition'])  # can't use base_args a.t.m. since it has the simulation/bcr-phylo args in it
    cmdfos = []
    print '  running %d combinations of: %s' % (len(valstrs), ' '.join(varnames))
    for icombo, vstrs in enumerate(valstrs):
        print '   %s' % ' '.join(vstrs)

        partition_fname = get_partition_fname(varnames, vstrs)
        if utils.output_exists(args, partition_fname, outlabel='partition', offset=8):
            continue

        action = 'partition'  # 'get-tree-metrics'  # maybe eventually make this a cmd line option, but this is fine for now
        cmd = './bin/partis %s --is-simu --infname %s --parameter-dir %s --plotdir %s --outfname %s' % (action, get_bcr_phylo_outfname(varnames, vstrs), get_parameter_dir(varnames, vstrs), get_partition_plotdir(varnames, vstrs), partition_fname)
        if action == 'partition':
            cmd += ' --n-final-clusters 1 --get-tree-metrics --n-procs %d' % args.n_procs
        cmd += ' --lb-tau %s' % vstrs[varnames.index('lb-tau')]
        cmd += ' --seed %s' % args.random_seed  # NOTE second/commented version this is actually wrong: vstrs[varnames.index('seed')]  # there isn't actually a reason for different seeds here (we want the different seeds when running bcr-phylo), but oh well, maybe it's a little clearer this way
        cmdfos += [{
            'cmd_str' : cmd,
            'outfname' : partition_fname,
            'logdir' : get_partition_outdir(varnames, vstrs),
        }]
        print '     %s %s' % (utils.color('red', 'run'), cmd)
    utils.run_cmds(cmdfos, debug='write', batch_system='slurm' if args.slurm else None)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('action', choices=['get-max-lbi', 'run-bcr-phylo', 'partition', 'plot'])
parser.add_argument('--carry-cap-list', default='250:1000:4000')
parser.add_argument('--n-sim-seqs-per-gen-list', default='50,75,80:200,250')
parser.add_argument('--obs-times-list', default='30,40,50:125,150')
parser.add_argument('--lb-tau-list', default='0.0005:0.001:0.002:0.003:0.005:0.008:0.012')
parser.add_argument('--n-tau-lengths-list', help='set either this or --n-generations-list')
parser.add_argument('--n-generations-list', default='4:5:6:7:8:9:10', help='set either this or --n-tau-lengths-list')
parser.add_argument('--seq-len', default=400, type=int)
parser.add_argument('--n-replicates', default=1, type=int)
parser.add_argument('--n-procs', default=1, type=int)
parser.add_argument('--random-seed', default=0, type=int, help='note that if --n-replicates is greater than 1, this is only the random seed of the first replicate')
parser.add_argument('--base-outdir', default='%s/partis/tree-metrics' % os.getenv('fs', default=os.getenv('HOME')))
parser.add_argument('--label', default='test')
parser.add_argument('--make-plots', action='store_true')
parser.add_argument('--overwrite', action='store_true')  # not really propagated to everything I think
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
if [args.n_tau_lengths_list, args.n_generations_list].count(None) != 1:
    raise Exception('have to set exactly one of --n-tau-lengths, --n-generations')

if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work')))

# ----------------------------------------------------------------------------------------
if args.action == 'get-max-lbi':
    calc_max_lbi(args)
elif args.action == 'run-bcr-phylo':
    run_bcr_phylo(args)
elif args.action == 'partition':
    partition(args)
elif args.action == 'plot':
    make_plots(args)
