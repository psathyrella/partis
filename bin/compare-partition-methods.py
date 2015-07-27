#!/usr/bin/env python
import os
import argparse
import glob
import sys
import numpy
import random
import re
from collections import OrderedDict
import time
import csv
from subprocess import check_call, Popen, check_output
sys.path.insert(1, './python')

from humans import humans
import plotting
from hist import Hist
import seqfileopener
import utils
from clusterplot import ClusterPlot
from glomerator import Glomerator

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')
parser.add_argument('--dataset', choices=['stanford', 'adaptive'], default='adaptive')
parser.add_argument('--only-run')  # colon-separated list of human,subset pairs to run, e.g. A,3:C,8
parser.add_argument('--mutation-multipliers', default='1')
parser.add_argument('--n-leaf-list', default='10')
parser.add_argument('--subset', type=int)
parser.add_argument('--n-subsets', type=int)
parser.add_argument('--istartstop')  # NOTE usual zero indexing
all_actions = ['cache-data-parameters', 'simulate', 'cache-simu-parameters', 'partition', 'naive-hamming-partition', 'vsearch-partition', 'run-viterbi', 'run-changeo', 'write-plots', 'compare-sample-sizes', 'compare-subsets']
parser.add_argument('--actions', required=True)  #default=':'.join(all_actions))
args = parser.parse_args()
args.only_run = utils.get_arg_list(args.only_run)
args.actions = utils.get_arg_list(args.actions)
args.mutation_multipliers = utils.get_arg_list(args.mutation_multipliers, intify=True)
args.n_leaf_list = utils.get_arg_list(args.n_leaf_list, intify=True)
args.istartstop = utils.get_arg_list(args.istartstop, intify=True)

assert args.subset is None or args.istartstop is None  # dosn't make sense to set both of them

legends = {'vollmers-0.9' : 'VJ CDR3 0.9',
           'partition partis' : 'full partis',
           'partition' : 'full partis',
           'naive-hamming-partition partis' : 'naive partis',
           'naive-hamming-partition' : 'naive partis',
           'vsearch-partition partis' : 'vsearch partis',
           'vsearch-partition' : 'vsearch partis',
           'changeo' : 'Change-O',
           '0.1-true-singletons' : '10% random singletons',
           '0.1-true-reassign' : '10% random reassign'
           }


changeorandomcrapstr = '_db-pass_parse-select_clone-pass.tab'

if args.dataset == 'stanford':
    datadir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers'
    files = [ datadir + '/' + f for f in os.listdir(datadir)]
elif args.dataset == 'adaptive':
    datadirs = [ '/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/' + h for h in humans['adaptive'] ]
    files = []
    for datadir in datadirs:
        files += [ datadir + '/' + fname for fname in os.listdir(datadir) if '-M_merged.tsv.bz2' in fname ]  # if you switch to naive (N), be careful 'cause A is split in pieces

# ----------------------------------------------------------------------------------------
def leafmutstr(n_leaves, mut_mult):
    return 'simu-' + str(n_leaves) + '-leaves-' + str(mut_mult) + '-mutate'

# ----------------------------------------------------------------------------------------
def get_simfname(label, n_leaves, mut_mult):
    return fsdir + '/' + label + '/' + leafmutstr(n_leaves, mut_mult) + '.csv'  # NOTE duplicate code

# ----------------------------------------------------------------------------------------
def generate_incorrect_partition(true_partition, n_misassigned, error_type, debug=False):
    new_partition = list(true_partition)
    if debug:
        print '  before', new_partition
    for _ in range(n_misassigned):
        iclust = random.randint(0, len(new_partition) - 1)  # choose a cluster from which to remove a sequence
        iseq = random.randint(0, len(new_partition[iclust]) - 1)  # and choose the sequence to remove from this cluster
        uid = new_partition[iclust][iseq]
        new_partition[iclust].remove(uid)  # remove it
        if [] in new_partition:
            new_partition.remove([])
        if error_type == 'singletons':  # put the sequence in a cluster by itself
            new_partition.append([uid, ])
            if debug:
                print '    %s: %d --> singleton' % (uid, iclust)
        elif error_type == 'reassign':  # choose a different cluster to add it to
            inewclust = iclust
            while inewclust == iclust:
                inewclust = random.randint(0, len(new_partition) - 1)
            new_partition[inewclust].append(uid)
            if debug:
                print '    %s: %d --> %d' % (uid, iclust, inewclust)
        else:
            raise Exception('%s not among %s' % (error_type, 'singletons, reassign'))
    if debug:
        print '  after', new_partition
    return new_partition

# ----------------------------------------------------------------------------------------
def write_latex_table(adj_mis):
    for mut_mult in args.mutation_multipliers:
        for n_leaves in args.n_leaf_list:
            if n_leaves == args.n_leaf_list[0]:
                print '\\textbf{multiplier} & \\textbf{program}  ',
            print ' & %d  ' % n_leaves,
        print '\\\\'
        print '\\hline'
        iname = 0
        for name in adj_mis[args.n_leaf_list[0]][args.mutation_multipliers[0]]:
            if name == 'vollmers-0.5':
                continue
            if iname == 0:
                print '\\multirow{3}{*}{$\\times %d$} & ' % mut_mult,
            else:
                print '& ',
            print '%25s' % legends.get(name, name),
            iname += 1
            for n_leaves in args.n_leaf_list:
                try:
                    val, err = adj_mis[n_leaves][mut_mult][name]
                    print '  &    %5.2f $\\pm$ %.2f' % (val, err),
                except TypeError:
                    print '  &    %5.2f' % adj_mis[n_leaves][mut_mult][name],
            print '\\\\'

# ----------------------------------------------------------------------------------------
def parse_vollmers(these_hists, these_adj_mis, simfname, outdir):
    vollmers_fname = simfname.replace('.csv', '-run-viterbi.csv')
    with open(vollmers_fname) as vfile:
        vreader = csv.DictReader(vfile)
        for line in vreader:
            vhist = plotting.get_cluster_size_hist(utils.get_partition_from_str(line['clusters']))
            histfname = outdir + '/hists/vollmers-'  + line['threshold'] + '.csv'
            vhist.write(histfname)
            these_hists['vollmers-' + line['threshold']] = vhist
            these_adj_mis['vollmers-' + line['threshold']] = float(line['adj_mi'])
            write_adj_mi(float(line['adj_mi']), outdir + '/adj_mis/' + os.path.basename(histfname))

            truehist = plotting.get_cluster_size_hist(utils.get_partition_from_str(line['true_clusters']))  # true partition is also written here, for lack of a better place (note that it's of course the same for all thresholds)
            truehist.write(outdir + '/hists/true.csv')  # will overwite itself a few times
            these_hists['true'] = truehist

# ----------------------------------------------------------------------------------------
def parse_changeo(these_hists, these_adj_mis, simfname, simfbase, outdir):

    input_info, reco_info = seqfileopener.get_seqfile_info(simfname, is_data=False)
    
    indir = fsdir.replace('/partis-dev/_output', '/changeo')
    # indir = indir.replace('changeo', 'changeo.bak')
    infname = indir + '/' + simfbase.replace('-', '_') + '/' + changeorandomcrapstr
    if args.subset is not None:
        infname = infname.replace(changeorandomcrapstr, 'subset-' + str(args.subset) + changeorandomcrapstr)
    if args.istartstop is not None:  # TODO not yet functional
        assert False
    
    id_clusters = {}  # map from cluster id to list of seq ids
    with open(infname) as chfile:
        reader = csv.DictReader(chfile, delimiter='\t')
        for line in reader:
            clid = line['CLONE']
            uid = line['SEQUENCE_ID']
            if clid not in id_clusters:
                id_clusters[clid] = []
            id_clusters[clid].append(uid)
    
    partition = [ids for ids in id_clusters.values()]
    these_hists['changeo'] = plotting.get_cluster_size_hist(partition)
    these_hists['changeo'].write(outdir + '/hists/changeo.csv')
    adj_mi_fname = infname.replace(changeorandomcrapstr, '-adj_mi.csv')
    check_call(['cp', adj_mi_fname, outdir + '/adj_mis/changeo.csv'])
    these_adj_mis['changeo'] = read_adj_mi(adj_mi_fname)

# ----------------------------------------------------------------------------------------
def parse_partis(action, these_hists, these_adj_mis, simfname, outdir):
    args.infnames = [simfname.replace('.csv', '-' + action + '.csv'), ]  # NOTE make sure not to add any args here that conflict with the real command line args
    args.is_data = False
    args.use_all_steps = False
    args.normalize_axes = []
    args.xbounds, args.adjmi_bounds, args.logprob_bounds = None, None, None
    cplot = ClusterPlot(args)
    cplot.tmp_cluster_size_hist.write(outdir + '/hists/' + action + '.csv')
    these_hists[action + ' partis'] = cplot.tmp_cluster_size_hist
    these_adj_mis[action + ' partis'] = cplot.adj_mi_at_max_logprob
    write_adj_mi(cplot.adj_mi_at_max_logprob, outdir + '/adj_mis/' + action + '.csv')

# ----------------------------------------------------------------------------------------
def write_adj_mi(adj_mi, fname):
    if not os.path.exists(os.path.dirname(fname)):
        os.makedirs(os.path.dirname(fname))
    with open(fname, 'w') as outfile:
        writer = csv.DictWriter(outfile, ['adj_mi'])
        writer.writerow({'adj_mi' : adj_mi})

# ----------------------------------------------------------------------------------------
def read_adj_mi(fname):
    """ read file with a single number """
    with open(fname) as infile:
        reader = csv.DictReader(infile, fieldnames=['adj_mi'])
        for line in reader:
            return float(line['adj_mi'])

# ----------------------------------------------------------------------------------------
def write_all_plot_csvs(label):
    hists, adj_mis = {}, {}
    for n_leaves in args.n_leaf_list:
        for mut_mult in args.mutation_multipliers:
            write_each_plot_csvs(label, n_leaves, mut_mult, hists, adj_mis)

    write_latex_table(adj_mis)

# ----------------------------------------------------------------------------------------
def write_each_plot_csvs(label, n_leaves, mut_mult, hists, adj_mis):
    if n_leaves not in hists:
        hists[n_leaves] = {}
        adj_mis[n_leaves] = {}
    if mut_mult not in hists[n_leaves]:
        hists[n_leaves][mut_mult] = OrderedDict()
        adj_mis[n_leaves][mut_mult] = OrderedDict()
    these_hists = hists[n_leaves][mut_mult]
    these_adj_mis = adj_mis[n_leaves][mut_mult]

    simfname = get_simfname(label, n_leaves, mut_mult)
    if args.subset is not None:  # TODO not yet functional
        simfname = simfname.replace(label + '/', label + '/subset-' + str(args.subset) + '/')
    if args.istartstop is not None:  # TODO not yet functional
        assert False
        # subsimfname = simfname.replace(label + '/', label + '/subset-' + str(args.subset) + '/')
    simfbase = leafmutstr(n_leaves, mut_mult)
    csvdir = os.path.dirname(simfname) + '/' + simfbase

    # then vollmers annotation (and true hists)
    parse_vollmers(these_hists, these_adj_mis, simfname, csvdir)

    # then changeo
    parse_changeo(these_hists, these_adj_mis, simfname, simfbase, csvdir)

    # first do partis stuff
    for ptype in ['vsearch-', 'naive-hamming-', '']:
        parse_partis(ptype + 'partition', these_hists, these_adj_mis, simfname, csvdir)

    plotdir = os.getenv('www') + '/partis/clustering/' + label
    plotting.plot_cluster_size_hists(plotdir + '/plots/' + simfbase + '.svg', these_hists, title='%d leaves, %dx mutation' % (n_leaves, mut_mult), legends=legends, xmax=n_leaves*3.01)
    check_call(['./bin/makeHtml', plotdir, '3', 'null', 'svg'])
    check_call(['./bin/permissify-www', plotdir])

    # for name, adj_mi in these_adj_mis.items():
    #     with open(os.path.dirname(simfname) + '/adj_mis/' + name.replace(' ', '_') + '.csv', 'w') as adjmifile:
    #         adjmifile.write(adj_mi + '\n')

# ----------------------------------------------------------------------------------------
def compare_all_subsets(label):
    hists, adj_mis = {}, {}
    for n_leaves in args.n_leaf_list:
        for mut_mult in args.mutation_multipliers:
            compare_each_subsets(label, n_leaves, mut_mult, hists, adj_mis)

    write_latex_table(adj_mis)

# ----------------------------------------------------------------------------------------
def compare_each_subsets(label, n_leaves, mut_mult, hists, adj_mis):
    if n_leaves not in hists:
        hists[n_leaves] = {}
        adj_mis[n_leaves] = {}
    if mut_mult not in hists[n_leaves]:
        hists[n_leaves][mut_mult] = OrderedDict()
        adj_mis[n_leaves][mut_mult] = OrderedDict()
    these_hists = hists[n_leaves][mut_mult]
    these_adj_mis = adj_mis[n_leaves][mut_mult]

    basedir = fsdir + '/' + label
    expected_methods = ['true', 'vollmers-0.9', 'changeo', 'vsearch-partition', 'naive-hamming-partition', 'partition']  # mostly so we can specify the order
    tmp_adj_mis = OrderedDict()
    for method in expected_methods:
        method_hists = []
        for isub in range(args.n_subsets):
            subdir = basedir + '/subset-' + str(isub)

            # hists
            histfname = subdir + '/' + leafmutstr(n_leaves, mut_mult) + '/hists/' + method + '.csv'
            method_hists.append(Hist(fname=histfname))
            if method == 'true':
                continue

            # adj mis
            fname = subdir + '/' + leafmutstr(n_leaves, mut_mult) + '/adj_mis/' + method + '.csv'
            if method not in tmp_adj_mis:
                tmp_adj_mis[method] = []
            tmp_adj_mis[method].append(read_adj_mi(fname))

        these_hists[method] = plotting.make_mean_hist(method_hists)

    plotdir = os.getenv('www') + '/partis/clustering/subsets/' + label
    plotting.plot_cluster_size_hists(plotdir + '/plots/' + leafmutstr(n_leaves, mut_mult) + '.svg', these_hists, title='%d leaves, %dx mutation' % (n_leaves, mut_mult), legends=legends, xmax=n_leaves*3.01)
    check_call(['./bin/makeHtml', plotdir, '3', 'null', 'svg'])
    check_call(['./bin/permissify-www', plotdir])
    for meth, vals in tmp_adj_mis.items():
        mean = numpy.mean(vals)
        std = numpy.std(vals)
        these_adj_mis[meth] = (mean, std)
        print '  %30s %.3f +/- %.3f' % (meth, mean, std)

# ----------------------------------------------------------------------------------------
def get_misassigned_adj_mis(simfname, misassign_fraction, nseq_list, error_type):
    input_info, reco_info = seqfileopener.get_seqfile_info(simfname, is_data=False)
    n_reps = 1
    uid_list = input_info.keys()
    new_partitions = {}
    for nseqs in nseq_list:
        for irep in range(n_reps):  # repeat <nreps> times
            istart = irep * nseqs
            istop = istart + nseqs
            uids = uid_list[istart : istop]
            true_partition = utils.get_true_clusters(uids, reco_info).values()
            n_misassigned = int(misassign_fraction * nseqs)
            new_partition = generate_incorrect_partition(true_partition, n_misassigned, error_type=error_type)
            # new_partition = generate_incorrect_partition(true_partition, n_misassigned, error_type='singletons')
            new_partitions[nseqs] = new_partition
    return {nseqs : utils.mutual_information(new_partitions[nseqs], reco_info) for nseqs in nseq_list}

# ----------------------------------------------------------------------------------------
def compare_sample_sizes(label, n_leaves, mut_mult):
    expected_methods = ['vsearch-partition', 'naive-hamming-partition', 'partition']  # mostly so we can specify the order
    basedir = fsdir + '/' + label
    nseq_list = []
    method_results = OrderedDict()
    for meth in expected_methods:
        method_results[meth] = {}
    for dirname in glob.glob(basedir + '/istartstop-*'):
        dumstr, istart, istop = dirname.split('/')[-1].split('-')
        istart, istop = int(istart), int(istop)
        nseqs = istop - istart
        # if nseqs > 900:
        #     continue
        nseq_list.append(nseqs)
        for meth in expected_methods:
            csvfname = dirname + '/' + leafmutstr(n_leaves, mut_mult) + '-' + meth + '.csv'
            args.infnames = [csvfname, ]  # NOTE make sure not to add any args here that conflict with the real command line args
            args.is_data = False
            args.use_all_steps = False
            args.normalize_axes = []
            args.xbounds, args.adjmi_bounds, args.logprob_bounds = None, None, None
            cplot = ClusterPlot(args)
            # cplot.tmp_cluster_size_hist.write(outdir + '/hists/' + action + '.csv')
            # these_hists[action + ' partis'] = cplot.tmp_cluster_size_hist
            method_results[meth][nseqs] = cplot.adj_mi_at_max_logprob

    nseq_list.sort()

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    fsize = 20
    mpl.rcParams.update({
        # 'font.size': fsize,
        'legend.fontsize': fsize,
        'axes.titlesize': fsize,
        # 'axes.labelsize': fsize,
        'xtick.labelsize': fsize,
        'ytick.labelsize': fsize,
        'axes.labelsize': fsize
    })
    fig, ax = plt.subplots()
    fig.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.95)

    plots = {}

    colors = {'vsearch-partition' : '#008b8b',
              'naive-hamming-partition' : '#ff8c00',
              'partition' : '#cc0000',
              '0.1-true-singletons' : '#006600',
              '0.1-true-reassign' : '#32cd32'}
    method_results['0.1-true-singletons'] = get_misassigned_adj_mis(get_simfname(label, n_leaves, mut_mult), 0.1, nseq_list, 'singletons')
    method_results['0.1-true-reassign'] = get_misassigned_adj_mis(get_simfname(label, n_leaves, mut_mult), 0.1, nseq_list, 'reassign')
    for meth in method_results:
        linewidth = 2
        linestyle = '-'
        if 'true' in meth:
            linewidth = 4
            linestyle = '--'
        plots[meth] = ax.plot(nseq_list, [method_results[meth][ns] for ns in nseq_list], linewidth=linewidth, label=legends.get(meth, meth), color=colors.get(meth, 'grey'), linestyle=linestyle, alpha=1)
        if 'true' not in meth:
            plt.scatter(nseq_list, [method_results[meth][ns] for ns in nseq_list], color=colors.get(meth, 'grey'), linestyle='-', alpha=1, s=[30 for _ in range(len(nseq_list))])

    legend = ax.legend(loc=(0.55, 0.1))
    sns.despine(trim=True, bottom=True)
    sns.set_style("ticks")
    plt.title('%d leaves, %dx mutation' % (n_leaves, mut_mult))
    plt.xlabel('sample size')
    plt.ylabel('adjusted MI')
    plt.subplots_adjust(bottom=0.14, left=0.14)
    ax.set_xscale('log')
    xmin, xmax = nseq_list[0], nseq_list[-1] + 10
    plt.xlim(xmin, xmax)
    plt.ylim(0, 1.08)
    potential_xticks = [5, 10, 25, 50, 100, 300, 1000]
    xticks = [xt for xt in potential_xticks if xt < xmax]
    plt.xticks(xticks, [str(xt) for xt in xticks])
    plotdir = os.getenv('www') + '/partis/clustering/sample-sizes/' + label
    outfname = plotdir + '/plots/sample-sizes.svg'
    if not os.path.exists(plotdir + '/plots'):
        os.makedirs(plotdir + '/plots')
    plt.savefig(outfname)
    check_call(['./bin/makeHtml', plotdir, '3', 'null', 'svg'])
    check_call(['./bin/permissify-www', plotdir])

    return

# ----------------------------------------------------------------------------------------
def execute(action, label, datafname, n_leaves=None, mut_mult=None):
    real_action = action
    if 'partition' in action:
        real_action = 'partition'
    cmd = './bin/run-driver.py --label ' + label + ' --action ' + real_action
    if n_leaves is not None:
        simfname = get_simfname(label, n_leaves, mut_mult)
        if args.subset is not None:
            ntot = int(check_output(['wc', '-l', simfname]).split()[0]) - 1
            subsimfname = simfname.replace(label + '/', label + '/subset-' + str(args.subset) + '/')
            if os.path.exists(subsimfname):
                print '      subset file exists %s' % subsimfname
            else:
                print '      subsetting %d / %d' % (args.subset, args.n_subsets)
                if not os.path.exists(os.path.dirname(subsimfname)):
                    os.makedirs(os.path.dirname(subsimfname))
                check_call('head -n1 ' + simfname + ' >' + subsimfname, shell=True)
                n_per_subset = int(float(ntot) / args.n_subsets)  # NOTE ignores remainders, i.e. last few sequences
                istart = args.subset * n_per_subset + 2  # NOTE sed indexing (one-indexed with inclusive bounds). Also note extra +1 to avoid header
                istop = istart + n_per_subset - 1
                check_call('sed -n \'' + str(istart) + ',' + str(istop) + ' p\' ' + simfname + '>>' + subsimfname, shell=True)
            simfname = subsimfname
        if args.istartstop is not None:
            subsimfname = simfname.replace(label + '/', label + '/istartstop-' + '-'.join([str(i) for i in args.istartstop]) + '/')
            if os.path.exists(subsimfname):
                print '      subset file exists %s' % subsimfname
            else:
                print '      subsetting %d seqs with indices %d --> %d' % (args.istartstop[1] - args.istartstop[0], args.istartstop[0], args.istartstop[1])
                if not os.path.exists(os.path.dirname(subsimfname)):
                    os.makedirs(os.path.dirname(subsimfname))
                check_call('head -n1 ' + simfname + ' >' + subsimfname, shell=True)
                check_call('sed -n \'' + str(args.istartstop[0] + 2) + ',' + str(args.istartstop[1] + 1) + ' p\' ' + simfname + '>>' + subsimfname, shell=True)  # NOTE conversion from standdard zero indexing to sed inclusive one-indexing (and +1 for header line)
            simfname = subsimfname
    extras = []

    def output_exists(outfname):
        if os.path.exists(outfname):
            print '                      partition output exists, skipping (%s)' % outfname
            return True
        else:
            return False

    if action == 'cache-data-parameters':
        if os.path.exists(fsdir + '/' + label + '/data'):
            print '                      data parametes exist in %s, skipping ' % (fsdir + '/' + label + '/data')
            return
        cmd += ' --datafname ' + datafname
        extras += ['--n-max-queries', + n_data_to_cache]
        n_procs = n_data_to_cache / 500
    elif action == 'simulate':
        if os.path.exists(simfname):
            print '                      simfile exists, skipping (%s)' % simfname
            return
        cmd += ' --simfname ' + simfname
        extras += ['--n-sim-events', int(float(n_sim_seqs) / n_leaves)]
        extras += ['--n-leaves', n_leaves, '--mutation-multiplier', mut_mult]
        n_procs = 10
    elif action == 'cache-simu-parameters':
        if os.path.exists(simfname.replace('.csv', '')):
            print '                      simu parameters exist in %s, skipping ' % simfname.replace('.csv', '')
            return
        cmd += ' --simfname ' + simfname
        n_procs = 20
    elif action == 'partition':
        if os.path.exists(simfname.replace('.csv', '-partition.csv')):
            print '                      partition output exists, skipping (%s)' % simfname.replace('.csv', '-partition.csv')
            return
        cmd += ' --simfname ' + simfname
        cmd += ' --outfname ' + simfname.replace('.csv', '-' + action + '.csv')
        extras += ['--n-max-queries', n_to_partition]
        n_procs = max(1, n_to_partition / 15)  # something like 15 seqs/process to start with
    elif action == 'naive-hamming-partition':
        if os.path.exists(simfname.replace('.csv', '-naive-hamming-partition.csv')):
            print '                      partition output exists, skipping (%s)' % simfname.replace('.csv', '-naive-hamming-partition.csv')
            return
        cmd += ' --simfname ' + simfname
        cmd += ' --outfname ' + simfname.replace('.csv', '-' + action + '.csv')
        extras += ['--n-max-queries', n_to_partition, '--auto-hamming-fraction-bounds']
        n_procs = max(1, n_to_partition / 30)
    elif action == 'vsearch-partition':
        outfname = simfname.replace('.csv', '-' + action + '.csv')
        if output_exists(outfname):
            return
        cmd += ' --simfname ' + simfname
        cmd += ' --outfname ' + outfname
        extras += ['--n-max-queries', n_to_partition, '--naive-vsearch']
        n_procs = max(1, n_to_partition / 100)  # only used for ighutil step
    elif action == 'run-viterbi':
        if os.path.exists(simfname.replace('.csv', '-run-viterbi.csv')):
            print '                      vollmers output exists, skipping (%s)' % simfname.replace('.csv', '-run-viterbi.csv')
            return
        extras += ['--annotation-clustering', 'vollmers', '--annotation-clustering-thresholds', '0.5:0.9']
        cmd += ' --simfname ' + simfname
        extras += ['--n-max-queries', n_to_partition]
        n_procs = max(1, n_to_partition / 50)
    elif action == 'run-changeo':
        changeodir = '/home/dralph/work/changeo/changeo'
        changeo_fsdir = '/fh/fast/matsen_e/dralph/work/changeo'  #.bak'
        _simfbase = leafmutstr(n_leaves, mut_mult).replace('-', '_')
        imgtdir = changeo_fsdir + '/' + _simfbase
        if os.path.isdir(imgtdir):
            print '                      already untar\'d into %s' % imgtdir
        else:
            tar_cmd = 'mkdir ' + imgtdir + ';'
            tar_cmd += ' tar Jxvf ' + imgtdir + '.txz --exclude=\'IMGT_HighV-QUEST_individual_files_folder/*\' -C ' + imgtdir
            check_call(tar_cmd, shell=True)

        if args.subset is not None:
            subset_dir = imgtdir + '/subset-' + str(args.subset)
            if not os.path.exists(subset_dir):
                os.makedirs(subset_dir)
                tsvfnames = glob.glob(imgtdir + '/*.txt')
                check_call(['cp', '-v', imgtdir + '/11_Parameters.txt', subset_dir + '/'])
                tsvfnames.remove(imgtdir + '/11_Parameters.txt')
                tsvfnames.remove(imgtdir + '/README.txt')
                input_info, reco_info = seqfileopener.get_seqfile_info(simfname, is_data=False)
                subset_ids = input_info.keys()
                utils.subset_files(subset_ids, tsvfnames, subset_dir)
            imgtdir = subset_dir

        def run(cmdstr):
            print 'RUN %s' % cmdstr
            check_call(cmdstr.split())
        
        resultfname = imgtdir + changeorandomcrapstr
        if os.path.exists(resultfname):
            print '                         changeo already finished'
            return

        check_call(['./bin/csv2fasta', simfname])
        check_call(['mv', simfname.replace('.csv', '.fa'), simfname.replace('.csv', '.fasta')])
        cmd = changeodir + '/MakeDb.py imgt -i ' + imgtdir + ' -s ' + simfname.replace('.csv', '.fasta')
        run(cmd)
        cmd = changeodir + '/ParseDb.py select -d ' + imgtdir + '_db-pass.tab -f FUNCTIONAL -u T'
        run(cmd)
        cmd = changeodir + '/DefineClones.py bygroup -d ' + imgtdir + '_db-pass_parse-select.tab --act first --model m1n --dist 7'
        run(cmd)

        # read changeo's output and toss it into a csv
        input_info, reco_info = seqfileopener.get_seqfile_info(simfname, is_data=False)
        id_clusters = {}  # map from cluster id to list of seq ids
        with open(resultfname) as chfile:
            reader = csv.DictReader(chfile, delimiter='\t')
            for line in reader:
                clid = line['CLONE']
                uid = line['SEQUENCE_ID']
                if clid not in id_clusters:
                    id_clusters[clid] = []
                id_clusters[clid].append(uid)
        
        partition = [ids for ids in id_clusters.values()]
        # these_hists['changeo'] = plotting.get_cluster_size_hist(partition)
        write_adj_mi(utils.mutual_information(partition, reco_info, debug=True), imgtdir + '-adj_mi.csv')
        return
    elif action == 'compare-sample-sizes':
        compare_sample_sizes(label, n_leaves, mut_mult)
        return
    elif action == 'compare-subsets':
        assert False
    else:
        raise Exception('bad action %s' % action)

    cmd +=  ' --plotdir ' + os.getenv('www') + '/partis'
    n_proc_str = str(n_procs)
    if n_procs > 10:
        n_proc_str += ':10'
        extras += ['--slurm', '--workdir', fsdir + '/_tmp/' + str(random.randint(0,99999))]
        
    cmd += ' --n-procs ' + n_proc_str
    # extras += ['--print-partitions', ]

    cmd += utils.get_extra_str(extras)
    print '   ' + cmd
    # return
    # check_call(cmd.split())
    # return
    logbase = fsdir + '/' + label + '/_logs/' + leafmutstr(n_leaves, mut_mult) + '-' + action
    if args.subset is not None:
        logbase = logbase.replace('_logs/', '_logs/subset-' + str(args.subset) + '/')
    if args.istartstop is not None:
        logbase = logbase.replace('_logs/', '_logs/istartstop-' + '-'.join([str(i) for i in args.istartstop]) + '/')
    if not os.path.exists(os.path.dirname(logbase)):
        os.makedirs(os.path.dirname(logbase))
    proc = Popen(cmd.split(), stdout=open(logbase + '.out', 'w'), stderr=open(logbase + '.err', 'w'))
    procs.append(proc)
    time.sleep(5)

# ----------------------------------------------------------------------------------------
n_to_partition = 1300
if args.istartstop is not None:
    n_to_partition = args.istartstop[1] - args.istartstop[0]
n_data_to_cache = 50000
# mutation_multipliers = ['1']  #['1', '4']
# n_leaf_list = [5]  #, 10, 25, 50]
n_sim_seqs = 10000
fsdir = '/fh/fast/matsen_e/' + os.getenv('USER') + '/work/partis-dev/_output'
procs = []

# ----------------------------------------------------------------------------------------
for datafname in files:
    if args.dataset == 'stanford':
	    human = os.path.basename(datafname).replace('_Lineages.fasta', '')
    elif args.dataset == 'adaptive':
	    human = re.findall('[ABC]', datafname)[0]
    print 'run', human
    label = human
    # label += '.bak'

    if args.only_run is not None and human not in args.only_run:
        continue

    performance_fname = fsdir + '/' + label + '/partition-performance.csv'
    performance_file_header = ['method', 'mut_mult', 'n_leaves', 'adj_mi', 'n_clusters', 'n_true_clusters']

    for action in args.actions:
        if action == 'write-plots' or action == 'compare-subsets':
            continue
        print '      ----> ', action
        if action == 'cache-data-parameters':
            execute(action, label, datafname)
            continue
        for n_leaves in args.n_leaf_list:
            print '  ----> ', n_leaves, ' leaves'
            for mut_mult in args.mutation_multipliers:
                print '         ----> mutate', mut_mult
                execute(action, label, datafname, n_leaves, mut_mult)
                # sys.exit()

    if 'write-plots' in args.actions:
        write_all_plot_csvs(label)
    if 'compare-subsets' in args.actions:
        compare_all_subsets(label)

    break
