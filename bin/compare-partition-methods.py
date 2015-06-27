#!/usr/bin/env python
import os
import argparse
import sys
import random
import re
import time
import csv
from subprocess import check_call, Popen
sys.path.insert(1, './python')

from humans import humans
import plotting
import utils
from clusterplot import ClusterPlot

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')
parser.add_argument('--dataset', choices=['stanford', 'adaptive'], default='adaptive')
parser.add_argument('--only-run')  # colon-separated list of human,subset pairs to run, e.g. A,3:C,8
all_actions = ['cache-data-parameters', 'simulate', 'cache-simu-parameters', 'partition', 'run-viterbi', 'write-plots']
parser.add_argument('--actions', required=True)  #default=':'.join(all_actions))
args = parser.parse_args()
args.only_run = utils.get_arg_list(args.only_run)
args.actions = utils.get_arg_list(args.actions)

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
def write_latex_table(adj_mis):
    for mut_mult in mutation_multipliers:
        print 'mutation x %s' % mut_mult
        for n_leaves in n_leaf_list:
            print '  &  %s  ' % n_leaves,
        print '\\\\'
        for name in adj_mis[n_leaf_list[0]][mutation_multipliers[0]]:
            print '%25s' % name,
            for n_leaves in n_leaf_list:
                print '  &    %5.2f' % adj_mis[n_leaves][mut_mult][name],
            print '\\\\'

# ----------------------------------------------------------------------------------------
def write_all_plot_csvs(label):
    hists, adj_mis = {}, {}
    for n_leaves in n_leaf_list:
        for mut_mult in mutation_multipliers:
            write_each_plot_csvs(label, n_leaves, mut_mult, hists, adj_mis)

# ----------------------------------------------------------------------------------------
def write_each_plot_csvs(label, n_leaves, mut_mult, hists, adj_mis):
    if n_leaves not in hists:
        hists[n_leaves] = {}
        adj_mis[n_leaves] = {}
    if mut_mult not in hists[n_leaves]:
        hists[n_leaves][mut_mult] = {}
        adj_mis[n_leaves][mut_mult] = {}

    simfname = fsdir + '/' + label + '/' + leafmutstr(n_leaves, mut_mult) + '.csv'  # NOTE duplicate code
    simfbase = leafmutstr(n_leaves, mut_mult)
    csvdir = os.path.dirname(simfname) + '/plots'

    # first do partis stuff
    args.infnames = [simfname.replace('.csv', '-partition.csv'), ]  # NOTE make sure not to add any args here that conflict with the real command line args
    args.is_data = False
    args.use_all_steps = False
    args.normalize_axes = []
    args.xbounds, args.adjmi_bounds, args.logprob_bounds = None, None, None
    cplot = ClusterPlot(args)
    cplot.tmp_cluster_size_hist.write(csvdir + '/' + simfbase + '-partis-hist.csv')
    hists[n_leaves][mut_mult]['partis'] = cplot.tmp_cluster_size_hist
    adj_mis[n_leaves][mut_mult]['partis'] = cplot.adj_mi_at_max_logprob

    # then vollmers annotation
    vollmers_fname = simfname.replace('.csv', '-run-viterbi.csv')
    with open(vollmers_fname) as vfile:
        vreader = csv.DictReader(vfile)
        for line in vreader:
            vhist = plotting.get_cluster_size_hist(utils.get_partition_from_str(line['clusters']))
            vhist.write(csvdir + '/' + simfbase + '-vollmers-' + line['threshold'] + '-hist.csv')
            hists[n_leaves][mut_mult]['vollmers-' + line['threshold']] = vhist
            adj_mis[n_leaves][mut_mult]['vollmers-' + line['threshold']] = float(line['adj_mi'])

            truehist = plotting.get_cluster_size_hist(utils.get_partition_from_str(line['true_clusters']))  # true partition is also written here, for lack of a better place (note that it's of course the same for all thresholds)
            truehist.write(csvdir + '/' + simfbase + '-true-hist.csv')  # will overwite itself a few times
            hists[n_leaves][mut_mult]['true'] = truehist

    plotdir = os.getenv('www') + '/partis/clustering/' + label
    plotting.plot_cluster_size_hists(plotdir + '/plots/' + simfbase + '.svg', hists[n_leaves][mut_mult], title='mean leaves %s, mutation x %s' % (n_leaves, mut_mult), xmax=n_leaves*3)
    check_call(['./bin/makeHtml', plotdir, '3', 'null', 'svg'])
    check_call(['./bin/permissify-www', plotdir])

# ----------------------------------------------------------------------------------------
def execute(action, label, datafname, n_leaves=None, mut_mult=None):
    cmd = './bin/run-driver.py --label ' + label + ' --action ' + action
    if n_leaves is not None:
        simfname = fsdir + '/' + label + '/' + leafmutstr(n_leaves, mut_mult) + '.csv'  # NOTE duplicate code
    extras = []

    if action == 'cache-data-parameters':
        if os.path.exists(fsdir + '/' + label + '/data'):
            print '                      data parametes exist in %s, skipping ' % os.path.dirname(simfname)
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
            print '                      simu parametes exist in %s, skipping ' % os.path.dirname(simfname)
            return
        cmd += ' --simfname ' + simfname
        n_procs = 20
    elif action == 'partition':
        if os.path.exists(simfname.replace('.csv', '-partition.csv')):
            print '                      partition output exists, skipping (%s)' % simfname.replace('.csv', '-partition.csv')
            return
        cmd += ' --simfname ' + simfname
        extras += ['--n-max-queries', n_to_partition]
        n_procs = n_to_partition / 15  # something like 15 seqs/process to start with
    elif action == 'run-viterbi':
        if os.path.exists(simfname.replace('.csv', '-run-viterbi.csv')):
            print '                      vollmers output exists, skipping (%s)' % simfname.replace('.csv', '-run-viterbi.csv')
            return
        extras += ['--annotation-clustering', 'vollmers', '--annotation-clustering-thresholds', '0.5:0.9']
        cmd += ' --simfname ' + simfname
        extras += ['--n-max-queries', n_to_partition]
        n_procs = n_to_partition / 50
    else:  # for plotting actions don't do anything here
        return

    cmd +=  ' --plotdir ' + os.getenv('www') + '/partis'
    n_proc_str = str(n_procs)
    if n_procs > 10:
        n_proc_str += ':10'
        extras += ['--slurm', '--workdir', fsdir + '/_tmp/' + str(random.randint(0,99999))]
        
    cmd += ' --n-procs ' + n_proc_str
    # extras += ['--print-partitions', ]

    cmd += utils.get_extra_str(extras)
    print '   ' + cmd
    # check_call(cmd.split())
    logbase = os.path.dirname(simfname) + '/_logs/' + os.path.basename(simfname).replace('.csv', '') + '-' + action
    proc = Popen(cmd.split(), stdout=open(logbase + '.out', 'w'), stderr=open(logbase + '.err', 'w'))
    procs.append(proc)
    # sys.exit()

# ----------------------------------------------------------------------------------------
n_to_partition = 1000
n_data_to_cache = 50000
mutation_multipliers = ['1', '2', '4']
n_leaf_list = [5, 5, 10, 25, 50]
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

    if args.only_run is not None and human not in args.only_run:
        continue

    performance_fname = fsdir + '/' + label + '/partition-performance.csv'
    performance_file_header = ['method', 'mut_mult', 'n_leaves', 'adj_mi', 'n_clusters', 'n_true_clusters']

    for action in args.actions:
        if action == 'write-plots':
            continue
        print '      ----> ', action
        if action == 'cache-data-parameters':
            execute(action, label, datafname)
            continue
        for n_leaves in n_leaf_list:
            print '  ----> ', n_leaves, ' leaves'
            for mut_mult in mutation_multipliers:
                print '         ----> mutate', mut_mult
                execute(action, label, datafname, n_leaves, mut_mult)

    if 'write-plots' in args.actions:
        write_all_plot_csvs(label)

    break
