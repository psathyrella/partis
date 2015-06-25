#!/usr/bin/env python
import os
import sys
import random
import argparse
import re
import time
import csv
from subprocess import check_call, Popen
sys.path.insert(1, './python')

from humans import humans
import plotting
import utils

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')
parser.add_argument('--dataset', choices=['stanford', 'adaptive'], default='adaptive')
parser.add_argument('--only-run')  # colon-separated list of human,subset pairs to run, e.g. A,3:C,8
all_actions = ['cache-data-parameters', 'simulate', 'cache-simu-parameters', 'partition', 'run-viterbi', 'plot']
parser.add_argument('--actions', default=':'.join(all_actions))
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
def execute(action, label, n_leaves, fname, mut_mult):
    simfname = fsdir + '/_output/' + label + '/simu-' + str(n_leaves) + '-leaves-' + str(mut_mult) + '-mutate.csv'
    cmd = './bin/run-driver.py --label ' + label + ' --action ' + action
    extras = []

    if action == 'cache-data-parameters':
        if os.path.exists(os.path.dirname(simfname) + '/data'):
            print '                      data parametes exist in %s, skipping ' % os.path.dirname(simfname)
            return
        cmd += ' --datafname ' + fname
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
        extras += ['--annotation-clustering', 'vollmers', '--annotation-clustering-thresholds', '0.5:0.7:0.8:0.9:0.95']
        cmd += ' --simfname ' + simfname
        extras += ['--n-max-queries', n_to_partition]
        n_procs = n_to_partition / 100
    elif action == 'plot':
        if n_leaves not in hists:
            hists[n_leaves] = {}
            adj_mis[n_leaves] = {}
        if mut_mult not in hists[n_leaves]:
            hists[n_leaves][mut_mult] = {}
            adj_mis[n_leaves][mut_mult] = {}

        # first do partis stuff
        from clusterplot import ClusterPlot
        # make sure not to add any args here that conflict with the real command line args
        args.infnames = [simfname.replace('.csv', '-partition.csv'), ]
        args.is_data = False
        args.use_all_steps = False
        args.normalize_axes = []
        args.xbounds, args.adjmi_bounds, args.logprob_bounds = None, None, None
        cplot = ClusterPlot(args)
        pwriter.writerow({'method' : 'partis',
                          'mut_mult' : mut_mult,
                          'n_leaves' : n_leaves,
                          'adj_mi' : cplot.adj_mi_at_max_logprob,
                          'n_clusters' : cplot.tmp_n_clusters,
                          'n_true_clusters' : cplot.tmp_n_true_clusters})
        cplot.tmp_cluster_size_hist.write(simfname.replace('.csv', '-partition-hist.csv'))
        hists[n_leaves][mut_mult]['partis'] = cplot.tmp_cluster_size_hist
        adj_mis[n_leaves][mut_mult]['partis'] = cplot.adj_mi_at_max_logprob

        # then vollmers annotation
        vollmers_fname = simfname.replace('.csv', '-run-viterbi.csv')
        with open(vollmers_fname) as vfile:
            vreader = csv.DictReader(vfile)
            for line in vreader:
                pwriter.writerow({'method' : 'vollmers-thresh-' + line['threshold'],
                                  'mut_mult' : mut_mult,
                                  'n_leaves' : n_leaves,
                                  'adj_mi' : line['adj_mi'],
                                  'n_clusters' : line['n_clusters'],
                                  'n_true_clusters' : line['n_true_clusters']})
                vhist = plotting.get_cluster_size_hist(utils.get_partition_from_str(line['clusters']))
                vhist.write(vollmers_fname.replace('.csv', '-' + line['threshold'] + '-hist.csv'))
                hists[n_leaves][mut_mult]['vollmers-' + line['threshold']] = vhist
                adj_mis[n_leaves][mut_mult]['vollmers-' + line['threshold']] = float(line['adj_mi'])
                # true partition is also written here, for lack of a better place (not that it's of course the same for all thresholds)
                truehist = plotting.get_cluster_size_hist(utils.get_partition_from_str(line['true_clusters']))
                truehist.write(vollmers_fname.replace('.csv', '-true-hist.csv'))  # will overwite itself a few times
                hists[n_leaves][mut_mult]['true'] = truehist

        plotting.plot_cluster_size_hists(os.path.dirname(simfname), hists[n_leaves][mut_mult])

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
    logbase = fsdir + '/_logs/' + os.path.basename(simfname).replace('.csv', '')
    proc = Popen(cmd.split(), stdout=open(logbase + '.out', 'w'), stderr=open(logbase + '.err', 'w'))
    procs.append(proc)
    # sys.exit()
        
# ----------------------------------------------------------------------------------------
n_to_partition = 1000
n_data_to_cache = 50000
mutation_multipliers = ['1', '2', '4']
n_leaf_list = [3, 5, 10, 25, 50]
n_sim_seqs = 10000
fsdir = '/fh/fast/matsen_e/' + os.getenv('USER') + '/work/partis-dev'
hists, adj_mis = {}, {}
procs = []
# ----------------------------------------------------------------------------------------
for fname in files:
    if args.dataset == 'stanford':
	    human = os.path.basename(fname).replace('_Lineages.fasta', '')
    elif args.dataset == 'adaptive':
	    human = re.findall('[ABC]', fname)[0]
    print 'run', human
    label = human

    if args.only_run is not None and human not in args.only_run:
        continue

    performance_fname = fsdir + '/_output/' + label + '/partition-performance.csv'
    performance_file_header = ['method', 'mut_mult', 'n_leaves', 'adj_mi', 'n_clusters', 'n_true_clusters']
    pfile =  open(performance_fname, 'a')
    pwriter = csv.DictWriter(pfile, performance_file_header)
    pwriter.writeheader()

    # procs = []
    for n_leaves in n_leaf_list:
        # if os.path.exists('_output/' + label + '/simu.csv'):
        #     print 'exists:', '_output/' + label + '/simu.csv'
        #     continue
        print '  ----> ', n_leaves, ' leaves'
        for action in args.actions:
            print '      ----> ', action
            for mut_mult in mutation_multipliers:
                print '         ----> mutate', mut_mult
                execute(action, label, n_leaves, fname, mut_mult)
                if action == 'cache-data-parameters':
                    continue

    
        # check_call(['limit_procs', 'python'])
        # time.sleep(10)

    if 'plot' in args.actions:
        for n_leaves in n_leaf_list:
            for mut_mult in mutation_multipliers:
                print '  %s leaves   %s mut_mult' % (n_leaves, mut_mult)
                for name, adj_mi in adj_mis[n_leaves][mut_mult].items():
                    print '    %20s  %5.2f' % (name, adj_mi)

    break


pfile.close()
