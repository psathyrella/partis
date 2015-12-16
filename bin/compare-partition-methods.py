#!/usr/bin/env python
import os
import argparse
import sys
import re
import csv
sys.path.insert(1, './python')
csv.field_size_limit(sys.maxsize)
from humans import humans
import utils
import compareutils


# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--fsdir', default='/fh/fast/matsen_e/' + os.getenv('USER') + '/work/partis-dev/_output')
parser.add_argument('--dataset', choices=['stanford', 'adaptive'], default='adaptive')
parser.add_argument('--mutation-multipliers', default='1')
parser.add_argument('--data', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--indels', action='store_true')
parser.add_argument('--lonely-leaves', action='store_true')
parser.add_argument('--extra-label-str')
parser.add_argument('--bak', action='store_true')
parser.add_argument('--count-distances', action='store_true')
parser.add_argument('--n-leaf-list', default='10')
parser.add_argument('--subset', type=int)
parser.add_argument('--n-to-partition', type=int, default=5000)
parser.add_argument('--n-data-to-cache', type=int, default=50000)
parser.add_argument('--n-sim-seqs', type=int, default=10000)
parser.add_argument('--n-subsets', type=int)
parser.add_argument('--istartstop')  # NOTE usual zero indexing
parser.add_argument('--startstoplist')  # list of istartstops for comparisons
parser.add_argument('--dont-normalize', action='store_true')
parser.add_argument('--logaxis', action='store_true')
parser.add_argument('--zoom', action='store_true')
parser.add_argument('--humans', default=None)  #'A')
all_actions = ['cache-data-parameters', 'simulate', 'cache-simu-parameters', 'partition', 'naive-hamming-partition', 'vsearch-partition', 'run-viterbi', 'run-changeo', 'run-mixcr', 'run-igscueal', 'write-plots', 'compare-sample-sizes', 'compare-subsets']
parser.add_argument('--actions', required=True)  #, choices=all_actions)  #default=':'.join(all_actions))
args = parser.parse_args()
args.actions = utils.get_arg_list(args.actions)
args.mutation_multipliers = utils.get_arg_list(args.mutation_multipliers, intify=True)
args.n_leaf_list = utils.get_arg_list(args.n_leaf_list, intify=True)
args.istartstop = utils.get_arg_list(args.istartstop, intify=True)
args.startstoplist = utils.get_arg_list(args.startstoplist)
args.humans = utils.get_arg_list(args.humans)

if 'cache-data-parameters' in args.actions:
    args.data = True

assert args.subset is None or args.istartstop is None  # dosn't make sense to set both of them

if args.subset is not None:
    if 'write-plots' not in args.actions:
        assert args.n_subsets == 10  # for all the subset plots, I split into ten subsets, then ended up only using the first thre of 'em, so you have to set n_subsets to 10 if you're running methods, but then to 3 when you're writing plots
    args.n_to_partition = 1300
if args.istartstop is not None:
    assert False  # I think I'm not actually using the istartstop stuff for anything... maybe? in any case I moved all the results to a bak/ subdir of args.fsdir
    args.n_to_partition = args.istartstop[1] - args.istartstop[0]

# ----------------------------------------------------------------------------------------

# # ----------------------------------------------------------------------------------------
# from clusterplot import ClusterPlot
# import plotting
# legends = {}
# legends['v-true'] = 'true (V indels)'
# legends['cdr3-true'] = 'true (CDR3 indels)'
# legends['v-indels'] = 'full partis (V indels)'
# legends['cdr3-indels'] = 'full partis (CDR3 indels)'
# args.is_data = False
# args.debug = False
# args.use_all_steps = False
# args.normalize_axes = []
# args.xbounds, args.adjmi_bounds, args.logprob_bounds = None, None, None

# input_info, reco_info = seqfileopener.get_seqfile_info('v-indels.csv', is_data=False)
# v_truehist = plotting.get_cluster_size_hist(utils.get_true_partition(reco_info))
# v_cpath = ClusterPath()
# v_cpath.readfile('v-indels-partitions-new.csv')

# input_info, reco_info = seqfileopener.get_seqfile_info('cdr3-indels.csv', is_data=False)
# cdr3_truehist = plotting.get_cluster_size_hist(utils.get_true_partition(reco_info))
# cdr3_cpath = ClusterPath()
# cdr3_cpath.readfile('cdr3-indels-partitions-new.csv')
# # print 'v-indels: %f' % v_cpath.adj_mi_at_max_logprob
# # print 'cdr3-indels: %f' % cdr3_cplot.adj_mi_at_max_logprob
# plotting.plot_cluster_size_hists(os.getenv('www') + '/partis/clustering/indel-performance.svg',
#                                  OrderedDict([
#                                      ['v-true', v_truehist],
#                                      ['v-indels', plotting.get_cluster_size_hist(v_cpath.partitions[v_cpath.i_best])],
#                                      ['cdr3-true', cdr3_truehist],
#                                      ['cdr3-indels', plotting.get_cluster_size_hist(cdr3_cpath.partitions[cdr3_cpath.i_best])]
#                                  ]),
#                                  title='%d leaves, %dx mutation, indels' % (10, 1), xmax=10*3.01)
# sys.exit()
# ----------------------------------------------------------------------------------------

if args.dataset == 'stanford':
    datadir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers'
    files = [datadir + '/' + f for f in os.listdir(datadir)]
elif args.dataset == 'adaptive':
    # files = []
    # datadirs = [ '/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/' + h for h in humans['adaptive'] ]
    # for datadir in datadirs:
    #     files += [ datadir + '/' + fname for fname in os.listdir(datadir) if '-M_merged.tsv.bz2' in fname ]  # if you switch to naive (N), be careful 'cause A is split in pieces
    files = [args.fsdir.replace('_output', 'data') + '/' + args.dataset + '/' + human + '/shuffled.csv' for human in humans[args.dataset]]

# ----------------------------------------------------------------------------------------
for datafname in files:
    if args.dataset == 'stanford':
        human = os.path.basename(datafname).replace('_Lineages.fasta', '')
    elif args.dataset == 'adaptive':
        human = re.findall('[ABC]', datafname)[0]

    if args.humans is not None and human not in args.humans:
        continue

    label = human
    if args.extra_label_str is not None:
        label += '-' + args.extra_label_str
    if args.bak:
        label += '.bak'

    print 'run', human
    n_leaves, mut_mult = None, None  # values if we're runing on data
    for action in args.actions:
        if action == 'write-plots' or action == 'compare-subsets':
            continue
        print '      ----> ', action
        if action == 'cache-data-parameters':
            compareutils.execute(args, action, datafname, label, n_leaves, mut_mult)
            continue
        for n_leaves in args.n_leaf_list:
            print '  ----> ', n_leaves, ' leaves'
            for mut_mult in args.mutation_multipliers:
                print '         ----> mutate', mut_mult
                compareutils.execute(args, action, datafname, label, n_leaves, mut_mult)
                # sys.exit()

    if 'write-plots' in args.actions:
        compareutils.write_all_plot_csvs(args, label)
    if 'compare-subsets' in args.actions:
        compareutils.compare_all_subsets(args, label)

    # break
