#!/usr/bin/env python
import os
import argparse
import sys
import re
import csv
from subprocess import Popen
sys.path.insert(1, './python')
csv.field_size_limit(sys.maxsize)
from humans import humans
import utils
import compareutils
import glob

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--fsdir', default='/fh/fast/matsen_e/' + os.getenv('USER') + '/work/partis-dev/_output')
parser.add_argument('--mutation-multipliers', default='1')
parser.add_argument('--data', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--expected-methods', default='vollmers-0.9:mixcr:changeo:vsearch-partition:naive-hamming-partition:partition')
parser.add_argument('--synthetic-partitions', default='misassign-0.10-singletons:misassign-0.10-reassign:misassign-0.90-singletons:misassign-0.90-reassign')
parser.add_argument('--indels', action='store_true')
parser.add_argument('--lonely-leaves', action='store_true')
parser.add_argument('--mimic', action='store_true')
parser.add_argument('--extra-label-str')
parser.add_argument('--bak', action='store_true')
parser.add_argument('--count-distances', action='store_true')
parser.add_argument('--n-leaf-list', default='10')
parser.add_argument('--hfrac-bound-list')
parser.add_argument('--subset', type=int)
parser.add_argument('--n-to-partition', type=int, default=5000)
parser.add_argument('--n-data-to-cache', type=int, default=100000)
parser.add_argument('--n-sim-seqs', type=int, default=10000)
parser.add_argument('--n-subsets', type=int)
parser.add_argument('--istartstop')  # NOTE usual zero indexing
parser.add_argument('--istartstoplist')  # list of istartstops for comparisons
parser.add_argument('--plot-mean-of-subsets', action='store_true')
# parser.add_argument('--dont-normalize', action='store_true')
# parser.add_argument('--logaxis', action='store_true')
# parser.add_argument('--zoom', action='store_true')
parser.add_argument('--humans', default=None)  #'A')
parser.add_argument('--no-similarity-matrices', action='store_true')
all_actions = ['cache-data-parameters', 'simulate', 'cache-simu-parameters', 'partition', 'naive-hamming-partition', 'vsearch-partition', 'run-viterbi', 'run-changeo', 'run-mixcr', 'run-igscueal', 'synthetic', 'write-plots', 'compare-subsets']
parser.add_argument('--actions', required=True)  #, choices=all_actions)  #default=':'.join(all_actions))
args = parser.parse_args()
args.actions = utils.get_arg_list(args.actions)
print 'TODO fix intify/floatify'
print 'TODO get it to stop subsetting the files every time through'
args.mutation_multipliers = utils.get_arg_list(args.mutation_multipliers, floatify=True)
args.n_leaf_list = utils.get_arg_list(args.n_leaf_list, intify=True)
args.istartstop = utils.get_arg_list(args.istartstop, intify=True)
args.istartstoplist = utils.get_arg_list(args.istartstoplist, intify=True, list_of_pairs=True)
args.humans = utils.get_arg_list(args.humans)
args.hfrac_bound_list = utils.get_arg_list(args.hfrac_bound_list, floatify=True, list_of_pairs=True)
args.expected_methods = utils.get_arg_list(args.expected_methods)
args.synthetic_partitions = utils.get_arg_list(args.synthetic_partitions)

if 'cache-data-parameters' in args.actions:
    args.data = True

assert args.subset is None or args.istartstop is None  # dosn't make sense to set both of them

if args.subset is not None:
    if 'write-plots' not in args.actions:
        assert args.n_subsets == 10  # for all the subset plots, I split into ten subsets, then ended up only using the first thre of 'em, so you have to set n_subsets to 10 if you're running methods, but then to 3 when you're writing plots
if args.istartstop is not None: 
    args.n_to_partition = args.istartstop[1] - args.istartstop[0]

# ----------------------------------------------------------------------------------------
# compareutils.FOOP()
# sys.exit()

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
procs = []
for human in args.humans:
    if human in humans['stanford']:
        datadir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers'
        datafname = glob.glob(datadir + '/*' + human + '*')[0]  # should throw an index error if length is less than one... but still, this is hackey
    elif human in humans['adaptive']:
        datafname = args.fsdir.replace('_output', 'data') + '/adaptive/' + human + '/shuffled.csv'
    else:
        assert False

    label = human
    if args.extra_label_str is not None:
        label += '-' + args.extra_label_str
    if args.bak:
        label += '.bak'

    print 'run', human
    n_leaves, mut_mult = None, None  # values if we're runing on data
    if not args.data:
        if args.hfrac_bound_list is None:
            parameterlist = [{'n_leaves' : nl, 'mut_mult' : mm, 'hfrac_bounds' : None} for nl in args.n_leaf_list for mm in args.mutation_multipliers]
        else:
            parameterlist = [{'n_leaves' : nl, 'mut_mult' : mm, 'hfrac_bounds' : hbs} for nl in args.n_leaf_list for mm in args.mutation_multipliers for hbs in args.hfrac_bound_list]

    for action in args.actions:
        if action == 'write-plots' or action == 'compare-subsets':
            continue
        print ' ', action
        if action == 'cache-data-parameters':
            compareutils.execute(args, action, datafname, label, n_leaves, mut_mult, procs)
            continue

        for params in parameterlist:
            compareutils.execute(args, action, datafname, label, params['n_leaves'], params['mut_mult'], procs, params['hfrac_bounds'])

    if 'write-plots' in args.actions:
        compareutils.write_all_plot_csvs(args, label, parameterlist, datafname)
    if 'compare-subsets' in args.actions:
        compareutils.compare_subsets(args, label)

if len(procs) > 0:
    exit_codes = [p.wait() for p in procs]
    print 'exit ', exit_codes
