#!/usr/bin/env python
import os
import argparse
import sys
import re
import csv
from subprocess import Popen
csv.field_size_limit(sys.maxsize)

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')

import humans
import utils
import compareutils
import glob


# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--gldir', default=partis_dir + '/data/germlines/human', help='Directory from which to read non-sample-specific information (e.g. germline genes)')
parser.add_argument('--fsdir', default='/fh/fast/matsen_e/' + os.getenv('USER') + '/work/partis-dev/_output/update-17')  # remove the update-17 bit for stuff from when I was actually preparing the paper
parser.add_argument('--mutation-multipliers', default='1')
parser.add_argument('--is-simu', action='store_true')
parser.add_argument('--print-metrics', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--expected-methods', default='vollmers-0.9:mixcr:changeo:vsearch-partition:naive-hamming-partition:partition')
parser.add_argument('--synthetic-partitions', default='distance-0.03:0.10-reassign:0.60-singletons') # 0.75-singletons
parser.add_argument('--indels', action='store_true')
parser.add_argument('--indel-location')
parser.add_argument('--lonely-leaves', action='store_true')
parser.add_argument('--mimic', action='store_true')
parser.add_argument('--box', action='store_true')
parser.add_argument('--zipf', action='store_true')
parser.add_argument('--extra-label-str')
parser.add_argument('--bak', action='store_true')
parser.add_argument('--count-distances', action='store_true')
parser.add_argument('--n-leaf-list', default='10')
parser.add_argument('--hfrac-bound-list')
parser.add_argument('--subset', type=int)
parser.add_argument('--n-max-queries')
parser.add_argument('--n-subsets', type=int)
parser.add_argument('--istartstop')  # NOTE usual zero indexing
parser.add_argument('--istartstoplist')  # list of istartstops for comparisons
parser.add_argument('--plot-mean-of-subsets', action='store_true')
parser.add_argument('--humans', required=True)  #'A')
parser.add_argument('--no-similarity-matrices', action='store_true')
parser.add_argument('--no-slurm', action='store_true')
parser.add_argument('--dry-run', action='store_true')
parser.add_argument('--seed-cluster-bounds', default='20:30')
parser.add_argument('--iseed')
parser.add_argument('--old-output-structure', action='store_true', help='output paths corresponding to clustering paper, i.e. with everything in /fh/fast/matsen_e/dralph')
all_actions = ['cache-parameters', 'simulate', 'vjcdr3-partition', 'partition', 'naive-hamming-partition', 'vsearch-partition', 'seed-partition', 'seed-naive-hamming-partition', 'run-viterbi', 'run-changeo', 'run-mixcr', 'run-igscueal', 'synthetic', 'write-plots', 'compare-subsets', 'annotate-seed-clusters']
parser.add_argument('--actions', required=True, choices=all_actions)  #default=':'.join(all_actions))
args = parser.parse_args()
args.actions = utils.get_arg_list(args.actions)
args.mutation_multipliers = utils.get_arg_list(args.mutation_multipliers, floatify=True)
args.n_leaf_list = utils.get_arg_list(args.n_leaf_list, floatify=True)
args.istartstop = utils.get_arg_list(args.istartstop, intify=True)
args.istartstoplist = utils.get_arg_list(args.istartstoplist, intify=True, list_of_pairs=True)
args.humans = utils.get_arg_list(args.humans)
args.hfrac_bound_list = utils.get_arg_list(args.hfrac_bound_list, floatify=True, list_of_pairs=True)
args.expected_methods = utils.get_arg_list(args.expected_methods)
args.synthetic_partitions = utils.get_arg_list(args.synthetic_partitions)
for isp in range(len(args.synthetic_partitions)):  # I really shouldn't have set it up this way
    args.synthetic_partitions[isp] = 'misassign-' + args.synthetic_partitions[isp]
args.seed_cluster_bounds = utils.get_arg_list(args.seed_cluster_bounds, intify=True)

assert args.subset is None or args.istartstop is None  # dosn't make sense to set both of them

if args.subset is not None:
    if 'write-plots' not in args.actions:
        assert args.n_subsets == 10  # for all the subset plots, I split into ten subsets, then ended up only using the first thre of 'em, so you have to set n_subsets to 10 if you're running methods, but then to 3 when you're writing plots
if args.istartstop is not None: 
    args.n_to_partition = args.istartstop[1] - args.istartstop[0]

if args.bak:
    args.fsdir = args.fsdir.replace('_output', '_output.bak')

if 'simulate' in args.actions:
    if args.n_max_queries is None:
        raise Exception('you have to specify a number of simulated sequences (use --n-max-queries)')
    if not args.is_simu:
        print 'Autoforcing --is-simu, since actions include \'simulate\''
        args.is_simu = True

# # ----------------------------------------------------------------------------------------
# compareutils.FOOP()
# sys.exit()

# ----------------------------------------------------------------------------------------
procs = []
for human in args.humans:

    if args.old_output_structure:
        if human in humans.humans['vollmers']:
            datadir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers'
            datafname = glob.glob(datadir + '/*' + human + '*')[0]  # should throw an index error if length is less than one... but still, this is hackey
        elif human in humans.humans['adaptive']:
            datafname = args.fsdir.replace('_output', 'data') + '/adaptive/' + human + '/shuffled.csv'
        else:
            assert False
    else:
        datafname = humans.get_datafname(human)
        if args.n_max_queries is None:
            if args.istartstop is not None:
                args.n_max_queries = args.istartstop[1] - args.istartstop[0]
            else:
                args.n_max_queries = humans.get_nseqs(human)

    label = human
    if args.extra_label_str is not None:
        label += '-' + args.extra_label_str

    print 'run', human
    n_leaves, mut_mult = None, None  # values if we're runing on data
    parameterlist = [{'n_leaves' : None, 'mut_mult' : None, 'hfrac_bounds' : None}]
    if args.is_simu:
        if args.hfrac_bound_list is None:
            parameterlist = [{'n_leaves' : nl, 'mut_mult' : mm, 'hfrac_bounds' : None} for nl in args.n_leaf_list for mm in args.mutation_multipliers]
        else:
            parameterlist = [{'n_leaves' : nl, 'mut_mult' : mm, 'hfrac_bounds' : hbs} for nl in args.n_leaf_list for mm in args.mutation_multipliers for hbs in args.hfrac_bound_list]

    for action in args.actions:
        if action == 'write-plots' or action == 'compare-subsets':
            continue
        print ' ', action
        if action == 'cache-parameters' and not args.is_simu:
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
