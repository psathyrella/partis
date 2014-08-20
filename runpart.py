#!/usr/bin/env python

import argparse
import sys

# from clusterer import Clusterer
# clust = Clusterer(10, debug=True)  #-300.0
# clust.cluster('/tmp/pairs.csv')
# sys.exit()

from partitiondriver import PartitionDriver

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--simfile')
# restrict to certain query seqs or recombination events
parser.add_argument('--queries')
parser.add_argument('--reco_ids')

parser.add_argument('--algorithm', default='viterbi', choices=['viterbi', 'forward'])
parser.add_argument('--n_max_per_region', type=int, default=5)
parser.add_argument('--n_best_events', type=int, default=5)
parser.add_argument('--n_total_queries', type=int, default=-1)
parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
parser.add_argument('--pair', action='store_true')
parser.add_argument('--no_clean', action='store_true')
parser.add_argument('--write_hmm_files', action='store_true')
parser.add_argument('--write_hmms_on_fly', action='store_true', default=False, help='write, on the fly, only the hmm files we need? Or look for cached versions?')
parser.add_argument('--skip_sw', action='store_true', default=False, help='look for cached output, or run smith-waterman from scratch?')
parser.add_argument('--human', default='A', choices=['A', 'B', 'C'])
parser.add_argument('--naivety', default='M', choices=['N', 'M'])
args = parser.parse_args()
if args.queries != None:
    args.queries = args.queries.strip().split(':')  # to allow ids with minus signs, need to add a space, which you then have to strip() off
if args.reco_ids != None:
    args.reco_ids = args.reco_ids.strip().split(':')

# ----------------------------------------------------------------------------------------
if not args.simfile:
    args.simfile = '/home/dralph/Dropbox/work/recombinator/output/' + args.human + '/' + args.naivety + '/simu.csv'
datadir = '/home/dralph/Dropbox/work/recombinator/data'
parter = PartitionDriver(datadir, args, default_v_right_length=89)

if args.write_hmm_files:
    parter.write_all_hmms(89)
    sys.exit()

parter.run()
