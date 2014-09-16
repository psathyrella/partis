#!/usr/bin/env python

import argparse
import sys

# from clusterer import Clusterer
# clust = Clusterer(10, debug=True)  #-300.0
# clust.cluster('/tmp/pairs.csv')
# sys.exit()

from partitiondriver import PartitionDriver

workdir = '/home/dralph/Dropbox/work'
datadir = workdir + '/partis/recombinator/data'
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
parser.add_argument('--v_right_length', type=int, default=89)
parser.add_argument('--pair', action='store_true')
parser.add_argument('--no_clean', action='store_true')
parser.add_argument('--write_all_hmm_files', action='store_true')
parser.add_argument('--write_hmm_files')
parser.add_argument('--human', default='A', choices=['A', 'B', 'C'])
parser.add_argument('--naivety', default='M', choices=['N', 'M'])
parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
args = parser.parse_args()
if args.queries != None:
    args.queries = args.queries.strip().split(':')  # to allow ids with minus signs, need to add a space, which you then have to strip() off
if args.reco_ids != None:
    args.reco_ids = args.reco_ids.strip().split(':')

# ----------------------------------------------------------------------------------------
if not args.simfile:
    args.simfile = workdir + '/partis/recombinator/output/' + args.human + '/' + args.naivety + '/simu.csv'
parter = PartitionDriver(datadir, args, default_v_right_length=89, stochhmm_dir=workdir + '/StochHMM')

if args.write_all_hmm_files:
    parter.write_all_hmms(args.v_right_length)
    sys.exit()
if args.write_hmm_files != None:
    args.write_hmm_files = args.write_hmm_files.split(':')
    parter.write_specified_hmms(args.write_hmm_files, args.v_right_length)
    sys.exit()

parter.run()
# IGHV3-23*04:IGHJ3*02_F:IGHJ4*02_F:IGHD3-22*01:IGHV3-64*04:IGHV3-72*01:IGHJ6*02_F:IGHJ2*01_F:IGHD6-19*01:IGHD4-23*01:IGHV1-18*01:IGHJ5*02_F:IGHD3-9*01:IGHD4-17*01:IGHD3-10*01:IGHV5-51*01
