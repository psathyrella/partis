#!/usr/bin/env python
import argparse
import sys
import os

from utils import utils
from partitiondriver import PartitionDriver

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
parser.add_argument('--no_clean', action='store_true')  # don't remove the various temp files
# actions:
parser.add_argument('--cache_parameters', action='store_true')  # cache parameter counts and hmm files in dir specified by <parameter_dir>
parser.add_argument('--point_estimate', action='store_true')
parser.add_argument('--parameter_dir')

parser.add_argument('--n_bases_skip', type=int, default=0)  # number of bases to skip on the left side of the sequence

parser.add_argument('--human', default='A', choices=['A', 'B', 'C'])
parser.add_argument('--naivety', default='M', choices=['N', 'M'])
parser.add_argument('--seqfile')  # input file
parser.add_argument('--queries')  # restrict to certain query seqs
parser.add_argument('--reco_ids')  # or recombination events
parser.add_argument('--is_data', action='store_true')
parser.add_argument('--n_total_queries', type=int, default=-1)  # stop after this many queries
parser.add_argument('--only_genes')  # skip all gene matches except for these when parsing the sw output
parser.add_argument('--pair', action='store_true')

parser.add_argument('--n_max_per_region', type=int, default=3)
parser.add_argument('--n_best_events', type=int, default=5)
parser.add_argument('--default_v_fuzz', type=int, default=2)  # TODO play around with these default fuzzes
parser.add_argument('--default_d_fuzz', type=int, default=2)
parser.add_argument('--v_right_length', type=int, default=89)  # length of v gene starting from junction with d (used mainly for writing hmms model files)
parser.add_argument('--datadir', default='./data')  # non-sample-specific information
parser.add_argument('--stochhmm_dir', default=os.getenv('HOME') + '/work/StochHMM')
parser.add_argument('--workdir', default='/tmp/' + os.getenv('USER') + '/hmms/' + str(os.getpid()))

args = parser.parse_args()
if args.queries != None:
    args.queries = args.queries.strip().split(':')  # to allow ids with minus signs, need to add a space, which you then have to strip() off
if args.reco_ids != None:
    args.reco_ids = args.reco_ids.strip().split(':')
if args.only_genes != None:
    args.only_genes = args.only_genes.strip().split(':')

utils.prep_dir(args.workdir)
# ----------------------------------------------------------------------------------------
parter = PartitionDriver(args)
if args.cache_parameters:
    parter.cache_parameters()
elif args.point_estimate:
    parter.point_estimate()
else:
    parter.partition()
