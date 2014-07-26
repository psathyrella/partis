#!/usr/bin/env python

import argparse
import sys
from partitiondriver import PartitionDriver    

# from clusterer import Clusterer
# clust = Clusterer()
# clust.cluster('pairwise-scores.csv')
# sys.exit()

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--simfile')  #, default='/home/dralph/Dropbox/work/recombinator/output/' + human + '/' + naivety + '/simu.csv')
parser.add_argument('--algorithm', default='viterbi', choices=['viterbi', 'forward'])
parser.add_argument('--n_max_per_region', type=int, default=5)
parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
parser.add_argument('--pair', action='store_true')
parser.add_argument('--write_hmms_on_fly', action='store_true', default=False, help='write, on the fly, only the hmm files we need? Or look for cached versions?')
parser.add_argument('--run_sw', action='store_true', default=True, help='run smith-waterman from scratch, or look for cached output?')
parser.add_argument('--human', default='A', choices=['A', 'B', 'D'])
parser.add_argument('--naivety', default='M', choices=['N', 'M'])
args = parser.parse_args()

# ----------------------------------------------------------------------------------------
human = 'A'
naivety = 'M'
datadir = '/home/dralph/Dropbox/work/recombinator/data'
parter = PartitionDriver(datadir, args, default_v_right_length=89)
# parter.write_all_hmms(89)
parter.run()
