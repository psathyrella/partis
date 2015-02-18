#!/usr/bin/env python
import sys
import argparse
from subprocess import check_call

parser = argparse.ArgumentParser()
parser.add_argument('--leaves', type=int, required=True)
args = parser.parse_args()

cmd = './bin/partis.py'

events = 50  #int(100./args.leaves)
cmdstr = cmd + ' --simulate --outfname _output/viterbi-cluster-' + str(args.leaves) + '/simu.csv --n-procs 10 --parameter-dir _output/every-10-A-subset-0/data/hmm --n-max-queries ' + str(events) + ' --n-leaves ' + str(args.leaves)
print cmdstr
check_call(cmdstr.split())

cmdstr = cmd + ' --run-algorithm viterbi --n-procs 10'
cmdstr += ' --seqfile _output/viterbi-cluster-' + str(args.leaves) + '/simu.csv'
cmdstr += ' --parameter-dir _output/every-10-A-subset-0/simu/hmm'
cmdstr += ' --no-plot --pants-seated-clustering'
print cmdstr
check_call(cmdstr.split())

#               
# leaves  |   true events | inferred events |    ratio
# ------- | --------------| --------------- | -------
#   5     |      100      |       260       |     0.38
#  10     |      100      |       449       |     0.22
#  20     |      100      |       774       |     0.13
#  50     |       20      |       233       |     0.085
# 100     |       20      |       511       |     0.04
# 200     |       50      |      2713       |     0.018
# 300     |       50      |      3058       |     0.016
# 500     |       50      |      6693       |     0.0075
# 600     |       50      |      6441       |     0.0078
# 700     |       50      |      6328       |     0.0079
