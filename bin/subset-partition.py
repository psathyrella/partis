#!/usr/bin/env python
import time
import sys
sys.path.insert(1, 'python/')
import argparse
from subprocess import check_call, Popen, PIPE

parser = argparse.ArgumentParser()
parser.add_argument('--action', required=True)
args = parser.parse_args()

n_subsets = 10
for isub in range(1, n_subsets):
    cmd = './bin/compare-partition-methods.py --actions ' + args.action + ' --subset ' + str(isub) + ' --n-subsets ' + str(n_subsets)
    cmd += ' --mutation-multipliers 1:4'
    cmd += ' --n-leaf-list 1:5:10:25:50'
    cmd += ' --overwrite'
    # cmd += ' --indels'
    # cmd += ' --lonely-leaves'
    print cmd
    Popen(cmd.split())
    # check_call(cmd.split())
    time.sleep(60)
    # sys.exit()
