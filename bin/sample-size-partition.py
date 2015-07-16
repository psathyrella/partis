#!/usr/bin/env python
import sys
import time
import os
from subprocess import Popen, PIPE, check_call
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument('--action', required=True)
args = parser.parse_args()

istart = 0
for n_queries in [4, 5, 8, 15, 20, 25, 40, 50, 60, 75, 100, 125, 150, 175, 200]:
# for n_queries in (10, 25, 100, 200, 300, 400, 500, 750, 1000, 2000):
    cmd = './bin/compare-partition-methods.py --actions ' + args.action
    istop = istart + n_queries
    # print '  %d queries from %d --> %d' % (n_queries, istart, istop)
    cmd += ' --istartstop ' + str(istart) + ':' + str(istop)
    print cmd
    check_call(cmd.split())
    # time.sleep(0.1)
    istart = istop
