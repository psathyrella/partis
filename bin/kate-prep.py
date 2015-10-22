#!/usr/bin/env python
import csv
import sys
from collections import OrderedDict
# from subprocess import PIPE, check_call, check_output, CalledProcessError
import argparse

# ----------------------------------------------------------------------------------------
# questions
#   1. email says 300 nt total, but looks like there's 299?

parser = argparse.ArgumentParser()
parser.add_argument('--infname', default='kate-test.txt')  #s_1_1_2109_qseq-randomized.txt')
parser.add_argument('--read', type=int, choices=[1, 2], required=True)
parser.add_argument('--datadir', default='data/kate')  # where to get primer and "sample sheet"
args = parser.parse_args()

# ----------------------------------------------------------------------------------------
# hardcoded crap

adapters = { 1 : 'CCGGGAGCTGCATGTGTCAGAGG',
             2 : 'GGGCTGGCAAGCCACGTTTGGTG'}

# name : (istart, istop)
readinfo = OrderedDict()
readinfo[1] = OrderedDict()
readinfo[1]['random_id'] = (0, 5)
readinfo[1]['barcode'] = (6, 13)
readinfo[1]['adapter_1'] = (14, 36)
readinfo[1]['c_primer'] = (38, None)  # variable endpoint
readinfo[1]['antisense_c'] = (None, 298)
readinfo[2] = OrderedDict()
readinfo[2]['barcode'] = (0, 7)
readinfo[2]['random_id'] = (8, 13)
readinfo[2]['adapter_2'] = (14, 36)
readinfo[2]['v_primer'] = (37, None)  # variable endpoint
readinfo[2]['sense_leader'] = (None, 298)

# ----------------------------------------------------------------------------------------
def process_sequence(seq):
    """ pull off the primer and whatnot """
    seqfo = {}
    for name, (istart, istop) in readinfo[args.read].items():
        if 'primer' in name:

# ----------------------------------------------------------------------------------------
with open(args.infname) as infile:
    for line in infile:
        lineinfo = line.split()
        fullseq = lineinfo[8]
        seqfo = process_sequence(fullseq)
        sys.exit()
