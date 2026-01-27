#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import collections
import numpy
import random
import sys
import os
import csv
from io import open

import partis.utils as utils
import partis.seqfileopener as seqfileopener

parser = argparse.ArgumentParser()
parser.add_argument('infile')
parser.add_argument('outfile')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--chimera-freq', default=1., type=float, help='fraction of sequences to make chimeric')
parser.add_argument('--min-chunk-len', default=15, type=int, help='require that each bit of the chimera is at least this long')
args = parser.parse_args()

input_info, _, _ = seqfileopener.read_sequence_file(args.infile, is_data=False)
if len(input_info) < 50:
    print('%s making chimeras with only %d sequences, and since we choose from among the existing sequence for templates this won\'t be very effective' % (utils.color('yellow', 'warning'), len(input_info)))

n_chimeric = 0
outfo = collections.OrderedDict()
for uid, seqfo in input_info.items():
    if args.debug:
        print(uid)

    if numpy.random.uniform(0, 1) > args.chimera_freq:  # no chimeras for this sequence
        if args.debug:
            print('        non-chimeric')
        continue

    break_point = random.randint(args.min_chunk_len, len(seqfo['seqs'][0]) - args.min_chunk_len)
    switch_uid = numpy.random.choice(input_info)
    switch_seq = input_info[switch_uid]['seqs'][0][ : break_point]

    if args.debug:
        print('    switching to %s at %d:' % (switch_uid, break_point))
        print('          %s' % switch_seq)
        print('          %s%s' % (' ' * len(switch_seq), seqfo['seqs'][0][break_point : ]))

    outfo[uid] = switch_seq + seqfo['seqs'][0][break_point : ]
    n_chimeric += 1

print('writing %d / %d chimeric sequences to %s' % (n_chimeric, len(input_info), args.outfile))
with open(args.outfile, 'w') as outfile:
    for uid, seq in outfo.items():
        outfile.write('>%s\n%s\n' % (uid, seq))
