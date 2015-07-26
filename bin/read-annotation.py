#!/usr/bin/env python
from subprocess import check_call, Popen
import sys
sys.path.insert(1, './python')
import os
import argparse
import csv

import utils

parser = argparse.ArgumentParser()
parser.add_argument('--infname', required=True)
parser.add_argument('--datadir', default='data/imgt')
# parser.add_argument('--is-data', action='store_true')
args = parser.parse_args()

germline_seqs = utils.read_germlines(args.datadir)
# if not args.is_data:
#     reco_info = 

with open(args.infname) as infile:
    reader = csv.DictReader(infile)
    for line in reader:
        utils.process_input_line(line,  # NOTE duplicates code in read_annotation_output() in partitiondriver
                                 splitargs=('unique_ids', 'seqs'),
                                 int_columns=('nth_best', 'v_5p_del', 'd_5p_del', 'cdr3_length', 'j_5p_del', 'j_3p_del', 'd_3p_del', 'v_3p_del'),
                                 float_columns=('logprob'))
        for iseq in range(len(line['seqs'])):
            tmpline = dict(line)
            tmpline['seq'] = line['seqs'][iseq]
            ack fack
        utils.print_reco_event(germline_seqs, line)
