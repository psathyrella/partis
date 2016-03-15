#!/usr/bin/env python
import sys
sys.path.insert(1, './python')
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import argparse

from seqfileopener import get_seqfile_info
import utils

parser = argparse.ArgumentParser()
parser.add_argument('infname')
parser.add_argument('--datadir', default='data/imgt')
# parser.add_argument('--simfname')
parser.add_argument('--is-data', action='store_true')
args = parser.parse_args()

glfo = utils.read_germline_set(args.datadir)

# reco_info = None
# if args.simfname is not None:
#     input_info, reco_info = get_seqfile_info(args.simfname, args.is_data, glfo=glfo)

with open(args.infname) as infile:
    reader = csv.DictReader(infile)
    for line in reader:
        utils.process_input_line(line)
        utils.add_implicit_info(glfo, line, multi_seq=True, existing_implicit_keys=('aligned_d_seqs', 'aligned_j_seqs', 'aligned_v_seqs', 'cdr3_length', 'naive_seq', 'in_frames', 'mutated_invariants', 'stops'))
        utils.print_reco_event(glfo['seqs'], line, print_uid=True)
