#!/usr/bin/env python
import sys
import os
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import argparse
current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
if not os.path.exists(current_script_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % current_script_dir
sys.path.insert(1, current_script_dir)

from seqfileopener import get_seqfile_info
import utils

parser = argparse.ArgumentParser()
parser.add_argument('infname')
parser.add_argument('--datadir', default='data/imgt')
parser.add_argument('--simfname')
parser.add_argument('--is-data', action='store_true')
parser.add_argument('--print-uid', action='store_true')
args = parser.parse_args()

glfo = utils.read_germline_set(args.datadir)

reco_info = None
if args.simfname is not None:
    _, reco_info = get_seqfile_info(args.simfname, args.is_data, glfo=glfo)

with open(args.infname) as infile:
    reader = csv.DictReader(infile)
    for line in reader:
        utils.process_input_line(line)
        if args.simfname is not None:
            utils.print_true_events(glfo, reco_info, line, print_uid=args.print_uid)
        utils.add_implicit_info(glfo, line, multi_seq=True, existing_implicit_keys=('aligned_d_seqs', 'aligned_j_seqs', 'aligned_v_seqs', 'cdr3_length', 'naive_seq', 'in_frames', 'mutated_invariants', 'stops'))
        print 'inferred:\n'
        utils.print_reco_event(glfo['seqs'], line, print_uid=args.print_uid)
