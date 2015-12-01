#!/usr/bin/env python
import sys
sys.path.insert(1, './python')
import csv
import argparse

from clusterpath import ClusterPath
from seqfileopener import get_seqfile_info
import utils

parser = argparse.ArgumentParser()
parser.add_argument('--infname', required=True)
parser.add_argument('--dont-abbreviate', action='store_true', help='Print full seq IDs (otherwise just prints an \'o\')')
parser.add_argument('--n-to-print', type=int, help='How many partitions to print (centered on the best partition)')
parser.add_argument('--datadir', default='data/imgt')
parser.add_argument('--simfname')
parser.add_argument('--is-data', action='store_true')
args = parser.parse_args()

germline_seqs = utils.read_germlines(args.datadir)
cyst_positions = utils.read_cyst_positions(args.datadir)
with open(args.datadir + '/j_tryp.csv') as csv_file:  # get location of <end> tryptophan in each j region
    tryp_reader = csv.reader(csv_file)
    tryp_positions = {row[0]:row[1] for row in tryp_reader}  # WARNING: this doesn't filter out the header line

reco_info = None
if args.simfname is not None:
    input_info, reco_info = get_seqfile_info(args.simfname, args.is_data, germline_seqs, cyst_positions, tryp_positions)

cp = ClusterPath()
cp.readfile(args.infname)
cp.print_partitions(abbreviate=(not args.dont_abbreviate), n_to_print=args.n_to_print, reco_info=reco_info)
