#!/usr/bin/env python
import sys
sys.path.insert(1, './python')
import argparse

from clusterpath import ClusterPath
from seqfileopener import get_seqfile_info
import utils

parser = argparse.ArgumentParser()
parser.add_argument('--infname', required=True)
parser.add_argument('--dont-abbreviate', action='store_true', help='Print full seq IDs (otherwise just prints an \'o\')')
parser.add_argument('--n-to-print', type=int, help='How many partitions to print (centered on the best partition)')
parser.add_argument('--simfname')
parser.add_argument('--is-data', action='store_true')
args = parser.parse_args()

self.germline_seqs = utils.read_germlines(args.datadir)
if args.simfname is not None:
    input_info, reco_info = get_seqfile_info(args.simfname, args.is_data, self.germline_seqs, self.cyst_positions, self.tryp_positions,
                                                       self.args.n_max_queries, self.args.queries, self.args.reco_ids)
cp = ClusterPath(-1)
cp.readfile(args.infname)
cp.print_partitions(one_line=True, abbreviate=(not args.dont_abbreviate), n_to_print=args.n_to_print)
