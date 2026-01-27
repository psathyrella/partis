#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import csv
import os
import sys
from io import open
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import argparse
import colored_traceback.always
import yaml
import json
import operator
import random
import numpy
from pathlib import Path

import partis.utils as utils

dstr = """
Extract heavy/light chain pairing info from fasta file <infname> and write it to yaml/json file <outfname>.
Should have the same effect as setting --guess-pairing-info when running bin/split-loci.py.
"""
parser = argparse.ArgumentParser(description=dstr,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)  # why tf isn't this printing the defaults?
parser.add_argument('infname')
parser.add_argument('outfname')
parser.add_argument('--droplet-id-separators', help=utils.did_help['seps'])
parser.add_argument('--droplet-id-indices', help=utils.did_help['indices'])
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--for-testing-n-max-queries', type=int, default=-1, help='only for testing, applied when reading initial fasta file, just in case it\'s huge and you want to run quickly without having to read the whole file')
parser.add_argument('--n-max-queries', type=int, default=-1, help='see partis help (although here it applies to droplets, not individual seqs)')
parser.add_argument('--n-random-queries', type=int, help='see partis help (although here it applies to droplets, not individual seqs)')
parser.add_argument('--input-metafname', help='json/yaml file with additional (beyond pairing info) input meta info (see partis help)')
parser.add_argument('--random-seed', type=int, default=1)
args = parser.parse_args()
random.seed(args.random_seed)
numpy.random.seed(args.random_seed)
args.droplet_id_indices = utils.get_arg_list(args.droplet_id_indices, intify=True)

if utils.output_exists(args, args.outfname, offset=4, debug=False):
    print('  extract-pairing-info.py output exists and --overwrite was not set, so not doing anything: %s' % args.outfname)
    sys.exit(0)

seqfos = utils.read_fastx(args.infname, n_max_queries=args.for_testing_n_max_queries)
if args.n_max_queries != -1 or args.n_random_queries is not None:
    seqfos = utils.subset_paired_queries(seqfos, args.droplet_id_separators, args.droplet_id_indices, n_max_queries=args.n_max_queries, n_random_queries=args.n_random_queries)
metafos = utils.extract_pairing_info(seqfos, droplet_id_separators=args.droplet_id_separators, droplet_id_indices=args.droplet_id_indices, input_metafname=args.input_metafname)

utils.mkdir(args.outfname, isfile=True)
utils.jsdump(args.outfname, metafos)
