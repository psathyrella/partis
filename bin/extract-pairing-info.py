#!/usr/bin/env python
import csv
import os
import sys
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import argparse
import colored_traceback.always
import yaml
import json
import operator

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils

dstr = """
Extract heavy/light chain pairing info from fasta file <infname> and write it to yaml/json file <outfname>.
"""
parser = argparse.ArgumentParser(description=dstr,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)  # why tf isn't this printing the defaults?
parser.add_argument('infname')
parser.add_argument('outfname')
parser.add_argument('--droplet-id-separators', help=utils.did_help['seps'])
parser.add_argument('--droplet-id-indices', help=utils.did_help['indices'])
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--n-max-queries', type=int, default=-1, help='Maximum number of query sequences to read from input file, starting from beginning of file')
parser.add_argument('--input-metafname', help='json/yaml file with additional (beyond pairing info) input meta info (see partis help)')
args = parser.parse_args()
args.droplet_id_indices = utils.get_arg_list(args.droplet_id_indices, intify=True)

if utils.output_exists(args, args.outfname, offset=4, debug=False):
    print '  extract-pairing-info.py output exists and --overwrite was not set, so not doing anything: %s' % args.outfname
    sys.exit(0)

seqfos = utils.read_fastx(args.infname, n_max_queries=args.n_max_queries)
metafos = utils.extract_pairing_info(seqfos, droplet_id_separators=args.droplet_id_separators, droplet_id_indices=args.droplet_id_indices, input_metafname=args.input_metafname)

utils.mkdir(args.outfname, isfile=True)
with open(args.outfname, 'w') as outfile:
    json.dump(metafos, outfile)
