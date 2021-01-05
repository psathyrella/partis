#!/usr/bin/env python
import csv
import os
import sys
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import argparse
import colored_traceback.always
import json
import numpy

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
from clusterpath import ClusterPath

parser = argparse.ArgumentParser()
parser.add_argument('--infiles', nargs='+')
parser.add_argument('--outfile', required=True)
parser.add_argument('--metafile', required=True)
parser.add_argument('--ig-or-tr', default='ig')
parser.add_argument('--mean-cells-per-droplet', type=float)
args = parser.parse_args()

# ----------------------------------------------------------------------------------------
# TODO this should go somewhere else, like maybe in recombinator, bin/partis, or utils
# ----------------------------------------------------------------------------------------

outfos, metafos = [], {}
# for ltmp in utils.sub_loci(args.ig_or_tr):
#     ifn = utils.get_single_entry([f for f in args.infiles if '-'+ltmp in f])
for ifn in args.infiles:
    ltmp = utils.getprefix(os.path.basename(ifn)).split('-')[0]
    assert ltmp in utils.loci
    _, annotation_list, _ = utils.read_output(ifn)
    for line in annotation_list:
        for uid, seq, pids in zip(line['unique_ids'], line['input_seqs'], line['paired-uids']):
            outfos.append({'name' : uid, 'seq' : seq})
            metafos[uid] = {'locus' : ltmp, 'paired-uids' : pids}

if args.mean_cells_per_droplet is not None:  # NOTE that this modifies the input meta info file, which is where partis will then get the paird uid info, but does *not* modify the 'paired-uids' key in the original simulation files
    n_droplets = int(0.5 * float(len(outfos)) / args.mean_cells_per_droplet)  # 0.5 is because <outfos> includes both heavy and light sequences
    droplet_ids = [[] for _ in range(n_droplets)]
    sfo_dict = {s['name'] : s for s in outfos}
    while len(sfo_dict) > 0:
        tid = next(iter(sfo_dict))
        idrop = numpy.random.choice(range(len(droplet_ids)))
        droplet_ids[idrop] += [tid] + metafos[tid]['paired-uids']  # add <tid> plus its paired ids to this drop
        for uid in [tid] + metafos[tid]['paired-uids']:
            sfo_dict[uid]['droplet-ids'] = droplet_ids[idrop]
            del sfo_dict[uid]
    for sfo in outfos:
        metafos[sfo['name']]['paired-uids'] = sfo['droplet-ids']

with open(args.outfile, 'w') as outfile:
    for sfo in outfos:
        outfile.write('>%s\n%s\n' % (sfo['name'], sfo['seq']))
with open(args.metafile, 'w') as mfile:
    json.dump(metafos, mfile)
