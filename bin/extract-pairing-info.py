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
import hist

parser = argparse.ArgumentParser()
parser.add_argument('infname')
parser.add_argument('outfname')
parser.add_argument('--droplet-id-separator', default='_', help='everything in the sequence id before this character is treated as the droplet id, e.g. for the default, the uid AAACGGGCAAGCGAGT-1_contig_2 has a droplet id of AAACGGGCAAGCGAGT-1')
args = parser.parse_args()

seqfos = utils.read_fastx(args.infname)
droplet_ids = {}
for sfo in seqfos:
    did = utils.get_droplet_id(sfo['name'])
    if did not in droplet_ids:
        droplet_ids[did] = []
    droplet_ids[did].append(sfo['name'])

print '  read %d sequences with %d droplet ids' % (len(seqfos), len(droplet_ids))
count_info = {}
for dlist in droplet_ids.values():
    if len(dlist) not in count_info:
        count_info[len(dlist)] = 0
    count_info[len(dlist)] += 1
print '    contigs per'
print '      droplet     count   fraction'
total = sum(count_info.values())
for size, count in sorted(count_info.items(), key=operator.itemgetter(0)):
    frac = count / float(total)
    print '       %2d        %5d     %s' % (size, count, ('%.3f'%frac) if frac > 0.001 else ('%.1e'%frac))

metafos = {}
for sfo in seqfos:
    if sfo['name'] in metafos:
        raise Exception('duplicate uid \'%s\' in %s' % (sfo['name'], args.infname))
    metafos[sfo['name']] = {'paired-uids' : [u for u in droplet_ids[utils.get_droplet_id(sfo['name'])] if u != sfo['name']]}

if not os.path.exists(os.path.dirname(args.outfname)):
    os.makedirs(os.path.dirname(args.outfname))
with open(args.outfname, 'w') as outfile:
    # yaml.dump(metafos, outfile, Dumper=yaml.CDumper)
    json.dump(metafos, outfile)
