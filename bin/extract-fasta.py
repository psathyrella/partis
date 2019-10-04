#!/usr/bin/env python
import csv
import os
import sys
import argparse
import colored_traceback.always

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
from clusterpath import ClusterPath

helpstr = """
Script to extract sequences from a partis output file and write them to a fasta file.
For details of partis output files, see the manual.
Example usage:
    ./bin/extract-fasta.py --input-file partis-output.yaml --fasta-output-file out.fa  # extact all sequences from best partition in <partis-output.yaml>
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('--input-file', required=True, help='partis output file')
parser.add_argument('--fasta-output-file', required=True, help='output fasta file name')
parser.add_argument('--partition-index', type=int, help='if set, use the partition at this index in the cluster path, rather than the default of using the best partition')
parser.add_argument('--seed-unique-id', help='if set, take sequences only from the cluster containing this seed sequence, rather than the default of taking all sequences from all clusters')
parser.add_argument('--cluster-index', type=int, help='if set, take sequences only from the cluster at this index in the partition, rather than the default of taking all sequences from all clusters')
parser.add_argument('--indel-reversed-seqs', action='store_true', help='if set, take sequences that have had any shm indels "reversed" (i.e. insertions are reversed, and deletions are replaced with the germline bases) rather than the default of using sequences from the original input file. Indel-reversed sequences can be convenient because they are by definition the same length as and aligned to the naive sequence.')
parser.add_argument('--glfo-dir', default=partis_dir + '/data/germlines/human', help='Directory with germline info. Only necessary for old-style csv output files.')
parser.add_argument('--locus', default='igh', help='only used for old-style csv output files')
args = parser.parse_args()

glfo = None
if utils.getsuffix(args.input_file) == '.csv':
    print '  reading deprecated csv format, so need to read germline info from somewhere else, using --glfo-dir %s, hopefully it works (you can set it to the proper thing with --glfo-dir)' % args.glfo_dir
    glfo = glutils.read_glfo(args.glfo_dir, locus=args.locus)

glfo, annotation_list, cpath = utils.read_output(args.input_file, glfo=glfo)

if cpath is None:
    clusters_to_use = [l['unique_ids'] for l in annotation_list]
    print '  no cluster path in input file, so just using all %d sequences (in %d clusters) in annotations' % (sum(len(c) for c in clusters_to_use), len(clusters_to_use))
else:
    ipartition = cpath.i_best if args.partition_index is None else args.partition_index
    print '  found %d clusters in %s' % (len(cpath.partitions[ipartition]), 'best partition' if args.partition_index is None else 'partition at index %d (of %d)' % (ipartition, len(cpath.partitions)))
    if args.cluster_index is None:
        clusters_to_use = cpath.partitions[ipartition]
        print '    taking all %d clusters' % len(clusters_to_use)
    else:
        clusters_to_use = [cpath.partitions[ipartition][args.cluster_index]]
        print '    taking cluster at index %d' % args.cluster_index
    if args.seed_unique_id is not None:
        clusters_to_use = [c for c in clusters_to_use if args.seed_unique_id in c]  # NOTE can result in more than one cluster with the seed sequence (e.g. if this file contains intermediate annotations from seed partitioning))
        print '    removing clusters not containing sequence \'%s\' (leaving %d)' % (args.seed_unique_id, len(clusters_to_use))

seqfos = []
annotations = {':'.join(adict['unique_ids']) : adict for adict in annotation_list}  # collect the annotations in a dictionary so they're easier to access
for cluster in clusters_to_use:
    if ':'.join(cluster) not in annotations:
        print '  %s cluster with size %d not in annotations, so skipping it' % (utils.color('red', 'warning'), len(cluster))
        continue
    cluster_annotation = annotations[':'.join(cluster)]
    seqfos += [{'name' : u, 'seq' : s} for u, s in zip(cluster_annotation['unique_ids'], cluster_annotation['seqs' if args.indel_reversed_seqs else 'input_seqs'])]

if not os.path.exists(os.path.dirname(os.path.abspath(args.fasta_output_file))):
    os.makedirs(os.path.dirname(os.path.abspath(args.fasta_output_file)))
print '  writing %d sequences to %s' % (len(seqfos), args.fasta_output_file)
with open(args.fasta_output_file, 'w') as outfile:
    for sfo in seqfos:
        outfile.write('>%s\n%s\n' % (sfo['name'], sfo['seq']))
