#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import csv
import os
import sys
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import argparse
import operator
import colored_traceback.always
import collections

import partis.utils as utils
import partis.glutils as glutils
from partis.clusterpath import ClusterPath

dstr = """
Add seqs from the fasta file --new-seq-file to an annotation from --partis-output-file.
Looks for a cluster in the best partition that has sequences in common with the fasta file (and crashes if there's more than one such cluster).
Writes a single modified annotation to --outfile.
"""
parser = argparse.ArgumentParser(description=dstr,
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)  # why tf isn't this printing the defaults?
parser.add_argument('--new-seq-file', required=True, help='fasta input file with seqs to be added to annotations + partitions in partis output yaml')
parser.add_argument('--partis-output-file', required=True, help='partis output file to which to add the seqs from --new-seq-file')
parser.add_argument('--partition-index', type=int, help='index of partition from which to take the clusters/annotations (if not set, uses the best partition)')
parser.add_argument('--glfo-dir', default=partis_dir + '/data/germlines/human', help='germline info directory. Only used if --partis-output-file is an old-style .csv, and this default dir may work if your output file doesn\'t have novel inferred genes. Otherwise, is the germline info dir from the partis inferred parameter directory corresponding to your output file --partis-output-file.')
parser.add_argument('--locus', default='igh')
parser.add_argument('--outfile', required=True, help='output partis yaml file')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--n-test-subset-seqs', type=int, help='take only the first N seqs from both the fasta file and the annotation in the partis output file (e.g. for testing when the family is huge)')
args = parser.parse_args()

new_seqfos = utils.read_fastx(args.new_seq_file, sanitize_seqs=True)
print('    read %d seqs from %s' % (len(new_seqfos), args.new_seq_file))

glfo = None
if utils.getsuffix(args.partis_output_file) == '.csv':
    print('    reading deprecated csv format, so need to read germline info from somewhere else, using --glfo-dir %s, hopefully it works' % args.glfo_dir)
    glfo = glutils.read_glfo(args.glfo_dir, locus=args.locus)

glfo, annotation_list, cpath = utils.read_output(args.partis_output_file, glfo=glfo, locus=args.locus)
if args.partition_index is not None:
    print('  using non-best partition index %d (best is %d)' % (args.partition_index, cpath.i_best))
partition = cpath.partitions[cpath.i_best if args.partition_index is None else args.partition_index]
print('    read partition with %d clusters from %s' % (len(partition), args.partis_output_file))

new_uids = set(sfo['name'] for sfo in new_seqfos)
clusters_with_overlap = []
for cluster in partition:
    overlap_uids = set(cluster) & new_uids
    if len(overlap_uids) > 0:
        clusters_with_overlap.append((cluster, overlap_uids))

if len(clusters_with_overlap) == 0:
    raise Exception('no clusters in partition have any overlap with sequences from fasta file')
elif len(clusters_with_overlap) > 1:
    # raise Exception('too many clusters %d in the partition overlaps with sequences from the fasta file' % len(clusters_with_overlap))
    clusters_with_overlap = sorted(clusters_with_overlap, key=lambda p: len(p[1]), reverse=True)
    ostrs = ['%d %d'%(len(c), len(o)) for c, o in clusters_with_overlap]
    print('  %s more than one cluster overlaps with sequences from fasta file, just taking first one (size overlap): %s,  %s' % (utils.color('yellow', 'warning'), utils.color('red', ostrs[0]), ',  '.join(ostrs[1:])))
old_cluster = clusters_with_overlap[0][0]

print('    adding %d fasta sequences to cluster of size %d (%d fasta sequences were already in cluster)' % (len(new_uids - set(old_cluster)), len(old_cluster), len(new_uids & set(old_cluster))))
sfos_to_add = [sfo for sfo in new_seqfos if sfo['name'] not in old_cluster]
annotation_dict = utils.get_annotation_dict(annotation_list)
annotation = annotation_dict[':'.join(old_cluster)]

if args.n_test_subset_seqs is not None:
    print('  taking only first %d seqs from fasta and annotation' % args.n_test_subset_seqs)
    utils.restrict_to_iseqs(annotation, list(range(args.n_test_subset_seqs)), glfo)
    sfos_to_add = sfos_to_add[:args.n_test_subset_seqs]
utils.add_seqs_to_line(annotation, sfos_to_add, glfo, debug=args.debug)

output_headers = list(set(annotation_list[0].keys()) | set(utils.annotation_headers))  # try to pick up any extra headers that were written to the file
utils.write_annotations(args.outfile, glfo, [annotation], output_headers)
