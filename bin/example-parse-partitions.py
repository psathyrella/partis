#!/usr/bin/env python
import csv
import os
import sys
import argparse

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
from clusterpath import ClusterPath

parser = argparse.ArgumentParser()
parser.add_argument('--partition-file', default=partis_dir + '/test/reference-results/partition-ref-simu.csv')
parser.add_argument('--annotation-file', default=partis_dir + '/test/reference-results/partition-ref-simu-cluster-annotations.csv')
parser.add_argument('--germline-info-dir', default=partis_dir + '/test/reference-results/test/parameters/simu/hmm/germline-sets')
parser.add_argument('--locus', default='igh')
args = parser.parse_args()

glfo = glutils.read_glfo(args.germline_info_dir, locus=args.locus)  # read germline info (it'll be nice to switch to yaml output files so we can store this in the same file as the rest of the output)

# get list of most likely partitions
cp = ClusterPath()
cp.readfile(args.partition_file)
print 'list of partitions:'
cp.print_partitions(abbreviate=True)  # print little 'o's instead of the full sequence ids

# get annotations for each cluster in the single most likely partition
annotations = {}
with open(args.annotation_file) as csvfile:
    reader = csv.DictReader(csvfile)
    for line in reader:  # each line corresponds to a single cluster from the most likely partition
        if line['v_gene'] == '':  # failed (i.e. couldn't find an annotation)
            continue
        utils.process_input_line(line)  # converts strings in the csv file to floats/ints/dicts/etc.
        utils.add_implicit_info(glfo, line)  # add stuff to <line> that's useful, but isn't written to the csv since it's redundant
        annotations[':'.join(line['unique_ids'])] = line  # can't hash a list, so convert back to the colon-separated-string format (no, this isn't the best way to do this)

# print annotations for the three biggest clusters in the most likely partition
most_likely_partition = cp.partitions[cp.i_best]  # a partition is represented as a list of lists of strings, with each string a sequence id
biggest_clusters = sorted(most_likely_partition, key=len, reverse=True)[:3]
print '\nannotations for the three biggest clusters:'
for cluster in biggest_clusters:
    cluster_annotation = annotations[':'.join(cluster)]
    utils.print_reco_event(cluster_annotation)
