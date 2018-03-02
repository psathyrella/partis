#!/usr/bin/env python
import csv
import sys
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import argparse
import numpy
import itertools
from Bio.Seq import Seq

partis_path = '.'  # edit this if you're not running from the main partis dir
sys.path.insert(1, partis_path + '/python')
import utils
import glutils
from clusterpath import ClusterPath

parser = argparse.ArgumentParser()
parser.add_argument('--infiles')
parser.add_argument('--labels')
parser.add_argument('--locus')
parser.add_argument('--param')
args = parser.parse_args()

args.infiles = utils.get_arg_list(args.infiles)
args.labels = utils.get_arg_list(args.labels)
assert len(args.infiles) == len(args.labels)
glfo = glutils.read_glfo(args.param + '/hmm/germline-sets', locus=args.locus)

# ----------------------------------------------------------------------------------------
def getkey(uid_list):
    return ':'.join(uid_list)

# ----------------------------------------------------------------------------------------
def read_annotations(fname):
    annotations = {}
    with open(fname.replace('.csv', '-cluster-annotations.csv')) as csvfile:
        reader = csv.DictReader(csvfile)
        for line in reader:  # there's a line for each cluster
            if line['v_gene'] == '':  # failed (i.e. couldn't find an annotation)
                continue
            utils.process_input_line(line)  # converts strings in the csv file to floats/ints/dicts/etc.
            utils.add_implicit_info(glfo, line)  # add stuff to <line> that's useful, isn't written to the csv since it's redundant
            # utils.print_reco_event(line)  # print ascii-art representation of the rearrangement event
            annotations[getkey(line['unique_ids'])] = line
    return annotations

# ----------------------------------------------------------------------------------------
def naive_cdr3(info):
    naiveseq, _ = utils.subset_sequences(info, iseq=0, restrict_to_region='cdr3')
    return naiveseq

# ----------------------------------------------------------------------------------------
def naive_hdist(line1, line2):
    if line1['cdr3_length'] != line2['cdr3_length']:
        return 9999999
    return utils.hamming_distance(naive_cdr3(line1), naive_cdr3(line2))

# ----------------------------------------------------------------------------------------
def cdr3_translation(info):
    naive_cdr3_seq = naive_cdr3(info)
    if len(naive_cdr3_seq) % 3 != 0:
        # print '  out of frame: adding %s' % ((3 - len(naive_cdr3_seq) % 3) * 'N')
        naive_cdr3_seq += (3 - len(naive_cdr3_seq) % 3) * 'N'
    return Seq(naive_cdr3_seq).translate()

# ----------------------------------------------------------------------------------------
min_cluster_size = 4
cpaths = [ClusterPath() for _ in range(len(args.infiles))]
for ifile in range(len(args.infiles)):
    cpaths[ifile].readfile(args.infiles[ifile])
partitions = [sorted(cp.partitions[cp.i_best], key=len, reverse=True) for cp in cpaths]
partitions = [[c for c in partition if len(c) > min_cluster_size] for partition in partitions]
annotations = [read_annotations(fn) for fn in args.infiles]


nearest_cluster_lists = {l1 : {l2 : [] for l2 in args.labels if l2 != l1} for l1 in args.labels}
for if1 in range(len(args.infiles)):
    label1 = args.labels[if1]
    for if2 in range(len(args.infiles)):
        if if1 == if2:
            continue
        label2 = args.labels[if2]
        print '  %s vs %s' % (label1, label2)
        for cluster1 in partitions[if1]:  # for each cluster in the first partition
            info1 = annotations[if1][getkey(cluster1)]
            print '      %3d    %s' % (len(cluster1), cdr3_translation(info1))
            def keyfcn(c2):
                return naive_hdist(info1, annotations[if2][getkey(c2)])
            sorted_clusters = sorted([c for c in partitions[if2] if keyfcn(c) != 9999999], key=keyfcn)  # make a list of the clusters in the other partition that's sorted by how similar their naive sequence are
            nearest_cluster_lists[label1][label2].append(sorted_clusters)
            # print '             %s' % ' '.join([str(len(c)) for c in sorted_clusters])
            for nclust in sorted_clusters:
                nclust_naive_cdr3 = cdr3_translation(annotations[if2][getkey(nclust)])
                print '        %3d  %s' % (len(nclust), utils.color_mutants(cdr3_translation(info1), nclust_naive_cdr3))

# sorted_clusters = sorted(best_partition, key=len, reverse=True)  # sort by size
# # sort by size
# def keyfunc(q):
#     return len(annotations[q]['unique_ids'])

# sorted_clusters = sorted(annotations, key=keyfunc, reverse=True)

# # # loop over ten biggest
# # for cluster in sorted_clusters[:10]:
# #     print len(annotations[cluster]['unique_ids'])

# # add more criteria
# def boolfunc(q):
#     if annotations[cluster]['cdr3_length'] < 50:
#         return False
#     if len(annotations[cluster]['unique_ids']) < 20:
#         return False
#     return True

# interesting_clusters = [cluster for cluster in sorted_clusters if boolfunc(cluster)]
# print '  found %d interesting clusters' % len(interesting_clusters)
# for cluster in interesting_clusters:
#     print_stuff(annotations[cluster])

# sys.exit()
