#!/usr/bin/env python
import csv
import sys
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import argparse
import numpy
import itertools
from Bio.Seq import Seq

# # example usage:
# fsd=/fh/fast/matsen_e/processed-data/partis
# ./bin/cross-check-partitions.py \
#     --locus igh \
#     --min-cluster-sizes 150:5 --max-cdr3-distance 5 \
#     --param $fsd/laura-mb/v17/Hs-LN-D-5RACE-IgG:$fsd/laura-mb-2/v17/BF520-g-M9 \
#     --labels laura-mb-D:laura-mb-2-M9 \
#     --infiles $fsd/laura-mb/v17/partitions/Hs-LN-D-5RACE-IgG-isub-2/partition.csv:$fsd/laura-mb-2/v17/partitions/BF520-g-M9-isub-2/partition.csv

partis_path = '.'  # edit this if you're not running from the main partis dir
sys.path.insert(1, partis_path + '/python')
import utils
import glutils
from clusterpath import ClusterPath

parser = argparse.ArgumentParser()
parser.add_argument('--infiles')
parser.add_argument('--labels')
parser.add_argument('--locus')
parser.add_argument('--parameter-dirs')
parser.add_argument('--min-cluster-sizes', default='4:2', help='first number is for outer cluster loop (first column), second is for inner loop (second column)')
parser.add_argument('--max-cdr3-distance', default=5, type=int, help='ignore clusters with a cdr3 that differs by more than this many nucleotides')
args = parser.parse_args()

args.infiles = utils.get_arg_list(args.infiles)
args.labels = utils.get_arg_list(args.labels)
args.parameter_dirs = utils.get_arg_list(args.parameter_dirs)
args.min_cluster_sizes = utils.get_arg_list(args.min_cluster_sizes, intify=True)
assert len(args.infiles) == len(args.labels)
if len(args.parameter_dirs) == 1:
    print '  note: using same glfo for all infiles'
    args.parameter_dirs = [args.parameter_dirs[0] for _ in args.labels]
assert len(args.parameter_dirs) == len(args.labels)

glfos = [glutils.read_glfo(pdir + '/hmm/germline-sets', locus=args.locus) for pdir in args.parameter_dirs]

# ----------------------------------------------------------------------------------------
def getkey(uid_list):
    return ':'.join(uid_list)

# ----------------------------------------------------------------------------------------
def read_annotations(fname, glfo):
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
def naive_hdist_or_none(line1, line2):
    if line1['cdr3_length'] != line2['cdr3_length']:
        return None
    hdist = utils.hamming_distance(naive_cdr3(line1), naive_cdr3(line2))
    if hdist > args.max_cdr3_distance:
        return None
    return hdist

# ----------------------------------------------------------------------------------------
def cdr3_translation(info):
    naive_cdr3_seq = naive_cdr3(info)
    if len(naive_cdr3_seq) % 3 != 0:
        # print '  out of frame: adding %s' % ((3 - len(naive_cdr3_seq) % 3) * 'N')
        naive_cdr3_seq += (3 - len(naive_cdr3_seq) % 3) * 'N'
    return Seq(naive_cdr3_seq).translate()

# ----------------------------------------------------------------------------------------
cpaths = [ClusterPath() for _ in range(len(args.infiles))]
for ifile in range(len(args.infiles)):
    cpaths[ifile].readfile(args.infiles[ifile])
partitions = [sorted(cp.partitions[cp.i_best], key=len, reverse=True) for cp in cpaths]
partitions = [[c for c in partition if len(c) > args.min_cluster_sizes[1]] for partition in partitions]
annotations = [read_annotations(args.infiles[ifn], glfos[ifn]) for ifn in range(len(args.infiles))]

nearest_cluster_lists = {l1 : {l2 : [] for l2 in args.labels if l2 != l1} for l1 in args.labels}
for if1 in range(len(args.infiles)):
    label1 = args.labels[if1]
    print '%s' % utils.color(None, label1)
    for if2 in range(len(args.infiles)):
        if if1 == if2:
            continue
        label2 = args.labels[if2]
        print '\n      %5s      %5s    cdr3' % ('', label2)
        print '   index size  index size  dist'
        for cluster1 in partitions[if1]:  # for each cluster in the first partition
            if len(cluster1) < args.min_cluster_sizes[0]:
                continue
            info1 = annotations[if1][getkey(cluster1)]
            def keyfcn(c2):
                return naive_hdist_or_none(info1, annotations[if2][getkey(c2)])
            sorted_clusters = sorted([c for c in partitions[if2] if keyfcn(c) is not None], key=keyfcn)  # make a list of the clusters in the other partition that's sorted by how similar their naive sequence are
            nearest_cluster_lists[label1][label2].append(sorted_clusters)

            size_index_str = '%4d %4d' % (partitions[if1].index(cluster1), len(cluster1))
            extra_str = ''
            if len(sorted_clusters) == 0:
                size_index_str = utils.color('yellow', size_index_str)
                extra_str = utils.color('yellow', '  x')
            print '     %s                   %-30s%s' % (size_index_str, cdr3_translation(info1), extra_str)
            for nclust in sorted_clusters:
                nclust_naive_cdr3 = cdr3_translation(annotations[if2][getkey(nclust)])
                hdist = naive_hdist_or_none(info1, annotations[if2][getkey(nclust)])
                print '               %4d %4d    %2s   %-30s' % (partitions[if2].index(nclust), len(nclust), '%d' % hdist if hdist > 0 else '',
                                                                 utils.color_mutants(cdr3_translation(info1), nclust_naive_cdr3, amino_acid=True))
