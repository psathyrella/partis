#!/usr/bin/env python
import csv
import os
import sys
import argparse
import colored_traceback.always
import glob
import subprocess
import operator

# if you move this script, you'll need to change this method of getting the imports
# partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
partis_dir = os.getcwd()
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
from clusterpath import ClusterPath

parser = argparse.ArgumentParser()
parser.add_argument('--basedir', default='/fh/fast/matsen_e/dralph/partis/mmshipley-bg505-naive-seeds-only')
parser.add_argument('--locus', required=True)
parser.add_argument('--prob-to-ignore', default=0.15)
args = parser.parse_args()

# ----------------------------------------------------------------------------------------
def read_linearham_output():
    lh_info = {}
    clusterdirs = glob.glob('%s/linearham/%s/cluster*' % (args.basedir, args.locus))  # /mcmciter10000_mcmcthin10_tuneiter5000_tunethin100_numrates4_seed0/burninfrac0.1_subsampfrac0.05/aa_naive_seqs.dnamap'
    for cdir in clusterdirs:
        input_seqfos = utils.read_fastx('%s/input_seqs.fasta' % cdir)
        input_uids = [sfo['name'] for sfo in input_seqfos if sfo['name'] != 'naive']
        outfnames = subprocess.check_output(['find', cdir, '-name', 'aa_naive_seqs.dnamap']).strip().split()
        if len(outfnames) == 0:
            print '  no linearham output for %s in %s' % (os.path.basename(cdir), cdir)
            continue
        elif len(outfnames) != 1:
            raise Exception('too many outfnames %s' % outfnames)
        clusterfo = []
        with open(outfnames[0]) as outfile:
            for line in outfile:
                if line[0] == '>':  # just skip these for now, we're just printing nucleotide level stuff, not aa
                    continue
                prob, naive_seq = line.strip().split(',')
                clusterfo.append((naive_seq, float(prob)))
        clusterfo = sorted(clusterfo, key=operator.itemgetter(1), reverse=True)  # it's sorted by aa naive seq in the file, and within that I think by nuc naive seq? Anyway, we need to make sure
        lh_info[':'.join(input_uids)] = clusterfo
    return lh_info

# ----------------------------------------------------------------------------------------
def print_lines(nseq_info, ref_seq, namestr, namecolor):
    assert nseq_info == sorted(nseq_info, key=operator.itemgetter(1), reverse=True)
    total_prob = 0.
    for naive_seq, prob in nseq_info:
        print '  %s   %5.2f  %s' % (utils.color_mutants(naive_seq if ref_seq is None else ref_seq, naive_seq), prob, utils.color(namecolor, namestr))
        if ref_seq is None:
            ref_seq = naive_seq
        if 1. - total_prob < args.prob_to_ignore:
            break
        total_prob += prob
    return ref_seq

glfo, annotation_list, cpath = utils.read_output('%s/%s/partition-with-alternative-annotations.yaml' % (args.basedir, args.locus))
lh_info = read_linearham_output()

# print annotations for the biggest cluster in the most likely partition
annotations = {':'.join(adict['unique_ids']) : adict for adict in annotation_list}  # collect the annotations in a dictionary so they're easier to access
most_likely_partition = cpath.partitions[cpath.i_best]  # a partition is represented as a list of lists of strings, with each string a sequence id
sorted_clusters = sorted(most_likely_partition, key=len, reverse=True)
for cluster in sorted_clusters:
    line = annotations[':'.join(cluster)]
    print ':'.join(line['unique_ids'])
    utils.print_reco_event(line, extra_str='  ')
    print ''

    lh_clusters = [(uidstr, cfo) for uidstr, cfo in lh_info.items() if set(uidstr.split(':')) & set(line['unique_ids'])]
    lh_naive_seqs = []
    if len(lh_clusters) == 0:
        print '  %s zero linearham clusters with any of these uids' % utils.color('red', 'error')
    elif len(lh_clusters) != 1:
        raise Exception('expected 1 linearham cluster but found %d' % len(lh_clusters))
    else:
        lh_uidstr, lh_naive_seqs = lh_clusters[0]
        if set(lh_uidstr.split(':')) != set(line['unique_ids']):
            print '  %s different uids\n       extra in linearham: %s\n   missing from linearham: %s' % (utils.color('red', 'error'), set(lh_uidstr.split(':')) - set(line['unique_ids']), set(line['unique_ids']) - set(lh_uidstr.split(':')))

    ref_seq = None
    if len(lh_naive_seqs) == 0:
        print '  %s no prob/naive_seq pairs in linearham for this cluster' % utils.color('red', 'error')

    ref_seq = print_lines(lh_naive_seqs, ref_seq, 'linearham', 'green')
    ref_seq = print_lines(sorted(line['alternative-annotations']['naive-seqs'], key=operator.itemgetter(1), reverse=True), ref_seq, 'partis', 'blue')

    print ''
