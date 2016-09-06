#!/usr/bin/env python
import numpy
import csv
import operator
import sys
import os

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils

def build_v_gene_set(glfo, introns):
    for d_gene, counts in introns.items():
        print utils.color_gene(d_gene)
        refseq = None
        n_ok = 0
        mutecounts = {}
        for seq in sorted(counts, key=len, reverse=True):
            count = counts[seq]
            if refseq is None:
                refseq = seq
            # print '    %3d   %150s' % (count, seq)
            partial_refseq = refseq[len(refseq) - len(seq):]
            if seq == partial_refseq:
                n_ok += 1
            else:
                # utils.color_mutants(partial_refseq, seq, print_result=True, extra_str='                ')
                n_mutes = utils.hamming_distance(partial_refseq, seq)
                if n_mutes not in mutecounts:
                    mutecounts[n_mutes] = 0
                mutecounts[n_mutes] += 1
        print '    %d / %d ok' % (n_ok, n_ok + sum(mutecounts.values())),
        if len(mutecounts) > 0:
            print '(mean of %.1f mutations among the other %d' % (numpy.average(mutecounts.keys(), weights=mutecounts.values()), sum(mutecounts.values())),
        print ''

# tmpglfo = glutils.read_glfo('tmp-germlines', 'h')
glfo = glutils.read_glfo('data/germlines/human', 'h')

infname = '/fh/fast/matsen_e/data/2016-06-02-katie/VMO_Memory-3/VMO_Memory-3.tsv'

n_failed, n_v_ok, n_total = 0, 0, 0
introns = {}

with open(infname) as infile:
    reader = csv.DictReader(infile, delimiter='\t')
    iline = 0
    introns = {}
    for line in reader:
        iline += 1
        if iline > 5000:
            break
        n_total += 1
        if line[utils.adaptive_headers['v_gene']] != 'unresolved':  # for the moment, skip seqs that have v genes
            n_v_ok += 1
            continue
        line = utils.convert_from_adaptive_headers(glfo, line, uid='%09d' % iline, only_dj_rearrangements=True)
        if line['failed']:
            n_failed += 1
            continue

        d_gene = line['d_gene']
        if d_gene not in introns:
            introns[d_gene] = {}
        iseq = glfo['seqs']['v'][line['v_gene']]
        if iseq not in introns[d_gene]:
            introns[d_gene][iseq] = 0
        introns[d_gene][iseq] += 1

build_v_gene_set(glfo, introns)

print 'processed %d (of these:   %d failed   %d v ok)' % (n_total, n_failed, n_v_ok)
