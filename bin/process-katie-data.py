#!/usr/bin/env python
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

# tmpglfo = glutils.read_glfo('tmp-germlines', 'h')
glfo = glutils.read_glfo('data/germlines/human', 'h')

infname = '/fh/fast/matsen_e/data/2016-06-02-katie/VMO_Memory-3/VMO_Memory-3.tsv'

n_failed = 0
introns = {}

with open(infname) as infile:
    reader = csv.DictReader(infile, delimiter='\t')
    iline = 0
    introns = {}
    for line in reader:
        iline += 1
        if iline > 1000:
            break
        if line[utils.adaptive_headers['v_gene']] != 'unresolved':  # for the moment, skip seqs that have v genes
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

for d_gene, counts in introns.items():
    print d_gene
    for seq, count in sorted(counts.items(), key=operator.itemgetter(1), reverse=True):
        print '    %3d   %s' % (count, seq)

print 'failed: %d' % n_failed
