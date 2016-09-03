#!/usr/bin/env python
import csv
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

with open(infname) as infile:
    reader = csv.DictReader(infile, delimiter='\t')
    iline = 0
    introns = {}
    for line in reader:
        iline += 1
        if iline > 10:
            break
        # print '%5s%5s%5s' % (line['vIndex'], line['dIndex'], line['jIndex'])
        # continue
        print ''
        line = utils.convert_from_adaptive_headers(glfo, line, uid='%09d' % iline)
        if line['v_gene'] is not None:  # for the moment, skip seqs that have v genes
            continue
        if line['d_gene'] is None or line['j_gene'] is None:
            continue

        # seq = line['nucleotide']
        # d_seq = glfo['seqs']['d'][line['dGeneName']]
        # dstart = seq.find(
