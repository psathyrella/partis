#!/usr/bin/env python
import os
import sys
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')
import utils
import glutils

glfo = glutils.read_glfo('data/germlines/human', 'h')
# for g, s in glfo['seqs']['v'].items():
#     print '%s  %3d' % (utils.color_gene(g, width=20), len(s) - glfo['cyst-positions'][g])
# sys.exit()
base, a1, a2 = '4-30-2', '05', '03'
gene1, gene2 = 'IGHV' + base + '*' + a1, 'IGHV' + base + '*' + a2
seq1 = glfo['seqs']['v'][gene1]
seq2 = glfo['seqs']['v'][gene2]
min_length = min(len(seq1), len(seq2))
utils.color_mutants(seq1[:min_length], seq2[:min_length], print_result=True, print_isnps=True)
if min_length < len(seq1):
    print 'extra for %s: %s' % (utils.color_gene(gene1), seq1[min_length:])
if min_length < len(seq2):
    print 'extra for %s: %s' % (utils.color_gene(gene2), seq2[min_length:])
