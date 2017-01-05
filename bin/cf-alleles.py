#!/usr/bin/env python
import argparse
import os
import sys

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils

parser = argparse.ArgumentParser()
parser.add_argument('--base', required=True)
parser.add_argument('--alleles')
parser.add_argument('--other-genes')
parser.add_argument('--region', default='v')
parser.add_argument('--chain', default='h')
parser.add_argument('--glfo-dir', default='data/germlines/human')
args = parser.parse_args()
glfo = glutils.read_glfo(args.glfo_dir, args.chain)
if args.alleles is None:
    args.alleles = [utils.allele(g) for g in glfo['seqs'][args.region] if args.base == utils.primary_version(g) + '-' + utils.sub_version(g)]
else:
    args.alleles = utils.get_arg_list(args.alleles)
args.other_genes = utils.get_arg_list(args.other_genes)

# for g, s in glfo['seqs']['v'].items():
#     print '%s  %3d' % (utils.color_gene(g, width=20), len(s) - glfo['cyst-positions'][g])
# sys.exit()

# base = '4-59'
# a1, a2 = '12', '01'
# gene1, gene2 = 'IGHV' + base + '*' + a1, 'IGHV' + base + '*' + a2

genes = ['IG' + args.chain.upper() + args.region.upper() + args.base + '*' + al for al in args.alleles]
if args.other_genes is not None:
    genes += args.other_genes

codon_positions = glfo[utils.conserved_codons[args.chain][args.region] + '-positions']

def print_str(gene, seq):
    return '%s   %s' % (utils.color_gene(gene, width=15), seq)

ref_gene = genes[0]
ref_seq = glfo['seqs'][args.region][ref_gene]
print print_str(ref_gene, utils.color_mutants(ref_seq, ref_seq, emphasis_positions=[codon_positions[ref_gene] + i for i in range(3)])), '   (reference)'

for igene in range(1, len(genes)):
    gene = genes[igene]
    seq = glfo['seqs'][args.region][gene]
    min_length = min(len(seq), len(ref_seq))
    colored_seq = utils.color_mutants(ref_seq[:min_length], seq[:min_length], print_isnps=True, emphasis_positions=[codon_positions[gene] + i for i in range(3)])
    print print_str(gene, colored_seq)
    if min_length < len(ref_seq) and igene == 0:
        print 'extra for %s: %s' % (utils.color_gene(ref_gene), ref_seq[min_length:])
    if min_length < len(seq):
        print 'extra for %s: %s' % (utils.color_gene(gene), seq[min_length:])
