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
parser.add_argument('--ref-allele', help='print this one first')
parser.add_argument('--other-genes')
parser.add_argument('--region', default='v')
parser.add_argument('--locus', default='igh', choices=utils.loci.keys())
parser.add_argument('--glfo-dir', default='data/germlines/human')
args = parser.parse_args()
glfo = glutils.read_glfo(args.glfo_dir, args.locus)
if args.alleles is None:
    args.alleles = [utils.allele(g) for g in glfo['seqs'][args.region] if args.base == utils.primary_version(g) + ('-' + utils.sub_version(g) if utils.sub_version(g) is not None else '')]
else:
    args.alleles = utils.get_arg_list(args.alleles)
if len(args.alleles) == 0:
    raise Exception('couldn\'t find any alleles for --base %s. Other choices:\n    %s' % (args.base, ' '.join(glfo['seqs'][args.region])))
args.other_genes = utils.get_arg_list(args.other_genes)

# for g, s in glfo['seqs']['v'].items():
#     print '%s  %3d' % (utils.color_gene(g, width=20), len(s) - glfo['cyst-positions'][g])
# sys.exit()

# base = '4-59'
# a1, a2 = '12', '01'
# gene1, gene2 = 'IGHV' + base + '*' + a1, 'IGHV' + base + '*' + a2

genes = [args.locus.upper() + args.region.upper() + args.base + '*' + al for al in args.alleles]
if args.other_genes is not None:
    genes += args.other_genes
seqstrs = ['' for _ in range(len(genes))]
snpstrs = ['' for _ in range(len(genes))]

gene_str_width = max([utils.len_excluding_colors(utils.color_gene(g)) for g in genes])
codon_positions = glfo[utils.conserved_codons[args.locus][args.region] + '-positions'] if args.region != 'd' else None
max_seq_len = max([len(glfo['seqs'][args.region][g]) for g in genes])

ref_gene = genes[0] if args.ref_allele is None else utils.rejoin_gene(args.locus, args.region, utils.primary_version(genes[0]), utils.sub_version(genes[0]), args.ref_allele)
if ref_gene != genes[0]:
    genes.remove(ref_gene)
    genes.insert(0, ref_gene)
ref_seq = glfo['seqs'][args.region][ref_gene]
ref_pos = codon_positions[ref_gene]

for igene in range(0, len(genes)):
    gene = genes[igene]
    seq = glfo['seqs'][args.region][gene]
    pos = codon_positions[gene]
    if pos < ref_pos:  # align the codon position in the case that this seq is shorter up to the codon
        seq = (ref_pos - pos) * 'N' + seq
        pos += (ref_pos - pos)

    right_pad_str = ''
    if len(seq) < max_seq_len:
        right_pad_str = (max_seq_len - len(seq)) * ' '

    min_length = min(len(seq), len(ref_seq))  # can only compute hamming distance out as far as the shorter of the two
    emph_positions = None if args.region == 'd' else [pos + i for i in range(3)]
    colored_seq, isnps = utils.color_mutants(ref_seq[:min_length], seq[:min_length], return_isnps=True, emphasis_positions=emph_positions)
    seqstrs[igene] += '%s%s' % (colored_seq, right_pad_str)
    if min_length < len(seq):
        seqstrs[igene] += utils.color('blue', seq[min_length:])
    if len(isnps) > 0:
        snpstrs[igene] = '%2d (%s)' % (len(isnps), ' '.join([str(i) for i in isnps]))

# ----------------------------------------------------------------------------------------
def print_str(gene, seqstr, snpstr):
    return '%s  %s  %s  %s' % (utils.color_gene(gene, width=gene_str_width), seqstr, utils.color_gene(gene, width=gene_str_width), snpstr)

for igene in range(len(genes)):
    print print_str(genes[igene], seqstrs[igene], snpstrs[igene])
