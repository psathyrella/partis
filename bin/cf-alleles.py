#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import os
import sys

import partis.utils as utils
import partis.glutils as glutils

parser = argparse.ArgumentParser()
parser.add_argument('--bases', required=True, help='colon-separated list of the bits before the stars, e.g. 1-18:2-2 (set to \'all\' to print entire germline set)')
parser.add_argument('--allele-numbers')
parser.add_argument('--ref-allele', help='print this one first')
parser.add_argument('--other-genes')
parser.add_argument('--region', default='v')
parser.add_argument('--locus', default='igh', choices=utils.loci)
parser.add_argument('--species', default='human')
parser.add_argument('--glfo-dir', help='default set below')
args = parser.parse_args()

if args.glfo_dir is None:
    args.glfo_dir = 'data/germlines/' + args.species

glfo = glutils.read_glfo(args.glfo_dir, args.locus)

# ----------------------------------------------------------------------------------------
def get_base(gene):
    basestr = utils.primary_version(gene)
    if utils.sub_version(gene) is not None:
        basestr += '-' + utils.sub_version(gene)
    return basestr

# ----------------------------------------------------------------------------------------
def get_genes(base, alleles=None):
    if alleles is None:  # take all of 'em
        alleles = [utils.allele(g) for g in glfo['seqs'][args.region] if base == get_base(g)]
    return [args.locus.upper() + args.region.upper() + base + '*' + al for al in alleles]

if args.bases == 'all':
    input_groupfcn = None  # lambda g: str(utils.primary_version(g) in ['4', '5'])  # this example puts all the 4 and 5 primary versions in one group, and everybody else in another
    glutils.print_glfo(glfo, only_region=(args.region if args.region != 'v' else None), input_groupfcn=input_groupfcn)  # not much point in doing only v, since it's the one that takes most of the time
    sys.exit(0)

args.bases = utils.get_arg_list(args.bases)
args.allele_numbers = utils.get_arg_list(args.allele_numbers)
genes = [g for base in args.bases for g in get_genes(base, args.allele_numbers)]
if len(genes) == 0:
    raise Exception('couldn\'t find any genes for the specified --bases %s\n  choices:\n    %s' % (' '.join(args.bases), ' '.join(sorted(set([get_base(g) for g in glfo['seqs'][args.region]])))))
args.other_genes = utils.get_arg_list(args.other_genes)
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
        seq = (ref_pos - pos) * '-' + seq
        pos += (ref_pos - pos)

    right_pad_str = ''  # i think i don't need this any more since i have the align option in color_mutants
    # if len(seq) < max_seq_len:
    #     right_pad_str = (max_seq_len - len(seq)) * ' '

    emph_positions = None if args.region == 'd' else [pos + i for i in range(3)]
    colored_seq, isnps = utils.color_mutants(ref_seq, seq, return_isnps=True, emphasis_positions=emph_positions, align=True)
    seqstrs[igene] += '%s%s' % (colored_seq, right_pad_str)
    if len(isnps) > 0:
        snpstrs[igene] = '%2d (%s)' % (len(isnps), ' '.join([str(i) for i in isnps]))

# ----------------------------------------------------------------------------------------
def print_str(gene, seqstr, snpstr):
    return '%s  %s  %s  %s' % (utils.color_gene(gene, width=gene_str_width), seqstr, utils.color_gene(gene, width=gene_str_width), snpstr)

for igene in range(len(genes)):
    print(print_str(genes[igene], seqstrs[igene], snpstrs[igene]))
