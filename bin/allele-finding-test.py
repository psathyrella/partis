#!/usr/bin/env python
import sys
import os
import random
from subprocess import check_call
sys.path.insert(1, './python')

import utils

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print 'RUN', cmd_str
    sys.stdout.flush()
    check_call(cmd_str.split())

outdir = '_tmp/allele-finder'
base_cmd = './bin/partis'

# ----------------------------------------------------------------------------------------
def join_gene_names(gene_name_str):
    return ':'.join([utils.sanitize_name(g) for g in gene_name_str.split(':')])

# ----------------------------------------------------------------------------------------
def get_label(existing_genes, new_allele):
    return '_existing_' + join_gene_names(existing_genes) + '_new_' + join_gene_names(new_allele)

# ----------------------------------------------------------------------------------------
def run_test(existing_v_genes, new_v_allele, dj_genes, seed=None):
    if seed is not None:
        random.seed(seed)

    label = 'test'  #get_label(existing_genes, new_allele)
    simfname = outdir + '/simu-' + label + '.csv'
    outpdir = outdir + '/simu-' + label
    plotdir = os.getenv('www') + '/partis/allele-finding/' + label

    existing_genes = existing_v_genes + ':' + dj_genes

    snps_to_add = {'IGHV3-71*01' : 4}
    utils.rewrite_germline_fasta('data/imgt', outdir + '/germlines-for-simulation', only_genes=existing_genes.split(':'), snps_to_add=snps_to_add, rename_snpd_genes=True)

    # simulate
    cmd_str = base_cmd + ' simulate --n-sim-events 1000 --n-procs 10 --simulate-partially-from-scratch --mutation-multiplier 0.5'
    cmd_str += ' --datadir ' + outdir + '/germlines-for-simulation'
    cmd_str += ' --outfname ' + simfname
    if seed is not None:
        cmd_str += ' --seed ' + str(seed)
    run(cmd_str)

    utils.rewrite_germline_fasta('data/imgt', outdir + '/germlines', only_genes=existing_genes.split(':'))

    # cache-parameters
    cmd_str = base_cmd + ' cache-parameters --infname ' + simfname + ' --n-procs 10 --find-new-alleles --only-smith-waterman'
    cmd_str += ' --datadir ' + outdir + '/germlines'
    cmd_str += ' --only-genes ' + existing_genes
    cmd_str += ' --parameter-dir ' + outpdir
    cmd_str += ' --plotdir ' + plotdir
    if seed is not None:
        cmd_str += ' --seed ' + str(seed)
    run(cmd_str)

# ----------------------------------------------------------------------------------------

seed = 1  # None
dj_genes = 'IGHD6-19*01:IGHJ4*02'
existing_v_genes = 'IGHV3-71*01:IGHV3-71*03'  # 1-18*01
new_v_allele = 'IGHV3-71*03' #1-18*04

run_test(existing_v_genes, new_v_allele, dj_genes, seed=seed)

# glfo = utils.read_germline_set('data/imgt')
# allelic_groups = utils.separate_into_allelic_groups(glfo['seqs'])
# for primary_version in allelic_groups['v']:
#     for sub_version in allelic_groups['v'][primary_version]:
#         if len(allelic_groups['v'][primary_version][sub_version]) == 1:
#             continue
#         print '    %15s   %15s   %s' % (primary_version, sub_version, allelic_groups['v'][primary_version][sub_version])

# print ''
# print glfo['seqs']['v']['IGHV3-30*01']
# for g in glfo['seqs']['v']:
#     if '3-30*' not in g:
#         continue
#     print utils.color_mutants(glfo['seqs']['v']['IGHV3-30*01'], glfo['seqs']['v'][g])

