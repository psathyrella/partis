#!/usr/bin/env python
import sys
import os
import random
from subprocess import check_call
sys.path.insert(1, './python')

import utils
import glutils

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print 'RUN', cmd_str
    sys.stdout.flush()
    check_call(cmd_str.split())

outdir = '_tmp/allele-finder'
base_cmd = './bin/partis'
chain = 'h'

# ----------------------------------------------------------------------------------------
def join_gene_names(gene_name_str):
    return ':'.join([utils.sanitize_name(g) for g in gene_name_str.split(':')])

# ----------------------------------------------------------------------------------------
def get_label(existing_genes, new_allele):
    return '_existing_' + join_gene_names(existing_genes) + '_new_' + join_gene_names(new_allele)

# ----------------------------------------------------------------------------------------
def run_test(simulation_v_genes, inference_v_genes, dj_genes, seed=None):
    if seed is not None:
        random.seed(seed)

    label = 'test'  #get_label(existing_genes, new_allele)
    simfname = outdir + '/simu-' + label + '.csv'
    outpdir = outdir + '/simu-' + label
    if os.getenv('www') is not None:
        plotdir = os.getenv('www') + '/partis/allele-finding/' + label
    else:
        plotdir = '_www/partis/allele-finding/' + label

    snps_to_add = [
        {'gene' : 'IGHV3-71*01', 'positions' : (35, )},
        {'gene' : 'IGHV3-71*01', 'positions' : (35, 50)},
        # {'gene' : 'IGHV3-71*01', 'positions' : (35, 45, 20, 50, 77)},
        # {'gene' : 'IGHV3-71*01', 'positions' : (35, 60, 50)},
        # {'gene' : 'IGHV1-18*01', 'positions' : (100, 101)},
        # {'gene' : 'IGHV1-18*01', 'positions' : (20, )}
    ]
    simulation_genes = simulation_v_genes + ':' + dj_genes
    sglfo = read_glfo('data/imgt', chain=chain, only_genes=simulation_genes.split(':'), debug=True)
    glutils.add_some_snps(snps_to_add, sglfo, debug=True)#, remove_template_genes=True)
    glutils.write_glfo(outdir + '/germlines-for-simulation', sglfo)

    # simulate
    cmd_str = base_cmd + ' simulate --n-sim-events 1000 --n-procs 10 --simulate-partially-from-scratch --mutation-multiplier 0.5'
    cmd_str += ' --initial-datadir ' + outdir + '/germlines-for-simulation'
    cmd_str += ' --outfname ' + simfname
    if seed is not None:
        cmd_str += ' --seed ' + str(seed)
    run(cmd_str)

    inference_genes = inference_v_genes + ':' + dj_genes
    iglfo = read_glfo('data/imgt', chain=chain, only_genes=inference_genes.split(':'), debug=True)
    glutils.write_glfo(outdir + '/germlines-for-inference', iglfo)

    # generate germline set and cache parameters
    cmd_str = base_cmd + ' cache-parameters --infname ' + simfname + ' --n-procs 10 --only-smith-waterman'
    cmd_str += ' --find-new-alleles --new-allele-fname ' + outdir + '/new-alleles.fa'
    # cmd_str += ' --generate-germline-set'
    cmd_str += '  --debug-new-allele-finding'
    cmd_str += ' --initial-datadir ' + outdir + '/germlines-for-inference'
    cmd_str += ' --parameter-dir ' + outpdir
    cmd_str += ' --plotdir ' + plotdir
    if seed is not None:
        cmd_str += ' --seed ' + str(seed)
    run(cmd_str)

# ----------------------------------------------------------------------------------------

seed = 1  # None
dj_genes = 'IGHD6-19*01:IGHJ4*02'
inference_v_genes = 'IGHV3-71*01' #:IGHV1-18*01'
simulation_v_genes = inference_v_genes  # + ':IGHV3-71*02:IGHV3-71*03'  #:IGHV1-18*01'

run_test(simulation_v_genes, inference_v_genes, dj_genes, seed=seed)
sys.exit()

allelic_groups = utils.separate_into_allelic_groups(glfo['seqs'])
for gene in allelic_groups['v']['1']['69']:
    print glfo['seqs']['v'][gene]
sys.exit()
for primary_version in allelic_groups['v']:
    for sub_version in allelic_groups['v'][primary_version]:
        if len(allelic_groups['v'][primary_version][sub_version]) == 1:
            continue
        print '    %15s   %15s   %s' % (primary_version, sub_version, allelic_groups['v'][primary_version][sub_version])

# print ''
# print glfo['seqs']['v']['IGHV3-30*01']
# for g in glfo['seqs']['v']:
#     if '3-30*' not in g:
#         continue
#     print utils.color_mutants(glfo['seqs']['v']['IGHV3-30*01'], glfo['seqs']['v'][g])
