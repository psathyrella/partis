#!/usr/bin/env python
import sys
import os
import glob
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

    # simulate
    if True:
        simevents = 8000
        cmd_str = base_cmd + ' simulate --n-sim-events ' + str(simevents) + ' --n-leaves 1 --constant-number-of-leaves --simulate-partially-from-scratch --mutation-multiplier 0.5 --outfname ' + simfname
        # cmd_str += ' --n-procs 10'
        cmd_str += ' --subsimproc --slurm --n-procs 30'

        if simulation_v_genes is not None:
            simulation_genes = simulation_v_genes + ':' + dj_genes
            sglfo = glutils.read_glfo('data/germlines/human', chain=chain, only_genes=simulation_genes.split(':'), debug=True)
            # snps_to_add = [
            #     {'gene' : 'IGHV4-59*01', 'positions' : (94, 30, 138, 13, 62, 205, 77, 237, 93, 218, 76, 65, 31, 120, 22, 216, 79, 56, 109)},
            # ]
            # glutils.add_some_snps(snps_to_add, sglfo, debug=True)
            glutils.write_glfo(outdir + '/germlines/simulation', sglfo)
            cmd_str += ' --initial-germline-dir ' + outdir + '/germlines/simulation'
        else:
            cmd_str += ' --generate-germline-set'
            cmd_str += ' --n-genes-per-region 1:5:3'
            cmd_str += ' --n-alleles-per-gene 2,3:1,2:1,2'

        if seed is not None:
            cmd_str += ' --seed ' + str(seed)
        run(cmd_str)

    # remove any old sw cache files
    sw_cachefiles = glob.glob(outpdir + '/sw-cache-*.csv')
    if len(sw_cachefiles) > 0:
        check_call(['rm', '-v'] + sw_cachefiles)

    # generate germline set and cache parameters
    cmd_str = base_cmd + ' cache-parameters --infname ' + simfname + ' --only-smith-waterman --debug-allele-finding'
    # cmd_str += ' --n-procs 10'
    cmd_str += ' --n-procs 30 --slurm'

    if inference_v_genes is not None:
        inference_genes = inference_v_genes + ':' + dj_genes
        iglfo = glutils.read_glfo('data/germlines/human', chain=chain, only_genes=inference_genes.split(':'), debug=True)
        glutils.write_glfo(outdir + '/germlines/inference', iglfo)
        cmd_str += ' --initial-germline-dir ' + outdir + '/germlines/inference'
        cmd_str += ' --find-new-alleles'  # --new-allele-fname ' + outdir + '/new-alleles.fa'
    else:
        cmd_str += ' --generate-germline-set'

    cmd_str += ' --parameter-dir ' + outpdir
    cmd_str += ' --plotdir ' + plotdir
    if seed is not None:
        cmd_str += ' --seed ' + str(seed)
    run(cmd_str)

# ----------------------------------------------------------------------------------------

seed = None  # 1
dj_genes = 'IGHD6-19*01:IGHJ4*02'
inference_v_genes = None #'IGHV4-30-2*05' #'IGHV1-18*01'
simulation_v_genes = None #inference_v_genes + ':IGHV4-30-2*03' # + ':IGHV9-99*02'  # 'IGHV4-59*01:IGHV4-59*04' # + IGHV1-18*01'

run_test(simulation_v_genes, inference_v_genes, dj_genes, seed=seed)

# allelic_groups = utils.separate_into_allelic_groups(glfo['seqs'])
