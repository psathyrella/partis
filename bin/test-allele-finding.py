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
        cmd_str = base_cmd + ' simulate --n-sim-events 500 --n-procs 10 --simulate-partially-from-scratch --mutation-multiplier 0.5 --outfname ' + simfname

        if simulation_v_genes is not None:
            simulation_genes = simulation_v_genes + ':' + dj_genes
            sglfo = glutils.read_glfo('data/germlines/human', chain=chain, only_genes=simulation_genes.split(':'), debug=True)
            # snps_to_add = [{'gene' : 'IGHV4-59*01', 'positions' : (50, 200)}, ]
            # glutils.add_some_snps(snps_to_add, sglfo, remove_template_genes=False, debug=True)
            prevalence_fname = outdir + '/v_gene-probs.csv'  # NOTE there's some infrastructure for coming up with this file name automatically in utils.py
            prevalence_counts = {}
            for g in sglfo['seqs']['v']:
                prevalence_counts[g] = 1
            glutils.write_allele_prevalence_file('v', prevalence_fname, sglfo, prevalence_counts)
            glutils.write_glfo(outdir + '/germlines/simulation', sglfo)
            cmd_str += ' --initial-germline-dir ' + outdir + '/germlines/simulation'
            # cmd_str += ' --allele-prevalence-fnames ' +  prevalence_fname + '::'
        else:
            cmd_str += ' --generate-germline-set'

        if seed is not None:
            cmd_str += ' --seed ' + str(seed)
        run(cmd_str)
        sys.exit()

    # remove any old sw cache files
    sw_cachefiles = glob.glob(outpdir + '/sw-cache-*.csv')
    if len(sw_cachefiles) > 0:
        check_call(['rm', '-v'] + sw_cachefiles)

    # generate germline set and cache parameters
    cmd_str = base_cmd + ' cache-parameters --infname ' + simfname + ' --n-procs 10 --only-smith-waterman --debug-allele-finding'

    if inference_v_genes is not None:
        inference_genes = inference_v_genes + ':' + dj_genes
        iglfo = glutils.read_glfo('data/germlines/human', chain=chain, only_genes=inference_genes.split(':'), debug=True)
        glutils.write_glfo(outdir + '/germlines/inference', iglfo)
        cmd_str += ' --initial-germline-dir ' + outdir + '/germlines/inference'
        # cmd_str += ' --find-new-alleles'  # --new-allele-fname ' + outdir + '/new-alleles.fa'
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
inference_v_genes = None  # 'IGHV9-99*01' #'IGHV1-18*01'
simulation_v_genes = None  # 'IGHV4-59*01:IGHV4-59*04' # + IGHV1-18*01'

run_test(simulation_v_genes, inference_v_genes, dj_genes, seed=seed)

# allelic_groups = utils.separate_into_allelic_groups(glfo['seqs'])
