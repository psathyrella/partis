#!/usr/bin/env python
import argparse
import time
import sys
import os
import glob

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
def run_test(args):
    print 'seed %d' % args.seed
    label = 'test'  #get_label(existing_genes, new_allele)
    simfname = outdir + '/simu-' + label + '.csv'
    outpdir = outdir + '/simu-' + label
    if os.getenv('www') is not None:
        plotdir = os.getenv('www') + '/partis/allele-finding/' + label
    else:
        plotdir = '_www/partis/allele-finding/' + label

    # simulate
    if not args.nosim:
        cmd_str = base_cmd + ' simulate --n-sim-events ' + str(args.n_sim_events) + ' --n-leaves 1 --constant-number-of-leaves --simulate-partially-from-scratch --mutation-multiplier 0.5 --outfname ' + simfname
        if args.slurm:
            cmd_str += ' --n-procs 30 --slurm'
        else:
            cmd_str += ' --n-procs 10'

        if args.gen_gset:
            cmd_str += ' --generate-germline-set'
            cmd_str += ' --n-genes-per-region 1:5:3'
            cmd_str += ' --n-alleles-per-gene 2,3:1,2:1,2'
        else:
            simulation_genes = args.sim_v_genes + ':' + args.dj_genes
            sglfo = glutils.read_glfo('data/germlines/human', chain=chain, only_genes=simulation_genes.split(':'), debug=True)
            # snps_to_add = [
            #     {'gene' : 'IGHV4-59*01', 'positions' : (94, 30, 138, 13, 62, 205, 77, 237, 93, 218, 76, 65, 31, 120, 22, 216, 79, 56, 109)},
            # ]
            # glutils.add_some_snps(snps_to_add, sglfo, debug=True)
            glutils.write_glfo(outdir + '/germlines/simulation', sglfo)
            cmd_str += ' --initial-germline-dir ' + outdir + '/germlines/simulation'

        if args.seed is not None:
            cmd_str += ' --seed ' + str(args.seed)
        run(cmd_str)

    # remove any old sw cache files
    sw_cachefiles = glob.glob(outpdir + '/sw-cache-*')
    if len(sw_cachefiles) > 0:
        check_call(['rm', '-rv'] + sw_cachefiles)

    # generate germline set and cache parameters
    cmd_str = base_cmd + ' cache-parameters --infname ' + simfname + ' --only-smith-waterman --debug-allele-finding'
    if args.slurm:
        cmd_str += ' --n-procs 30 --slurm'
    else:
        cmd_str += ' --n-procs 10'

    if args.gen_gset:
        cmd_str += ' --generate-germline-set'
    else:
        inference_genes = args.inf_v_genes + ':' + args.dj_genes
        iglfo = glutils.read_glfo('data/germlines/human', chain=chain, only_genes=inference_genes.split(':'), debug=True)
        glutils.write_glfo(outdir + '/germlines/inference', iglfo)
        cmd_str += ' --initial-germline-dir ' + outdir + '/germlines/inference'
        cmd_str += ' --find-new-alleles'  # --new-allele-fname ' + outdir + '/new-alleles.fa'

    cmd_str += ' --parameter-dir ' + outpdir
    cmd_str += ' --plotdir ' + plotdir
    if args.seed is not None:
        cmd_str += ' --seed ' + str(args.seed)
    run(cmd_str)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--nosim', action='store_true')
parser.add_argument('--n-sim-events', type=int, default=2000)
parser.add_argument('--seed', type=int, default=int(time.time()))
parser.add_argument('--gen-gset', action='store_true')
parser.add_argument('--dj-genes', default='IGHD6-19*01:IGHJ4*02')
parser.add_argument('--sim-v-genes', default='IGHV4-39*01')
parser.add_argument('--inf-v-genes', default='IGHV4-39*06')
parser.add_argument('--slurm', action='store_true')
args = parser.parse_args()
if not args.gen_gset is None:
    args.sim_v_genes += ':' + args.inf_v_genes

run_test(args)
