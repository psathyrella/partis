#!/usr/bin/env python
import copy
import random
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
    simfname = args.outdir + '/simu-' + label + '.csv'
    outpdir = args.outdir + '/simu-' + label
    if os.getenv('www') is not None:
        plotdir = os.getenv('www') + '/partis/allele-finding/' + label
    else:
        plotdir = '_www/partis/allele-finding/' + label

    # simulate
    if not args.nosim:
        cmd_str = base_cmd + ' simulate --n-sim-events ' + str(args.n_sim_events) + ' --n-leaves 1 --constant-number-of-leaves --simulate-partially-from-scratch --outfname ' + simfname
        cmd_str += ' --mutation-multiplier ' + str(args.mut_mult)

        cmd_str += ' --n-procs ' + str(args.n_procs)
        if args.slurm:
            cmd_str += ' --slurm --subsimproc'

        if args.gen_gset:
            cmd_str += ' --generate-germline-set'
            cmd_str += ' --n-genes-per-region 1:5:3'
            cmd_str += ' --n-alleles-per-gene 2,3:1,2:1,2'
        else:
            simulation_genes = ':'.join(args.sim_v_genes + args.dj_genes)
            sglfo = glutils.read_glfo('data/germlines/human', chain=chain, only_genes=simulation_genes.split(':'), debug=True)

            if args.snp_positions is not None:
                snps_to_add = [{'gene' : g, 'positions' : args.snp_positions} for g in args.sim_v_genes]
                glutils.add_some_snps(snps_to_add, sglfo, debug=True)

            glutils.write_glfo(args.outdir + '/germlines/simulation', sglfo)
            cmd_str += ' --initial-germline-dir ' + args.outdir + '/germlines/simulation'

        if args.seed is not None:
            cmd_str += ' --seed ' + str(args.seed)
        run(cmd_str)

    # remove any old sw cache files
    sw_cachefiles = glob.glob(outpdir + '/sw-cache-*.csv')
    if len(sw_cachefiles) > 0:
        for cachefname in sw_cachefiles:
            check_call(['rm', '-v', cachefname])
            sw_cache_gldir = cachefname.replace('.csv', '-glfo')
            if os.path.exists(sw_cache_gldir):  # if stuff fails halfway through, you can get one but not the other
                glutils.remove_glfo_files(sw_cache_gldir, chain)
                # os.rmdir(sw_cache_gldir)

    # generate germline set and cache parameters
    cmd_str = base_cmd + ' cache-parameters --infname ' + simfname + ' --only-smith-waterman --debug-allele-finding'
    # cmd_str = 'python -m cProfile -s tottime -o prof.out ' + cmd_str
    cmd_str += ' --n-procs ' + str(args.n_procs)
    if args.slurm:
        cmd_str += ' --slurm'

    if args.gen_gset:
        cmd_str += ' --generate-germline-set'
    else:
        inference_genes = ':'.join(args.inf_v_genes + args.dj_genes)
        iglfo = glutils.read_glfo('data/germlines/human', chain=chain, only_genes=inference_genes.split(':'), debug=True)
        glutils.write_glfo(args.outdir + '/germlines/inference', iglfo)
        cmd_str += ' --initial-germline-dir ' + args.outdir + '/germlines/inference'
        cmd_str += ' --find-new-alleles'  # --new-allele-fname ' + args.outdir + '/new-alleles.fa'

    cmd_str += ' --parameter-dir ' + outpdir
    cmd_str += ' --plotdir ' + plotdir
    if args.seed is not None:
        cmd_str += ' --seed ' + str(args.seed)
    run(cmd_str)

# ----------------------------------------------------------------------------------------
def comprehensive_test(args):
    def cmd_str(iproc):
        clist = copy.deepcopy(sys.argv)
        utils.remove_from_arglist(clist, '--comprehensive')
        utils.remove_from_arglist(clist, '--n-tests', has_arg=True)
        utils.replace_in_arglist(clist, '--outdir', args.outdir + '/' + str(iproc))
        utils.replace_in_arglist(clist, '--seed', str(args.seed + iproc))
        # clist.append('--slurm')
        return ' '.join(clist)

    cmdfos = [{'cmd_str' : cmd_str(iproc),
               'workdir' : args.workdir + '/' + str(iproc),
               'logdir' : args.outdir + '/' + str(iproc),
               'outfname' : args.outdir + '/' + str(iproc)}
        for iproc in range(args.n_tests)]
    for iproc in range(args.n_tests):
        if os.path.exists(cmdfos[iproc]['outfname']):
            check_call(['rm', '-r', cmdfos[iproc]['outfname']])
    print '  look for logs in %s' % args.outdir
    utils.run_cmds(cmdfos, debug='write')

# ----------------------------------------------------------------------------------------
fsdir = '/fh/fast/matsen_e/dralph'

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--nosim', action='store_true')
parser.add_argument('--n-sim-events', type=int, default=200)
parser.add_argument('--n-procs', type=int, default=2)
parser.add_argument('--seed', type=int, default=int(time.time()))
parser.add_argument('--gen-gset', action='store_true')
parser.add_argument('--dj-genes', default='IGHD6-19*01:IGHJ4*02', help='.')
parser.add_argument('--sim-v-genes', default='IGHV4-39*01:IGHV4-39*06', help='.')
parser.add_argument('--inf-v-genes', default='IGHV4-39*01', help='.')
parser.add_argument('--snp-positions')
parser.add_argument('--mut-mult', type=float, default=0.5)
parser.add_argument('--slurm', action='store_true')
parser.add_argument('--outdir', default=fsdir + '/partis/allele-finder')
parser.add_argument('--workdir', default=fsdir + '/_tmp/hmms/' + str(random.randint(0, 999999)))
parser.add_argument('--comprehensive', action='store_true')
parser.add_argument('--n-tests', type=int, default=3)
args = parser.parse_args()
args.dj_genes = utils.get_arg_list(args.dj_genes)
args.sim_v_genes = utils.get_arg_list(args.sim_v_genes)
args.inf_v_genes = utils.get_arg_list(args.inf_v_genes)
args.snp_positions = utils.get_arg_list(args.snp_positions, intify=True)

if args.comprehensive:
    comprehensive_test(args)
else:
    run_test(args)
