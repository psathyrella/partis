#!/usr/bin/env python
import numpy
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
        cmd_str = base_cmd + ' simulate --n-sim-events ' + str(args.n_sim_events) + ' --n-leaves ' + str(args.n_leaves) + ' --constant-number-of-leaves --rearrange-from-scratch --outfname ' + simfname
        cmd_str += ' --mutation-multiplier ' + str(args.mut_mult)

        cmd_str += ' --n-procs ' + str(args.n_procs)
        if args.slurm:
            cmd_str += ' --batch-system slurm --subsimproc'

        if args.gen_gset:
            cmd_str += ' --generate-germline-set'
            cmd_str += ' --n-genes-per-region 1:5:3'
            cmd_str += ' --n-alleles-per-gene 2,3:1,2:1,2'
        else:
            simulation_genes = ':'.join(args.sim_v_genes + args.dj_genes)
            sglfo = glutils.read_glfo('data/germlines/human', chain=chain, only_genes=simulation_genes.split(':'))

            added_snp_names = None
            if args.snp_positions is not None:
                snps_to_add = [{'gene' : args.sim_v_genes[ig], 'positions' : args.snp_positions[ig]} for ig in range(len(args.sim_v_genes))]
                added_snp_names = glutils.add_some_snps(snps_to_add, sglfo, debug=True, remove_template_genes=args.remove_template_genes)

            if args.allele_prevalence_freqs is not None:
                if len(args.allele_prevalence_freqs) != len(sglfo['seqs']['v']):
                    raise Exception('--allele-prevalence-freqs not the right length')
                gene_list = sorted(sglfo['seqs']['v']) if added_snp_names is None else list(set(args.sim_v_genes)) + added_snp_names
                prevalence_freqs = {'v' : {g : f for g, f in zip(gene_list, args.allele_prevalence_freqs)}, 'd' : {}, 'j' : {}}
                glutils.write_allele_prevalence_freqs(prevalence_freqs, args.workdir + '/allele-prevalence-freqs.csv')
                cmd_str += ' --allele-prevalence-fname ' + args.workdir + '/allele-prevalence-freqs.csv'

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
    cmd_str = base_cmd + ' cache-parameters --infname ' + simfname + ' --only-smith-waterman --debug-allele-finding --always-find-new-alleles --n-max-allele-finding-iterations 2'
    # cmd_str = 'python -m cProfile -s tottime -o prof.out ' + cmd_str
    cmd_str += ' --n-procs ' + str(args.n_procs)
    if args.slurm:
        cmd_str += ' --batch-system slurm'

    if args.gen_gset:
        cmd_str += ' --find-new-alleles'
    else:
        inference_genes = ':'.join(args.inf_v_genes + args.dj_genes)
        iglfo = glutils.read_glfo('data/germlines/human', chain=chain, only_genes=inference_genes.split(':'), debug=True)
        glutils.write_glfo(args.outdir + '/germlines/inference', iglfo)
        cmd_str += ' --initial-germline-dir ' + args.outdir + '/germlines/inference'
        cmd_str += ' --find-new-alleles --dont-remove-unlikely-alleles'  # --new-allele-fname ' + args.outdir + '/new-alleles.fa'
        # cmd_str += ' --n-max-snps 12'

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
parser.add_argument('--n-sim-events', type=int, default=20)
parser.add_argument('--n-leaves', type=int, default=1)
parser.add_argument('--n-procs', type=int, default=2)
parser.add_argument('--seed', type=int, default=int(time.time()))
parser.add_argument('--gen-gset', action='store_true')
parser.add_argument('--dj-genes', default='IGHD6-19*01:IGHJ4*02', help='.')
parser.add_argument('--sim-v-genes', default='IGHV4-39*01:IGHV4-39*06', help='.')
parser.add_argument('--inf-v-genes', default='IGHV4-39*01', help='.')
parser.add_argument('--snp-positions')
parser.add_argument('--remove-template-genes', action='store_true')
parser.add_argument('--mut-mult', type=float, default=0.5)
parser.add_argument('--slurm', action='store_true')
parser.add_argument('--outdir', default=fsdir + '/partis/allele-finder')
parser.add_argument('--workdir', default=fsdir + '/_tmp/hmms/' + str(random.randint(0, 999999)))
parser.add_argument('--comprehensive', action='store_true')
parser.add_argument('--n-tests', type=int, default=3)
parser.add_argument('--allele-prevalence-freqs')
args = parser.parse_args()
args.dj_genes = utils.get_arg_list(args.dj_genes)
args.sim_v_genes = utils.get_arg_list(args.sim_v_genes)
args.inf_v_genes = utils.get_arg_list(args.inf_v_genes)
args.snp_positions = utils.get_arg_list(args.snp_positions)
args.allele_prevalence_freqs = utils.get_arg_list(args.allele_prevalence_freqs, floatify=True)
if args.snp_positions is not None:
    args.snp_positions = [[int(p) for p in pos_str.split(',')] for pos_str in args.snp_positions]
    assert len(args.snp_positions) == len(args.sim_v_genes)
    # args.snp_positions = {args.sim_v_genes[ig] : args.snp_positions[ig] for ig in range(len(args.sim_v_genes))}

if args.seed is not None:
    random.seed(args.seed)
    numpy.random.seed(args.seed)

if args.comprehensive:
    comprehensive_test(args)
else:
    run_test(args)
