#!/usr/bin/env python
from collections import OrderedDict
import os
import random
import argparse
import sys
import subprocess
sys.path.insert(1, './python')

import utils
import glutils
from hist import Hist

fsdir = '/fh/fast/matsen_e'
alfdir = fsdir + '/dralph/partis/allele-finder'
locus = 'igh'
region = 'v'

# # ----------------------------------------------------------------------------------------
# sys.path.insert(1, './datascripts')
# import heads
# label = 'vz'
# ptype = 'sw'
# subject = None  # 'GMC'
# studies = [
#     'kate-qrs-2016-09-09',
#     'laura-mb-2016-12-22',
#     'chaim-donor-45-2016-08-04',
#     'adaptive-billion-read-2016-04-07',
#     'vollmers-2016-04-08',
#     # 'jason-mg-2017-02-01',
#     # 'jason-influenza-2017-02-03',
# ]

# merged_names, merged_dirs = [], []
# for study in studies:
#     names, dirs = [], []
#     metafo = heads.read_metadata(study)
#     print study
#     for dset in metafo:
#         if subject is not None and metafo[dset]['subject'] != subject:
#             continue
#         bdir = fsdir + '/processed-data/partis/' + study + '/' + label + '/' + dset
#         if metafo[dset]['timepoint'] == 'merged':
#             continue
#         if dset == 'Hs-LN3-5RACE-IgG':  # bad one
#             continue
#         if not os.path.exists(bdir):
#             print '    %s missing' % dset
#             continue
#         names.append(metafo[dset]['shorthand'])
#         dirs.append(bdir + '/plots/' + ptype + '/mute-freqs/overall')
#     outdir = fsdir + '/dralph/partis/tmp/lots-of-mfreqs/' + study
#     if subject is not None:
#         outdir += '/' + subject
#     subprocess.check_call(['./bin/compare-plotdirs.py', '--outdir', outdir, '--plotdirs', ':'.join(dirs), '--names', ':'.join(names), '--normalize'])
#     merged_names += names
#     merged_dirs += dirs
# # subprocess.check_call(['./bin/compare-plotdirs.py', '--outdir', fsdir + '/dralph/partis/tmp/lots-of-mfreqs/merged', '--plotdirs', ':'.join(merged_dirs), '--names', ':'.join(merged_names), '--normalize'])
# sys.exit()

# # ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print '%s %s' % (utils.color('red', 'run'), cmd_str)
    sys.stdout.flush()
    # subprocess.Popen(cmd_str.split())
    subprocess.check_call(cmd_str.split())

# ----------------------------------------------------------------------------------------
def get_outdir(baseoutdir, n_events, varname, varval):
    outdir = baseoutdir
    outdir += '/' + varname + '-' + str(varval)
    return outdir + '/n-events-' + str(n_events).replace('000', 'k')

# ----------------------------------------------------------------------------------------
def get_single_performance(outdir, debug=False):
    sglfo = glutils.read_glfo(outdir + '/germlines/simulation', locus=locus)
    iglfo = glutils.read_glfo(outdir + '/simu-test/sw/germline-sets', locus=locus)
    missing_alleles = set(sglfo['seqs'][region]) - set(iglfo['seqs'][region])
    spurious_alleles = set(iglfo['seqs'][region]) - set(sglfo['seqs'][region])
    if debug:
        if len(missing_alleles) > 0:
            print '    %2d  missing %s' % (len(missing_alleles), ' '.join([utils.color_gene(g) for g in missing_alleles]))
        if len(spurious_alleles) > 0:
            print '    %2d spurious %s' % (len(spurious_alleles), ' '.join([utils.color_gene(g) for g in spurious_alleles]))
        if len(missing_alleles) == 0 and len(spurious_alleles) == 0:
            print '    none missing'
    return {'missing' : len(missing_alleles), 'spurious' : len(spurious_alleles)}

# ----------------------------------------------------------------------------------------
def plot_test(args, baseoutdir, varvals, debug=False):
    import plotting
    plot_types = ['missing', 'spurious']

    def get_performance(varname, varval):
        perf_vals = {pt : [] for pt in plot_types}
        for iproc in range(args.n_tests):
            single_vals = get_single_performance(get_outdir(baseoutdir, n_events, varname, varval) + '/' + str(iproc))
            for ptype in plot_types:
                perf_vals[ptype].append(single_vals[ptype])
        for ptype in plot_types:
            assert len(set(perf_vals[ptype]) - set([0, 1])) == 0  # should only be zeroes and/or ones
        return perf_vals

    plotvals = []
    for varval in varvals:
        print '%s %s' % (args.action, varval)
        plotvals.append({pt : {k : [] for k in ['xvals', 'ycounts', 'ytotals']} for pt in plot_types})
        for n_events in args.n_event_list:
            if debug:
                print '    %d' % n_events
            perf_vals = get_performance(varname=args.action, varval=varval)
            for ptype in plot_types:
                count = perf_vals[ptype].count(1)
                plotvals[-1][ptype]['xvals'].append(n_events)
                plotvals[-1][ptype]['ycounts'].append(count)
                plotvals[-1][ptype]['ytotals'].append(len(perf_vals[ptype]))
    for ptype in plot_types:
        plotting.plot_gl_inference_fractions(baseoutdir, ptype, [pv[ptype] for pv in plotvals], labels=varvals, xlabel='sample size', ylabel='fraction %s' % ptype)

# ----------------------------------------------------------------------------------------
def get_base_cmd(args, n_events):
    cmd = './bin/test-allele-finding.py'
    cmd += ' --n-procs 5 --n-tests ' + str(args.n_tests) + ' --slurm'
    cmd += ' --inf-v-genes ' + args.v_genes[0]
    cmd += ' --n-sim-events ' + str(n_events)
    return cmd

# ----------------------------------------------------------------------------------------
def run_mfreq_test(args, baseoutdir):
    for mfreq in args.mfreqs:
        for n_events in args.n_event_list:
            cmd = get_base_cmd(args, n_events)
            cmd += ' --sim-v-genes ' + args.v_genes[0]
            cmd += ' --mut-mult ' + str(mfreq)
            cmd += ' --nsnp-list 1'
            cmd += ' --outdir ' + get_outdir(baseoutdir, n_events, args.action, mfreq)
            run(cmd)

# ----------------------------------------------------------------------------------------
def run_nsnp_test(args, baseoutdir):
    for nsnp in args.nsnp_list:
        for n_events in args.n_event_list:
            cmd = get_base_cmd(args, n_events)
            nsnpstr = nsnp
            sim_v_genes = [args.v_genes[0]]
            if args.action == 'multi-nsnp':
                nsnpstr = ':'.join([str(n) for n in nsnp])
                sim_v_genes *= len(nsnp)
            cmd += ' --sim-v-genes ' + ':'.join(sim_v_genes)
            cmd += ' --nsnp-list ' + nsnpstr
            cmd += ' --outdir ' + get_outdir(baseoutdir, n_events, args.action, nsnpstr)
            run(cmd)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['mfreq', 'nsnp', 'multi-nsnp'])
parser.add_argument('--nsnp-list', default='1:2:3')
parser.add_argument('--mfreqs', default='0.1:1:2')
parser.add_argument('--n-event-list', default='5000')
parser.add_argument('--n-tests', type=int, default=5)
parser.add_argument('--plot', action='store_true')
parser.add_argument('--v-genes', default='IGHV4-39*01')
parser.add_argument('--label')
args = parser.parse_args()

if args.action == 'nsnp':
    args.nsnp_list = utils.get_arg_list(args.nsnp_list, intify=True)
elif args.action == 'multi-nsnp':  # list of nsnps for each test, e.g. '1,1:2,2' runs two tests: 1) two new alleles, each with one snp and 2) two new alleles each with 2 snps
    args.nsnp_list = [[int(n) for n in gstr.split(',')] for gstr in utils.get_arg_list(args.nsnp_list)]
args.mfreqs = utils.get_arg_list(args.mfreqs, floatify=True)
args.n_event_list = utils.get_arg_list(args.n_event_list, intify=True)
args.v_genes = utils.get_arg_list(args.v_genes)

original_glfo = glutils.read_glfo('data/germlines/human', locus=locus)

# ----------------------------------------------------------------------------------------
baseoutdir = alfdir
if args.label is not None:
    baseoutdir += '/' + args.label
baseoutdir += '/' + args.action

if args.action == 'mfreq':
    if args.plot:
        plot_test(args, baseoutdir, varvals=args.mfreqs)
    else:
        run_mfreq_test(args, baseoutdir)
elif args.action == 'nsnp' or args.action == 'multi-nsnp':
    if args.plot:
        plot_test(args, baseoutdir, varvals=args.nsnp_list)
    else:
        run_nsnp_test(args, baseoutdir)
