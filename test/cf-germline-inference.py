#!/usr/bin/env python
import copy
import numpy
from collections import OrderedDict
import os
import random
import argparse
import sys
import subprocess
import colored_traceback.always
sys.path.insert(1, './python')

import utils
import glutils
from hist import Hist
sys.path.insert(1, './datascripts')
import heads

sim_locus = 'igh'
region = 'v'

legend_titles = {
    'mfreq' : 'mutation',
    'nsnp' : 'N SNPs',
    'multi-nsnp' : 'N SNPs',
    'prevalence' : 'prevalence',
    'n-leaves' : 'mean N leaves',
}

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
#         bdir = utils.fsdir() + '/processed-data/partis/' + study + '/' + label + '/' + dset
#         if metafo[dset]['timepoint'] == 'merged':
#             continue
#         if dset == 'Hs-LN3-5RACE-IgG':  # bad one
#             continue
#         if not os.path.exists(bdir):
#             print '    %s missing' % dset
#             continue
#         names.append(metafo[dset]['shorthand'])
#         dirs.append(bdir + '/plots/' + ptype + '/mute-freqs/overall')
#     outdir = utils.fsdir() + '/partis/tmp/lots-of-mfreqs/' + study
#     if subject is not None:
#         outdir += '/' + subject
#     subprocess.check_call(['./bin/compare-plotdirs.py', '--outdir', outdir, '--plotdirs', ':'.join(dirs), '--names', ':'.join(names), '--normalize'])
#     merged_names += names
#     merged_dirs += dirs
# # subprocess.check_call(['./bin/compare-plotdirs.py', '--outdir', utils.fsdir() + '/partis/tmp/lots-of-mfreqs/merged', '--plotdirs', ':'.join(merged_dirs), '--names', ':'.join(merged_names), '--normalize'])
# sys.exit()

# # ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
def varvalstr(name, val):
    if name == 'multi-nsnp':
        valstr = ':'.join([str(v) for v in val])
    elif name == 'alcluster':
        valstr = '-'.join(['n%s-%s' % (k, v) for k, v in val.items()])
    elif name == 'gls-gen':
        valstr = 'simu'
    else:
        valstr = str(val)
    return valstr

# ----------------------------------------------------------------------------------------
def legend_str(args, val):
    if args.action == 'mfreq':
        # lstr = '%.1fx' % val
        if val < 0.15:
            return 'low'
        elif val > 0.5 and val < 1.5:
            return 'typical'
        elif val > 1.5:
            return 'high'
        else:
            assert False
    elif args.action == 'nsnp':
        lstr = '%d' % val
    elif args.action == 'multi-nsnp':
        lstr = '%s' % '+'.join([str(v) for v in val])
    elif args.action == 'prevalence':
        lstr = '%d%%' % (100*val)
    elif args.action == 'n-leaves':
        lstr = '%.1f' % val
    elif args.action == 'weibull':
        lstr = '%.1f' % val
    elif args.action == 'alcluster':
        lstr = 'er...'
    else:
        assert False
    return lstr

# ----------------------------------------------------------------------------------------
def get_outdir(args, baseoutdir, varname, varval, n_events=None):
    outdir = baseoutdir
    if args.action == 'gls-gen':
        outdir += '/' + varvalstr(varname, varval) + '/' + args.gls_gen_difficulty
    if args.action == 'data':
        outdir += '/' + varvalstr(varname, varval)
    else:
        outdir += '/' + varname + '-' + varvalstr(varname, varval) + '/n-events-' + str(n_events)
    return outdir

# ----------------------------------------------------------------------------------------
def get_single_performance(outdir, method, debug=False):
    sglfo = glutils.read_glfo(outdir + '/germlines/simulation', locus=sim_locus)
    iglfo = glutils.read_glfo(outdir + '/' + method + '/sw/germline-sets', locus=sim_locus)
    glutils.synchronize_glfos(ref_glfo=sglfo, new_glfo=iglfo, region=region)
    missing_alleles = set(sglfo['seqs'][region]) - set(iglfo['seqs'][region])
    spurious_alleles = set(iglfo['seqs'][region]) - set(sglfo['seqs'][region])
    if debug:
        if len(missing_alleles) > 0:
            print '    %2d  missing %s' % (len(missing_alleles), ' '.join([utils.color_gene(g) for g in missing_alleles]))
        if len(spurious_alleles) > 0:
            print '    %2d spurious %s' % (len(spurious_alleles), ' '.join([utils.color_gene(g) for g in spurious_alleles]))
        if len(missing_alleles) == 0 and len(spurious_alleles) == 0:
            print '    none missing'
    return {
        'missing' : len(missing_alleles),
        'spurious' : len(spurious_alleles),
        'total' : len([g for g in sglfo['seqs'][region] if '+' in g]),  # anybody with a '+' should be a new allele
    }

# ----------------------------------------------------------------------------------------
def get_gls_fname(outdir, method, locus, sim_truth=False, data=False):  # NOTE duplicates/depends on code in test-allele-finding.py
    if data:
        if method == 'partis' or method == 'full':
            outdir += '/hmm/germline-sets'  # NOTE this is inside the datascripts output dir, also NOTE doesn't use <method> (since we only have partis for a method a.t.m., although could use --label or --extra-str to differentiate)
        else:
            outdir += '/' + method
    elif sim_truth:
        outdir += '/germlines/simulation'
    elif method == 'partis' or method == 'full':
        outdir += '/' + method + '/sw/germline-sets'
    elif 'tigger' in method or method == 'igdiscover':
        outdir += '/' + method
    else:
        assert False
    return glutils.get_fname(outdir, locus, region)

# ----------------------------------------------------------------------------------------
def make_gls_tree_plot(args, plotdir, plotname, glsfnames, glslabels, locus):
    # ete3 requires its own python version, so we run as a subprocess
    cmdstr = 'export PATH=%s:$PATH && xvfb-run -a ./bin/plot-gl-set-trees.py' % args.ete_path
    cmdstr += ' --plotdir ' + plotdir
    cmdstr += ' --plotname ' + plotname
    cmdstr += ' --glsfnames ' + ':'.join(glsfnames)
    cmdstr += ' --glslabels ' + ':'.join(glslabels)
    cmdstr += ' --locus ' + locus
    if args.plotcache:
        cmdstr += ' --use-cache'
    utils.simplerun(cmdstr, shell=True, debug=False)

# ----------------------------------------------------------------------------------------
def get_gls_gen_plots(args, baseoutdir, method):
    varname = args.action
    varval = 'simu'

    for iproc in range(args.iteststart, args.n_tests):
        outdir = get_outdir(args, baseoutdir, varname, varval, n_events=args.gls_gen_events) + '/' + str(iproc)
        print '%-2d                            %s' % (iproc, outdir)
        simfname = get_gls_fname(outdir, method=None, locus=sim_locus, sim_truth=True)
        inffname = get_gls_fname(outdir, method, sim_locus)
        make_gls_tree_plot(args, outdir + '/' + method + '/gls-gen-plots', varvalstr(varname, varval), glsfnames=[simfname, inffname], glslabels=['sim', 'inf'])

# ----------------------------------------------------------------------------------------
def get_data_plots(args, baseoutdir, method):
    for var in args.varvals:
        study, dset = var.split('/')
        mfo = heads.read_metadata(study)[dset]
        data_outdir = heads.get_datadir(study, 'processed', extra_str='gls-gen-paper-' + args.label) + '/' + dset
        outdir = get_outdir(args, baseoutdir, varname='data', varval=study + '/' + dset)  # for data, only the plots go here, since datascripts puts its output somewhere else
        print '%-15s  %10s    %s' % (study, dset, outdir)
        make_gls_tree_plot(args, outdir + '/' + method + '/gls-gen-plots', study + '-' + dset, glsfnames=[get_gls_fname(data_outdir, method, locus=mfo['locus'], data=True)], glslabels=['data'], locus=mfo['locus'])

# ----------------------------------------------------------------------------------------
def get_data_pair_plots(args, baseoutdir, method, study, dsets):
    mfo = heads.read_metadata(study)[dsets[0]]
    assert heads.read_metadata(study)[dsets[1]]['locus'] == mfo['locus']
    data_outdirs = [heads.get_datadir(study, 'processed', extra_str='gls-gen-paper-' + args.label) + '/' + ds for ds in dsets]
    outdir = get_outdir(args, baseoutdir, varname='data', varval=study + '/' + '-vs-'.join(dsets))  # for data, only the plots go here, since datascripts puts its output somewhere else
    make_gls_tree_plot(args, outdir + '/' + method + '/gls-gen-plots', study + '-' + '-vs-'.join(dsets), glsfnames=[get_gls_fname(dout, method, locus=mfo['locus'], data=True) for dout in data_outdirs], glslabels=dsets, locus=mfo['locus'])

# ----------------------------------------------------------------------------------------
def plot_single_test(args, baseoutdir, method):
    import plotting
    plot_types = ['missing', 'spurious']

    def get_performance(varname, varval):
        perf_vals = {pt : [] for pt in plot_types + ['total']}
        for iproc in range(args.iteststart, args.n_tests):
            single_vals = get_single_performance(get_outdir(args, baseoutdir, varname, varval, n_events=n_events) + '/' + str(iproc), method=method)
            for ptype in plot_types + ['total']:
                perf_vals[ptype].append(single_vals[ptype])
        return perf_vals

    plotvals = []
    for varval in args.varvals:
        print '%s %s' % (args.action, varvalstr(args.action, varval))
        plotvals.append({pt : {k : [] for k in ['xvals', 'ycounts', 'ytotals']} for pt in plot_types})
        for n_events in args.n_event_list:
            perf_vals = get_performance(varname=args.action, varval=varval)
            print '  %d' % n_events
            print '    iproc    %s' % ' '.join([str(i) for i in range(args.iteststart, args.n_tests)])
            print '    missing  %s' % ' '.join([str(v) for v in perf_vals['missing']]).replace('0', ' ')
            print '    spurious %s' % ' '.join([str(v) for v in perf_vals['spurious']]).replace('0', ' ')
            for ptype in plot_types:
                count = sum(perf_vals[ptype])
                plotvals[-1][ptype]['xvals'].append(n_events)
                plotvals[-1][ptype]['ycounts'].append(count)
                plotvals[-1][ptype]['ytotals'].append(sum(perf_vals['total']))
    for ptype in plot_types:
        plotting.plot_gl_inference_fractions(baseoutdir, ptype, [pv[ptype] for pv in plotvals], labels=[legend_str(args, v) for v in args.varvals], xlabel='sample size', ylabel='fraction %s' % ptype, leg_title=legend_titles.get(args.action, None), title=ptype + ' alleles')

# ----------------------------------------------------------------------------------------
def plot_tests(args, baseoutdir, method):
    if args.action == 'gls-gen':
        get_gls_gen_plots(args, baseoutdir, method)
    elif args.action == 'data':
        get_data_plots(args, baseoutdir, method)
        all_dsets = [v.split('/')[1] for v in args.varvals]
        for study in data_pairs:
            for dp in data_pairs[study]:
                ds_1, ds_2 = dp
                if ds_1 in all_dsets and ds_2 in all_dsets:
                    print study, ds_1, ds_2
                    get_data_pair_plots(args, baseoutdir, method, study, [ds_1, ds_2])
    else:
        plot_single_test(args, baseoutdir, method)

# ----------------------------------------------------------------------------------------
def get_base_cmd(args, n_events, method):
    cmd = './bin/test-allele-finding.py'
    cmd += ' --n-procs ' + str(args.n_procs_per_test) + ' --n-tests ' + str(args.n_tests)
    if args.iteststart != 0:
        cmd += ' --iteststart ' + str(args.iteststart)
    cmd += ' --methods ' + method
    cmd += ' --n-sim-events ' + str(n_events)
    if not args.no_slurm:
        cmd += ' --slurm'
    if args.action != 'gls-gen':
        cmd += ' --inf-v-genes ' + args.v_genes[0]
    return cmd

# ----------------------------------------------------------------------------------------
def run_single_test(args, baseoutdir, val, n_events, method):
    cmd = get_base_cmd(args, n_events, method)
    outdir = get_outdir(args, baseoutdir, args.action, val, n_events=n_events)
    sim_v_genes = [args.v_genes[0]]
    mut_mult = None
    nsnpstr, nindelstr = '1', ''
    if args.action == 'mfreq':
        mut_mult = val
    elif args.action == 'nsnp':
        nsnpstr = str(val)
    elif args.action == 'multi-nsnp':
        nsnpstr = ':'.join([str(n) for n in val])
        sim_v_genes *= len(val)
    elif args.action == 'prevalence':
        cmd += ' --allele-prevalence-freqs ' + str(1. - val) + ':' + str(val)  # i.e. previously-known allele has 1 - p, and new allele has p
    elif args.action == 'n-leaves':
        cmd += ' --n-leaves ' + str(val)  # NOTE default of 1 (for other tests) is set in test-allele-finding.py
        cmd += ' --n-leaf-distribution geometric'
        cmd += ' --n-max-queries ' + str(n_events)  # i.e. we simulate <n_events> rearrangement events, but then only use <n_events> sequences for inference
    elif args.action == 'weibull':
        cmd += ' --n-leaves 5'  # NOTE default of 1 (for other tests) is set in test-allele-finding.py
        cmd += ' --n-leaf-distribution geometric'
        cmd += ' --n-max-queries ' + str(n_events)  # i.e. we simulate <n_events> rearrangement events, but then only use <n_events> sequences for inference
    elif args.action == 'alcluster':
        nsnpstr = val['snp']
        nindelstr = val['indel']
        sim_v_genes *= len(val['snp'].split(':'))
    elif args.action == 'gls-gen':
        nsnpstr = '1:2:3:4:10:10'
        nindelstr = '0:0:0:0:3:3'
        if args.gls_gen_difficulty == 'easy':
            genes_per_region_str = '20:5:3'
            n_sim_alleles_per_gene_str = '1,2:1,2:1,2'
            min_allele_prevalence_freq = 0.1
            mut_mult = 0.3
        elif args.gls_gen_difficulty == 'hard':
            genes_per_region_str = '25:5:3'
            n_sim_alleles_per_gene_str = '1,2,3:1,2:1,2'
            min_allele_prevalence_freq = 0.04
            mut_mult = 1.
        else:
            assert False
        cmd += ' --n-genes-per-region ' + genes_per_region_str
        cmd += ' --n-sim-alleles-per-gene ' + n_sim_alleles_per_gene_str
        cmd += ' --min-allele-prevalence-freq ' + str(min_allele_prevalence_freq)
        cmd += ' --gls-gen'
    else:
        assert False

    if mut_mult is not None:
        cmd += ' --mut-mult ' + str(mut_mult)

    if args.action != 'gls-gen':
        cmd += ' --sim-v-genes ' + ':'.join(sim_v_genes)
    if '--nosim' not in cmd:
        if nsnpstr != '':
            cmd += ' --nsnp-list ' + nsnpstr
        if nindelstr != '':
            cmd += ' --nindel-list ' + nindelstr
    cmd += ' --outdir ' + outdir
    utils.simplerun(cmd, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def run_data(args, baseoutdir, study, dset, method):
    cmd = './datascripts/run.py cache-parameters'
    cmd += ' --study ' + study
    cmd += ' --samples ' + dset
    assert args.label is not None  # it's got a default now, so it shouldn't anymore be None
    cmd += ' --extra-str gls-gen-paper-' + args.label
    if args.no_slurm:
        cmd += ' --no-slurm'
    cmd += ' --n-procs ' + str(args.n_procs_per_test)
    if args.n_random_queries is not None:
        assert method == 'partis'  # I don't think it works for any others a.t.m.
        cmd += ' --n-random-queries ' + str(args.n_random_queries)
    if args.check:
        cmd += ' --check'
    if method != 'partis':
        cmd += ' --other-method ' + method

    utils.simplerun(cmd, dryrun=args.dry_run)

# ----------------------------------------------------------------------------------------
def run_tests(args, baseoutdir, method):
    if args.action == 'gls-gen':
        n_events = args.gls_gen_events
        val = 'simu'
        run_single_test(args, baseoutdir, val, n_events, method)
    elif args.action == 'data':
        for var in args.varvals:
            study, dset = var.split('/')
            run_data(args, baseoutdir, study, dset, method)
    else:
        for val in args.varvals:
            for n_events in args.n_event_list:
                run_single_test(args, baseoutdir, val, n_events, method)

# ----------------------------------------------------------------------------------------
default_varvals = {
    'mfreq' : '0.1:1.0:2.0',
    'nsnp' : '1:2:3:4',
    'multi-nsnp' : '1,1:1,3:2,3',
    'prevalence' : '0.1:0.2:0.3',
    'n-leaves' : '1.5:3:10',
    'weibull' : '0.3:0.5:1.3',
    'alcluster' : [
        # {'snp' : 25, 'indel' : 3},
        {'snp' : '75:100:125', 'indel' : '3:4:5'},
    ],
    'gls-gen' : None,
    'data' : {
        'jason-mg' : ['HD07-igh', 'HD07-igk', 'HD07-igl', 'AR03-igh', 'AR03-igk', 'AR03-igl'],
        # 'kate-qrs' : ['1g', '4g'],
        # 'kate-qrs' : ['1g', '1k', '1l', '4g', '4k', '4l', '2k', '2l', '3k', '3l'],
        # 'jason-influenza' : ['FV-igh-m8d', 'FV-igh-p7d', 'FV-igh-p28d'],
    }
}
data_pairs = {'kate-qrs' : [
    ['1g', '4g'],
    ['1k', '4k'],
    ['1l', '4l'],
    ]}
default_varvals['data'] = ':'.join([study + '/' + heads.full_dataset(heads.read_metadata(study), dset) for study in default_varvals['data'] for dset in default_varvals['data'][study]])
for study in data_pairs:
    for idp in range(len(data_pairs[study])):
        data_pairs[study][idp] = [heads.full_dataset(heads.read_metadata(study), ds) for ds in data_pairs[study][idp]]
# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['mfreq', 'nsnp', 'multi-nsnp', 'prevalence', 'n-leaves', 'weibull', 'alcluster', 'gls-gen', 'data'])
parser.add_argument('--methods', default='partis') # not using <choices> 'cause it's harder since it's a list
parser.add_argument('--v-genes', default='IGHV4-39*01')
parser.add_argument('--varvals')
parser.add_argument('--n-event-list', default='1000:2000:4000:8000')  # NOTE modified later for multi-nsnp also NOTE not used for gen-gset
parser.add_argument('--gls-gen-events', type=int, default=300000)
parser.add_argument('--gls-gen-difficulty', choices=['easy', 'hard'])
parser.add_argument('--n-random-queries', type=int)
parser.add_argument('--n-tests', type=int, default=10)
parser.add_argument('--iteststart', type=int, default=0)
parser.add_argument('--n-procs-per-test', type=int, default=5)
parser.add_argument('--plot', action='store_true')
parser.add_argument('--no-slurm', action='store_true')
parser.add_argument('--plotcache', action='store_true')
parser.add_argument('--check', action='store_true')
parser.add_argument('--dry-run', action='store_true')
parser.add_argument('--label', default='xxx')
parser.add_argument('--ete-path', default='/home/' + os.getenv('USER') + '/anaconda_ete/bin')
args = parser.parse_args()

args.methods = utils.get_arg_list(args.methods)
args.v_genes = utils.get_arg_list(args.v_genes)
args.n_event_list = utils.get_arg_list(args.n_event_list, intify=True)

# ----------------------------------------------------------------------------------------
alfdir = utils.fsdir() + '/partis/allele-finder'
baseoutdir = alfdir
if args.label is not None:
    baseoutdir += '/' + args.label
baseoutdir += '/' + args.action

if args.varvals is None:
    args.varvals = default_varvals[args.action]
kwargs = {}
if args.action == 'mfreq' or args.action == 'prevalence' or args.action == 'n-leaves' or args.action == 'weibull':
    kwargs['floatify'] = True
if args.action == 'nsnp':
    kwargs['intify'] = True
if args.action != 'alcluster':  # could also do this for data i think, if i remove that line up there ^
    args.varvals = utils.get_arg_list(args.varvals, **kwargs)
if args.action == 'multi-nsnp':
    args.varvals = [[int(n) for n in gstr.split(',')] for gstr in args.varvals]  # list of nsnps for each test, e.g. '1,1:2,2' runs two tests: 1) two new alleles, each with one snp and 2) two new alleles each with 2 snps
    factor = numpy.median([(len(nl) + 1) / 2. for nl in args.varvals])  # i.e. the ratio of (how many alleles we'll be dividing the events among), to (how many we'd be dividing them among for the other [single-nsnp] tests)
    args.n_event_list = [int(factor * n) for n in args.n_event_list]

for method in args.methods:
    if args.plot:
        if method == 'simu':
            continue
        plot_tests(args, baseoutdir, method)
    else:
        run_tests(args, baseoutdir, method)
