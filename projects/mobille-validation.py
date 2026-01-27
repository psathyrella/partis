#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import sys
import csv
from io import open
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import argparse
import colored_traceback.always
import glob
import itertools
import json
import collections
import numpy
import math

import partis.utils as utils
import partis.glutils as glutils
import partis.hutils as hutils

mdir = 'packages/MobiLLe/Data/Simulated_datasets'
base_odir = '/fh/fast/matsen_e/dralph/partis/mobille-validation'

# ----------------------------------------------------------------------------------------
def get_true_ptn(spval, stype):
    tcfn = '%s/True_cluster_by_simulator/%s_%s_true_cluster.txt' % (mdir, spval, stype)
    true_partition = []
    with open(tcfn) as tcfile:
        reader = csv.DictReader(tcfile, delimiter=str('\t'), fieldnames=['iclust', 'uids'])
        for line in reader:
            cluster = [u.strip() for u in line['uids'].split()]
            true_partition.append(cluster)
    return true_partition
# ----------------------------------------------------------------------------------------
def fst_infn(spval, stype):
    return '%s/Fasta/%s_%s_simulated.fasta' % (mdir, spval, stype)
# ----------------------------------------------------------------------------------------
def imgt_odir():  # imgt output files, from the mobille repo
    return '%s/IMGT_highvquest_output' % mdir
# ----------------------------------------------------------------------------------------
def bodir(spval, stype):
    return '%s/%s/%s-%s' % (base_odir, args.version, spval, stype)
# ----------------------------------------------------------------------------------------
def bmdir(spval, stype, mthd, iseed=None):
    return '%s/%s%s' % (bodir(spval, stype), mthd, '' if args.n_random_seeds is None or iseed is None else '/seed-%d'%iseed)
# ----------------------------------------------------------------------------------------
def paramdir(spval, stype, iseed=None):
    return '%s/parameters' % bmdir(spval, stype, 'partis', iseed=iseed)
# ----------------------------------------------------------------------------------------
def ptnfn(spval, stype, mthd, iseed=None):
    return '%s/partition.yaml' % bmdir(spval, stype, mthd, iseed=iseed)
# ----------------------------------------------------------------------------------------
def pltdir(simu=False):
    return '%s/%s/%splots' % (base_odir, args.version, 'simu-' if simu else '')
# ----------------------------------------------------------------------------------------
def mtrfn(spval, stype, mthd):
    return '%s/metrics.yaml' % bmdir(spval, stype, mthd)

# ----------------------------------------------------------------------------------------
def write_metrics(spval, stype, mthd, debug=False):
    # ----------------------------------------------------------------------------------------
    def addval(mvals, mtype, mname, val):
        if mtype not in mvals:
            mvals[mtype] = {}
        if mname not in mvals[mtype]:
            mvals[mtype][mname] = []
        mvals[mtype][mname].append(val)
    # ----------------------------------------------------------------------------------------
    ofn = mtrfn(spval, stype, mthd)
    if utils.output_exists(args, ofn, 'metric'):
        return
    tru_ptn = get_true_ptn(spval, stype)
    mvals = {}
    for iseed in range(1 if args.n_random_seeds is None else args.n_random_seeds):
        _, _, cpath = utils.read_output(ptnfn(spval, stype, mthd, iseed=iseed))
        inf_ptn = cpath.best()
        # import partis.plotting as plotting
        # print plotting.plot_cluster_similarity_matrix(pltdir()+'/csim-plots', 'csim-matrix-%s-%s%s'%(stype, spval, '' if args.n_random_seeds is None else '-seed-%d'%iseed), 'true', tru_ptn, mthd, inf_ptn, 100, debug=True) #, debug=True)
        vdict = {}
        for mtstr in ['pairwise', 'closeness']:
            for mname, mval in utils.pairwise_cluster_metrics(mtstr, inf_ptn, tru_ptn).items():
                addval(mvals, 'mobille-%s'%mtstr, mname, mval)
        ccfs = utils.per_seq_correct_cluster_fractions(inf_ptn, tru_ptn) #, debug=True)
        addval(mvals, 'partis', 'purity', ccfs[0])
        addval(mvals, 'partis', 'completeness', ccfs[1])
        addval(mvals, 'partis', 'f1', 2 * ccfs[0] * ccfs[1] / float(sum(ccfs)))  # same as scipy.stats.hmean([precis, recall])
    for mname in mvals:
        for mtype, mvlist in mvals[mname].items():
            mvals[mname][mtype] = {'mean' : numpy.mean(mvlist), 'err' : numpy.std(mvlist, ddof=1) / math.sqrt(len(mvlist)), 'vals' : mvlist}
            if debug:
                print('    %17s %12s %.4f +/- %.4f   %s' % (mname, mtype, mvals[mname][mtype]['mean'], mvals[mname][mtype]['err'], ' '.join('%.4f'%v for v in mvlist)))
    utils.mkdir(ofn, isfile=True)
    print('  writing metrics to %s' % ofn)
    utils.jsdump(ofn, mvals)

# ----------------------------------------------------------------------------------------
def run_method(mthd, spval, stype, iseed=0):
    ofn = ptnfn(spval, stype, mthd, iseed=iseed)
    if utils.output_exists(args, ofn, mthd):
        return

    if mthd == 'partis':
        for action in ['cache-parameters', 'partition']:
            cmd = './bin/partis %s --infname %s --parameter-dir %s --debug-allele-finding --random-seed %d' % (action, fst_infn(spval, stype), paramdir(spval, stype, iseed=iseed), iseed)
            if action == 'partition':
                cmd += ' --outfname %s' % ofn
            # else:
            #     cmd += ' --plotdir %s/plots' % os.path.dirname(ofn)
            utils.simplerun(cmd, logfname='%s/%s.log'%(os.path.dirname(ofn), action)) #, dryrun=True)
    else:
        if mthd == 'mobille':
            cmd = './test/mobille-igblast-run.py mobille --inpath %s --outdir %s --base-imgt-outdir %s --single-chain --id-str %s_%s' % (fst_infn(spval, stype), os.path.dirname(ofn), imgt_odir(), spval, stype)
        elif mthd == 'scoper':
            cmd = './test/scoper-run.py --indir %s --outdir %s --single-chain' % (paramdir(spval, stype, iseed=iseed), os.path.dirname(ofn))
        else:
            assert False
        utils.simplerun(cmd, logfname='%s/%s.log'%(os.path.dirname(ofn), mthd)) #, dryrun=True)
# ----------------------------------------------------------------------------------------
def lzv(lvstr):
    assert lvstr[:3] == 'l00' and len(lvstr) == 5
    return float(int(lvstr[3:5])) / 100

# ----------------------------------------------------------------------------------------
def make_plots(swarm=False, debug=False):
    # ----------------------------------------------------------------------------------------
    def plot_metric_type(mtr_type):
        # ----------------------------------------------------------------------------------------
        def read_files(mthd):
            plotvals = {m : {t : collections.OrderedDict([[v, None] for v in spvals]) for t in stypes} for m in args.methods}
            for stype in stypes:
                for spval in spvals:
                    if not os.path.exists(mtrfn(spval, stype, mthd)):
                        continue
                    for mthd in args.methods:
                        with open(mtrfn(spval, stype, mthd)) as mfile:
                            vals = json.load(mfile)
                            plotvals[mthd][stype][spval] = vals[mtr_type]['f1']['vals'] if swarm else [vals[mtr_type]['f1'][k] for k in ('mean', 'err')]  # NOTE mtr_type is 'partis' or 'mobille', i.e. has the same values as args.methods
            return plotvals
        # ----------------------------------------------------------------------------------------
        for ist, stype in enumerate(stypes):
            fig, ax = plotting.mpl_init()
            for mthd in args.methods:
                plotvals = read_files(mthd)
                if debug:
                    print('  %s' % mthd)
                xvals, yvals = list(zip(*list(plotvals[mthd][stype].items())))
                if debug:
                    if ist==0:
                        print('  %-18s %2s %s' % (mtr_type, '', '  '.join('%5s'%lzv(v) for v in xvals)))
                    # print '    %8s  %8s  %s' % (stype, mthd, '  '.join('%.3f'%v for v in yvals))
                xvals = [lzv(v) for v in xvals]
                if not swarm:
                    yvals, yerrs = list(zip(*yvals))
                if args.n_random_seeds is None:
                    ax.plot(xvals, yvals, label=mthd, alpha=0.6, linewidth=3, markersize=13, marker='.', color=method_colors.get(mthd))
                else:
                    if swarm:
                        import seaborn as sns
                        sns.swarmplot(data=yvals)
                        # sns.boxplot(data=yvals, boxprops={'alpha' : 0.6}, showmeans=True, meanprops={'marker' : 'o', 'markerfacecolor' : 'white', 'markeredgecolor' : 'black'})
                    else:
                        ax.errorbar(xvals, yvals, yerr=yerrs, label=mthd, alpha=0.6, linewidth=3, markersize=13, marker='.')
            if swarm:
                ax.plot((0, 4) if swarm else (xvals[0], xvals[-1]), (1, 1), linewidth=1.5, alpha=0.5, color='grey') #, linestyle='--') #, label='1/seq len')
            fn = plotting.mpl_finish(ax, pltdir(), 'f1-%s-%s'%(stype, 'f1' if mtr_type=='partis' else mtr_type.split('-')[1]), ylabel=metric_labels.get(mtr_type, mtr_type), title='%sclonal'%stype,
                                     xlabel='lambda 0', xticks=None if swarm else xvals, xticklabels=['%.2f'%v for v in xvals], ybounds=(0, 1.05), leg_loc=(0.1, 0.2))
            fnames[0 if 'pairwise' in mtr_type else (2 if mtr_type=='partis' else 1)].append(fn)
    # ----------------------------------------------------------------------------------------
    import partis.plotting as plotting
    utils.prep_dir(pltdir(), wildlings=['*.csv', '*.svg'])
    fnames = [[], [], []]
    for mtstr in ['pairwise', 'closeness']:
        plot_metric_type('mobille-%s'%mtstr)
    plot_metric_type('partis')
    plotting.make_html(pltdir(), fnames=fnames)

# ----------------------------------------------------------------------------------------
def plot_simulation():
    import partis.plotting as plotting
    utils.prep_dir(pltdir(simu=True), wildlings=['*.csv', '*.svg'])
    for stype in stypes:
        # csize_hists = {}
        for spval in spvals:
            tru_ptn = get_true_ptn(spval, stype)
            bubfos = [{'id' : str(i), 'radius' : len(c), 'texts' : [{'tstr' : str(len(c)), 'fsize' : 12, 'tcol' : 'black'}], 'fracs' : None} for i, c in enumerate(tru_ptn)]
            fn = plotting.bubble_plot('%s-%s-bubbles'%(stype, spval), pltdir(simu=True), bubfos, title='%s: %s'%(stype, lzv(spval)))
            # cslist = [len(c) for c in tru_ptn]
            # csize_hists[spval] = hutils.make_hist_from_list_of_values(cslist, 'int', lzv(spval), is_log_x=True)
        # fn = plotting.plot_cluster_size_hists(pltdir(simu=True), '%s-cluster-size'%stype, csize_hists, title='%sclonal'%stype)
    plotting.make_html(pltdir(simu=True), n_columns=4) #, fnames=fnames)

# ----------------------------------------------------------------------------------------
# ./projects/mobille-validation.py --version v1 # --actions simplot  --methods mobille:scoper:partis
# ./projects/mobille-validation.py --version trand --actions plot --n-random-seeds 10
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='run:write:plot')  # also: simplot
parser.add_argument('--methods', default='partis:mobille:scoper')
parser.add_argument('--version', default='test')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--n-random-seeds', type=int, help='number of replicates with different random seeds to run')
args = parser.parse_args()
args.actions = utils.get_arg_list(args.actions)
args.methods = utils.get_arg_list(args.methods)

spvals = ['l00%d'%c for c in [16, 26, 36, 46]]
stypes = ['mono', 'oligo', 'poly']
metric_labels = {'mobille-pairwise' : 'pairwise f1', 'mobille-closeness' : 'closeness f1', 'partis' : 'partis f1'}
method_colors = {'mobille' : '#debd36', 'partis' : '#c532b8', 'scoper' : '#62d156'}

for stype in stypes:
    for spval in spvals:
        # if stype != 'mono' or spval != 'l0036':
        #     continue
        for mthd in args.methods:
            if 'run' in args.actions:
                for iseed in range(1 if args.n_random_seeds is None else args.n_random_seeds):
                    run_method(mthd, spval, stype, iseed=iseed)
            if 'write' in args.actions:
                write_metrics(spval, stype, mthd)
if 'plot' in args.actions:
    make_plots(swarm=args.n_random_seeds is not None)
if 'simplot' in args.actions:
    plot_simulation()
