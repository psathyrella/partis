#!/usr/bin/env python
import sys
import csv
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

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/projects', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
import hutils

mdir = 'packages/MobiLLe/Data/Simulated_datasets'
base_odir = '/fh/fast/matsen_e/dralph/partis/mobille-validation'

# ----------------------------------------------------------------------------------------
def get_true_ptn(spval, stype):
    tcfn = '%s/True_cluster_by_simulator/%s_%s_true_cluster.txt' % (mdir, spval, stype)
    true_partition = []
    with open(tcfn) as tcfile:
        reader = csv.DictReader(tcfile, delimiter='\t', fieldnames=['iclust', 'uids'])
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
def mb_metrics(mtstr, inf_ptn, tru_ptn, debug=False):
    # ----------------------------------------------------------------------------------------
    def id_dict(ptn):
        reco_info = utils.build_dummy_reco_info(ptn)  # not actually reco info unless it's the true partition
        return {uid : reco_info[uid]['reco_id'] for cluster in ptn for uid in cluster}  # speed optimization
    # ----------------------------------------------------------------------------------------
    # inf_ptn, tru_ptn = [['c'], ['a', 'b', 'e'], ['d', 'g', 'f']], [['a', 'b', 'c'], ['d', 'e', 'f', 'g']]  # example from paper, should be pairwise: (0.666666, 0.44444444, ?), closeness: (0.857, 0.6, ?) (see note below, they calculate it wrong)
    # inf_ptn, tru_ptn = [['a'], ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n']], [['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n']]
    utils.check_intersection_and_complement(inf_ptn, tru_ptn, a_label='true', b_label='inferred')
    if mtstr == 'pairwise':
        tp, fp, fn, n_tot = 0, 0, 0, 0
        tru_ids, inf_ids = [id_dict(ptn) for ptn in [tru_ptn, inf_ptn]]
        for u1, u2 in itertools.combinations(set(u for c in tru_ptn for u in c), 2):
            is_tru_clonal, is_inf_clonal = [tids[u1] == tids[u2] for tids in [tru_ids, inf_ids]]
            n_tot += 1
            if is_tru_clonal and is_inf_clonal:
                tp += 1
            elif is_tru_clonal:
                fn += 1
            elif is_inf_clonal:
                fp += 1
            else:  # singletons
                pass
    elif mtstr == 'closeness':
        tp, fp, fn, n_tot = set(), set(), set(), set()
        if debug:
            print '    infcl   trucl      tp   fp     fn'
        # for infcl, trucl in itertools.product(inf_ptn, tru_ptn):  # NOTE this is *not* what they mean, see next line
        for infcl in inf_ptn:  # NOTE this is apparently what they mean by "we first identified the best correspondence between inferred clonal lineages and correct clonal assignments" but note you'd get a *different* answer if you looped over tru_ptn
            trucl = sorted(tru_ptn, key=lambda c: len(set(c) & set(infcl)), reverse=True)[0]
            infset, truset = set(infcl), set(trucl)
            if len(infset & truset) == 0:
                continue
            tp |= infset & truset  # OMFG their example figure is wrong, it (correctly) shows 6 entries, but then when they calculate recall they switch it to 5
            fp |= infset - truset
            fn |= truset - infset
            n_tot |= truset
            if debug:
                print '  %20s   %20s   %20s   %20s   %20s' % (infcl, trucl, infset & truset, infset - truset, truset - infset)
        tp, fp, fn = [len(s) for s in [tp, fp, fn]]
    else:
        assert False
    precis = tp / float(tp + fp)
    recall = tp / float(tp + fn)
    return {'precision' : precis, 'recall' : recall, 'f1' : 2 * precis * recall / float(precis + recall)}  # same as scipy.stats.hmean([precis, recall])

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
        # import plotting
        # print plotting.plot_cluster_similarity_matrix(pltdir()+'/csim-plots', 'csim-matrix-%s-%s%s'%(stype, spval, '' if args.n_random_seeds is None else '-seed-%d'%iseed), 'true', tru_ptn, mthd, inf_ptn, 100, debug=True) #, debug=True)
        vdict = {}
        for mtstr in ['pairwise', 'closeness']:
            for mname, mval in mb_metrics(mtstr, inf_ptn, tru_ptn).items():
                addval(mvals, 'mobille-%s'%mtstr, mname, mval)
        ccfs = utils.per_seq_correct_cluster_fractions(inf_ptn, tru_ptn) #, debug=True)
        addval(mvals, 'partis', 'purity', ccfs[0])
        addval(mvals, 'partis', 'completeness', ccfs[1])
        addval(mvals, 'partis', 'f1', 2 * ccfs[0] * ccfs[1] / float(sum(ccfs)))  # same as scipy.stats.hmean([precis, recall])
    for mname in mvals:
        for mtype, mvlist in mvals[mname].items():
            mvals[mname][mtype] = {'mean' : numpy.mean(mvlist), 'err' : numpy.std(mvlist, ddof=1) / math.sqrt(len(mvlist)), 'vals' : mvlist}
            if debug:
                print '    %17s %12s %.4f +/- %.4f   %s' % (mname, mtype, mvals[mname][mtype]['mean'], mvals[mname][mtype]['err'], ' '.join('%.4f'%v for v in mvlist))
    utils.mkdir(ofn, isfile=True)
    print '  writing metrics to %s' % ofn
    with open(ofn, 'w') as mfile:
        json.dump(mvals, mfile)

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
                    print '  %s' % mthd
                xvals, yvals = zip(*plotvals[mthd][stype].items())
                if debug:
                    if ist==0:
                        print '  %-18s %2s %s' % (mtr_type, '', '  '.join('%5s'%lzv(v) for v in xvals))
                    # print '    %8s  %8s  %s' % (stype, mthd, '  '.join('%.3f'%v for v in yvals))
                xvals = [lzv(v) for v in xvals]
                if not swarm:
                    yvals, yerrs = zip(*yvals)
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
    import plotting
    utils.prep_dir(pltdir(), wildlings=['*.csv', '*.svg'])
    fnames = [[], [], []]
    for mtstr in ['pairwise', 'closeness']:
        plot_metric_type('mobille-%s'%mtstr)
    plot_metric_type('partis')
    plotting.make_html(pltdir(), fnames=fnames)

# ----------------------------------------------------------------------------------------
def plot_simulation():
    import plotting
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
