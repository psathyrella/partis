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

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/projects', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils

mdir = 'packages/MobiLLe/Data/Simulated_datasets'
base_odir = '/fh/fast/matsen_e/dralph/partis/mobille-validation'
vsn = 'v0'

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
def bodir(spval, stype):
    return '%s/%s/%s-%s' % (base_odir, vsn, spval, stype)
# ----------------------------------------------------------------------------------------
def bmdir(spval, stype, mthd):
    return '%s/%s' % (bodir(spval, stype), mthd)
# ----------------------------------------------------------------------------------------
def paramdir(spval, stype):
    return '%s/parameters' % bmdir(spval, stype, 'partis')
# ----------------------------------------------------------------------------------------
def ptnfn(spval, stype, mthd):
    return '%s/partition.yaml' % bmdir(spval, stype, mthd)
# ----------------------------------------------------------------------------------------
def pltdir():
    return '%s/%s/plots' % (base_odir, vsn)
# ----------------------------------------------------------------------------------------
def mtrfn(spval, stype):
    return '%s/metrics.yaml' % bodir(spval, stype)

# ----------------------------------------------------------------------------------------
def mb_metrics(mtype, inf_ptn, tru_ptn, debug=False):
    # ----------------------------------------------------------------------------------------
    def id_dict(ptn):
        reco_info = utils.build_dummy_reco_info(ptn)  # not actually reco info unless it's the true partition
        return {uid : reco_info[uid]['reco_id'] for cluster in ptn for uid in cluster}  # speed optimization
    # ----------------------------------------------------------------------------------------
    # inf_ptn, tru_ptn = [['c'], ['a', 'b', 'e'], ['d', 'g', 'f']], [['a', 'b', 'c'], ['d', 'e', 'f', 'g']]  # example from paper, should be pairwise: (0.666666, 0.44444444, ?), closeness: (0.857, 0.6, ?) (see note below, they calculate it wrong)
    # inf_ptn, tru_ptn = [['a'], ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n']], [['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n']]
    utils.check_intersection_and_complement(inf_ptn, tru_ptn, a_label='true', b_label='inferred')
    if mtype == 'pairwise':
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
    elif mtype == 'closeness':
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
def write_metrics(spval, stype, mthd):
    ofn = mtrfn(spval, stype)
    if utils.output_exists(args, ofn, 'metric'):
        return
    tru_ptn = get_true_ptn(spval, stype)
    _, _, cpath = utils.read_output(ptnfn(spval, stype, mthd))
    inf_ptn = cpath.best()
    # import plotting
    # print plotting.plot_cluster_similarity_matrix(pltdir(spval, stype, mthd), 'csim-matrix', 'true', true_partition, 'partis', inf_ptn, 30) #, debug=True)
    vals = {}
    for mtype in ['pairwise', 'closeness']:
        vals['mobille-%s'%mtype] = mb_metrics(mtype, inf_ptn, tru_ptn)
    ccfs = utils.per_seq_correct_cluster_fractions(inf_ptn, tru_ptn) #, debug=True)
    vals['partis'] = {'purity' : ccfs[0], 'completeness' : ccfs[1]}
    vals['partis']['f1'] = 2 * ccfs[0] * ccfs[1] / float(sum(ccfs))  # same as scipy.stats.hmean([precis, recall])
    utils.mkdir(ofn, isfile=True)
    with open(ofn, 'w') as mfile:
        json.dump(vals, mfile)

# ----------------------------------------------------------------------------------------
def run_partis(spval, stype):
    ofn = ptnfn(spval, stype, mthd)
    if utils.output_exists(args, ofn, 'partis'):
        return
    for action in ['cache-parameters', 'partition']:
        cmd = './bin/partis %s --infname %s --parameter-dir %s' % (action, '%s/Fasta/%s_%s_simulated.fasta' % (mdir, spval, stype), paramdir(spval, stype))
        if action == 'partition':
            cmd += ' --outfname %s' % ofn
        utils.simplerun(cmd, logfname=utils.replace_suffix(ofn, '.log')) #, dryrun=True)

# ----------------------------------------------------------------------------------------
def make_plots():
    # ----------------------------------------------------------------------------------------
    def plot_metric_type(mtr_type):
        # ----------------------------------------------------------------------------------------
        def read_files():
            plotvals = {m : {t : collections.OrderedDict([[v, None] for v in spvals]) for t in stypes} for m in args.methods}
            for stype in stypes:
                for spval in spvals:
                    if not os.path.exists(mtrfn(spval, stype)):
                        continue
                    for mthd in args.methods:
                        with open(mtrfn(spval, stype)) as mfile:
                            vals = json.load(mfile)
                            plotvals[mthd][stype][spval] = vals[mtr_type]['f1']  # NOTE mtr_type is 'partis' or 'mobille', i.e. has the same values as args.methods
            return plotvals
        # ----------------------------------------------------------------------------------------
        def lzv(lvstr):
            assert lvstr[:3] == 'l00' and len(lvstr) == 5
            return float(int(lvstr[3:5])) / 100
        # ----------------------------------------------------------------------------------------
        plotvals = read_files()
        fig, ax = plotting.mpl_init()
        for stype in ['mono']: #stypes:
            for mthd in args.methods:
                xvals, yvals = zip(*plotvals[mthd][stype].items())
                xvals = [lzv(v) for v in xvals]
                ax.plot(xvals, yvals, label=mthd, alpha=0.6, linewidth=3, markersize=13, marker='.')
                ax.plot((xvals[0], xvals[-1]), (1, 1), linewidth=1.5, alpha=0.5, color='grey') #, linestyle='--') #, label='1/seq len')
                print plotting.mpl_finish(ax, pltdir(), 'f1-%s'%stype, ylabel=metric_labels.get(mtr_type, mtr_type), xlabel='lambda 0', xticks=xvals, ybounds=(0, 1.05)) #, xticklabels=['%.3f'%v for v in xvals]) #, log=log, xticks=xticks, xticklabels=xticks, leg_loc=(0.1, 0.2), xbounds=(xticks[0], xticks[-1]), ybounds=(ymin, 1.01), title=title, xlabel=xlabel, ylabel='metric value')
    # ----------------------------------------------------------------------------------------
    import plotting
    plot_metric_type('mobille-pairwise')

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--actions', default='run:write:plot')
parser.add_argument('--methods', default='partis')
parser.add_argument('--overwrite', action='store_true')
args = parser.parse_args()
args.actions = utils.get_arg_list(args.actions)
args.methods = utils.get_arg_list(args.methods)

spvals = ['l00%d'%c for c in [16, 26, 36, 46]]
stypes = ['mono', 'oligo', 'poly']
metric_labels = {'mobille-pairwise' : 'pairwise f1', 'mobille-closeness' : 'closeness f1'}

for stype in stypes:
    for spval in spvals:
        # if stype != 'mono' or spval != 'l0036':
        #     continue
        if 'run' in args.actions:
            run_partis(spval, stype)
        if 'write' in args.actions:
            write_metrics(spval, stype, 'partis')
if 'plot' in args.actions:
    make_plots()
