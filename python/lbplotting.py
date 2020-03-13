import copy
import pickle
from scipy import stats
import os
import sys
import numpy
import yaml
import json
import time
import warnings
import collections
import itertools
import math

import utils
from hist import Hist
import treeutils
import plotting

# ----------------------------------------------------------------------------------------
# this name is terrible, but it's complicated and I can't think of a better one
def lb_metric_axis_cfg(metric_method, final_plots=False):  # x axis variables against which we plot each lb metric (well, they're the x axis on the scatter plots, not the ptile plots)
    base_cfg = collections.OrderedDict([('lbi', [('affinity', 'affinity')]),
                                        ('lbr', [('n-ancestor', 'N ancestors')]),  # , ('branch-length', 'branch length')])  # turning off branch length at least for now (for run time reasons)
    ])
    if metric_method is None:
       return base_cfg.items()
    base_cfg['dtr'] = [('affinity', 'affinity'),] if final_plots else [('affinity', 'affinity'), ('n-ancestor', 'N ancestors')]  # hack hack hack
    if metric_method in base_cfg:
        return [(m, cfg) for m, cfg in base_cfg.items() if m == metric_method]
    else:  # shm, delta-lbi, cons-dist-*, etc
        xv = 'n-ancestor' if metric_method == 'delta-lbi' else 'affinity'  # also hack hack hack
        return [[metric_method, [(xv, xv.replace('n-a', 'N a'))]]]

# ----------------------------------------------------------------------------------------
def single_lbma_cfg_vars(metric_method, final_plots=False):
    assert metric_method is not None  # only use this for single metrics
    return lb_metric_axis_cfg(metric_method, final_plots=final_plots)[0][1]  # [0] gets you the first metric's cfg (there's ony one because we specified a single <metric_method>, [1] gets you the (var, label)

# ----------------------------------------------------------------------------------------
def add_use_relative_affy_stuff(cfg_list, include_relative_affy_plots=False):  # NOTE acts on sub lists of the return value of the above (i.e. the .values())
    cfg_list = [[s, l, False] for s, l in cfg_list]
    if include_relative_affy_plots:
        for ic, (s, l, u) in enumerate(cfg_list):
            if s == 'affinity':  # add it just after the existing (non-relative) 'affinity'
                cfg_list.insert(ic + 1, ['affinity', 'affinity', True])
                break
    return cfg_list
# ----------------------------------------------------------------------------------------
def rel_affy_str(var, use_relative_affy):
    if var == 'affinity' and use_relative_affy:
        return 'relative-'
    else:
        return ''

def meanmaxfcns(): return (('mean', lambda line, plotvals: numpy.mean(plotvals)), ('max', lambda line, plotvals: max(plotvals)))
def mean_of_top_quintile(vals):  # yeah, yeah could name it xtile and have another parameter, but maybe I won't actually need to change it
    frac = 0.2  # i.e. top quintile
    n_to_take = int(frac * len(vals))  # NOTE don't use numpy.percentile(), since affinity is fairly discrete-valued, which causes bad stuff (e.g. you don't take anywhere near the number of cells that you were trying to)
    return numpy.mean(sorted(vals)[len(vals) - n_to_take:])
mean_max_metrics = ['lbi', 'lbr', 'shm', 'lbi-cons-aa', 'lbi-cons-nuc'] + treeutils.dtr_metrics
cluster_summary_cfg = collections.OrderedDict()
for k in mean_max_metrics:
    cluster_summary_cfg[k] = meanmaxfcns()
cluster_summary_cfg['affinity'] = (('top-quintile', lambda line, plotvals: mean_of_top_quintile(plotvals)), )
cluster_summary_cfg['fay-wu-h'] = (('fay-wu-h', lambda line, plotvals: -utils.fay_wu_h(line)), )
cluster_summary_cfg['cons-dist-nuc'] = (('cons-seq-shm-nuc', lambda line, plotvals: treeutils.lb_cons_seq_shm(line, aa=False)), )  # NOTE the cluster_summary_cfg key doesn't really make sense here any more (used to be 'consensus'), but i'm not really using these cluster summary things any more, since it turns out to not work very well to choose entire families
cluster_summary_cfg['cons-dist-aa'] = (('cons-seq-shm-aa', lambda line, plotvals: treeutils.lb_cons_seq_shm(line, aa=True)), )  # NOTE the cluster_summary_cfg key doesn't really make sense here any more (used to be 'consensus'), but i'm not really using these cluster summary things any more, since it turns out to not work very well to choose entire families
cluster_summary_cfg['is_leaf'] = (('x-dummy-x', lambda line, plotvals: None), )  # just to keep things from breaking, doesn't actually get used
def get_lbscatteraxes(lb_metric):
    return ['affinity', lb_metric]
def getptvar(xvar): return xvar if xvar == 'affinity' else 'n-ancestor'  # NOTE if I start using 'branch-length' again, that'd also need to be included in the latter case
def ungetptvar(xvar): return xvar if xvar == 'affinity' else 'delta-affinity'  # ok this name sucks, and these two functions are anyway shitty hacks
def ungetptlabel(xvar): return xvar if xvar == 'affinity' else 'affinity change'

per_seq_metrics = ['lbi', 'lbr', 'shm', 'cons-dist-nuc', 'cons-dist-aa', 'delta-lbi', 'lbi-cons-aa', 'lbi-cons-nuc'] + treeutils.dtr_metrics
# per_clust_metrics = ('lbi', 'lbr', 'shm', 'fay-wu-h', 'cons-dist-nuc')  # don't need this atm since it's just all of them (note that 'cons-dist-nuc' doesn't really make sense here, see cluster_summary_cfg)
mtitle_cfg = {'per-seq' : {'cons-dist-nuc' : '- nuc distance to cons seq', 'cons-dist-aa' : '- AA distance to cons seq', 'shm' : '- N mutations', 'delta-lbi' : 'change in lb index', 'z-score-err' : 'z score diff (lb - affy)', 'edge-dist' : 'min root/tip dist',
                           'affinity-ptile' : 'affinity percentile', 'lbi-ptile' : 'lbi percentile', 'lbr-ptile' : 'lbr percentile', 'within-families-affinity-dtr' : 'in-family dtr', 'within-families-delta-affinity-dtr' : 'in-family dtr', 'among-families-affinity-dtr' : 'among-families dtr', 'among-families-delta-affinity-dtr' : 'among-families dtr'},
              'per-cluster' : {'fay-wu-h' : '- Fay-Wu H', 'cons-seq-shm-nuc' : 'N mutations in cons seq', 'shm' : '- N mutations', 'affinity' : 'top quintile affinity'}}
mtitle_shorts = {'cons-dist-aa' : 'AA cons dist', 'cons-dist-nuc' : 'nuc cons dist'}
def mtitlestr(pchoice, lbm, short=False, max_len=13):
    if pchoice == 'per-cluster' and 'cons-dist' in lbm:
        lbm = cluster_summary_cfg[lbm][0][0]  # hack hack hack
    mtstr = mtitle_cfg[pchoice].get(lbm, treeutils.lb_metrics.get(lbm, lbm))
    if short and len(mtstr) > max_len:
        tmpstr = mtitle_shorts.get(lbm, lbm)
        if len(tmpstr) < len(mtstr):
            mtstr = tmpstr
    return mtstr
# ----------------------------------------------------------------------------------------
metric_for_target_distance_labels = {
    'aa' : 'AA hamming',
    'nuc' : 'nuc',
    'aa-sim-ascii' : 'AA ascii',
    'aa-sim-blosum' : 'AA BLOSUM',
}
cdist_keys = ['cons-dist-'+s for s in ['nuc', 'aa']]

# ----------------------------------------------------------------------------------------
def add_fn(fnames, fn=None, init=False, n_per_row=4):
    if fnames is None:
        if init:
            fnames = []
        else:
            return
    if len(fnames) == 0 or len(fnames[-1]) >= n_per_row:
        fnames.append([])
    if fn is not None:
        fnames[-1].append(fn)
    return fnames

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_selection_hists(histfname, plotdir, plotname, plot_all=False, n_plots=7, title='', xlabel=''):
    import joypy
    # ----------------------------------------------------------------------------------------
    def plot_this_time(otime, numpyhists):
        if plot_all:
            return True
        if otime == 0:
            return False
        if otime in (len(numpyhists),):
            return True
        if otime % max(1, int(len(numpyhists) / float(n_plots))) == 0:
            return True
        return False
    # ----------------------------------------------------------------------------------------
    def get_hists(hfname):
        with open(hfname) as runstatfile:
            numpyhists = pickle.load(runstatfile)
        xmin, xmax = None, None
        hists, ylabels, xtralabels = [], [], []
        for ihist in range(len(numpyhists)):
            nphist = numpyhists[ihist]  # numpy.hist is two arrays: [0] is bin counts, [1] is bin x values (not sure if low, high, or centers)
            obs_time = ihist  #  + 1  # I *think* it's right without the 1 (although I guess it's really a little arbitrary)
            if not plot_this_time(obs_time, numpyhists):
                continue
            if nphist is None:  # time points at which we didn't sample
                hists.append(None)
                ylabels.append(None)
                xtralabels.append(None)
                continue
            bin_contents, bin_edges = nphist
            assert len(bin_contents) == len(bin_edges) - 1
            # print ' '.join('%5.1f' % c for c in bin_contents)
            # print ' '.join('%5.1f' % c for c in bin_edges)
            hist = Hist(len(bin_edges) - 1, bin_edges[0], bin_edges[-1])
            for ibin in range(len(bin_edges) - 1):  # nphist indexing, not Hist indexing
                lo_edge = bin_edges[ibin]
                hi_edge = bin_edges[ibin + 1]
                xmin = lo_edge if xmin is None else min(xmin, lo_edge)
                xmax = hi_edge if xmax is None else max(xmax, hi_edge)
                bin_center = (hi_edge + lo_edge) / 2.
                for _ in range(bin_contents[ibin]):
                    hist.fill(bin_center)
            hists.append(hist)
            ylabels.append('%d' % obs_time)
            xtralabels.append('(%.1f, %.0f)' % (hist.get_mean(), hist.integral(include_overflows=True)))

        hists, ylabels, xtralabels = zip(*[(h, yl, xl) for h, yl, xl in zip(hists, ylabels, xtralabels) if h is not None])  # remove None time points
        return hists, ylabels, xtralabels, xmin, xmax

    # ----------------------------------------------------------------------------------------
    all_hists, all_ylabels, all_xtralabels, xmin, xmax = get_hists(histfname)  # these xmin, xmax are the actual (ORd) bounds of the histograms (whereas below we also get the ranges that around filled)
    if sum(h.integral(include_overflows=True) for h in all_hists) == 0:
        print '  %s no/empty hists in %s' % (utils.color('yellow', 'warning'), histfname)
        return
    jpdata = []
    for hist in all_hists:
        jpdata.append([x for x, y in zip(hist.get_bin_centers(), hist.bin_contents) for _ in range(int(y)) if x > xmin and x < xmax])  # NOTE this is repeating the 'for _ in range()' in the fcn above, but that's because I used to be actually using the Hist()s, and maybe I will again

    fbin_xmins, fbin_xmaxs = zip(*[h.get_filled_bin_xbounds(extra_pads=2) for h in all_hists])
    xmin_filled, xmax_filled = min(fbin_xmins), max(fbin_xmaxs)

    pre_fig, pre_ax = plotting.mpl_init()  # not sure to what extent these really get used after joypy is done with things
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # i don't know why it has to warn me that it's clearing the fig/ax I'm passing in, and I don't know how else to stop it
        fig, axes = joypy.joyplot(jpdata, labels=all_ylabels, fade=True, hist=True, overlap=0.5, ax=pre_ax, x_range=(xmin_filled, xmax_filled), bins=int(xmax_filled - xmin_filled), xlabelsize=15) #, ylabelsize=15)
    # xtextpos = 0.85 * (xmax_filled - xmin_filled) + xmin_filled  # this is from before I found transform=ax.transAxes, but I don't want to remove it yet
    fsize = 15
    for ax, lab in zip(axes, all_xtralabels):
        ax.text(0.85, 0.2, lab, fontsize=fsize, transform=ax.transAxes)
    fig.text(0.03, 0.9, 'generation', fontsize=fsize)
    fig.text(0.8, 0.87, '(mean, N cells)', fontsize=fsize)
    # NOTE do *not* set your own x ticks/labels in the next line, since they'll be in the wrong place (i.e. not the same as where joypy puts them) (also note, the stupid y labels don't work, but setting them on the joyplot axes also doesn't work)
    plotting.mpl_finish(pre_ax, plotdir, plotname, title=title, xlabel=xlabel) #, ylabel='generation') #, leg_loc=(0.7, 0.45)) #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_kd_vals(plotdir, event):
    kd_changes = []
    dtree = treeutils.get_dendro_tree(treestr=event['tree'])
    for node in dtree.preorder_internal_node_iter():
        if node.taxon.label not in event['unique_ids']:
            continue
        inode = event['unique_ids'].index(node.taxon.label)
        node_affinity = event['affinities'][inode]
        for child in node.child_nodes():
            if child.taxon.label not in event['unique_ids']:
                continue
            ichild = event['unique_ids'].index(child.taxon.label)
            child_affinity = event['affinities'][ichild]
            kd_changes.append(1./child_affinity - 1./node_affinity)

    if len(kd_changes) > 0:
        hist = Hist(30, min(kd_changes), max(kd_changes))
        for val in kd_changes:
            hist.fill(val)
        fig, ax = plotting.mpl_init()
        hist.mpl_plot(ax, square_bins=True, errors=False)  #remove_empty_bins=True)
        plotname = 'kd-changes'
        plotting.mpl_finish(ax, plotdir,  plotname, xlabel='parent-child kd change', ylabel='branches', log='y') #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

    plotvals = {'shm' : [], 'kd_vals' : []}
    for iseq, uid in enumerate(event['unique_ids']):
        plotvals['shm'].append(event['n_mutations'][iseq])
        plotvals['kd_vals'].append(1. / event['affinities'][iseq])
    # new_cmap = plotting.truncate_colormap(plt.cm.Blues, 0, 1)
    # ax.hexbin(kd_changes, shms, gridsize=25, cmap=plt.cm.Blues) #, info['ccf_under'][meth], label='clonal fraction', color='#cc0000', linewidth=4)
    fig, ax = plotting.mpl_init()
    ax.scatter(plotvals['kd_vals'], plotvals['shm'], alpha=0.4)
    plotname = 'kd-vs-shm'
    plotting.mpl_finish(ax, plotdir, plotname, xlabel='Kd', ylabel='N mutations') #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_target_attraction(plotdir, event):  # plots of which sequences are going toward which targets
    fig, ax = plotting.mpl_init()

    # affinity vs stuff:
    # xvals = [1. / af for line in mutated_events for af in line['affinities']]
    # yvals = [nm for line in mutated_events for nm in line['n_mutations']]

    # # min distance to target:
    # yvals = [hd for line in mutated_events for hd in get_min_target_hdists(line['input_seqs'], line['target_seqs'])]
    # ax.scatter(xvals, yvals, alpha=0.65)

    hist = Hist(len(event['target_seqs']), -0.5, len(event['target_seqs']) - 0.5, value_list=event['nearest_target_indices'])
    hist.mpl_plot(ax, alpha=0.7, ignore_overflows=True)

    plotname = 'nearest-target-identities'
    plotting.mpl_finish(ax, plotdir, plotname, xlabel='index (identity) of nearest target sequence', ylabel='counts') #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_simulation(outdir, event, extrastr, metric_for_target_distance_label):
    utils.prep_dir(outdir + '/plots', wildlings=['*.csv', '*.svg'])

    plot_bcr_phylo_kd_vals(outdir + '/plots', event)
    plot_bcr_phylo_target_attraction(outdir + '/plots', event)

    plot_bcr_phylo_selection_hists('%s/%s_min_aa_target_hdists.p' % (outdir, extrastr), outdir + '/plots', 'min-aa-target-all-cells', title='all cells', xlabel='%s distance to nearest target seq' % metric_for_target_distance_label)
    plot_bcr_phylo_selection_hists('%s/%s_sampled_min_aa_target_hdists.p' % (outdir, extrastr), outdir + '/plots', 'min-aa-target-sampled-cells', plot_all=True, title='sampled cells (excluding ancestor sampling)', xlabel='%s distance to nearest target seq' % metric_for_target_distance_label)
    plot_bcr_phylo_selection_hists('%s/%s_n_mutated_nuc_hdists.p' % (outdir, extrastr), outdir + '/plots', 'n-mutated-nuc-all-cells', title='SHM all cells', xlabel='N nucleotide mutations to naive')

    plotting.make_html(outdir + '/plots')

# ----------------------------------------------------------------------------------------
def get_tree_from_line(line, is_true_line):
    if is_true_line:
        return line['tree']
    if 'tree-info' not in line:  # if 'tree-info' is missing, it should be because it's a small cluster in data that we skipped when calculating lb values
        return None
    return line['tree-info']['lb']['tree']

# ----------------------------------------------------------------------------------------
def make_lb_scatter_plots(xvar, baseplotdir, lb_metric, lines_to_use, fnames=None, is_true_line=False, colorvar=None, only_overall=False, only_iclust=False, add_uids=False, yvar=None, choose_among_families=False,
                          add_jitter=False, min_ptile=80., iclust_fnames=None, use_relative_affy=False):  # <is_true_line> is there because we want the true and inferred lines to keep their trees in different places, because the true line just has the one, true, tree, while the inferred line could have a number of them (yes, this means I maybe should have called it the 'true-tree' or something)
    cdist_pt_keys = [s+'-ptile' for s in cdist_keys]
    if yvar is None:
        yvar = lb_metric
    add_fn(fnames, init=True)

    choice_str = 'among-families' if choose_among_families else 'within-families'
    xlabel = mtitlestr('per-seq', xvar).replace('- N', 'N')
    ylabel = mtitlestr('per-seq', yvar)
    if choose_among_families:
        assert '-ptile' in xvar or '-ptile' in yvar
        if xvar == 'affinity-ptile':
            all_affinities = sorted([a for l in lines_to_use for a in l['affinities']])  # I don't need to sort, but it seems like the stats.percentileofscore() calls might be faster if I do
            affinity_ptiles = {u : stats.percentileofscore(all_affinities, l['affinities'][i], kind='weak') for l in lines_to_use for i, u in enumerate(l['unique_ids'])}
        if yvar == '%s-ptile'%lb_metric:
            all_lbvals = sorted([l['tree-info']['lb'][lb_metric][u] for l in lines_to_use for u in l['unique_ids']])  # have to remove the ones that aren't in <line> (and sort isn't necessary, but is maybe faster)
            lb_ptiles = {u : stats.percentileofscore(all_lbvals, l['tree-info']['lb'][lb_metric][u], kind='weak') for l in lines_to_use for u in l['unique_ids']}
        elif yvar in cdist_pt_keys:
            assert False  # needs finishing
            # all_cdistvals = []
            # for ltmp in lines_to_use:
            #     cseq = XXX needs updating XXX utils.cons_seq_of_line(ltmp, aa=XXX)
            #     for utmp, stmp in zip(ltmp['unique_ids'], ltmp['seqs']):
            #         all_cdistvals.append(treeutils.XXX update (use add_cons_dists(), or just access the keys in <line>) XXX-lb_cons_dist(cseq, stmp, aa=XXX))
            # all_cdistvals = sorted(all_cdistvals)
            # XXX lb_ptiles = {u : stats.percentileofscore(all_cdistvals, l['tree-info']['lb'][lb_metric][u], kind='weak') for l in lines_to_use for u in l['unique_ids']}
    if '-ptile' in xvar:
        xlabel = '%s %s' % (choice_str, xlabel)
    if '-ptile' in yvar:
        ylabel = '%s %s' % (choice_str, ylabel)

    vtypes = [xvar, yvar]
    if add_uids: vtypes.append('uids')
    if colorvar is not None: vtypes.append(colorvar)
    plotdir = '%s/%s/%s%s-vs-%s%s' % (baseplotdir, lb_metric, rel_affy_str(yvar, use_relative_affy), yvar, rel_affy_str(xvar, use_relative_affy), xvar)
    if '-ptile' in xvar or '-ptile' in yvar:
        plotdir += '-%s' % choice_str
    utils.prep_dir(plotdir, wildlings='*.svg')
    plotvals = {x : [] for x in vtypes}
    basetitle = '%s %s vs %s' % ('true' if is_true_line else 'inferred', mtitlestr('per-seq', yvar, short=True), mtitlestr('per-seq', xvar, short=True).replace('- N', 'N'))  # here 'shm' the plain number of mutations, not 'shm' the non-lb metric, so we have to fiddle with the label in mtitle_cfg
    scatter_kwargs = {'xvar' : xvar, 'xlabel' : xlabel, 'colorvar' : colorvar, 'leg_loc' : (0.55, 0.75), 'log' : 'y' if yvar == 'lbr' else ''}
    if use_relative_affy:
        scatter_kwargs['warn_text'] = 'relative affinity'
    sorted_lines = sorted(lines_to_use, key=lambda l: len(l['unique_ids']), reverse=True)
    for iclust, line in enumerate(sorted_lines):  # get depth/n_mutations for each node
        iclust_plotvals = {x : [] for x in vtypes}
        lbfo = line['tree-info']['lb'][lb_metric]
        if len(set(lbfo) - set(line['unique_ids'])) > 0 or len(set(line['unique_ids']) - set(lbfo)) > 0:
            if 'ptile' in xvar or 'ptile' in yvar:  # if it's not a ptile, we just don't plot it unless we have a value for both x and y. But for ptiles I don't really want to think through right now how to handle it
                raise Exception('  %s different uids in <lbfo> and <line>, and xvar or yvar is a ptile: %3d extra uids in lbfo, %3d extra in line, %3d in common. This (the code below) needs to be checked.' % (utils.color('yellow', 'warning'), len(set(lbfo) - set(line['unique_ids'])), len(set(line['unique_ids']) - set(lbfo)), len(set(line['unique_ids']) & set(lbfo))))
        if colorvar in ['is_leaf', 'edge-dist'] or xvar == 'edge-dist':
            dtree = treeutils.get_dendro_tree(treestr=get_tree_from_line(line, is_true_line))
            if dtree is None:
                continue

        # TODO I should really combine the xvar and yvar stuff here
        if xvar == 'shm':
            def xvalfcn(i): return line['n_mutations'][i]
        elif xvar == 'affinity':
            def xvalfcn(i): return line['affinities'][i]
        elif xvar in cdist_keys:
            treeutils.add_cons_dists(line, aa='aa-' in xvar)
            tkey = xvar.replace('cons-dist-', 'cons_dists_')
            def xvalfcn(i): return -line[tkey][i]
        elif xvar == 'edge-dist':
            def xvalfcn(i): return treeutils.edge_dist_fcn(dtree, line['unique_ids'][i])
        elif xvar == 'affinity-ptile':
            if not choose_among_families:
                affinity_ptiles = {u : stats.percentileofscore(line['affinities'], line['affinities'][i], kind='weak') for i, u in enumerate(line['unique_ids'])}
            def xvalfcn(i): return affinity_ptiles[line['unique_ids'][i]]
        else:
            assert False

        if yvar == lb_metric:
            def yvalfcn(i): return lbfo.get(line['unique_ids'][i], None)
        elif yvar == 'z-score-err':  # prediction error, i.e. difference between lb z-score and affinity z-score (this isn't really good for much: see note below)
            # z score difference: (this isn't actually what we want, since 1) we don't care if the best lbi was +3 sigma while best affy was +7 sigma, but only if the ranking is right and 2) we don't care about all the cells that are below the median affy)
            zscores = {lb_metric : utils.get_z_scores([lbfo[u] for u in line['unique_ids'] if u in lbfo]),
                       'affinity' : utils.get_z_scores(line['affinities'])}
            def yvalfcn(i): return zscores[lb_metric][i] - zscores['affinity'][i]
        elif yvar == '%s-ptile'%lb_metric:
            lbvals = [lbfo[u] for u in line['unique_ids'] if u in lbfo]  # have to remove the ones that aren't in <line> (and vice versa)
            if not choose_among_families:
                lb_ptiles = {u : stats.percentileofscore(lbvals, lbfo[u], kind='weak') for u in line['unique_ids'] if u in lbfo}
            def yvalfcn(i): return lb_ptiles[line['unique_ids'][i]]
        elif yvar in cdist_pt_keys:
            treeutils.add_cons_dists(line, aa='aa-' in yvar)
            tkey = yvar.replace('-ptile', '').replace('cons-dist-', 'cons_dists_')
            cvals = [-v for v in line[tkey]]
            if not choose_among_families:
                cdist_ptiles = {u : stats.percentileofscore(cvals, cvals[i], kind='weak') for i, u in enumerate(line['unique_ids'])}
            def yvalfcn(i): return cdist_ptiles[line['unique_ids'][i]]
        else:
            assert False

        for iseq, uid in enumerate(line['unique_ids']):
            xval = xvalfcn(iseq)
            if xvar == 'affinity-ptile' and '-ptile' in yvar and xvalfcn(iseq) < min_ptile:  # and yvalfcn(iseq) < min_ptile:  the number of cells with high lbi but low affinity (last, commented criterion) is just too small to bother plotting -- all our errors come from the other direction
                continue
            yval = yvalfcn(iseq)
            if yvar == 'lbr' and yval == 0:
                continue
            if None in (xval, yval):
                continue
            iclust_plotvals[xvar].append(xval)
            iclust_plotvals[yvar].append(yval)
            if colorvar is not None:
                if colorvar == 'is_leaf':
                    node = dtree.find_node_with_taxon_label(uid)
                    colorval = node.is_leaf() if node is not None else None
                elif colorvar == 'affinity':
                    colorval = line['affinities'][iseq] if 'affinities' in line else None
                elif colorvar == 'edge-dist':
                    colorval = treeutils.edge_dist_fcn(dtree, line['unique_ids'][iseq])
                else:
                    assert False
                iclust_plotvals[colorvar].append(colorval)  # I think any uid in <line> should be in the tree, but may as well handle the case where it isn't
            if add_uids:
                iclust_plotvals['uids'].append(uid)  # use to add None here instead of <uid> if this node didn't have an affinity value, but that seems unnecessary, I can worry about uid config options later when I actually use the uid dots for something
        if add_jitter: #xvar == 'affinity-ptile' and '-ptile' in yvar:
            def jitter(frac=0.02):
                # delta = max(3, max(iclust_plotvals[xvar]) - min(iclust_plotvals[xvar]))
                delta = 100. - min_ptile if '-ptile' in xvar else max(iclust_plotvals[xvar]) - min(iclust_plotvals[xvar])
                return numpy.random.uniform(-frac * delta, frac * delta)
            iclust_plotvals[xvar] = [x + jitter() for x in iclust_plotvals[xvar]]
            if 'jitter' not in scatter_kwargs['xlabel']:
                scatter_kwargs['xlabel'] += ' (+jitter)'
        if not only_overall:
            title = '%s (%d observed, %d total)' % (basetitle, len(line['unique_ids']), len(lbfo))
            tmpplotname = '%s%s-vs-%s%s-iclust-%d' % (rel_affy_str(yvar, use_relative_affy), yvar, rel_affy_str(yvar, use_relative_affy), xvar, iclust)
            fn = plot_2d_scatter(tmpplotname, plotdir, iclust_plotvals, yvar, ylabel, title, **scatter_kwargs)
            if iclust_fnames is not None and iclust < iclust_fnames:
                fnames[-1].append(fn)
        assert len(set([len(plotvals[vt]) for vt in plotvals])) == 1  # make sure all of them are the same length
        for vtype in [vt for vt in plotvals if vt != 'uids']:
            plotvals[vtype] += iclust_plotvals[vtype]

    # uncomment this to only plot high affinity/lb values (mostly to speed up plotting)
    # if xvar == 'affinity' and len(plotvals[xvar]) > 250:
    #     print 'yep'
    #     min_ptile = 75
    #     min_xval, min_yval = [numpy.percentile(plotvals[v], min_ptile) for v in [xvar, yvar]]
    #     new_plotvals = {x : [] for x in vtypes}
    #     for ival in range(len(plotvals[xvar])):
    #         if plotvals[xvar][ival] < min_xval and plotvals[yvar][ival] < min_yval:
    #             continue
    #         for vt in plotvals:
    #             new_plotvals[vt].append(plotvals[vt][ival])
    #     plotvals = new_plotvals

    if not only_iclust:
        tmpplotname = '%s%s-vs-%s%s-all-clusters' % (rel_affy_str(yvar, use_relative_affy), yvar, rel_affy_str(xvar, use_relative_affy), xvar)
        fn = plot_2d_scatter(tmpplotname, plotdir, plotvals, yvar, ylabel, '%s (all clusters)' % basetitle, **scatter_kwargs)
        add_fn(fnames, fn=fn)
    if iclust_fnames is not None:
        fnames.append([])

# ----------------------------------------------------------------------------------------
def plot_lb_distributions(lb_metric, baseplotdir, lines_to_use, is_true_line=False, fnames=None, only_overall=False, affy_key='affinities', iclust_fnames=None):
    def make_hist(plotvals, n_total, n_skipped, iclust=None, affinities=None):
        if len(plotvals) == 0:
            return
        hist = Hist(30, min(plotvals), max(plotvals), value_list=plotvals)
        fig, ax = plotting.mpl_init()
        hist.mpl_plot(ax) #, square_bins=True, errors=False)
        # fig.text(0.7, 0.8, 'mean %.3f' % numpy.mean(plotvals), fontsize=15)
        # fig.text(0.7, 0.75, 'max %.3f' % max(plotvals), fontsize=15)
        # if affinities is not None:
        #     fig.text(0.38, 0.88, 'mean/max affinity: %.4f/%.4f' % (numpy.mean(affinities), max(affinities)), fontsize=15)
        plotname = '%s-%s' % (lb_metric, str(iclust) if iclust is not None else 'all-clusters')
        leafskipstr = ', skipped %d leaves' % n_skipped if n_skipped > 0 else ''  # ok they're not necessarily leaves, but almost all of them are leaves (and not really sure how a non-leaf could get zero, but some of them seem to)
        title = '%s %s: %d entries%s (%s)' % ('true' if is_true_line else 'inferred', mtitlestr('per-seq', lb_metric, short=True), n_total, leafskipstr, 'all clusters' if iclust is None else 'iclust %d'%iclust)
        fn = plotting.mpl_finish(ax, plotdir, plotname, xlabel=mtitlestr('per-seq', lb_metric), log='y', ylabel='counts', title=title)
        if iclust is None or (iclust_fnames is not None and iclust < iclust_fnames):
            add_fn(fnames, fn=fn)

    sorted_lines = sorted([l for l in lines_to_use if 'tree-info' in l], key=lambda l: len(l['unique_ids']), reverse=True)  # if 'tree-info' is missing, it should be because it's a small cluster we skipped when calculating lb values

    plotvals = []
    n_total_skipped_leaves = 0
    plotdir = '%s/%s/distributions' % (baseplotdir, lb_metric)
    utils.prep_dir(plotdir, wildlings=['*.svg'])
    for iclust, line in enumerate(sorted_lines):
        lbfo = line['tree-info']['lb'][lb_metric]  # NOTE this contains any ancestor nodes that the phylogenetic program has in the tree, but that aren't in <line['unique_ids']>
        if is_true_line:
            iclust_plotvals = [lbfo[u] for u in line['unique_ids'] if u in lbfo]  # for the true plots, we *don't* want to include any inferred ancestor nodes that weren't actually sampled, since we don't have affinity info for them, and it'd make it look like there's a mismatch between these distributions and the lb vs affinity plots (which won't have them)
        else:
            iclust_plotvals = lbfo.values()  # whereas for real data, we want to include the inferred ancestor nodes for which we don't have sequences (although I guess in most cases where we're really interested in them, we would've used a phylogenetics program that also inferred their sequences, so they'd presumably have been added to <line['unique_ids']>)
        cluster_size = len(iclust_plotvals)  # i.e. including leaves
        if lb_metric == 'lbr':
            iclust_plotvals = [v for v in iclust_plotvals if v > 0.]  # don't plot the leaf values, they just make the plot unreadable
        if not only_overall:
            affinities = line[affy_key] if affy_key in line else None
            make_hist(iclust_plotvals, cluster_size, cluster_size - len(iclust_plotvals), iclust=iclust, affinities=affinities)
        plotvals += iclust_plotvals
        n_total_skipped_leaves += cluster_size - len(iclust_plotvals)
    make_hist(plotvals, len(plotvals) + n_total_skipped_leaves, n_total_skipped_leaves)
    if iclust_fnames is not None:
        fnames.append([])

# ----------------------------------------------------------------------------------------
def make_lb_affinity_joyplots(plotdir, lines, lb_metric, fnames=None, n_clusters_per_joy_plot=25, n_max_joy_plots=25, n_plots_per_row=4):
    partition = utils.get_partition_from_annotation_list(lines)
    annotation_dict = {':'.join(l['unique_ids']) : l for l in lines}
    sorted_clusters = sorted(partition, key=lambda c: mean_of_top_quintile(annotation_dict[':'.join(c)]['affinities']), reverse=True)
    sorted_cluster_groups = [sorted_clusters[i : i + n_clusters_per_joy_plot] for i in range(0, len(sorted_clusters), n_clusters_per_joy_plot)]
    repertoire_size = sum([len(c) for c in sorted_clusters])
    max_affinity = max([a for c in sorted_clusters for a in annotation_dict[':'.join(c)]['affinities']])  # it's nice to keep track of the max values over the whole repertoire so all plots can have the same max values
    max_lb_val = max([annotation_dict[':'.join(c)]['tree-info']['lb'][lb_metric][u] for c in sorted_clusters for u in c])  # NOTE don't use all the values in the dict in 'tree-info', since non-sampled sequences (i.e. usually intermediate ancestors) are in there
    if max_lb_val == 0.:  # at least atm, this means this is lbr on a family with no common ancestor sampling
        return
    # print 'divided repertoire of size %d with %d clusters into %d cluster groups' % (repertoire_size, len(sorted_clusters), len(sorted_cluster_groups))
    iclustergroup = 0
    for subclusters in sorted_cluster_groups:
        if iclustergroup > n_max_joy_plots:
            continue
        title = 'affinity and %s (%d / %d)' % (mtitlestr('per-seq', lb_metric), iclustergroup + 1, len(sorted_cluster_groups))  # NOTE it's important that this denominator is still right even when we don't make plots for all the clusters (which it is, now)
        fn = plotting.make_single_joyplot(subclusters, annotation_dict, repertoire_size, plotdir, '%s-affinity-joyplot-%d' % (lb_metric, iclustergroup), x1key='affinities', x1label='affinity', x2key=lb_metric, x2label=mtitlestr('per-seq', lb_metric),
                                          global_max_vals={'affinities' : max_affinity, lb_metric : max_lb_val}, title=title)  # note that we can't really add cluster_indices> like we do in partitionplotter.py, since (i think?) the only place there's per-cluster plots we'd want to correspond to is in the bcr phylo simulation dir, which has indices unrelated to anything we're sorting by here, and which we can't reconstruct
        add_fn(fnames, fn=fn)
        iclustergroup += 1

# ----------------------------------------------------------------------------------------
def plot_2d_scatter(plotname, plotdir, plotvals, yvar, ylabel, title, xvar='affinity', xlabel='affinity', colorvar=None, log='', leg_loc=None, warn_text=None, markersize=15, stats=None):
    leafcolors = {'leaf' : '#006600', 'internal' : '#2b65ec'}  # green, blue
    if len(plotvals[xvar]) == 0:
        # print '    no %s vs affy info' % yvar
        return '%s/%s.svg' % (plotdir, plotname)
    fig, ax = plotting.mpl_init()
    alpha = 0.4
    if colorvar is None:
        ax.scatter(plotvals[xvar], plotvals[yvar], alpha=0.4)
    else:
        if colorvar == 'is_leaf':
            colorfcn = lambda x: leafcolors['leaf' if x else 'internal']
        else:
            smap = plotting.get_normalized_scalar_map([v for v in plotvals[colorvar] if v is not None], 'viridis')
            colorfcn = lambda x: 'grey' if x is None else plotting.get_smap_color(smap, None, val=x)
            alpha = 0.65
        for x, y, cval in zip(plotvals[xvar], plotvals[yvar], plotvals[colorvar]):  # we used to do the leaf/internal plots as two scatter() calls, which might be faster? but I think what really takes the time is writing the svgs, so whatever
            ax.plot([x], [y], color=colorfcn(cval), marker='.', markersize=markersize, alpha=alpha)
    if 'uids' in plotvals:
        for xval, yval, uid in zip(plotvals[xvar], plotvals[yvar], plotvals['uids']):  # note: two ways to signal not to do this: sometimes we have 'uids' in the dict, but don't fill it (so the zip() gives an empty list), but sometimes we populate 'uids' with None values
            if uid is None:
                continue
            ax.plot([xval], [yval], color='red', marker='.', markersize=markersize)
            ax.text(xval, yval, uid, color='red', fontsize=8)

    if warn_text is not None:
        ax.text(0.6 * ax.get_xlim()[1], 0.75 * ax.get_ylim()[1], warn_text, fontsize=30, fontweight='bold', color='red')
    xmin, xmax = [mfcn(plotvals[xvar]) for mfcn in [min, max]]
    ymin, ymax = [mfcn(plotvals[yvar]) for mfcn in [min, max]]
    xbounds = xmin - 0.02 * (xmax - xmin), xmax + 0.02 * (xmax - xmin)
    if 'y' in log:
        ybounds = 0.75 * ymin, 1.3 * ymax
    else:
        ybounds = ymin - 0.03 * (ymax - ymin), ymax + 0.08 * (ymax - ymin)
    if yvar in ['shm']+cdist_keys:
        ax.plot([xmin, xmax], [0, 0], linewidth=1, alpha=0.7, color='grey')
    leg_title, leg_prop, leg_iter = None, None, []
    if colorvar is not None:
        leg_loc = [0.1 if xvar in cdist_keys+['affinity'] else 0.7, 0.65]  # I think this is sometimes overriding the one that's passed in
        if yvar == 'cons-dist-aa':
            leg_loc = [0.75, 0.15]
        leg_prop = {'size' : 12}
        if colorvar == 'is_leaf':
            leg_iter = [(leafcolors[l], l) for l in ['leaf', 'internal']]
        elif plotvals[colorvar].count(None) == len(plotvals[colorvar]):  # no points have color values
            pass
        elif colorvar in ['affinity', 'edge-dist']:
            leg_title = mtitlestr('per-seq', colorvar)
            cmin, cmax = [mfcn(v for v in plotvals[colorvar] if v is not None) for mfcn in [min, max]]  # NOTE very similar to add_legend() in bin/plot-lb-tree.py
            if cmin != cmax:
                n_entries = 4
                max_diff = (cmax - cmin + utils.eps) / float(n_entries - 1)
                leg_iter = [(colorfcn(v), '%.3f'%v) for v in list(numpy.arange(cmin, cmax + 2*utils.eps, max_diff))]  # first value is exactly <cmin>, last value is exactly <cmax> (eps is to keep it from missing the last one)
        else:
            assert False
        for tcol, tstr in leg_iter:
            ax.plot([], [], color=tcol, label=tstr, marker='.', markersize=markersize, linewidth=0, alpha=alpha)

    if stats is not None:
        if stats == 'correlation':
            fig.text(0.7, 0.3, 'r = %.3f' % numpy.corrcoef(plotvals[xvar], plotvals[yvar])[0, 1], fontsize=20, fontweight='bold') #, color='red')
    fn = plotting.mpl_finish(ax, plotdir, plotname, title=title, xlabel=xlabel, ylabel=ylabel, xbounds=xbounds, ybounds=ybounds, log=log, leg_loc=leg_loc, leg_title=leg_title, leg_prop=leg_prop)
    return fn

# ----------------------------------------------------------------------------------------
def get_ptile_vals(lb_metric, plotvals, xvar, xlabel, ptile_range_tuple=(50., 100., 2.), dbgstr=None, affy_key='affinities', debug=False):
    # NOTE xvar and xlabel refer to the x axis on the scatter plot from which we make this ptile plot (i.e. are affinity, N ancestors, or branch length). On this ptile plot it's the y axis. (I tried calling it something else, but it was more confusing)
    xia = xvar == 'affinity'
    xkey = 'mean_%s_ptiles' % xvar
    tmp_ptvals = {'lb_ptiles' : [], xkey : [], 'perfect_vals' : []}  # , 'reshuffled_vals' : []}
    if len(plotvals[xvar]) == 0:
        return tmp_ptvals
    if debug:
        print '    getting ptile vals%s' % ('' if dbgstr is None else (' for %s' % utils.color('blue', dbgstr)))
        print '            %-12s N      mean     %s   |  perfect   perfect' % (lb_metric, 'mean   ' if xia else '')
        print '    ptile   threshold  taken    %s%-s %s|  N taken  mean %s'  % (utils.color('red', 'rel-') if xia and 'relative' in affy_key else '', 'affy' if xia else xlabel, '   affy ptile ' if xia else '', 'ptile' if xia else xlabel)
    max_lb_val = max(plotvals[lb_metric])
    sorted_xvals = sorted(plotvals[xvar], reverse=xia)
    if xia:
        np_arr_sorted_xvals = numpy.asarray(sorted_xvals)  # the stats calls in the next two lines make this conversion, and it's way too slow to do for every x
        corr_ptile_vals = [stats.percentileofscore(np_arr_sorted_xvals, x, kind='weak') for x in plotvals[xvar]]  # x (e.g. affinity) percentile of each plot val (could speed this up by only using the best half of each list (since ptiles are starting at 50))
        perf_ptile_vals = sorted(corr_ptile_vals, reverse=True)  # x (e.g. affinity) percentile or each plot val, *after sorting by x* (e.g. affinity)
    for percentile in numpy.arange(*ptile_range_tuple):
        lb_ptile_val = min(max_lb_val, numpy.percentile(plotvals[lb_metric], percentile))  # lb value corresponding to <percentile> (the min() is because the numpy call has precision issues that sometimes give you a value (very very slightly) larger than any of the actual values in the list)
        final_xvar_vals = [pt for lb, pt in zip(plotvals[lb_metric], corr_ptile_vals if xia else plotvals[xvar]) if lb >= lb_ptile_val]  # percentiles (if xia, else plain xvals [i.e. N ancestors or branch length]) corresponding to lb greater than <lb_ptile_val> (i.e. the ptiles for the x vals that you'd get if you took all the lb values greater than that)
        tmp_ptvals['lb_ptiles'].append(float(percentile))  # stupid numpy-specific float classes (I only care because I write it to a yaml file below)
        assert len(final_xvar_vals) > 0  # this shouldn't need to be here, but I'm removing the old block that handled the length-zero case (because it had a bad bug), and i want to be absolutely certain it doesn't come up any more. (it was necessary because the above line [and below in dbg] were '>' rather than '>=')
        tmp_ptvals[xkey].append(float(numpy.mean(final_xvar_vals)))

        # make a "perfect" line using the actual x values, as opposed to just a straight line (this accounts better for, e.g. the case where the top N affinities are all the same)
        n_to_take = len(final_xvar_vals)  # this used to be (in general) different than the number we took above, hence the weirdness/duplication (could probably clean up at this point)
        tmp_ptvals['perfect_vals'].append(float(numpy.mean((perf_ptile_vals if xia else sorted_xvals)[:n_to_take])))

        if debug:
            if xia:  # now we have to get these if we want to print them, since we no longer calculate them otherwise
                corr_xvals = [xv for lb, xv in zip(plotvals[lb_metric], plotvals[xvar]) if lb >= lb_ptile_val]  # x vals corresponding to lb greater than <lb_ptile_val> (i.e. the x vals that you'd get if you took all the lb values greater than that)
            v1str = ('%8.4f' % numpy.mean(corr_xvals if xia else final_xvar_vals)) if xia else ''
            f1str = '5.0f' if xia else '6.2f'
            f2str = '5.0f' if xia else ('8.2f' if xvar == 'n-ancestor' else '8.6f')
            print ('   %5.0f   %6.2f     %4d   %s  %'+f1str+'       |  %4d      %-'+f2str) % (percentile, lb_ptile_val, len(final_xvar_vals), v1str, tmp_ptvals[xkey][-1], n_to_take, tmp_ptvals['perfect_vals'][-1])
        # old way of adding a 'no correlation' line:
        # # add a horizontal line at 50 to show what it'd look like if there was no correlation (this is really wasteful... although it does have a satisfying wiggle to it. Now using a plain flat line [below])
        # shuffled_lb_vals = copy.deepcopy(plotvals[lb_metric])
        # random.shuffle(shuffled_lb_vals)
        # NON_corresponding_xvals = [affy for lb, affy in zip(shuffled_lb_vals, plotvals[xvar]) if lb > lb_ptile_val]
        # NON_corr_affy_ptiles = [stats.percentileofscore(plotvals[xvar], caffy, kind='weak') for caffy in NON_corresponding_xvals]
        # tmp_ptvals['reshuffled_vals'].append(numpy.mean(NON_corr_affy_ptiles))
    return tmp_ptvals

# ----------------------------------------------------------------------------------------
def get_mean_ptile_vals(n_clusters, ptile_vals, xvar, debug=False):  # NOTE kind of duplicates code in cf-tree-metrics.py (well except there we're averaging the *difference* between us and perfect
    non_empty_iclusts = [iclust for iclust in range(n_clusters) if 'iclust-%d'%iclust in ptile_vals and len(ptile_vals['iclust-%d'%iclust]['lb_ptiles']) > 0]
    if debug:
        if len(non_empty_iclusts) < n_clusters:
            print '  removed %d empty iclusts' % (n_clusters - len(non_empty_iclusts))
        print_var = 'lb_ptiles' # 'mean_%s_ptiles'%xvar
        fstr = '%5.1f' if ('affinity' in print_var or print_var == 'lb_ptiles') else '%4.2f'
        for iclust in non_empty_iclusts:
            print '    %3d   %s' % (iclust, '  '.join([fstr%v for v in ptile_vals['iclust-%d'%iclust][print_var]]))
    outvals = {k : [] for k in ptile_vals['iclust-%d'%non_empty_iclusts[0]]}
    for ival in range(len(ptile_vals['iclust-%d'%non_empty_iclusts[0]]['lb_ptiles'])):
        for tkey in outvals:
            tmpvals = [ptile_vals['iclust-%d'%iclust][tkey][ival] for iclust in non_empty_iclusts]
            if tkey == 'lb_ptiles':  # they should all be the same
                assert len(set(tmpvals)) == 1
                oval = tmpvals[0]
            else:
                oval = numpy.mean(tmpvals)
            outvals[tkey].append(oval)
    if debug:
        print '      --> %s' % '  '.join([fstr%v for v in outvals[print_var]])
    return outvals

# ----------------------------------------------------------------------------------------
def make_ptile_plot(tmp_ptvals, xvar, plotdir, plotname, xlabel=None, ylabel='?', title=None, fnames=None, ptile_range_tuple=(50., 100., 1.), true_inf_str='?', n_clusters=None, iclust=None, within_cluster_average=False, xlist=None, affy_key='affinities'):
    if len(tmp_ptvals['lb_ptiles']) == 0:
        return

    fig, ax = plotting.mpl_init()
    xia = xvar == 'affinity'
    if xlist is None or len(xlist) == 0:  # this won't actually be right for the latter case, but I don't want to set two different defaults, if it's empty we don't gaf
        xmean = 50
    else:
        xmean = numpy.mean(xlist)  # NOTE this is mean of "xvar", which is the x axis on the scatter plot, but here it's the y axis on the ptile plot
    xkey = 'mean_%s_ptiles' % xvar
    if xlabel is None:
        xlabel = xvar

    ax.plot(tmp_ptvals['lb_ptiles'], tmp_ptvals[xkey], linewidth=3, alpha=0.7)

    # lines corresponding to no correlation and perfect correlation to guide the eye
    bad_args = ((ax.get_xlim(), (xmean, xmean)), {'linewidth' : 3, 'alpha' : 0.7, 'color' : 'darkred', 'linestyle' : '--', 'label' : 'no correlation'})
    perf_args = ((tmp_ptvals['lb_ptiles'], tmp_ptvals['perfect_vals']), {'linewidth' : 3, 'alpha' : 0.7, 'color' : 'darkgreen', 'linestyle' : '--', 'label' : 'perfect correlation'})
    for (args, kwargs) in (perf_args, bad_args) if xia else (bad_args, perf_args):  # shenanigans are so their top/bottom ordering matches the actual lines
        ax.plot(*args, **kwargs)

    if xia:
        # ax.plot(ax.get_xlim(), [50 + 0.5 * x for x in ax.get_xlim()], linewidth=3, alpha=0.7, color='darkgreen', linestyle='--', label='perfect correlation')  # straight line
        # ax.plot(tmp_ptvals['lb_ptiles'], tmp_ptvals['reshuffled_vals'], linewidth=3, alpha=0.7, color='darkred', linestyle='--', label='no correlation')  # reshuffled vals
        ybounds = (45, 105)
        leg_loc = (0.5, 0.2)
        xlabel = xlabel.replace('affinity', 'affinities')
        ptile_ylabel = 'mean percentile of\nresulting %s' % xlabel
    else:
        ymax = max([xmean] + tmp_ptvals[xkey] + tmp_ptvals['perfect_vals'])
        ybounds = (-0.02*ymax, 1.1*ymax)
        leg_loc = (0.5, 0.6)
        ptile_ylabel = 'mean %s\nsince affinity increase' % xlabel

    if 'relative' in affy_key:
        fig.text(0.5, 0.82, 'relative %s' % xvar, fontsize=15, color='red', fontweight='bold')

    if n_clusters is not None:
        if within_cluster_average:
            fig.text(0.25, 0.88, 'within-cluster average over %d families' % n_clusters, fontsize=12, fontweight='bold')  # , color='red'
        else:
            fig.text(0.37, 0.88, 'choosing among %d families' % n_clusters, fontsize=12, fontweight='bold')  # , color='red'
    fn = plotting.mpl_finish(ax, plotdir, plotname, xbounds=ptile_range_tuple, ybounds=ybounds, leg_loc=leg_loc,
                             title='%s %s' % (ungetptlabel(xvar), '' if iclust is None else ', iclust %d'%iclust),
                             xlabel='%s threshold (percentile)' % ylabel, ylabel=ptile_ylabel, adjust={'left' : 0.21}, legend_fontsize=14)
    add_fn(fnames, fn=fn)

# ----------------------------------------------------------------------------------------
def plot_lb_vs_affinity(baseplotdir, lines, lb_metric, ptile_range_tuple=(50., 100., 1.), is_true_line=False, affy_key='affinities', only_csv=False, fnames=None, add_uids=False, colorvar='is_leaf', max_scatter_plot_size=2500, max_iclust_plots=10, make_distribution_plots=False, debug=False):
    # ----------------------------------------------------------------------------------------
    def get_plotvals(line):
        plotvals = {vt : [] for vt in vtypes + ['uids']}
        if colorvar is not None and colorvar == 'is_leaf':
            dtree = treeutils.get_dendro_tree(treestr=get_tree_from_line(line, is_true_line))  # keeping this here to remind myself how to get the tree if I need it
        if affy_key not in line:
            return plotvals
        for uid, affy in [(u, a) for u, a in zip(line['unique_ids'], line[affy_key]) if a is not None]:
            plotvals['affinity'].append(affy)
            if lb_metric in per_seq_metrics:
                plotvals[lb_metric].append(line['tree-info']['lb'][lb_metric][uid])  # NOTE there's lots of entries in the lb info that aren't observed (i.e. aren't in line['unique_ids'])
            if add_uids:
                plotvals['uids'].append(uid)
            if colorvar is not None and colorvar == 'is_leaf':
                node = dtree.find_node_with_taxon_label(uid)
                plotvals['is_leaf'].append(node.is_leaf() if node is not None else None)
        return plotvals
    # ----------------------------------------------------------------------------------------
    def getplotdir(extrastr=''):
        return '%s/%s/%s-vs-%saffinity%s' % (baseplotdir, lb_metric, lb_metric, rel_affy_str('affinity', use_relative_affy='relative' in affy_key), extrastr)
    # ----------------------------------------------------------------------------------------
    def icstr(iclust):
        return '-all-clusters' if iclust is None else '-iclust-%d' % iclust
    # ----------------------------------------------------------------------------------------
    def tmpstrs(iclust, vspstuff):
        lbstr, affystr, clstr = lb_metric, 'affinity', icstr(iclust)
        affystr = '%s%s' % (rel_affy_str('affinity', use_relative_affy='relative' in affy_key), affystr)
        pchoice = 'per-seq' if vspstuff is None else 'per-cluster'
        xlabel, ylabel = affystr, mtitlestr(pchoice, lb_metric)
        title = '%s on %s tree' % (mtitlestr(pchoice, lb_metric, short=True), true_inf_str)
        if vspstuff is not None:
            assert iclust is None
            lbstr = '%s-%s' % (vspstuff[lb_metric], lbstr)
            affystr = '%s-%s' % (vspstuff['affinity'], affystr)
            clstr = '-per-cluster'
            title += ' (per family)'
            xlabel = '%s %s' % (vspstuff['affinity'], xlabel)
            if lb_metric in mean_max_metrics:
                ylabel = '%s%s%s%s' % (vspstuff[lb_metric], ' ' if lb_metric in treeutils.lb_metrics else '(', ylabel, '' if lb_metric in treeutils.lb_metrics else ')')
        else:
            if iclust is None:
                title += ' (%d families together)' % len(lines)
        return lbstr, affystr, clstr, xlabel, ylabel, title
    # ----------------------------------------------------------------------------------------
    def tmpxlabel(iclust, vspstuff):
        _, _, _, xlabel, _, _ = tmpstrs(iclust, vspstuff)
        return xlabel
    # ----------------------------------------------------------------------------------------
    def tmpylabel(iclust, vspstuff):
        _, _, _, _, ylabel, _ = tmpstrs(iclust, vspstuff)
        return ylabel
    # ----------------------------------------------------------------------------------------
    def make_scatter_plot(plotvals, iclust=None, vspstuff=None):
        lbstr, affystr, clstr, xlabel, ylabel, title = tmpstrs(iclust, vspstuff)
        plotname = '%s-vs-%s-%s-tree%s' % (lbstr, affystr, true_inf_str, clstr)
        fn = plot_2d_scatter(plotname, getplotdir(), plotvals, lb_metric, ylabel, title, xlabel=xlabel, colorvar=colorvar if vspstuff is None else None)
        if iclust is None and vspstuff is None:
            add_fn(fnames, fn=fn)
    # ----------------------------------------------------------------------------------------
    def ptile_plotname(iclust=None, vspstuff=None, extra_str=None):
        lbstr, affystr, clstr, _, _, _ = tmpstrs(iclust, vspstuff)
        return '%s-vs-%s-%s-tree-ptiles%s%s' % (lbstr, affystr, true_inf_str, clstr, '' if extra_str is None else '-'+extra_str )
    # ----------------------------------------------------------------------------------------
    def getcorr(xvals, yvals):
        return numpy.corrcoef(xvals, yvals)[0, 1]
    # ----------------------------------------------------------------------------------------
    def getcorrkey(xstr, ystr):
        return '-vs-'.join([xstr, ystr])

    # ----------------------------------------------------------------------------------------
    add_fn(fnames, init=True)
    true_inf_str = 'true' if is_true_line else 'inferred'
    vtypes = get_lbscatteraxes(lb_metric)  # NOTE this puts relative affinity under the (plain) affinity key, which is kind of bad maybe i think probably UPDATE nah i think it's better
    if colorvar is not None:
        vtypes.append(colorvar)

    if not only_csv:
        if sum(len(l['unique_ids']) for l in lines) < max_scatter_plot_size:
            make_lb_scatter_plots('affinity', baseplotdir, lb_metric, lines, fnames=fnames, is_true_line=is_true_line, colorvar='edge-dist' if is_true_line else None,
                                  only_overall='among-families' in lb_metric, only_iclust='within-families' in lb_metric, add_jitter=is_true_line, use_relative_affy='relative' in affy_key)  # there's some code duplication between these two fcns, but oh well
            if make_distribution_plots:
                plot_lb_distributions(lb_metric, baseplotdir, lines, fnames=fnames, is_true_line=is_true_line, only_overall=True)
        else:  # ok this is hackey
            utils.prep_dir(getplotdir(), wildlings=['*.svg', '*.yaml'])
    for estr in ['-ptiles']:  # previous line does a prep_dir() call as well
        utils.prep_dir(getplotdir(estr), wildlings=['*.svg', '*.yaml'])

    per_seq_plotvals = {vt : [] for vt in vtypes}  # plot values for choosing single seqs/cells (only among all clusters, since the iclust ones don't need to kept outside the cluster loop)
    per_clust_plotvals = {vt : {sn : [] for sn, _ in cluster_summary_cfg[vt]} for vt in vtypes}  # each cluster plotted as one point using a summary over its cells (e.g. max, mean) for affinity and lb
    pt_vals = {'per-seq' : {}, 'per-cluster' : {}}  # 'per-seq': choosing single cells, 'per-cluster': choosing clusters; with subkeys in the former both for choosing sequences only within each cluster ('iclust-N', used later in cf-tree-metrics.py to average over all clusters in all processes) and for choosing sequences among all clusters together ('all-clusters')
    # correlation_vals = {'per-seq' : {}, 'per-cluster' : {}}
    if debug:
        print '                 %s   ' % ''.join([('  %-12s'%vt) for vt in vtypes[:2] for _ in cluster_summary_cfg[vt]])
        print '  iclust   size %s' % ''.join(('  %-12s'%st) for vt in vtypes[:2] for st, _ in cluster_summary_cfg[vt])
    for iclust, line in enumerate(lines):
        if debug:
            print '  %3d    %4d' % (iclust, len(line['unique_ids'])),
        iclust_plotvals = get_plotvals(line)  # if it's not in <per_seq_metrics> we still need the affinity values
        if lb_metric in per_seq_metrics:
            if iclust_plotvals[lb_metric].count(0.) == len(iclust_plotvals[lb_metric]):  # i.e. (atm) lbr on family that's only leaves (it would be nice to have a more sensible way to do this, but I guess it's not really a big deal since I think we're done sampling only leaves)
                if debug: print ''
                continue
            for vt in vtypes:
                per_seq_plotvals[vt] += iclust_plotvals[vt]
        for vt in vtypes[:2]:
            for sname, sfcn in cluster_summary_cfg[vt]:
                per_clust_plotvals[vt][sname].append(sfcn(line, iclust_plotvals[vt]))
                if debug:
                    print '%12.3f' % per_clust_plotvals[vt][sname][-1],
        if debug:
            print ''
        if lb_metric not in per_seq_metrics or 'among-families' in lb_metric:
            continue
        iclust_ptile_vals = get_ptile_vals(lb_metric, iclust_plotvals, 'affinity', 'affinity', dbgstr='iclust %d'%iclust, affy_key=affy_key, debug=debug)
        pt_vals['per-seq']['iclust-%d'%iclust] = iclust_ptile_vals
        # correlation_vals['per-seq']['iclust-%d'%iclust] = {getcorrkey(*vtypes[:2]) : getcorr(*[iclust_plotvals[vt] for vt in vtypes[:2]])}
        if not only_csv and len(iclust_plotvals['affinity']) > 0 and iclust < max_iclust_plots:
            # make_scatter_plot(iclust_plotvals, iclust=iclust)  # making these with make_lb_scatter_plots() now
            make_ptile_plot(iclust_ptile_vals, 'affinity', getplotdir('-ptiles'), ptile_plotname(iclust=iclust),
                            ylabel=tmpylabel(iclust, None), title=mtitlestr('per-seq', lb_metric, short=True), true_inf_str=true_inf_str, iclust=iclust, affy_key=affy_key)

    if lb_metric in per_seq_metrics and 'within-families' not in lb_metric:
        if per_seq_plotvals[lb_metric].count(0.) == len(per_seq_plotvals[lb_metric]):
            return
        # correlation_vals['per-seq']['all-clusters'] = {getcorrkey(*vtypes[:2]) : getcorr(*[per_seq_plotvals[vt] for vt in vtypes[:2]])}
        pt_vals['per-seq']['all-clusters'] = get_ptile_vals(lb_metric, per_seq_plotvals, 'affinity', 'affinity', affy_key=affy_key, debug=debug)  # choosing single cells from from every cell from every cluster together
        if not only_csv and len(per_seq_plotvals[lb_metric]) > 0:
            # make_scatter_plot(per_seq_plotvals)  # making these with make_lb_scatter_plots() now
            make_ptile_plot(pt_vals['per-seq']['all-clusters'], 'affinity', getplotdir('-ptiles'), ptile_plotname(),
                            ylabel=tmpylabel(None, None), title=mtitlestr('per-seq', lb_metric, short=True), fnames=fnames, true_inf_str=true_inf_str, n_clusters=len(lines), affy_key=affy_key)

    if lb_metric in per_seq_metrics and 'among-families' not in lb_metric and not only_csv:
        pt_vals['per-seq']['within-families-mean'] = get_mean_ptile_vals(len(lines), pt_vals['per-seq'], 'affinity')  # within-family mean
        make_ptile_plot(pt_vals['per-seq']['within-families-mean'], 'affinity', getplotdir('-ptiles'), ptile_plotname(iclust=None, extra_str='within-cluster-average'),
                        ylabel=tmpylabel(None, None), title=mtitlestr('per-seq', lb_metric, short=True), fnames=fnames, true_inf_str=true_inf_str, n_clusters=len(lines), within_cluster_average=True, affy_key=affy_key)

    # turning of the per-cluster plots for now, since it doesn't really make sense to choose families (and the plotting is the only really slow part)
    # # per-cluster plots
    # if 'dtr' not in lb_metric:
    #     # TODO combine with dtr training variable stuff
    #     for sn1, _ in cluster_summary_cfg[vtypes[0]]:  # I tried really hard to work out a way to get this in one (cleaner) loop
    #         for sn2, _ in cluster_summary_cfg[vtypes[1]]:
    #             vspairs = zip(vtypes[:2], (sn1, sn2))  # assign this (sn1, st2) combo to lb and affinity based on their order in <vtypes> (although now that we're using a double loop this is even weirder)
    #             vspdict = {v : s for v, s in vspairs}  # need to also access it by key
    #             tmpvals = {vt : per_clust_plotvals[vt][sn] for vt, sn in vspairs}  # e.g. 'affinity' : <max affinity value list>, 'lbi' : <mean lbi value list>
    #             tkey = getcorrkey('affinity%s' % vspdict['affinity'], '%s-%s' % (vspdict[lb_metric], lb_metric))  # can't use <vtypes> because of the stupid <affy_key_str> UPDATE now I've removed <affy_key_str> I should be able to UPDATE re-added affy_key/relative affinity stuff, but now I don't really care about this per-cluster stuff, but I think i don't need to change anything here? but not really sure
    #             # correlation_vals['per-cluster'][tkey] = getcorr(tmpvals['affinity'], tmpvals[lb_metric])
    #             tmp_ptile_vals = get_ptile_vals(lb_metric, tmpvals, 'affinity', 'affinity', affy_key=affy_key, debug=debug)
    #             pt_vals['per-cluster'][tkey] = tmp_ptile_vals
    #             if not only_csv:
    #                 make_scatter_plot(tmpvals, vspstuff=vspdict)
    #                 make_ptile_plot(tmp_ptile_vals, 'affinity', getplotdir('-ptiles'), ptile_plotname(vspstuff=vspdict),
    #                                 xlabel=tmpxlabel(None, vspdict), ylabel=tmpylabel(None, vspdict), title=mtitlestr('per-cluster', lb_metric, short=True), true_inf_str=true_inf_str, affy_key=affy_key)

    with open('%s/%s.yaml' % (getplotdir('-ptiles'), ptile_plotname()), 'w') as yfile:
        yamlfo = {'percentiles' : pt_vals} #, 'correlations' : correlation_vals}
        json.dump(yamlfo, yfile)

# ----------------------------------------------------------------------------------------
def plot_lb_vs_ancestral_delta_affinity(baseplotdir, lines, lb_metric, ptile_range_tuple=(50., 100., 1.), is_true_line=False, only_csv=False, fnames=None, max_scatter_plot_size=2500, max_iclust_plots=10, debug=False):
    # plot lb[ir] vs both number of ancestors and branch length to nearest affinity decrease (well, decrease as you move upwards in the tree/backwards in time)
    # ----------------------------------------------------------------------------------------
    def check_affinity_changes(affinity_changes):
        affinity_changes = sorted(affinity_changes)
        if debug:
            print '    checking affinity changes for negative values and unexpected variation: %s' % ' '.join(['%.4f' % a for a in affinity_changes])  # well, the variation isn't really unexpected, but it's worth keeping in mind
        if len(affinity_changes) == 0:
            if debug:
                print '      %s empty affinity changes list' % utils.color('yellow', 'note')
            return
        if any(a < 0. for a in affinity_changes):
            print '  %s negative affinity changes in %s' % (utils.color('red', 'error'), ' '.join(['%.4f' % a for a in affinity_changes]))
        max_diff = affinity_changes[-1] - affinity_changes[0]
        # if abs(max_diff) / numpy.mean(affinity_changes) > 0.2:  # this is almost always true, which is fine, and I don't really plan on doing anything to change it soon (it would be nice to at some point use a performance metric gives us credit for differential prediction of different affinity change magnitudes, but oh well)
        #     print'      %s not all affinity increases were the same size (min: %.4f   max: %.4f   abs(diff) / mean: %.4f' % (utils.color('yellow', 'warning'), affinity_changes[0], affinity_changes[-1], abs(max_diff) / numpy.mean(affinity_changes))
    # ----------------------------------------------------------------------------------------
    def get_plotvals(line, xvar):
        plotvals = {vt : [] for vt in [lb_metric, xvar]}  # , 'uids']}
        dtree = treeutils.get_dendro_tree(treestr=line['tree'])
        affinity_changes = []
        for uid in line['unique_ids']:
            node = dtree.find_node_with_taxon_label(uid)
            if node is dtree.seed_node:  # root doesn't have any ancestors
                continue
            lbval = line['tree-info']['lb'][lb_metric][uid]  # NOTE there's lots of entries in the lb info that aren't observed (i.e. aren't in line['unique_ids'])
            if lb_metric == 'lbr' and lbval == 0:  # lbr equals 0 should really be treated as None/missing
                continue
            n_steps, branch_len = treeutils.get_n_ancestors_to_affy_change(node, dtree, line, affinity_changes=affinity_changes, also_return_branch_len=True, debug=debug)  # also adds to <affinity_changes>
            if n_steps is None:
                continue
            plotvals[xvar].append(n_steps if xvar == 'n-ancestor' else branch_len)
            plotvals[lb_metric].append(lbval)
            # plotvals['uids'].append(uid)
        check_affinity_changes(affinity_changes)
        return plotvals

    # ----------------------------------------------------------------------------------------
    def getplotdir(xvar, extrastr=''):
        return '%s/%s-vs-%s%s' % (baseplotdir, lb_metric, xvar, extrastr)
    # ----------------------------------------------------------------------------------------
    def icstr(iclust):
        return '-all-clusters' if iclust is None else '-iclust-%d' % iclust
    # ----------------------------------------------------------------------------------------
    def make_scatter_plot(plotvals, xvar, iclust=None, tfns=None):
        if len(plotvals[xvar]) == 0 or len(plotvals[xvar]) > max_scatter_plot_size:
            return
        title = '%s on %s tree%s' % (mtitlestr('per-seq', lb_metric, short=True), true_inf_str, (' (%d families together)' % len(lines)) if iclust is None else ' (cluster %d)'%iclust)
        fn = plot_2d_scatter('%s-vs-%s-%s-tree%s' % (lb_metric, xvar, true_inf_str, icstr(iclust)), getplotdir(xvar), plotvals, lb_metric, mtitlestr('per-seq', lb_metric), title, xvar=xvar, xlabel='%s since affinity increase' % xlabel, log='y' if lb_metric == 'lbr' else '')
        add_fn(tfns, fn=fn)
    # ----------------------------------------------------------------------------------------
    def ptile_plotname(xvar, iclust, extra_str=None):
        return '%s-vs-%s-%s-tree-ptiles%s%s' % (lb_metric, xvar, true_inf_str, icstr(iclust), '' if extra_str is None else '-'+extra_str)

    # ----------------------------------------------------------------------------------------
    add_fn(fnames, init=True)  # oh, wait, I think this won't accomplish anything if init is actually needed
    true_inf_str = 'true' if is_true_line else 'inferred'
    xvar_list = collections.OrderedDict([(xvar, xlabel) for metric, cfglist in lb_metric_axis_cfg('lbr') for xvar, xlabel in cfglist])
    for xvar, estr in itertools.product(xvar_list, ['', '-ptiles']):
        utils.prep_dir(getplotdir(xvar, extrastr=estr), wildlings=['*.svg', '*.yaml'])
    if debug:
        print 'finding ancestors with most recent affinity increases'
    iclust_fnames = add_fn(None, init=True)
    for xvar, xlabel in xvar_list.items():
        per_seq_plotvals = {vt : [] for vt in [lb_metric, xvar]}  # , 'uids']}
        # not yet implemented: per_clust_plotvals = {st : {vt : [] for vt in vtypes} for st in cluster_summary_fcns}  # each cluster plotted as one point using a summary over its cells (max or mean) for affinity and lb
        pt_vals = {'per-seq' : {}, 'per-cluster' : {}}  # 'per-seq': choosing single cells, 'per-cluster': choosing clusters; with subkeys in the former both for choosing sequences only within each cluster ('iclust-N', used later in cf-tree-metrics.py to average over all clusters in all processes) and for choosing sequences among all clusters together ('all-clusters')
        # not yet implemented: correlation_vals = {'per-seq' : {}, 'per-cluster' : {}}
        for iclust, line in enumerate(lines):
            if debug:
                if iclust == 0: print ' %s' % utils.color('green', xvar)
                print '  %s' % utils.color('blue', 'iclust %d' % iclust)
                print '         node        ancestors  distance   affinity (%sX: change for chosen ancestor, %s: reached root without finding lower-affinity ancestor)' % (utils.color('red', '+'), utils.color('green', 'x'))
            iclust_plotvals = get_plotvals(line, xvar)
            for vtype in per_seq_plotvals:
                per_seq_plotvals[vtype] += iclust_plotvals[vtype]
            if 'among-families' in lb_metric:
                continue
            iclust_ptile_vals = get_ptile_vals(lb_metric, iclust_plotvals, xvar, xlabel, dbgstr='iclust %d'%iclust, debug=debug)
            pt_vals['per-seq']['iclust-%d'%iclust] = iclust_ptile_vals
            if not only_csv and iclust < max_iclust_plots:
                make_scatter_plot(iclust_plotvals, xvar, iclust=iclust, tfns=iclust_fnames if iclust < 3 else None)
                make_ptile_plot(iclust_ptile_vals, xvar, getplotdir(xvar, extrastr='-ptiles'), ptile_plotname(xvar, iclust), xlist=iclust_plotvals[xvar],
                                xlabel=xlabel, ylabel=mtitlestr('per-seq', lb_metric), true_inf_str=true_inf_str, iclust=iclust)
        fnames += iclust_fnames

        if 'among-families' not in lb_metric and not only_csv:
            pt_vals['per-seq']['within-families-mean'] = get_mean_ptile_vals(len(lines), pt_vals['per-seq'], xvar)
            make_ptile_plot(pt_vals['per-seq']['within-families-mean'], xvar, getplotdir(xvar, extrastr='-ptiles'), ptile_plotname(xvar, None, extra_str='within-cluster-average'),
                            xlabel=xlabel, ylabel=mtitlestr('per-seq', lb_metric), fnames=fnames, true_inf_str=true_inf_str, n_clusters=len(lines), within_cluster_average=True, xlist=per_seq_plotvals[xvar])

        if 'within-families' not in lb_metric:
            pt_vals['per-seq']['all-clusters'] = get_ptile_vals(lb_metric, per_seq_plotvals, xvar, xlabel, dbgstr='all clusters', debug=debug)  # "averaged" might be a better name than "all", but that's longer
            if not only_csv:
                # make_scatter_plot(per_seq_plotvals, xvar, tfns=fnames)
                make_ptile_plot(pt_vals['per-seq']['all-clusters'], xvar, getplotdir(xvar, extrastr='-ptiles'), ptile_plotname(xvar, None), xlist=per_seq_plotvals[xvar],
                                xlabel=xlabel, ylabel=mtitlestr('per-seq', lb_metric), fnames=fnames, true_inf_str=true_inf_str, n_clusters=len(lines))

        with open('%s/%s.yaml' % (getplotdir(xvar, extrastr='-ptiles'), ptile_plotname(xvar, None)), 'w') as yfile:
            yamlfo = {'percentiles' : pt_vals}
            json.dump(yamlfo, yfile)  # not adding the new correlation keys atm (like in the lb vs affinity fcn)

# ----------------------------------------------------------------------------------------
def plot_true_vs_inferred_lb(plotdir, true_lines, inf_lines, lb_metric, fnames=None, debug=False):
    plotvals = {val_type : {uid : l['tree-info']['lb'][lb_metric][uid] for l in lines for uid in l['unique_ids'] if uid in l['tree-info']['lb'][lb_metric]}  # NOTE there's lots of entries in the lb info that aren't observed (i.e. aren't in line['unique_ids'])
                for val_type, lines in (('true', true_lines), ('inf', inf_lines))}
    common_uids = set(plotvals['true']) & set(plotvals['inf'])  # there should/may be a bunch of internal nodes in the simulation lines but not in the inferred lines, but otherwise they should have the same uids
    # tmpvals = sorted([(u, plotvals['true'][u], plotvals['inf'][u]) for u in common_uids], key=lambda x: x[2] / x[1] if x[1] != 0 else 0)
    # for x in tmpvals:
    #     print ' %12s  %8.5f  %8.5f  %8.5f' % (x[0], x[2] / x[1] if x[1] != 0 else 0, x[1], x[2])
    # sys.exit()
    plotvals = {val_type : [plotvals[val_type][uid] for uid in common_uids] for val_type in plotvals}
    plotname = '%s-true-vs-inferred' % lb_metric
    mtstr = mtitlestr('per-seq', lb_metric)
    fn = plot_2d_scatter(plotname, plotdir, plotvals, 'inf', '%s on inferred tree' % mtstr, 'true vs inferred %s' % mtstr, xvar='true', xlabel='%s on true tree' % mtstr, stats='correlation')
    add_fn(fnames, fn=fn)

# ----------------------------------------------------------------------------------------
def plot_cons_seq_accuracy(baseplotdir, lines, n_total_bin_size=10000, fnames=None, debug=False):  # n_total_bin_size: merge together everybody with family size that's this close to each other
    csinfo = {}
    for line in lines:
        tmpfo = utils.get_cons_seq_accuracy_vs_n_sampled_seqs(line)
        n_total = len(line['unique_ids'])
        if n_total_bin_size is not None:
            n_total = n_total - n_total % n_total_bin_size
        if n_total not in csinfo:
            csinfo[n_total] = []
        csinfo[n_total].append(tmpfo)

    ctypes = ['nuc', 'aa']
    for ctype in ctypes:
        fig, ax = plotting.mpl_init()
        for n_total, subfos in csinfo.items():
            all_n_sampleds = [sfo[ctype]['n_sampled'] for sfo in subfos]
            # assert len(set(tuple(l) for l in all_n_sampleds)) == 1
            n_sample_list = sorted(all_n_sampleds, key=len)[0]  # they should all be similar lengths, but some families can be a bit smaller, so will be missing the last few values, so just ignore those last few values
            all_y_vals = [[sfo[ctype]['hdists'][i_n_sampled] for sfo in subfos] for i_n_sampled in range(len(n_sample_list))]
            plotvals = {'xvals' : n_sample_list,
                        'yvals' : [numpy.mean(all_y_vals[i_n_sampled]) for i_n_sampled in range(len(n_sample_list))],
                        'yerrs' : [numpy.std(all_y_vals[i_n_sampled], ddof=1) / math.sqrt(len(all_y_vals[i_n_sampled])) for i_n_sampled in range(len(n_sample_list))],
            }
            label = str(n_total) if len(csinfo) > 1 else '%d-%d' % (min(len(l['unique_ids']) for l in lines), max(len(l['unique_ids']) for l in lines))
            ax.errorbar(plotvals['xvals'], plotvals['yvals'], yerr=plotvals['yerrs'], label=label, alpha=0.7, linewidth=4, markersize=15, marker='.')
        metric = 'cons-dist-' + ctype
        def ctypetitle(ct):
            return ct.upper() if ct == 'aa' else ct + '.'
        fn = plotting.mpl_finish(ax, baseplotdir + '/' + metric, '%s-accuracy' % metric,
                                 xticks=n_sample_list, xticklabels=[str(x) for x in n_sample_list],
                                 xlabel='N sampled', ylabel='hamming dist. to\nfull-family cons seq', title='%s cons. seq accuracy'%ctypetitle(ctype),
                                 leg_title='N total', leg_loc=[0.7, 0.6], log='x')
        add_fn(fnames, fn=fn)

# ----------------------------------------------------------------------------------------
def get_lb_tree_cmd(treestr, outfname, lb_metric, affy_key, ete_path, subworkdir, metafo=None, tree_style=None, queries_to_include=None):
    treefname = '%s/tree.nwk' % subworkdir
    metafname = '%s/meta.yaml' % subworkdir
    if not os.path.exists(subworkdir):
        os.makedirs(subworkdir)
    with open(treefname, 'w') as treefile:
        treefile.write(treestr)
    cmdstr = '%s/bin/plot-lb-tree.py --treefname %s' % (utils.get_partis_dir(), treefname)
    if metafo is not None:
        with open(metafname, 'w') as metafile:
            yaml.dump(metafo, metafile)
        cmdstr += ' --metafname %s' % metafname
    if queries_to_include is not None:
        cmdstr += ' --queries-to-include %s' % ':'.join(queries_to_include)
    cmdstr += ' --outfname %s' % outfname
    cmdstr += ' --lb-metric %s' % lb_metric
    cmdstr += ' --affy-key %s' % utils.reversed_input_metafile_keys[affy_key]
    # cmdstr += ' --lb-tau %f' % lb_tau
    if lb_metric == 'lbr':
        cmdstr += ' --log-lbr'
    if tree_style is not None:
        cmdstr += ' --tree-style %s' % tree_style
    cmdstr, _ = utils.run_ete_script(cmdstr, ete_path, return_for_cmdfos=True, tmpdir=subworkdir, extra_str='        ')

    return {'cmd_str' : cmdstr, 'workdir' : subworkdir, 'outfname' : outfname, 'workfnames' : [treefname, metafname]}

# ----------------------------------------------------------------------------------------
def plot_lb_trees(metric_methods, baseplotdir, lines, ete_path, base_workdir, is_true_line=False, tree_style=None, queries_to_include=None):
    workdir = '%s/ete3-plots' % base_workdir
    plotdir = baseplotdir + '/trees'
    utils.prep_dir(plotdir, wildlings='*.svg')

    if not os.path.exists(workdir):
        os.makedirs(workdir)
    cmdfos = []
    for lb_metric in metric_methods:
        for iclust, line in enumerate(lines):  # note that <min_selection_metric_cluster_size> was already applied in treeutils
            treestr = get_tree_from_line(line, is_true_line)
            affy_key = 'affinities'  # turning off possibility of using relative affinity for now
            metafo = copy.deepcopy(line['tree-info']['lb'])  # NOTE there's lots of entries in the lb info that aren't observed (i.e. aren't in line['unique_ids'])
            if affy_key in line:  # either 'affinities' or 'relative_affinities'
                metafo[utils.reversed_input_metafile_keys[affy_key]] = {uid : affy for uid, affy in zip(line['unique_ids'], line[affy_key])}
            outfname = '%s/%s-tree-iclust-%d%s.svg' % (plotdir, lb_metric, iclust, '-relative' if 'relative' in affy_key else '')
            cmdfos += [get_lb_tree_cmd(treestr, outfname, lb_metric, affy_key, ete_path, '%s/sub-%d' % (workdir, len(cmdfos)), metafo=metafo, tree_style=tree_style, queries_to_include=queries_to_include)]

    if len(cmdfos) > 0:
        start = time.time()
        utils.run_cmds(cmdfos, clean_on_success=True, shell=True, n_max_procs=10, proc_limit_str='plot-lb-tree.py')  # I'm not sure what the max number of procs is, but with 21 it's crashing with some of them not able to connect to the X server, and I don't see a big benefit to running them all at once anyways
        print '    made %d ete tree plots (%.1fs)' % (len(cmdfos), time.time() - start)

    os.rmdir(workdir)

# ----------------------------------------------------------------------------------------
def plot_per_mutation_lonr(plotdir, lines_to_use, reco_info):
    fig, ax = plotting.mpl_init()

    plotvals = {'lonr' : [], 'affinity_change' : []}
    for line in lines_to_use:
        true_affinities = {uid : reco_info[uid]['affinities'][0] for uid in line['unique_ids']}
        nodefos = line['tree-info']['lonr']['nodes']
        for lfo in line['tree-info']['lonr']['values']:
            if lfo['parent'] not in true_affinities:
                print '    %s parent \'%s\' not in true affinities, skipping lonr values' % (utils.color('red', 'warning'), lfo['parent'])
                continue
            if lfo['child'] not in true_affinities:
                print '    %s child \'%s\' not in true affinities, skipping lonr values' % (utils.color('red', 'warning'), lfo['child'])
                continue

            plotvals['lonr'].append(lfo['lonr'])
            plotvals['affinity_change'].append(true_affinities[lfo['child']] - true_affinities[lfo['parent']])

    ax.scatter(plotvals['affinity_change'], plotvals['lonr'], alpha=0.7) #, info['ccf_under'][meth], label='clonal fraction', color='#cc0000', linewidth=4)
    plotname = 'lonr-per-mut-vs-affinity'
    plotting.mpl_finish(ax, plotdir, plotname, xlabel='change in affinity', ylabel='LONR') #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

# ----------------------------------------------------------------------------------------
def plot_aggregate_lonr(plotdir, lines_to_use, reco_info, debug=False):
    fig, ax = plotting.mpl_init()
    plotvals = {'S' : [], 'NS' : []}
    for line in lines_to_use:
        for lfo in line['tree-info']['lonr']['values']:
            if lfo['synonymous']:
                plotvals['S'].append(lfo['lonr'])
            else:
                plotvals['NS'].append(lfo['lonr'])
    # ax.plot(plotvals['S'], label='S', linewidth=3, alpha=0.7)
    # ax.plot(plotvals['NS'], label='NS', linewidth=3, alpha=0.7)
    xmin, xmax = [mfcn([x for mtlist in plotvals.values() for x in mtlist]) for mfcn in (min, max)]
    hists = {mt : Hist(30, xmin, xmax, value_list=plotvals[mt], title=mt, xtitle='LONR', ytitle='mutations') for mt in plotvals}
    plotname = 'lonr-ns-vs-s'

    lonr_score = hists['NS'].get_mean() - hists['S'].get_mean()
    draw_no_root(hists['NS'], more_hists=[hists['S']], plotname=plotname, plotdir=plotdir, alphas=[0.7, 0.7], plottitle='NS - S: %.2f' % lonr_score, errors=True, remove_empty_bins=True)

    # for mt, hist in hists.items():
    #     hist.mpl_plot(ax, label=mt, remove_empty_bins=True)
    # plotting.mpl_finish(ax, plotdir, plotname, xlabel='LONR', ylabel='mutations') #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)
