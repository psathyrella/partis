from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import copy
import pickle
from scipy import stats # this is fucking slow
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
import yaml
from io import open
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from . import utils
from .hist import Hist
from . import treeutils
from . import plotting
from . import hutils
from . import mds
from . import plotconfig

# ----------------------------------------------------------------------------------------
# this name is terrible, but it's complicated and I can't think of a better one
def lb_metric_axis_cfg(metric_method, final_plots=False):  # x axis variables against which we plot each lb metric (well, they're the x axis on the scatter plots, not the ptile plots)
    base_cfg = collections.OrderedDict([('lbi', [('affinity', 'affinity')]),
                                        ('lbr', [('n-ancestor', 'N ancestors'), ('daffy', 'affinity change')]), #, ('branch-length', 'branch length')]),  # turning off branch length at least for now (for run time reasons)
                                        ('lbf', [('n-ancestor', 'N ancestors')]),  # , ('branch-length', 'branch length')])  # turning off branch length at least for now (for run time reasons)
    ])
    if metric_method is None:
       return list(base_cfg.items())
    base_cfg['dtr'] = [('affinity', 'affinity'),] if final_plots else [('affinity', 'affinity'), ('n-ancestor', 'N ancestors')]  # hack hack hack
    if metric_method in base_cfg:
        return [(m, cfg) for m, cfg in base_cfg.items() if m == metric_method]
    else:  # shm, delta-lbi, cons-dist-*, etc
        xv = 'n-ancestor' if metric_method in ['delta-lbi', 'aa-lbr', 'aa-lbf'] else 'affinity'  # also hack hack hack
        return [[metric_method, [(xv, xv.replace('n-a', 'N a'))]]]

# ----------------------------------------------------------------------------------------
def ptile_range_tuple(xvar):
    if xvar in ['daffy', 'n-ancestor']:
        return (90., 100., 0.2)
    else:
        return (50., 100., 2.)

# ----------------------------------------------------------------------------------------
def single_lbma_cfg_vars(metric_method, final_plots=False):
    assert metric_method is not None  # only use this for single metrics
    return lb_metric_axis_cfg(metric_method, final_plots=final_plots)[0][1]  # [0] gets you the first metric's cfg (there's ony one because we specified a single <metric_method>, [1] gets you the (var, label)

# ----------------------------------------------------------------------------------------
def add_lbma_cfg_permutations(cfg_list, include_relative_affy_plots=False, make_distr_hists=False):  # NOTE acts on sub lists of the return value of the above (i.e. the .values())
    cfg_list = [[s, l, False, False] for s, l in cfg_list]
    if include_relative_affy_plots:
        for ic, (s, l, u) in enumerate(cfg_list):
            if s == 'affinity':  # add it just after the existing (non-relative) 'affinity'
                cfg_list.insert(ic + 1, ['affinity', 'affinity', True, False])
                break
    if make_distr_hists:
        assert not include_relative_affy_plots  # if i want to turn relative affinity plots back on, i'd just need to add both permutations here, but i don't care right now
        cfg_list.append(['affinity', 'affinity', False, True])
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
    n_to_take = max(1, int(frac * len(vals)))  # NOTE don't use numpy.percentile(), since affinity is fairly discrete-valued, which causes bad stuff (e.g. you don't take anywhere near the number of cells that you were trying to)
    return numpy.mean(sorted(vals)[len(vals) - n_to_take:])
mean_max_metrics = ['lbi', 'lbr', 'lbf', 'aa-lbi', 'aa-lbr', 'aa-lbf', 'shm', 'shm-aa', 'cons-lbi']
mean_max_metrics += treeutils.dtr_metrics
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

per_seq_metrics = ['lbi', 'lbr', 'lbf', 'aa-lbi', 'aa-lbr', 'lbf', 'shm', 'shm-aa', 'cons-dist-nuc', 'cons-dist-aa', 'delta-lbi', 'cons-lbi']
per_seq_metrics += treeutils.dtr_metrics
# per_clust_metrics = ('lbi', 'lbr', 'shm', 'fay-wu-h', 'cons-dist-nuc')  # don't need this atm since it's just all of them (note that 'cons-dist-nuc' doesn't really make sense here, see cluster_summary_cfg)
mtitle_cfg = {'per-seq' : {'cons-dist-nuc' : '- nuc distance to cons seq',
                           'cons-dist-aa' : '- AA distance to cons seq',
                           'delta-lbi' : 'change in lb index',
                           'z-score-err' : 'z score diff (lb - affy)',
                           'edge-dist' : 'root/tip dist',
                           'affinity-ptile' : 'affinity percentile',
                           'lbi-ptile' : 'lbi percentile',
                           'lbr-ptile' : 'lbr percentile',
                           'within-families-affinity-dtr' : 'in-family dtr',
                           'within-families-delta-affinity-dtr' : 'in-family dtr',
                           'among-families-affinity-dtr' : 'among-families dtr',
                           'among-families-delta-affinity-dtr' : 'among-families dtr'},
              'per-cluster' : {'fay-wu-h' : '- Fay-Wu H',
                               'cons-seq-shm-nuc' : 'N mutations in cons seq',
                               'affinity' : 'top quintile affinity'}}
mtitle_cfg['per-seq'].update(treeutils.legtexts)
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
def all_clust_str(n_clusters):
    return '%d cluster%s' % (n_clusters, utils.plural(n_clusters))

# ----------------------------------------------------------------------------------------
def add_fn(fnames, fn=None, init=False, n_per_row=4, new_row=False):
    if init: assert fnames is None  # otherwise there's no effect
    if fnames is None:
        if init:
            fnames = []
        else:
            return
    if len(fnames) == 0 or len(fnames[-1]) >= n_per_row or new_row:
        fnames.append([])
    if fn is not None:
        fnames[-1].append(fn)
    return fnames

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_selection_hists(histfname, plotdir, plotname, plot_all=False, n_plots=7, title='', xlabel='', fnames=None):
    # NOTE would be nice to use make_single_joyplot() in plot_bcr_phylo_selection_hists()
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
        with open(hfname, 'rb') as runstatfile:
            numpyhists = pickle.load(runstatfile, encoding='latin1')
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
            assert len(bin_contents) == len(bin_edges) - 1  # one time this failed because they were both length zero, but i'm not bothering to fix it since it shouldn't really happen again
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
        print('  %s no/empty hists in %s' % (utils.color('yellow', 'warning'), histfname))
        return
    jpdata = []
    for hist in all_hists:
        jpdata.append([x for x, y in zip(hist.get_bin_centers(), hist.bin_contents) for _ in range(int(y)) if x > xmin and x < xmax])  # NOTE this is repeating the 'for _ in range()' in the fcn above, but that's because I used to be actually using the Hist()s, and maybe I will again

    fbin_xmins, fbin_xmaxs = zip(*[h.get_filled_bin_xbounds(extra_pads=2) for h in all_hists])
    xmin_filled, xmax_filled = min(fbin_xmins), max(fbin_xmaxs)

    pre_fig, pre_ax = plotting.mpl_init()  # not sure to what extent these really get used after joypy is done with things
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # i don't know why it has to warn me that it's clearing the fig/ax I'm passing in, and I don't know how else to stop it
        fig, axes = joypy.joyplot(jpdata, labels=all_ylabels, fade=True, hist=True, overlap=1., ax=pre_ax, x_range=(xmin_filled, xmax_filled), bins=int(xmax_filled - xmin_filled), xlabelsize=15) #, ylabelsize=15)
    # xtextpos = 0.85 * (xmax_filled - xmin_filled) + xmin_filled  # this is from before I found transform=ax.transAxes, but I don't want to remove it yet
    fsize = 15
    for ax, lab in zip(axes, all_xtralabels):
        ax.text(0.85, 0.2, lab, fontsize=fsize, transform=ax.transAxes)
    fig.text(0.03, 0.9, 'generation', fontsize=fsize)
    fig.text(0.77, 0.9, '(mean, N cells)', fontsize=fsize)
    # NOTE do *not* set your own x ticks/labels in the next line, since they'll be in the wrong place (i.e. not the same as where joypy puts them) (also note, the stupid y labels don't work, but setting them on the joyplot axes also doesn't work)
    fn = plotting.mpl_finish(pre_ax, plotdir, plotname, title=title, xlabel=xlabel) #, ylabel='generation') #, leg_loc=(0.7, 0.45)) #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)
    add_fn(fnames, fn=fn)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_kd_vals(plotdir, event, legstr='', tdlabel='', fnames=None):
    # ----------------------------------------------------------------------------------------
    def scatter_plot(xvar, xlabel, yvar='affinity', ylabel='1/Kd', jitter=''):
        # # new_cmap = plotting.truncate_colormap(plt.cm.Blues, 0, 1)
        # # ax.hexbin(kd_changes, shms, gridsize=25, cmap=plt.cm.Blues) #, info['ccf_under'][meth], label='clonal fraction', color='#cc0000', linewidth=4)
        plotname = '%s-vs-%s' % (xvar, yvar)
        if legstr != '':
            plotname += '-' + legstr
        fn = plot_2d_scatter(plotname, plotdir, plotvals, yvar, ylabel, '%s sampled cells'%legstr, xvar=xvar, xlabel=xlabel, jitter=jitter, jitter_fracs={v : 0.05 if 'aa' in v else 0.01 for v in (xvar, yvar)})
        add_fn(fnames, fn=fn)
    # ----------------------------------------------------------------------------------------
    # scatter plot of kd/affinity vs shm
    plotvals = {'shm-nuc' : [], 'shm-aa' : [], 'target-distance' : [], 'lambda' : [], 'affinity' : []}
    for iseq, uid in enumerate(event['unique_ids']):
        plotvals['shm-nuc'].append(event['n_mutations'][iseq])
        plotvals['shm-aa'].append(utils.antnval(event, 'shm-aa', iseq))
        plotvals['target-distance'].append(utils.antnval(event, 'min_target_distances', iseq))
        plotvals['lambda'].append(utils.antnval(event, 'lambdas', iseq))
        plotvals['affinity'].append(event['affinities'][iseq])
    scatter_plot('shm-nuc', 'N nuc mutations', jitter='x')
    scatter_plot('shm-aa', 'N AA mutations', jitter='x')
    scatter_plot('target-distance', '%s distance to nearest target seq' % tdlabel) #, jitter=True)
    scatter_plot('target-distance', '%s distance to nearest target seq' % tdlabel, yvar='shm-aa', ylabel='N AA mutations', jitter='xy')
    scatter_plot('lambda', 'lambda')

    # make kd changes histogram
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
        if legstr != '':
            plotname += '-' + legstr
        fn = plotting.mpl_finish(ax, plotdir,  plotname, xlabel='parent-child kd change', ylabel='branches', log='y', title='%s sampled cells'%legstr) #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)
        add_fn(fnames, fn=fn, new_row=True)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_target_attraction(plotdir, event, legstr='', fnames=None):  # plots of which sequences are going toward which targets
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
    if legstr != '':
        plotname += '-' + legstr
    fn = plotting.mpl_finish(ax, plotdir, plotname, xlabel='index (identity) of nearest target sequence', ylabel='counts', title=legstr) #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)
    add_fn(fnames, fn=fn) #, new_row=True)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_mds(plotdir, line, legstr='', fnames=None):
    # labels = {u : i for u, i in zip(line['unique_ids'], line['nearest_target_indices'])}
    labels = {u : int(i) for u, i in zip(line['unique_ids'], line['min_target_distances'])}  # NOTE truncating to int here
    color_scale_vals = {u : i for u, i in zip(line['unique_ids'], line['min_target_distances'])}
    seqfos = utils.seqfos_from_line(line, aa=True)
    qtis = {}
    for itg, tseq in enumerate(line['target_seqs']):
        tid = 'target-%d'%itg
        seqfos.append({'name' : tid, 'seq' : utils.ltranslate(tseq)})
        qtis[tid] = tid
        # color_scale_vals[tid] = 0  # it's better to just not have it there since it often screws up the scale (e.g. if it's far from the other seqs)
        # labels[tid] = 0
    plotname = 'mds'
    if legstr != '':
        plotname += '-' + legstr
    # color_scale_vals = None
    labels = None
    leg_title = 'target dist.'
    mds.run_bios2mds(2, None, seqfos, plotdir + '/work', 1, aligned=True, remove_duplicates=True, plotdir=plotdir, plotname=plotname, queries_to_include=qtis, title='', leg_title=leg_title, cmapstr='viridis', color_scale_vals=color_scale_vals, labels=labels)
    os.rmdir(plotdir + '/work')
    add_fn(fnames, fn='%s/%s.svg'%(plotdir, plotname), new_row=True)
    add_fn(fnames, fn='%s/mds-legend.svg'%plotdir)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_simulation(plotdir, outdir, events, extrastr, metric_for_target_distance_label, lpair=None, plot_mds=False):
    # ----------------------------------------------------------------------------------------
    def pltevt(evt, legstr=''):
        plot_bcr_phylo_kd_vals(plotdir, evt, legstr=legstr, tdlabel=metric_for_target_distance_label, fnames=fnames)
        plot_bcr_phylo_target_attraction(plotdir, evt, legstr=legstr, fnames=fnames)
        if plot_mds:
            plot_bcr_phylo_mds(plotdir, evt, fnames=fnames)  # eh, doesn't really show anything useful (e.g. distance in plot doesn't relate well at all to target distance)
    # ----------------------------------------------------------------------------------------
    utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'], subdirs=['work'] if plot_mds else None)
    fnames = add_fn(None, init=True)

    legstr = '' if lpair is None else '+'.join(lpair)+' '
    plot_bcr_phylo_selection_hists('%s/%s_n_mutated_nuc_hdists.p' % (outdir, extrastr), plotdir, 'n-mutated-nuc-all-cells', title=legstr+'nuc SHM all cells', xlabel='N nucleotide mutations to naive', fnames=fnames)
    plot_bcr_phylo_selection_hists('%s/%s_n_mutated_aa_hdists.p' % (outdir, extrastr), plotdir, 'n-mutated-aa-all-cells', title=legstr+'AA SHM all cells', xlabel='N AA mutations to naive', fnames=fnames)
    plot_bcr_phylo_selection_hists('%s/%s_min_aa_target_hdists.p' % (outdir, extrastr), plotdir, 'min-aa-target-all-cells', title=legstr+'target dist. all cells', xlabel='%s distance to nearest target seq' % metric_for_target_distance_label, fnames=fnames)
    plot_bcr_phylo_selection_hists('%s/%s_sampled_min_aa_target_hdists.p' % (outdir, extrastr), plotdir, 'min-aa-target-sampled-cells', plot_all=True, title=legstr+'target dist. sampled cells (excluding ancestor sampling)', xlabel='%s distance to nearest target seq' % metric_for_target_distance_label, fnames=fnames)

    if lpair is None:
        assert len(events) == 1
        pltevt(events[0])
    else:
        assert len(events) == 2
        for ltmp, evt in zip(lpair, events):
            pltevt(evt, legstr=ltmp)

    plotting.make_html(plotdir, fnames=fnames)

# ----------------------------------------------------------------------------------------
def plot_subtree_purity(plotdir, base_plotname, dtree, antn, meta_key, meta_emph_formats=None, only_csv=False, swarm_plots=False, max_size=0):
    # ----------------------------------------------------------------------------------------
    def vlabel(var):
        return plotconfig.xtitles.get('subtree-purity-'+var)
    # ----------------------------------------------------------------------------------------
    def make_scatter_plot(yvar):
        add_jitter = True
        szvals = [sub_stat['size'] for tstats in st_stats.values() for sub_stat in tstats]
        xmax = int(max(max_size, max(szvals)))
        djit = 0.2 #0.01 * (xmax - min(szvals))
        fig, ax = plotting.mpl_init()
        for im, (mval, tstats) in enumerate(st_stats.items()):
            jitval = ((im+1) // 2) * djit * (-1)**im if add_jitter else 0
            for sub_stat in tstats:
                ax.scatter([sub_stat['size'] + jitval], [sub_stat[yvar]], facecolor=mcolors[mval], alpha=0.6) #, s=10)
        if xmax <= 5:
            xticks = [1, 2, 3, 4, 5]
        else:
            xticks = [1, 3, 5] + list(range(10, xmax-1, 5)) + [xmax]
        fn = plotting.mpl_finish(ax, plotdir, '%s-size-vs-%s'%(base_plotname, yvar), title='', xlabel=vlabel('size')+(' (+offset)' if add_jitter else ''), ylabel=vlabel(yvar), xticks=xticks, xbounds=[0, xmax+1])  # +1 is to account for jitter (well, might be nice to have it even without jitter)
        fnames[1].append(os.path.basename(fn).replace('.svg', ''))
    # ----------------------------------------------------------------------------------------
    st_nodes, st_stats = treeutils.find_pure_subtrees(dtree, antn, meta_key)
    if st_nodes is None:
        return [[]]
    all_emph_vals = set(st_stats)  # a lot of this stuff duplicates plotting.stack_meta_hists()
    all_emph_vals, emph_colors = plotting.meta_emph_init(meta_key, formats=meta_emph_formats, all_emph_vals=all_emph_vals)
    mcolors = {v : c for v, c in emph_colors}
    hkeys = ['size', 'mean-ancestor-distance', 'mean-root-depth']
    int_keys = ['size']

    # first do hists
    hist_lists, hist_colors = {tk : [] for tk in hkeys}, []
    plotvals = {tk : {} for tk in hkeys}
    for mval in sorted(st_stats):
        tstats = st_stats[mval]
        for tkey in hkeys:
            vlist = [s[tkey] for s in tstats]
            if tkey in int_keys:
                thist = Hist(init_int_bins=True, value_list=vlist, xtitle=vlabel(tkey), title=str(mval))
            else:
                n_bins = 10
                # xbins = hutils.autobins(vlist, n_bins)
                xbins, n_bins = hutils.auto_volume_bins(vlist, n_bins)
                thist = Hist(n_bins=n_bins, xmin=xbins[0], xmax=xbins[-1], xbins=xbins, value_list=vlist, xtitle=vlabel(tkey), title=str(mval))
            hist_lists[tkey].append(thist)
            plotvals[tkey][mval] = vlist
        hist_colors.append(mcolors[mval])
    fnames = [[], [], []]
    for tkey, hlist in hist_lists.items():  # would be nice to use stack_meta_hists() also for the hists here
        plotname = '%s-%s'%(base_plotname, tkey)
        plotting.draw_no_root(None, more_hists=hlist, plotname=plotname, colors=list(hist_colors), plotdir=plotdir, write_csv=True, only_csv=only_csv, shift_overflows=True,
                              leg_title='%s'%meta_key.rstrip('s'), alphas=[0.6 for _ in hlist], linewidths=[5, 3, 2], markersizes=[15, 10, 8], errors=True, remove_empty_bins=True, log='y' if tkey=='size' else '', plottitle='') #, normalize=True)
        fnames[0].append(plotname)

        if swarm_plots:
            plotting.stack_meta_hists(plotname, plotdir, meta_key, plotvals[tkey], formats=meta_emph_formats, no_hist=True, swarm_plots=True, xtitle=vlabel(tkey))
            fnames[-1].append('swarm-'+plotname)

    # then 2d plots
    for yvar in ['mean-ancestor-distance', 'mean-root-depth']:
        make_scatter_plot(yvar)

    fn = plotting.make_meta_info_legend(plotdir, 'subtree-size', meta_key.rstrip('s'), emph_colors, all_emph_vals, meta_emph_formats=meta_emph_formats, alpha=0.7)
    fnames[1].append(fn)

    return fnames

# ----------------------------------------------------------------------------------------
# NOTE partially duplicates partitionplotter.get_treefos()
def get_tree_in_line(line, is_true_line, aa=False):  # NOTE unlike treeutils.get_trees_for_annotations() this just looks to see if there's a tree already in the line
    if is_true_line and not aa:  # this is weird and ugly
        return line['tree']
    if 'tree-info' not in line:  # if 'tree-info' is missing, it should be because it's a small cluster in data that we skipped when calculating lb values
        return None
    # print '    %s tree' % utils.color('yellow', 'aa' if aa else 'nuc')
    return line['tree-info']['lb'].get('%stree'%('aa-' if aa else ''))  # also returns None if this key is missing (e.g. if we've only calculated cons-dist-aa)

# ----------------------------------------------------------------------------------------
# NOTE that this isn't symmetric wrt x/y vars -- some combos require one var to be x, and the other y (otherwise it'll crash cause it can't figure out how to calculate the values)
def make_lb_scatter_plots(xvar, baseplotdir, lb_metric, lines_to_use, fnames=None, is_true_line=False, colorvar=None, only_overall=False, only_iclust=False, add_uids=False, yvar=None, choose_among_families=False,
                          add_jitter=False, min_ptile=80., n_iclust_plot_fnames=None, use_relative_affy=False, queries_to_include=None, meta_info_to_emphasize=None, meta_emph_formats=None, add_stats=None, xlabel=None, ylabel=None, title_str=''):  # <is_true_line> is there because we want the true and inferred lines to keep their trees in different places, because the true line just has the one, true, tree, while the inferred line could have a number of them (yes, this means I maybe should have called it the 'true-tree' or something)
    # ----------------------------------------------------------------------------------------
    def add_warn(tstr, targs):
        if 'warn_text' in targs:
            targs['warn_text'] += '\n'
        else:
            targs['warn_text'] = ''
        targs['warn_text'] += tstr
    # ----------------------------------------------------------------------------------------
    if colorvar is not None and colorvar in [xvar, yvar]:
        raise Exception('colorvar \'%s\' can\'t be one of xvar, yvar (%s %s)' % (colorvar, xvar, yvar))
    if only_overall and only_iclust:
        print('  %s only_overall and only_iclust are both set, so not plotting anything' % utils.color('yellow', 'warning'))
        return
    if queries_to_include is not None:
        add_uids = True
    if meta_info_to_emphasize is not None:
        add_uids = True
        meta_emph_key, meta_emph_val = list(meta_info_to_emphasize.items())[0]
    cdist_pt_keys = [s+'-ptile' for s in cdist_keys]
    if yvar is None:
        yvar = lb_metric

    choice_str = 'among-families' if choose_among_families else 'within-families'
    if xlabel is None:
        xlabel = mtitlestr('per-seq', xvar).replace('- N', 'N')
    if ylabel is None:
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
    basetitle = '' # '%s %s vs %s%s ' % ('true' if is_true_line else 'inferred', mtitlestr('per-seq', yvar, short=True), mtitlestr('per-seq', xvar, short=True).replace('- N', 'N'), title_str)  # here 'shm' the plain number of mutations, not 'shm' the non-lb metric, so we have to fiddle with the label in mtitle_cfg
    scatter_kwargs = {'xvar' : xvar, 'xlabel' : xlabel, 'colorvar' : colorvar, 'leg_loc' : (0.55, 0.75), 'log' : 'y' if 'lbr' in yvar else '', 'stats' : add_stats}
    if use_relative_affy:
        add_warn('relative affinity', scatter_kwargs)
    sorted_lines = sorted(lines_to_use, key=lambda l: len(l['unique_ids']), reverse=True)
    n_total_null = {k : 0 for k in [xvar, yvar]}
    for iclust, line in enumerate(sorted_lines):  # get depth/n_mutations for each node
        iclust_plotvals = {x : [] for x in vtypes}
        lbfo = line['tree-info']['lb'][lb_metric]
        if len(set(lbfo) - set(line['unique_ids'])) > 0 or len(set(line['unique_ids']) - set(lbfo)) > 0:
            if 'ptile' in xvar or 'ptile' in yvar:  # if it's not a ptile, we just don't plot it unless we have a value for both x and y. But for ptiles I don't really want to think through right now how to handle it
                raise Exception('  %s different uids in <lbfo> and <line>, and xvar or yvar is a ptile: %3d extra uids in lbfo, %3d extra in line, %3d in common. This (the code below) needs to be checked.' % (utils.color('yellow', 'warning'), len(set(lbfo) - set(line['unique_ids'])), len(set(line['unique_ids']) - set(lbfo)), len(set(line['unique_ids']) & set(lbfo))))
        if colorvar in ['is_leaf', 'edge-dist'] or xvar == 'edge-dist':
            dtree = treeutils.get_dendro_tree(treestr=get_tree_in_line(line, is_true_line, aa='aa-lb' in lb_metric))
            if dtree is None:
                continue

        # TODO I should really combine the xvar and yvar stuff here
        if xvar == 'shm':
            def xvalfcn(i): return line['n_mutations'][i]
        elif xvar == 'affinity':
            def xvalfcn(i): return line['affinities'][i]
        elif xvar in cdist_keys:
            treeutils.add_cons_dists(line, aa='aa' in xvar)
            tkey = xvar.replace('cons-dist-', 'cons_dists_')
            def xvalfcn(i): return -line[tkey][i]
        elif xvar in ['aa-lbi', 'aa-lbr']:
            def xvalfcn(i): return line['tree-info']['lb'][xvar].get(line['unique_ids'][i], None)
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
            treeutils.add_cons_dists(line, aa='aa' in yvar)
            tkey = yvar.replace('-ptile', '').replace('cons-dist-', 'cons_dists_')
            cvals = [-v for v in line[tkey]]
            if not choose_among_families:
                cdist_ptiles = {u : stats.percentileofscore(cvals, cvals[i], kind='weak') for i, u in enumerate(line['unique_ids'])}
            def yvalfcn(i): return cdist_ptiles[line['unique_ids'][i]]
        else:
            assert False

        n_null = {k : 0 for k in [xvar, yvar]}
        for iseq, uid in enumerate(line['unique_ids']):
            xval = xvalfcn(iseq)
            # if xvar == 'affinity-ptile' and '-ptile' in yvar and xvalfcn(iseq) < min_ptile:  # and yvalfcn(iseq) < min_ptile:  the number of cells with high lbi but low affinity (last, commented criterion) is just too small to bother plotting -- all our errors come from the other direction
            #     continue
            yval = yvalfcn(iseq)
            if ('lbr' in yvar or 'lbf' in yvar) and yval == 0:
                continue
            if None in (xval, yval):
                for k, v in [(xvar, xval), (yvar, yval)]:
                    if v is None:
                        n_null[k] += 1
                        n_total_null[k] += 1
                continue
            iclust_plotvals[xvar].append(xval)
            iclust_plotvals[yvar].append(yval)
            if colorvar is not None: # and colorvar not in [xvar, yvar]:  # if its one of [xval, yval] we don't need to do anything
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
                if meta_info_to_emphasize is None:
                    iclust_plotvals['uids'].append(uid if queries_to_include is None or uid in queries_to_include else None)  # use to add None here instead of <uid> if this node didn't have an affinity value, but that seems unnecessary, I can worry about uid config options later when I actually use the uid dots for something
                else:
                    # estr = utils.meta_emph_str(meta_emph_key, meta_emph_val, formats=meta_emph_formats)  # use this if you want the emph value to show up, but for now i like having the uid
                    uval = None
                    if utils.meta_info_equal(meta_emph_key, meta_emph_val, utils.antnval(line, meta_emph_key, iseq), formats=meta_emph_formats):
                        uval = utils.antnval(line, 'alternate-uids', iseq, use_default=True, default_val=uid)
                    iclust_plotvals['uids'].append(uval)
        if len(iclust_plotvals[xvar]) == 0:
            continue
        iskargs = copy.deepcopy(scatter_kwargs)
        if add_jitter: #xvar == 'affinity-ptile' and '-ptile' in yvar:
            iclust_plotvals[xvar] = plotting.add_jitter(iclust_plotvals[xvar], delta=100. - min_ptile if '-ptile' in xvar else None, frac=0.02)
            iskargs['xlabel'] += ' (+jitter)'
        if not only_overall:
            title = '%siclust %d: %d observed, %d total' % (basetitle, iclust, len(line['unique_ids']), len(lbfo))
            tmpplotname = '%s%s-vs-%s%s-iclust-%d' % (rel_affy_str(yvar, use_relative_affy), yvar, rel_affy_str(yvar, use_relative_affy), xvar, iclust)
            for vname in [v for v in [xvar, yvar] if n_null[v] > 0]:
                add_warn('%s for %d / %d' % (vname, len(line['unique_ids']) - n_null[vname], len(line['unique_ids'])), iskargs)
            fn = plot_2d_scatter(tmpplotname, plotdir, iclust_plotvals, yvar, ylabel, title, **iskargs)
            if n_iclust_plot_fnames is not None and iclust < n_iclust_plot_fnames:
                add_fn(fnames, fn=fn)
        assert len(set([len(plotvals[vt]) for vt in plotvals])) == 1  # make sure all of them are the same length
        for vtype in plotvals:
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

    if not only_iclust and len(sorted_lines) > 1:
        tmpplotname = '%s%s-vs-%s%s-all-clusters' % (rel_affy_str(yvar, use_relative_affy), yvar, rel_affy_str(xvar, use_relative_affy), xvar)
        for vname in [v for v in [xvar, yvar] if n_total_null[v] > 0]:
            n_tot = sum(len(l['unique_ids']) for l in sorted_lines)
            add_warn('%s for %d / %d' % (vname, n_tot - n_total_null[vname], n_tot), scatter_kwargs)
        if add_jitter:
            scatter_kwargs['xlabel'] += ' (+jitter)'
        fn = plot_2d_scatter(tmpplotname, plotdir, plotvals, yvar, ylabel, '%s%s' % (basetitle, all_clust_str(len(sorted_lines))), **scatter_kwargs)
        add_fn(fnames, fn=fn)

# ----------------------------------------------------------------------------------------
def plot_lb_distributions(lb_metric, baseplotdir, lines_to_use, is_true_line=False, fnames=None, only_overall=False, n_iclust_plot_fnames=None, stats=''):
    def make_hist(plotvals, n_total, n_skipped, iclust=None):
        if len(plotvals) == 0:
            return
        if lb_metric in ['cons-dist-aa', 'cons-dist-nuc', 'shm', 'shm-aa']:
            hist = Hist(value_list=plotvals, init_int_bins=True)
        else:
            xmin, xmax = hutils.get_expanded_bounds([min(plotvals), max(plotvals)], 0)
            hist = Hist(30, xmin, xmax, value_list=plotvals)
        texts = []
        if 'mean' in stats:
            texts.append((0.7, 0.8, 'mean %.3f' % numpy.mean(plotvals)))
        if 'max' in stats:
            texts.append((0.7, 0.75, 'max %.3f' % max(plotvals)))
        plotname = '%s-%s' % (lb_metric, str(iclust) if iclust is not None else 'all-clusters')
        leafskipstr = ', skipped %d leaves' % n_skipped if n_skipped > 0 else ''  # ok they're not necessarily leaves, but almost all of them are leaves (and not really sure how a non-leaf could get zero, but some of them seem to)
        title = '%s %s: %d entries%s (%s)' % ('true' if is_true_line else 'inferred', mtitlestr('per-seq', lb_metric, short=True), n_total, leafskipstr, all_clust_str(len(sorted_lines)) if iclust is None else 'iclust %d'%iclust)
        fn = hist.fullplot(plotdir, plotname, fargs={'xlabel' : mtitlestr('per-seq', lb_metric), 'log' : 'y', 'ylabel' : 'counts', 'title' : title}, texts=texts)
        if iclust is None or (n_iclust_plot_fnames is not None and iclust < n_iclust_plot_fnames):
            add_fn(fnames, fn=fn)

    sorted_lines = sorted([l for l in lines_to_use if 'tree-info' in l], key=lambda l: len(l['unique_ids']), reverse=True)  # if 'tree-info' is missing, it should be because it's a small cluster we skipped when calculating lb values

    plotvals = []
    n_total_skipped_leaves = 0
    plotdir = '%s/%s/distributions' % (baseplotdir, lb_metric)
    utils.prep_dir(plotdir, wildlings=['*.svg', '*.csv'])
    for iclust, line in enumerate(sorted_lines):
        lbfo = line['tree-info']['lb'][lb_metric]  # NOTE this contains any ancestor nodes that the phylogenetic program has in the tree, but that aren't in <line['unique_ids']>
        if is_true_line:
            iclust_plotvals = [lbfo[u] for u in line['unique_ids'] if u in lbfo]  # for the true plots, we *don't* want to include any inferred ancestor nodes that weren't actually sampled, since we don't have affinity info for them, and it'd make it look like there's a mismatch between these distributions and the lb vs affinity plots (which won't have them)
        else:
            iclust_plotvals = list(lbfo.values())  # whereas for real data, we want to include the inferred ancestor nodes for which we don't have sequences (although I guess in most cases where we're really interested in them, we would've used a phylogenetics program that also inferred their sequences, so they'd presumably have been added to <line['unique_ids']>)
        cluster_size = len(iclust_plotvals)  # i.e. including leaves
        if 'lbr' in lb_metric or 'lbf' in lb_metric:
            iclust_plotvals = [v for v in iclust_plotvals if v > 0.]  # don't plot the leaf values, they just make the plot unreadable
        if not only_overall and len(sorted_lines) > 1:
            make_hist(iclust_plotvals, cluster_size, cluster_size - len(iclust_plotvals), iclust=iclust)
        plotvals += iclust_plotvals
        n_total_skipped_leaves += cluster_size - len(iclust_plotvals)
    make_hist(plotvals, len(plotvals) + n_total_skipped_leaves, n_total_skipped_leaves)
    if n_iclust_plot_fnames is not None:
        fnames.append([])

# ----------------------------------------------------------------------------------------
def make_lb_affinity_joyplots(plotdir, lines, lb_metric, fnames=None, n_clusters_per_joy_plot=25, n_max_joy_plots=25, n_plots_per_row=4):
    partition = utils.get_partition_from_annotation_list(lines)
    annotation_dict = utils.get_annotation_dict(lines)
    affy_lists = {k : [a for a in l['affinities'] if a is not None] for k, l in annotation_dict.items()}  # just skip None-type affinites, which i think is what we want? (note that the single joyplot fcn gets these lists separately -- these are just used for max vals over all clusters)
    if sum(len(l) for l in affy_lists.values()) == 0:
        return
    sorted_clusters = sorted(partition, key=lambda c: mean_of_top_quintile(affy_lists[':'.join(c)]), reverse=True)
    sorted_cluster_groups = [sorted_clusters[i : i + n_clusters_per_joy_plot] for i in range(0, len(sorted_clusters), n_clusters_per_joy_plot)]
    repertoire_size = sum([len(c) for c in sorted_clusters])
    max_affinity = max([a for alist in affy_lists.values() for a in alist])  # it's nice to keep track of the max values over the whole repertoire so all plots can have the same max values
    def lbv(c, u): return annotation_dict[':'.join(c)]['tree-info']['lb'][lb_metric].get(u)  # see note about affy_lists
    all_lb_vals = [lbv(c, u) for c in sorted_clusters for u in c if lbv(c, u) is not None]  # NOTE don't use all the values in the dict in 'tree-info', since non-sampled sequences (i.e. usually intermediate ancestors) are in there
    if len(all_lb_vals)==0 or 'lb' in lb_metric and max(all_lb_vals)==0.:  # at least atm, this means this is lbr on a family with no common ancestor sampling
        return
    if len(set(all_lb_vals)) == 1:
        return
    # print 'divided repertoire of size %d with %d clusters into %d cluster groups' % (repertoire_size, len(sorted_clusters), len(sorted_cluster_groups))
    iclustergroup = 0
    for subclusters in sorted_cluster_groups:
        if iclustergroup > n_max_joy_plots:
            continue
        title = 'affinity and %s (%d / %d)' % (mtitlestr('per-seq', lb_metric), iclustergroup + 1, len(sorted_cluster_groups))  # NOTE it's important that this denominator is still right even when we don't make plots for all the clusters (which it is, now)
        fn = plotting.make_single_joyplot(subclusters, annotation_dict, repertoire_size, plotdir, '%s-affinity-joyplot-%d' % (lb_metric, iclustergroup), x1key='affinities', x1label='affinity', x2key=lb_metric, x2label=mtitlestr('per-seq', lb_metric),
                                          global_max_vals={'affinities' : max_affinity, lb_metric : max(all_lb_vals)}, title=title, remove_none_vals=True, sortlabel='mean top quintile aff.')  # note that we can't really add <cluster_indices> like we do in partitionplotter.py, since (i think?) the only place there's per-cluster plots we'd want to correspond to is in the bcr phylo simulation dir, which has indices unrelated to anything we're sorting by here, and which we can't reconstruct
        add_fn(fnames, fn=fn)
        iclustergroup += 1

# ----------------------------------------------------------------------------------------
def plot_2d_scatter(plotname, plotdir, plotvals, yvar, ylabel, title, xvar='affinity', xlabel='affinity', colorvar=None, log='', leg_loc=None, warn_text=None, markersize=15, stats=None, jitter='', jitter_fracs=None):
    leafcolors = {'leaf' : '#006600', 'internal' : '#2b65ec'}  # green, blue
    chosecolors = {'chosen' : '#990012', 'nope' : '#006600'}  # red, green
    if len(plotvals[xvar]) == 0:
        # print '    no %s vs affy info' % yvar
        return '%s/%s.svg' % (plotdir, plotname)
    xvals, yvals = [plotvals[v] for v in (xvar, yvar)]
    if jitter_fracs is None:
        jitter_fracs = {}
    if 'x' in jitter:
        xvals = plotting.add_jitter(xvals, frac=jitter_fracs.get(xvar))
        xlabel += ' (+jitter)'
    if 'y' in jitter:
        yvals = plotting.add_jitter(yvals, frac=jitter_fracs.get(xvar))
        ylabel += ' (+jitter)'
    fig, ax = plotting.mpl_init()
    alpha = 0.4
    if colorvar is None:
        ax.scatter(xvals, yvals, alpha=0.4)
    else:
        if colorvar == 'is_leaf':
            colorfcn = lambda x: leafcolors['leaf' if x else 'internal']
        if colorvar == 'chosen':
            colorfcn = lambda x: chosecolors[x]
        else:
            smap = plotting.get_normalized_scalar_map([v for v in plotvals[colorvar] if v is not None], 'viridis')
            colorfcn = lambda x: 'grey' if x is None else plotting.get_smap_color(smap, None, val=x)
            alpha = 0.65
        for x, y, cval in zip(xvals, yvals, plotvals[colorvar]):  # we used to do the leaf/internal plots as two scatter() calls, which might be faster? but I think what really takes the time is writing the svgs, so whatever
            ax.plot([x], [y], color=colorfcn(cval), marker='.', markersize=markersize, alpha=alpha)
    if 'uids' in plotvals:
        for xval, yval, uid in zip(xvals, yvals, plotvals['uids']):  # note: two ways to signal not to do this: sometimes we have 'uids' in the dict, but don't fill it (so the zip() gives an empty list), but sometimes we populate 'uids' with None values
            if uid is None:
                continue
            # ax.plot([xval], [yval], color='red', marker='.', markersize=markersize)
            ax.text(xval, yval, uid, color='red', fontsize=8)

    if warn_text is not None:
        fig.text(0.6 if len(warn_text) < 15 else 0.4, 0.75 if len(warn_text) < 15 else 0.82, warn_text, fontsize=30 if len(warn_text) < 15 else 15, fontweight='bold', color='red')
    xmin, xmax = [mfcn([x for x in xvals if x is not None]) for mfcn in [min, max]]
    ymin, ymax = [mfcn([y for y in yvals if y is not None]) for mfcn in [min, max]]
    # if xvar == 'cons-dist-aa':
    #     xmin = max(-15, xmin)
    xbounds = xmin - 0.02 * (xmax - xmin), xmax + 0.02 * (xmax - xmin)
    if 'y' in log:
        ybounds = 0.75 * ymin, 1.3 * ymax
    else:
        ybounds = ymin - 0.03 * (ymax - ymin), ymax + 0.08 * (ymax - ymin)
    if yvar in ['shm']+cdist_keys:
        ax.plot([xmin, xmax], [0, 0], linewidth=1, alpha=0.7, color='grey')
    leg_title, leg_prop, leg_iter = None, None, []
    if colorvar is not None:
        leg_loc = [0.05 if xvar in cdist_keys+['affinity'] else 0.7, 0.55]  # I think this is sometimes overriding the one that's passed in
        if yvar == 'cons-dist-aa':
            leg_loc = [0.75, 0.15]
        leg_prop = {'size' : 12}
        if colorvar == 'is_leaf':
            leg_iter = [(leafcolors[l], l) for l in ['leaf', 'internal']]
        elif colorvar == 'chosen':
            leg_iter = [(chosecolors[l], l) for l in ['chosen', 'nope']]
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
            pcorr = numpy.corrcoef(xvals, yvals)[0, 1]
            if set([xvar, yvar]) == set(['lbi', 'aa-lbi']) and pcorr > 0.85:
                print('        %s correlation between lbi and aa-lbi is suspiciously high %.3f, which suggests that there weren\'t enough inferred ancestral sequences to rescale the nuc tree to amino acids, i.e. aa-lbi may have in effect basically been calculated on the nuc tree' % (utils.color('yellow', 'warning'), pcorr))
            sx, sy = (0.25, 0.8) if (colorvar is None or 'cons-dist' in yvar) else (leg_loc[0] + 0.15, leg_loc[1] + 0.28)
            fig.text(sx, sy, 'r = %.3f' % pcorr, fontsize=20, fontweight='bold') #, color='red')
    fn = plotting.mpl_finish(ax, plotdir, plotname, title=title, xlabel=xlabel, ylabel=ylabel, xbounds=xbounds, ybounds=ybounds, log=log, leg_loc=leg_loc, leg_title=leg_title, leg_prop=leg_prop)
    return fn

# ----------------------------------------------------------------------------------------
def get_actual_ptval(plotvals, pvr, ptile):
    npval = numpy.percentile(plotvals[pvr], ptile)  # numpy sometimes returns a value in between actual values in the list, or due to precision issues even slightly larger or smaller than any value in the list
    nearvals = sorted(set(plotvals[pvr]), key=lambda x: abs(x - npval))  # sort by nearness to <npval>
    return nearvals[0]  # return the nearest one

# ----------------------------------------------------------------------------------------
def get_ptile_vals(lb_metric, plotvals, xvar, xlabel, dbgstr=None, use_relative_affy=False, return_distr_hists=False, max_bin_width=0.2, min_bins=30, debug=False):
    # NOTE xvar and xlabel refer to the x axis on the scatter plot from which we make this ptile plot (i.e. are affinity, N ancestors, or branch length). On this ptile plot it's the y axis. (I tried calling it something else, but it was more confusing)
    if xvar == 'n-ancestor':
        plotvals = copy.deepcopy(plotvals)
        plotvals[xvar] = [abs(v) for v in plotvals[xvar]]  # TODO this maybe isn't what we want, see note in ptile plot fcn below
    xia = xvar == 'affinity'
    xkey = 'mean_%s_ptiles' % xvar
    tmp_ptvals = {'lb_ptiles' : [], xkey : [], 'perfect_vals' : []}  # , 'reshuffled_vals' : []}
    if len(plotvals[xvar]) == 0:
        return tmp_ptvals
    if debug:
        print('              getting ptile vals%s' % ('' if dbgstr is None else (' for %s' % utils.color('blue', dbgstr))))
        print('                      %-12s N      mean     %s     |  perfect   perfect' % (lb_metric, 'mean   ' if xia else ''))
        print('              ptile   threshold  taken    %s%-s %s  |  N taken  mean %s'  % (utils.color('red', 'rel-') if xia and use_relative_affy else '', xlabel.replace('affinity', 'affy'), '   affy ptile ' if xia else '', 'ptile' if xia else xlabel))
    max_lb_val = max(plotvals[lb_metric])
    sorted_xvals = sorted(plotvals[xvar], reverse=xia or xvar == 'daffy')
    if xia:
        np_arr_sorted_xvals = numpy.asarray(sorted_xvals)  # the stats calls in the next two lines make this conversion, and it's way too slow to do for every x
        corr_ptile_vals = [stats.percentileofscore(np_arr_sorted_xvals, x, kind='weak') for x in plotvals[xvar]]  # x (e.g. affinity) percentile of each plot val (could speed this up by only using the best half of each list (since ptiles are starting at 50))
        perf_ptile_vals = sorted(corr_ptile_vals, reverse=True)  # x (e.g. affinity) percentile or each plot val, *after sorting by x* (e.g. affinity)
    for percentile in numpy.arange(*ptile_range_tuple(xvar)):
        lb_ptile_val = get_actual_ptval(plotvals, lb_metric, percentile)  # lb value corresponding to <percentile> (the min() is because the numpy call has precision issues that sometimes give you a value (very very slightly) larger than any of the actual values in the list)
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
            f1str = '5.0f' if xia else '8.4f'
            f2str = '5.0f' if xia else ('8.2f' if xvar == 'n-ancestor' else '8.6f')
            print(('             %5.0f   %6.2f     %4d   %s  %'+f1str+'       |  %4d      %-'+f2str) % (percentile, lb_ptile_val, len(final_xvar_vals), v1str, tmp_ptvals[xkey][-1], n_to_take, tmp_ptvals['perfect_vals'][-1]))
        # old way of adding a 'no correlation' line:
        # # add a horizontal line at 50 to show what it'd look like if there was no correlation (this is really wasteful... although it does have a satisfying wiggle to it. Now using a plain flat line [below])
        # shuffled_lb_vals = copy.deepcopy(plotvals[lb_metric])
        # random.shuffle(shuffled_lb_vals)
        # NON_corresponding_xvals = [affy for lb, affy in zip(shuffled_lb_vals, plotvals[xvar]) if lb > lb_ptile_val]
        # NON_corr_affy_ptiles = [stats.percentileofscore(plotvals[xvar], caffy, kind='weak') for caffy in NON_corresponding_xvals]
        # tmp_ptvals['reshuffled_vals'].append(numpy.mean(NON_corr_affy_ptiles))

    if return_distr_hists:
        # ----------------------------------------------------------------------------------------
        def gethist(htype, ptbound, xmin, xmax, n_bins, ptval):
            assert htype in ['lo', 'hi']
            def vfcn(a): return a <= ptval if htype=='lo' else a >= ptval
            lblist, xvlist = zip(*[(lb, xv) for lb, xv in zip(plotvals[lb_metric], plotvals[xvar]) if vfcn(xv)])
            tstr = '%s %.0f%%\nmean %s' % ('bottom' if htype=='lo' else 'top', ptbound if htype=='lo' else 100. - ptbound, utils.round_to_n_digits(numpy.mean(xvlist), 2))
            if debug:
                print('  %s %.0f%%: affy %s= %.4f  %s' % (htype, ptbound, '<' if htype=='lo' else '>', ptval, sorted(lblist)))
            return Hist(n_bins=n_bins, xmin=xmin, xmax=xmax, title=tstr, value_list=lblist)
        # ----------------------------------------------------------------------------------------
        lower_ptile, upper_ptile = 20, 80
        lo_affy_ptile_val = get_actual_ptval(plotvals, xvar, lower_ptile)
        hi_affy_ptile_val = get_actual_ptval(plotvals, xvar, upper_ptile)
        if lb_metric in treeutils.int_metrics:
            xmin, xmax = min(plotvals[lb_metric]) - 0.5, max(plotvals[lb_metric]) + 0.5
            n_bins = xmax - xmin
            # if n_bins > 1.5 * min_bins:  # min_bins is the min number of bins for float vars, and here we're using it as the max bins for int vars, but... whatever
            #     n_bins = int(n_bins / 2.)
        else:
            xmin, xmax = 0., 1.01 * max(plotvals[lb_metric])  # this is the low edge of the overflow bin, so needs to be a bit above the biggest value
            n_bins = max(min_bins, int((xmax - xmin) / float(max_bin_width)))  # for super large lbr values like 100 you need way more bins
        distr_hists = [gethist(htp, ptb, xmin, xmax, n_bins, ptv) for htp, ptb, ptv in zip(['lo', 'hi'], [lower_ptile, upper_ptile], [lo_affy_ptile_val, hi_affy_ptile_val])]
        return tmp_ptvals, distr_hists
    else:
        return tmp_ptvals

# ----------------------------------------------------------------------------------------
# take a dict <ptile_vals> with a list of values for each cluster (e.g. 'lb_ptiles', 'perfect_vals' and 'mean_affinity_ptiles' where each list is over percentile values e.g. from 50 to 100), and return the same length list for each variable but averaged over clusters
def get_mean_ptile_vals(n_clusters, ptile_vals, xvar, debug=False):  # NOTE kind of duplicates code in scanplot.py (well except there we're averaging the *difference* between us and perfect
    non_empty_iclusts = [iclust for iclust in range(n_clusters) if 'iclust-%d'%iclust in ptile_vals and len(ptile_vals['iclust-%d'%iclust]['lb_ptiles']) > 0]
    if len(non_empty_iclusts) == 0:
        return {}  # not really sure this is right, adding it long after writing this
    if debug:
        if len(non_empty_iclusts) < n_clusters:
            print('  removed %d empty iclusts' % (n_clusters - len(non_empty_iclusts)))
        print_var = 'lb_ptiles' # 'mean_%s_ptiles'%xvar
        fstr = '%5.1f' if ('affinity' in print_var or print_var == 'lb_ptiles') else '%4.2f'
        for iclust in non_empty_iclusts:
            print('    %3d   %s' % (iclust, '  '.join([fstr%v for v in ptile_vals['iclust-%d'%iclust][print_var]])))
    outvals = {k : [] for k in ptile_vals['iclust-%d'%non_empty_iclusts[0]]}
    n_vals = len(ptile_vals['iclust-%d'%non_empty_iclusts[0]]['lb_ptiles'])
    for ival in range(n_vals):
        for tkey in outvals:
            tmpvals = [ptile_vals['iclust-%d'%iclust][tkey][ival] for iclust in non_empty_iclusts]  # list of values for variable <tkey>, one for each non-empty cluster
            if tkey == 'lb_ptiles':  # they should all be the same since it's the actual percentile value
                assert len(set(tmpvals)) == 1
                oval = tmpvals[0]
            else:
                oval = numpy.mean(tmpvals)  # average over non-empty iclusts
            outvals[tkey].append(oval)
    if debug:
        print('      --> %s' % '  '.join([fstr%v for v in outvals[print_var]]))
    return outvals

# ----------------------------------------------------------------------------------------
def make_ptile_plot(tmp_ptvals, xvar, plotdir, plotname, xlabel=None, ylabel='?', title=None, fnames=None, true_inf_str='?', n_clusters=None, iclust=None, within_cluster_average=False, xlist=None, use_relative_affy=False, title_str='', n_per_row=4):
    if 'lb_ptiles' not in tmp_ptvals or len(tmp_ptvals['lb_ptiles']) == 0:
        return

    if xvar == 'n-ancestor':  # TODO this probably doesn't really make sense since it gives credit for being off by one in the wrong direction (kind of -- i guess you could argue it's ok to be off either direction, and in any case you're off in the 'right' direction vastly more frequently than in the 'wrong' one so it probably doesn't matter)
        xlist = [abs(v) for v in xlist]

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
    line_args = [
        ((ax.get_xlim(), (xmean, xmean)), {'linewidth' : 3, 'alpha' : 0.7, 'color' : 'darkred', 'linestyle' : '--', 'label' : 'no correlation'}),
        ((tmp_ptvals['lb_ptiles'], tmp_ptvals['perfect_vals']), {'linewidth' : 3, 'alpha' : 0.7, 'color' : 'darkgreen', 'linestyle' : '--', 'label' : 'perfect correlation'}),
    ]
    if xvar == 'daffy':
        line_args.append(((ax.get_xlim(), (0, 0)), {'linewidth' : 1.5, 'alpha' : 0.7, 'color' : 'grey', 'linestyle' : '--', 'label' : 'zero'}))
    if xia:  # so their top/bottom ordering matches the actual lines
        line_args.reverse()
    for (args, kwargs) in line_args:
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
        if xvar == 'daffy':
            ymt = min([xmean] + tmp_ptvals[xkey])
            ymin = (ymt-0.05*abs(ymax-ymt))
        else:
            ymin = -0.02*ymax
        ybounds = (ymin, 1.1*ymax)
        leg_loc = (0.05, 0.62) if xvar=='daffy' else (0.5, 0.6)
        ptile_ylabel = 'mean %s\n%s' % (xlabel, 'from parent' if xvar=='daffy' else 'since affinity increase')

    if use_relative_affy:
        fig.text(0.5, 0.82, 'relative %s' % xvar, fontsize=15, color='red', fontweight='bold')

    if n_clusters is not None:
        if within_cluster_average:
            fig.text(0.25, 0.84, 'within-cluster average over %d families' % n_clusters, fontsize=12, fontweight='bold')  # , color='red'
        else:
            fig.text(0.37, 0.84, 'choosing among %d families' % n_clusters, fontsize=12, fontweight='bold')  # , color='red'
    fn = plotting.mpl_finish(ax, plotdir, plotname, xbounds=ptile_range_tuple(xvar), ybounds=ybounds, leg_loc=leg_loc,
                             title='%s%s%s' % (ungetptlabel(xvar), '' if iclust is None else ', iclust %d'%iclust, title_str),
                             xlabel='%s threshold (percentile)' % ylabel, ylabel=ptile_ylabel, adjust={'left' : 0.21}, legend_fontsize=14)
    add_fn(fnames, fn=fn, n_per_row=n_per_row)

# ----------------------------------------------------------------------------------------
def make_lb_vs_affinity_slice_plots(baseplotdir, lines, lb_metric, is_true_line=False, only_csv=False, fnames=None, separate_rows=False, use_quantile=False, paired=False, n_bin_cfg_fname=None, make_slice_ptile_plots=False, debug=False):
    from . import scanplot
    # make_slice_ptile_plots = True
    # debug = True
    n_bin_cfg = {'default' : 5}
    if n_bin_cfg_fname is not None:
        with open(n_bin_cfg_fname) as cfile:
            n_bin_cfg = yaml.load(cfile, Loader=yaml.CLoader)
    # ----------------------------------------------------------------------------------------
    def modify_lines(slvar, slbounds):
        # ----------------------------------------------------------------------------------------
        def null_affy(line, iseq):
            slval = utils.antnval(line, slvar, iseq) #utils.shm_aa(line, iseq=iseq) if slvar=='shm-aa' else treeutils.smvals(line, slvar, iseq=iseq)
            return slval < slbounds[0] or slval >= slbounds[1]  # include lower bound, exclude upper (for non-null) to match python slice conventions
        # ----------------------------------------------------------------------------------------
        def affy_val(line, iseq):
            if line['affinities'][iseq] is None or null_affy(line, iseq):
                return None
            return line['affinities'][iseq]
        # ----------------------------------------------------------------------------------------
        for antn in lines:
            antn['affinities'] = [affy_val(antn, iseq) for iseq in range(len(antn['unique_ids']))]
    # ----------------------------------------------------------------------------------------
    def slice_var(slvar, bincfg, int_bins):
        # ----------------------------------------------------------------------------------------
        def get_plot_vals(ix):
            slpd, no_plot = ['%s/%s-%s-slice-%d' % (baseplotdir, lb_metric, slvar, ix), False] if make_slice_ptile_plots else [None, True]
            ymfo = plot_lb_vs_affinity(slpd, lines, lb_metric, is_true_line=is_true_line, no_plot=no_plot, title_str=', slice %d'%ix, debug=False)  # this makes the distr hists, which is kind of wasteful, but turning them off is more trouble than it's worth
            if ymfo is None:
                return None, None
            icvals = []
            for iclust in range(len(lines)):
                if all(a is None for a in lines[iclust]['affinities']):
                    imissing[ix].append(iclust)
                    continue
                if len(set(lines[iclust]['affinities']) - set([None])) == 1:  # if all the affinity values are the same we can't (well, shouldn't, and don't any more) make the ptile plots
                    isame[ix].append(iclust)
                    continue
                diff_vals = scanplot.get_ptile_diff_vals(ymfo, iclust=iclust, spec_corr=not use_quantile) #, single_value=True)
                if len(diff_vals) > 0:
                    mdiff = numpy.mean(diff_vals)  # average over top 25% of ptiles (well, default is 25%)
                    icvals.append(mdiff)
                else:
                    imissing[ix].append(iclust)
            # mean: average over clusters (compare to 'for pvkey, ofvals in plotvals.items():' and 'for tau in tvd_keys:' bit in scaplot.make_plots())
            # err: standard error on mean (for standard deviation, comment out denominator)
            return numpy.mean(icvals), float('nan') if len(icvals)==1 else numpy.std(icvals, ddof=1) / math.sqrt(len(icvals))
        # ----------------------------------------------------------------------------------------
        def nonnullcnt():  # NOTE 'affinities' in <lines> are constantly getting modified
            return len([a for l in lines for a in l['affinities'] if a is not None])
        # ----------------------------------------------------------------------------------------
        def print_dbg(xv, mval, err, n_non_now):
            def fstr(a, w=7, j=''): return (utils.wfmt(utils.round_to_n_digits(a, 2), w, fmt='s', jfmt=j) if a is not None else utils.color('blue', '-', width=w))
            avals = sorted(set(a for l in lines for a in l['affinities'] if a is not None))
            def vstr(v, e=False):
                w = 3 if use_quantile else 5
                if v is None:
                    return w * ' '
                nd = 1 if use_quantile else 2
                if e: nd -= 1
                return utils.wfmt(v, w, fmt='.%df'%nd)
            print('    %s   %s  %s  %s %7.2f%7.2f   %s      %4d      %s%s%s     %s' % (utils.color('purple', str(ix), width=2), ('%7d' if int_bins else '%7.4f')%xv, vstr(mval), vstr(err, e=True), xbins[ix], xbins[ix+1],
                                                                                       utils.color('red' if n_non_now==0 else None, str(n_non_now), width=5),  # n non null
                                                                                       len(avals), fstr(None if len(avals)==0 else min(avals)), fstr(None if len(avals)==0 else max(avals)) if len(avals)>1 else fstr(None),
                                                                                       '' if len(affy_val_set)>1 or len(affy_val_set)==0 else utils.color('red', '   all affinities the same'), '  '.join([fstr(a) for a in avals]) if len(avals) < 10 else ''))
        # ----------------------------------------------------------------------------------------
        all_vals = [utils.antnval(l, slvar, i) for l in lines for i in range(len(l['unique_ids']))]
        all_vals = sorted(v for v in all_vals if v is not None)
        if len(set(all_vals)) == 1:
            print('    all %s values the same (%f), so not making slice plots' % (slvar, list(set(all_vals))[0]))
            return
        if isinstance(bincfg, list):  # explicit list of bin low edges
            xbins = [x for x in bincfg]  # maybe we'll modify it, so safer to copy
            n_bins = len(xbins) - 1
            if debug:
                print('    using bins specified in cfg:')
            if xbins != sorted(xbins):
                raise Exception('xbins not sorted: %s' % xbins)
            hutils.auto_bin_expand(all_vals, xbins, int_bins=int_bins, debug=debug)
        else:
            n_bins = bincfg
            # xbins = hutils.autobins(all_vals, n_bins)  # note that this fcn pads a bit to avoid over/underflows
            xbins, n_bins = hutils.auto_volume_bins(all_vals, n_bins, int_bins=int_bins, min_xdist=2, debug=debug)  # NOTE may reset n_bins (mostly will return fewer if int_bins is set)
        assert len(xbins) == n_bins + 1
        if all_vals[0] < xbins[0] or all_vals[-1] > xbins[-1]:
            print('  %s %s values (%.3f or %.3f) outside of lowest/highest x bins %.3f %.3f' % (utils.wrnstr(), slvar, all_vals[0], all_vals[-1], xbins[0], xbins[-1]))
        if debug:
            n_tot, n_non_null_cnt = sum(len(l['affinities']) for l in lines), 0
            initial_n_null = n_tot - nonnullcnt()
            print('   %s with %d bins from %.4f to %.4f (%d/%d initial null affinities)' % (utils.color('blue', slvar), n_bins, xbins[0], xbins[-1], initial_n_null, n_tot))
            hstr = 'qtp' if use_quantile else 's-corr'
            print('    %s   x     %s (+/-)    xmin   xmax  N-not-null N-distinct    min   max' % (utils.color('purple', 'slice'), hstr))
        affy_ranges = []
        imissing, isame = [[] for _ in range(n_bins)], [[] for _ in range(n_bins)]
        slhist = Hist(n_bins=n_bins, xbins=xbins, xmin=xbins[0], xmax=xbins[-1])
        for ix in range(n_bins):  # NOTE len(xbins) is one more than n_bins, for the low edge of the overflow bin
            modify_lines(slvar, (xbins[ix], xbins[ix+1]))
            xv, mval, err = 0.5 * (xbins[ix] + xbins[ix+1]), None, None
            affy_val_set = set(a for l in lines for a in l['affinities']) - set([None])
            affy_ranges.append([min(affy_val_set), max(affy_val_set)] if len(affy_val_set) > 0 else [0, 0])
            if len(affy_val_set) > 1:  # any(a is not None for l in lines for a in l['affinities']) and len(:
                mval, err = get_plot_vals(ix)
            if mval is None:
                slhist.set_ibin(slhist.find_bin(xv), 9999, 0)  # this is a hackey way to get it to not display (since it's out of bounds) NOTE <ix> is *not* the same as ibin
            else:
                slhist.set_ibin(slhist.find_bin(xv), mval, err)  # NOTE <ix> is *not* the same as ibin
            if debug:
                n_non_now = nonnullcnt()
                print_dbg(xv, mval, err, n_non_now)
                n_non_null_cnt += n_non_now
            for antn, o_affies in zip(lines, original_affinities):  # gotta reset them each time
                antn['affinities'] = o_affies
        if debug:
            print('      found %d / %d set to non-null' % (n_non_null_cnt, n_tot))
            assert n_non_null_cnt == n_tot - initial_n_null  # make sure everyone was set to non null exactly once
            assert initial_n_null == n_tot - nonnullcnt()  # and that we properly reset everybody to their initial values
        for ix in range(n_bins):
            if len(imissing[ix]) > 0:
                print('        slice %d: missing %d / %d iclusts: %s' % (ix, len(imissing[ix]), len(lines), ' '.join(str(i) for i in imissing[ix])))
            if len(isame[ix]) > 0:
                print('        slice %d: all affinity values the same for %d / %d iclusts: %s' % (ix, len(isame[ix]), len(lines), ' '.join(str(i) for i in isame[ix])))

        if only_csv:
            return

        fig, ax = plotting.mpl_init()
        slhist.mpl_plot(ax, square_bins=True, errors=True, no_vertical_bin_lines=True, color='#006600', linewidth=6)
        ax2 = ax.twinx()
        for ix, (amin, amax) in enumerate(affy_ranges):
            ax2.fill_between([xbins[ix], xbins[ix+1]], [amin, amin], [amax, amax], alpha=0.2, color='#2b65ec')
        ax2.set_yticks(ax2.get_yticks())  # NOTE for some fucking reason the ticks are misplaced if you don't reset them first like this
        ax2.set_yticklabels([str(utils.round_to_n_digits(a, 2)) for a in ax2.get_yticks()])
        # ax2.set_ylabel('affinity range')  # this doesn't show up, maybe just cause it's pushed off the right edge?
        fig.text(0.84, 0.85, 'affinity\nrange', fontsize=15, fontweight='bold', color='#2b65ec', alpha=0.6)
        plotname = '%s-%s-slice' % (lb_metric, slvar)
        legdict = plotting.legends
        legdict.update(treeutils.legtexts)
        ybounds = (0, 50) if use_quantile else (min(0, slhist.get_minimum() - 0.005), 1)
        ax.set_ylim(ybounds[0], ybounds[1])
        ax.set_ylabel('quantile to perfect' if use_quantile else 'specificity correlation', fontweight='bold', color='#006600', alpha=0.6)
        ax.set_xlabel(legdict.get(slvar, slvar))
        fn = plotting.mpl_finish(ax, baseplotdir + '/slices', plotname, title=mtitlestr('per-seq', lb_metric), right_y_axis=True) #, xlabel=legdict.get(slvar, slvar))  # , ybounds=ybounds ylabel=
        add_fn(fnames, fn=fn)
    # ----------------------------------------------------------------------------------------
    if debug:
        print('  %s' % utils.color('green', lb_metric))
    add_fn(fnames, new_row=True)
    original_affinities = [l['affinities'] for l in lines]
    if len(set(a for alist in original_affinities for a in alist)) == 1:
        print('    all original affinity values the same (%f), so not making slice plots' % list(set(a for alist in original_affinities for a in alist))[0])
        return
    int_vars = ['shm', 'shm-aa']
    slvars = ['affinities'] + int_vars
    if is_true_line and any('min_target_distances' in l for l in lines):
        slvars.append('min_target_distances')
    for slvar in slvars:
        slice_var(slvar, n_bin_cfg.get(slvar, n_bin_cfg['default']), slvar in int_vars)

    add_fn(fnames, new_row=True)

# ----------------------------------------------------------------------------------------
def plot_lb_vs_affinity(baseplotdir, lines, lb_metric, is_true_line=False, use_relative_affy=False, affy_label=None, only_csv=False, no_plot=False, fnames=None, separate_rows=False, add_uids=False,
                        colorvar=None, max_scatter_plot_size=5000, max_iclust_plots=10, title_str='', debug=False):
    # ----------------------------------------------------------------------------------------
    def get_plotvals(line):
        plotvals = {vt : [] for vt in vtypes + ['uids']}
        if colorvar is not None and colorvar == 'is_leaf':
            dtree = treeutils.get_dendro_tree(treestr=get_tree_in_line(line, is_true_line, aa='aa-lb' in lb_metric))  # keeping this here to remind myself how to get the tree if I need it
            if dtree is None:
                print('    %s asked for is_leaf colorvar, but couldn\'t find a tree, so will crash below (probably only specified to calculate cons-dist-aa)'  % utils.wrnstr())
        if 'affinities' not in line:
            if debug: print('  no affy key')
            return plotvals
        for uid, affy in [(u, a) for u, a in zip(line['unique_ids'], line['affinities']) if a is not None and u in line['tree-info']['lb'][lb_metric]]:
            plotvals['affinity'].append(affy)
            if lb_metric in per_seq_metrics:
                plotvals[lb_metric].append(line['tree-info']['lb'][lb_metric][uid])  # NOTE there's lots of entries in the lb info that aren't observed (i.e. aren't in line['unique_ids'])
            if add_uids:
                plotvals['uids'].append(uid)
            if colorvar is not None and colorvar == 'is_leaf':
                node = dtree.find_node_with_taxon_label(uid)
                plotvals['is_leaf'].append(node.is_leaf() if node is not None else None)
        if debug:  # this is a little confusing cause the min/max here look similar to the top-quintile/max whatever from per-cluster stuff, but i might remove the latter at some point
            mstrs = ['', '']
            if len(plotvals['affinity']) > 0:
                mvals = [mfcn(plotvals['affinity']) for mfcn in (min, max)]
                mstrs = ['%6.3f'%v for v in mvals]
                if len(set(mvals)) == 1:
                    mstrs = [utils.color('red', s) for s in mstrs]
            print('    %3d   %s %s' % (len(plotvals['affinity']), mstrs[0], mstrs[1]), end=' ')
        return plotvals
    # ----------------------------------------------------------------------------------------
    def getplotdir(extrastr=''):
        return '%s/%s/%s-vs-%saffinity%s' % (baseplotdir, lb_metric, lb_metric, rel_affy_str('affinity', use_relative_affy=use_relative_affy), extrastr)
    # ----------------------------------------------------------------------------------------
    def icstr(iclust):
        return '-all-clusters' if iclust is None else '-iclust-%d' % iclust
    # ----------------------------------------------------------------------------------------
    def tmpstrs(iclust, vspstuff):
        lbstr, affystr, clstr = lb_metric, 'affinity', icstr(iclust)
        affystr = '%s%s' % (rel_affy_str('affinity', use_relative_affy=use_relative_affy), affystr)
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
    def make_perf_distr_plot(dhists, iclust=None, tfns=None):  # NOTE duplicates a lot of stuff in the vs delta affinity fcn below
        if dhists is None:
            return
        title = '%s on %s tree%s' % (mtitlestr('per-seq', lb_metric, short=True), true_inf_str, (' (%d families together)' % len(lines)) if iclust is None else ' (cluster %d)'%iclust)
        plotname = '%s-vs-%s-%s-tree%s' % (lb_metric, 'affinity', true_inf_str, icstr(iclust))
        tpdir = getplotdir(extrastr='-perf-distr')
        normalize, colors = False, ['#006600', 'darkred']
        if len(dhists) > len(colors):
            normalize = True
            colors = ['#006600', 'royalblue', 'darkorange', 'darkred']
        plotting.draw_no_root(dhists[0], more_hists=dhists[1:], plotdir=tpdir, plotname=plotname, xtitle=mtitlestr('per-seq', lb_metric), plottitle=title, log='y' if iclust is None and 'lb' in lb_metric else '',  # NOTE don't normalize (and if you do, you have to deepcopy them first)
                              errors=True, alphas=[0.7 for _ in range(len(dhists))], colors=colors, linewidths=[5, 3, 2], leg_title='affinity', translegend=(0, -0.1), ytitle='freq.' if normalize else 'counts', normalize=normalize) #, markersizes=[0, 5, 11]) #, linestyles=['-', '-', '-.']) #'']) #, remove_empty_bins=True), '#2b65ec'
        add_fn(tfns, fn='%s/%s.svg'%(tpdir, plotname))

    # ----------------------------------------------------------------------------------------
    if no_plot:
        only_csv = True
    n_rows = 4
    scfns, wptfns, aptfns, dhfns = fnames, fnames, fnames, fnames
    if separate_rows:  # indices of the rows to which to append each type of plot (if not set, all go on the same row [unless row width exceeds the limit, set elsewhere])
        if len(fnames) < n_rows: fnames += [[] for _ in range(n_rows)]
        scfns, wptfns, aptfns, dhfns = [[fnames[i]] for i in range(len(fnames) - n_rows, len(fnames))]  # rows: 1. scatters, 2. within-family ptiles, 3. among-family ptiles 4. perf distributions
    true_inf_str = 'true' if is_true_line else 'inferred'
    vtypes = get_lbscatteraxes(lb_metric)  # NOTE this puts relative affinity under the (plain) affinity key, which is kind of bad maybe i think probably UPDATE nah i think it's better
    if colorvar is not None:
        vtypes.append(colorvar)

    if not only_csv:
        if sum(len(l['unique_ids']) for l in lines) < max_scatter_plot_size:
            make_lb_scatter_plots('affinity', baseplotdir, lb_metric, lines, fnames=scfns, n_iclust_plot_fnames=1 if len(lines)==1 else None, is_true_line=is_true_line, colorvar='edge-dist' if is_true_line else None,
                                  only_overall='among-families' in lb_metric, only_iclust='within-families' in lb_metric or len(lines)==1, add_jitter=is_true_line, use_relative_affy=use_relative_affy, xlabel=affy_label, title_str=title_str) #, add_stats='correlation')  # there's some code duplication between these two fcns, but oh well
        else:  # ok this is hackey
            print('    too many seqs %d >= %d, not making scatter plots' % (sum(len(l['unique_ids']) for l in lines), max_scatter_plot_size))
            utils.prep_dir(getplotdir(), wildlings=['*.svg', '*.yaml'])
    if not no_plot:
        for estr in ['-ptiles', '-perf-distr']:  # previous line does a prep_dir() call as well
            utils.prep_dir(getplotdir(estr), wildlings=['*.svg', '*.yaml'])

    per_seq_plotvals = {vt : [] for vt in vtypes}  # plot values for choosing single seqs/cells (only among all clusters, since the iclust ones don't need to be kept outside the cluster loop)
    per_clust_plotvals = {vt : {sn : [] for sn, _ in cluster_summary_cfg[vt]} for vt in vtypes}  # each cluster plotted as one point using a summary over its cells (e.g. max, mean) for affinity and lb
    pt_vals = {'per-seq' : {}, 'per-cluster' : {}}  # 'per-seq': choosing single cells, 'per-cluster': choosing clusters; with subkeys in the former both for choosing sequences only within each cluster ('iclust-N', used later in cf-tree-metrics.py to average over all clusters in all processes) and for choosing sequences among all clusters together ('all-clusters')
    distr_hists = {}
    # correlation_vals = {'per-seq' : {}, 'per-cluster' : {}}
    if debug:
        print('                     N with    affinity      %s   ' % ''.join([('  %-12s'%vt) for vt in vtypes[:2] for _ in cluster_summary_cfg[vt]]))
        print('      %s   size   affy    min    max    %s' % (utils.color('blue', 'iclust'), ''.join(('  %-12s'%st) for vt in vtypes[:2] for st, _ in cluster_summary_cfg[vt])))
    for iclust, line in enumerate(lines):
        if debug:
            print('      %s    %4d' % (utils.color('blue', str(iclust), width=3), len(line['unique_ids'])), end=' ')
        iclust_plotvals = get_plotvals(line)  # if it's not in <per_seq_metrics> we still need the affinity values
        all_the_same = [k for k in [lb_metric, 'affinity'] if len(set(iclust_plotvals[k])) < 2]
        if lb_metric in per_seq_metrics:
            if len(all_the_same) > 0:
                if debug: print('    all %s values the same' % utils.color('yellow', ', '.join(all_the_same)))
                continue  # this only stops them from being added to the all-clusters-together values
            for vt in vtypes:
                per_seq_plotvals[vt] += iclust_plotvals[vt]
        for vt in vtypes[:2]:
            if len(all_the_same) > 0:
                continue
            if len(iclust_plotvals[vt]) == 0:
                raise Exception('maybe need to add %s to per_seq_metrics?' % lb_metric)
            for sname, sfcn in cluster_summary_cfg[vt]:
                per_clust_plotvals[vt][sname].append(sfcn(line, iclust_plotvals[vt]))
                if debug:
                    print('%12.3f' % per_clust_plotvals[vt][sname][-1], end=' ')
        if debug:
            print('')
        if lb_metric not in per_seq_metrics or 'among-families' in lb_metric or len(all_the_same) > 0:
            continue
        iclust_ptile_vals, iclust_distr_hists = get_ptile_vals(lb_metric, iclust_plotvals, 'affinity', 'affinity', dbgstr='iclust %d'%iclust, use_relative_affy=use_relative_affy, return_distr_hists=True, debug=debug)  # NOTE we don't use the distr hists if this is for slice plots, but it's more trouble than it's worth to turn them off
        pt_vals['per-seq']['iclust-%d'%iclust] = iclust_ptile_vals
        distr_hists['iclust-%d'%iclust] = iclust_distr_hists
        # correlation_vals['per-seq']['iclust-%d'%iclust] = {getcorrkey(*vtypes[:2]) : getcorr(*[iclust_plotvals[vt] for vt in vtypes[:2]])}
        if not only_csv and len(iclust_plotvals['affinity']) > 0 and iclust < max_iclust_plots:
            make_perf_distr_plot(iclust_distr_hists, iclust=iclust, tfns=dhfns if len(lines)==1 else None) #iclust_fnames if iclust < 3 else None)
            make_ptile_plot(iclust_ptile_vals, 'affinity', getplotdir('-ptiles'), ptile_plotname(iclust=iclust),
                            ylabel=tmpylabel(iclust, None), title=mtitlestr('per-seq', lb_metric, short=True), true_inf_str=true_inf_str, iclust=iclust, use_relative_affy=use_relative_affy, title_str=title_str)

    if lb_metric in per_seq_metrics and 'within-families' not in lb_metric:
        if per_seq_plotvals[lb_metric].count(0.) == len(per_seq_plotvals[lb_metric]):
            return
        # correlation_vals['per-seq']['all-clusters'] = {getcorrkey(*vtypes[:2]) : getcorr(*[per_seq_plotvals[vt] for vt in vtypes[:2]])}
        pt_vals['per-seq']['all-clusters'], distr_hists['all-clusters'] = get_ptile_vals(lb_metric, per_seq_plotvals, 'affinity', 'affinity', dbgstr='all clusters together', use_relative_affy=use_relative_affy, return_distr_hists=True, debug=debug)  # choosing single cells from among all cells with all clusters together
        if len(lines) > 1 and not only_csv and len(per_seq_plotvals[lb_metric]) > 0:
            make_perf_distr_plot(distr_hists['all-clusters'], tfns=dhfns)
            make_ptile_plot(pt_vals['per-seq']['all-clusters'], 'affinity', getplotdir('-ptiles'), ptile_plotname(),
                            ylabel=tmpylabel(None, None), title=mtitlestr('per-seq', lb_metric, short=True), fnames=aptfns, true_inf_str=true_inf_str, n_clusters=len(lines), use_relative_affy=use_relative_affy, title_str=title_str)

    if lb_metric in per_seq_metrics and 'among-families' not in lb_metric:
        pt_vals['per-seq']['within-families-mean'] = get_mean_ptile_vals(len(lines), pt_vals['per-seq'], 'affinity')  # within-family mean
        if not only_csv:
            make_ptile_plot(pt_vals['per-seq']['within-families-mean'], 'affinity', getplotdir('-ptiles'), ptile_plotname(iclust=None, extra_str='within-cluster-average'),
                            ylabel=tmpylabel(None, None), title=mtitlestr('per-seq', lb_metric, short=True), fnames=wptfns, true_inf_str=true_inf_str, n_clusters=len(lines), within_cluster_average=True, use_relative_affy=use_relative_affy, title_str=title_str)

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

    yamlfo = {'percentiles' : pt_vals} #, 'correlations' : correlation_vals}
    if not no_plot:
        yamlfo['distr-hists'] = {}
        n_missing = len([hlist for hlist in distr_hists.values() if hlist is None])
        if n_missing > 0:
            print('        %s missing %d / %d h lists for distr hists on %s tree (probably weren\'t any affinity increases in the tree)' % (utils.color('yellow', 'warning'), n_missing, len(distr_hists), 'true' if is_true_line else 'inf'))
        if distr_hists is not None:  # not sure if this is right, but it had hlist here which was definitely wrong
            yamlfo['distr-hists'] = {k : {h.title : h.getdict() for h in hlist} for k, hlist in distr_hists.items()}
        utils.jsdump('%s/%s.yaml' % (getplotdir('-ptiles'), ptile_plotname()), yamlfo)

    return yamlfo

# ----------------------------------------------------------------------------------------
def plot_lb_vs_ancestral_delta_affinity(baseplotdir, lines, lb_metric, is_true_line=False, only_csv=False, fnames=None, separate_rows=False, max_scatter_plot_size=3000, max_iclust_plots=10,
                                        only_look_upwards=False, only_distr_plots=False, n_per_row=4, debug=False):
    # plot lb[ir] vs both number of ancestors and branch length to nearest affinity decrease (well, decrease as you move upwards in the tree/backwards in time)
    # ----------------------------------------------------------------------------------------
    def init_pvals():
        return {vt : [] for vt in [lb_metric, xvar]} #, 'daffy']}  # , 'uids']}
    # ----------------------------------------------------------------------------------------
    def get_plotvals(line, xvar, iclust):
        plotvals = init_pvals()
        dtree = treeutils.get_dendro_tree(treestr=get_tree_in_line(line, is_true_line))  # i don't think there's any reason for this to use the aa tree? , aa='aa-lb' in lb_metric))
        affy_increasing_edges, all_affy_changes = treeutils.find_affy_increases(dtree, line)
        if debug and iclust == 0:
            if debug > 1:
                print(utils.pad_lines(treeutils.get_ascii_tree(dendro_tree=dtree, width=250)))
            def e_affy(e): return utils.per_seq_val(line, 'affinities', e.head_node.taxon.label) - utils.per_seq_val(line, 'affinities', e.tail_node.taxon.label)
            print('    %d edges with affinity increases: %s' % (len(affy_increasing_edges), ' '.join('%.4f'%e_affy(e) for e in affy_increasing_edges)))
            print('                     and child nodes: %s' % ' '.join(e.head_node.taxon.label for e in affy_increasing_edges))
            print('                          ancestors                           chosen edge')
            print('         node    %7s   (desc.)  steps distance  affinity  affy change    (%s: reached root/leaves without finding lower-affinity ancestor)' % (lb_metric, utils.color('yellow', '?')))
        for uid in line['unique_ids']:
            node = dtree.find_node_with_taxon_label(uid)
            if node is dtree.seed_node:  # root doesn't have any ancestors
                continue
            if uid not in line['tree-info']['lb'][lb_metric]:
                continue
            lbval = line['tree-info']['lb'][lb_metric][uid]  # NOTE there's lots of entries in the lb info that aren't observed (i.e. aren't in line['unique_ids'])
            if ('lbr' in lb_metric or 'lbf' in lb_metric) and lbval == 0:  # lbr equals 0 should really be treated as None/missing
                continue
            n_steps, branch_len = treeutils.get_min_steps_to_affy_increase(affy_increasing_edges, node, dtree, line, also_return_branch_len=True, lbval=line['tree-info']['lb'][lb_metric][uid], only_look_upwards=only_look_upwards, debug=debug)
            if n_steps is None:
                continue
            if all_affy_changes.get(uid) is None:
                continue
            if xvar == 'n-ancestor':
                plotvals[xvar].append(n_steps)
            elif xvar == 'branch-length':
                plotvals[xvar].append(branch_len)
            elif xvar == 'daffy':
                plotvals[xvar].append(all_affy_changes[uid])
            else:
                assert False
            plotvals[lb_metric].append(lbval)
            # plotvals['uids'].append(uid)
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
        xlstr = 'affinity change' if xvar=='daffy' else '%s to nearest affinity increase' % xlabel
        fn = plot_2d_scatter('%s-vs-%s-%s-tree%s' % (lb_metric, xvar, true_inf_str, icstr(iclust)), getplotdir(xvar), plotvals, lb_metric, mtitlestr('per-seq', lb_metric), title, xvar=xvar, xlabel=xlstr, log='' if 'lbr' in lb_metric else '', stats='correlation' if xvar=='daffy' else None)
        add_fn(tfns, fn=fn, n_per_row=n_per_row)
    # ----------------------------------------------------------------------------------------
    def make_perf_distr_plot(dhists, xvar, iclust=None, tfns=None):
        if dhists is None:
            return
        title = '%s on %s tree%s' % (mtitlestr('per-seq', lb_metric, short=True), true_inf_str, (' (%d families together)' % len(lines)) if iclust is None else ' (cluster %d)'%iclust)
        plotname = '%s-vs-%s-%s-tree%s' % (lb_metric, xvar, true_inf_str, icstr(iclust))
        tpdir = getplotdir(xvar, extrastr='-perf-distr')
        normalize, colors = False, ['#006600', 'darkred']
        if xvar == 'n-ancestor':
            leg_title = 'N steps to\naff. increase'
        else:
            leg_title = 'affinity'
            colors = ['grey', 'red', 'blue']
        if len(dhists) > len(colors):
            colors = ['#006600', 'royalblue', 'darkorange', 'darkred']
        plotting.draw_no_root(dhists[0], more_hists=dhists[1:], plotdir=tpdir, plotname=plotname, xtitle=mtitlestr('per-seq', lb_metric), plottitle=title, log='y' if iclust is None else '',  # NOTE don't normalize (and if you do, you have to deepcopy them first)
                              errors=True, alphas=[0.7 for _ in range(len(dhists))], colors=colors, linewidths=[5, 3, 2], leg_title=leg_title, translegend=(0, -0.1), ytitle='freq.' if normalize else 'counts', normalize=normalize) #, markersizes=[0, 5, 11]) #, linestyles=['-', '-', '-.']) #'']) #, remove_empty_bins=True), '#2b65ec'
        add_fn(tfns, fn='%s/%s.svg'%(tpdir, plotname), n_per_row=n_per_row)
    # ----------------------------------------------------------------------------------------
    def get_distr_hists(plotvals, xvar, max_bin_width=1., min_bins=30, extra_hists=False, iclust=None):
        # ----------------------------------------------------------------------------------------
        def gethist(tstr, vfcn):  # NOTE similarity to fcn above
            vlist = [v for v, n in zip(plotvals[lb_metric], plotvals[xvar]) if vfcn(n)]
            return Hist(n_bins=n_bins, xmin=xmin, xmax=xmax, title=tstr, value_list=vlist)
        # ----------------------------------------------------------------------------------------
        if len(plotvals[xvar]) == 0:
            return None
        # xmin, xmax = [tfac * mfcn(plotvals[lb_metric]) for mfcn, tfac in zip((min, max), (0.9, 1.1))]
        xmin, xmax = 0., 1.01 * max(plotvals[lb_metric])  # this is the low edge of the overflow bin, so needs to be a bit above the biggest value
        n_bins = max(min_bins, int((xmax - xmin) / float(max_bin_width)))  # for super large lbr values like 100 you need way more bins
        if xvar == 'n-ancestor':
            zero_hist = gethist('zero', lambda n: n == 0)  # lb values for nodes that are immediately below affy-increasing branch
            # <0 is a bit further right than >0, and abs(v)==1 is a bit further right than >1, but these differences are all small compared to the difference to ==0, which is much further right (this is in quick[ish] tests, not doing full parameter scans)
            if extra_hists:
                m_one_hist = gethist('-1', lambda n: n == -1)
                p_one_hist = gethist('+1', lambda n: n == 1)
                g_one_hist = gethist('|n|>1', lambda n: abs(n) > 1)
                dhists = [g_one_hist, m_one_hist, p_one_hist, zero_hist]
            else:
                other_hist = gethist('not 0', lambda n: abs(n) > 0)  # not perfect
                dhists = [other_hist, zero_hist]
        else:
            const_hist = gethist('constant', lambda a: a == 0)
            incr_hist = gethist('increase', lambda a: a > 0)
            decr_hist = gethist('decrease', lambda a: a < 0)
            # const_hist = gethist('|da|<=0.3', lambda a: abs(a) <= 0.3)
            # incr_hist = gethist('da>0.3', lambda a: a > 0.3)
            # decr_hist = gethist('da<-0.3', lambda a: a < -0.3)
            dhists = [const_hist, decr_hist, incr_hist]
        return dhists
    # ----------------------------------------------------------------------------------------
    def ptile_plotname(xvar, iclust, extra_str=None):
        return '%s-vs-%s-%s-tree-ptiles%s%s' % (lb_metric, xvar, true_inf_str, icstr(iclust), '' if extra_str is None else '-'+extra_str)

    # ----------------------------------------------------------------------------------------
    n_rows = 5
    dhfns, scfns, nsfns, wptfns, aptfns = fnames, fnames, fnames, fnames, fnames
    if separate_rows:  # indices of the rows to which to append each type of plot (if not set, all go on the same row [unless row width exceeds the limit, set elsewhere])
        if len(fnames) < n_rows: fnames += [[] for _ in range(n_rows)]
        dhfns, scfns, nsfns, wptfns, aptfns = [[fnames[i]] for i in range(len(fnames) - n_rows, len(fnames))]
    true_inf_str = 'true' if is_true_line else 'inferred'
    xvar_list = collections.OrderedDict([(xvar, xlabel) for metric, cfglist in lb_metric_axis_cfg('lbr') for xvar, xlabel in cfglist])
    for xvar, estr in itertools.product(xvar_list, ['', '-ptiles', '-perf-distr']):
        utils.prep_dir(getplotdir(xvar, extrastr=estr), wildlings=['*.svg', '*.yaml'])
    if debug:
        print('%s finding ancestors with most recent affinity increases' % utils.color('blue', lb_metric))
    iclust_fnames = add_fn(None, init=True)
    for xvar, xlabel in xvar_list.items():
        per_seq_plotvals = init_pvals()
        pt_vals = {'per-seq' : {}, 'per-cluster' : {}}  # 'per-seq': choosing single cells, 'per-cluster': choosing clusters; with subkeys in the former both for choosing sequences only within each cluster ('iclust-N', used later in cf-tree-metrics.py to average over all clusters in all processes) and for choosing sequences among all clusters together ('all-clusters')
        distr_hists = {}
        for iclust, line in enumerate(lines):
            if debug:
                if iclust == 0: print(' %s' % utils.color('green', xvar))
                print('  %s' % utils.color('blue', 'iclust %d' % iclust))
            iclust_plotvals = get_plotvals(line, xvar, iclust)
            for vtype in per_seq_plotvals:
                per_seq_plotvals[vtype] += iclust_plotvals[vtype]
            if 'among-families' in lb_metric:
                continue
            iclust_ptile_vals = get_ptile_vals(lb_metric, iclust_plotvals, xvar, xlabel, dbgstr='iclust %d'%iclust, debug=debug)
            pt_vals['per-seq']['iclust-%d'%iclust] = iclust_ptile_vals
            distr_hists['iclust-%d'%iclust] = get_distr_hists(iclust_plotvals, xvar, iclust=iclust)
            if not only_csv and iclust < max_iclust_plots:
                make_perf_distr_plot(distr_hists['iclust-%d'%iclust], xvar, iclust=iclust, tfns=dhfns if len(lines)==1 else None) #iclust_fnames if iclust < 3 else None)
                if not only_distr_plots:
                    make_scatter_plot(iclust_plotvals, xvar, iclust=iclust, tfns=nsfns if len(lines)==1 else None) #iclust_fnames if iclust < 3 else None)
                    make_ptile_plot(iclust_ptile_vals, xvar, getplotdir(xvar, extrastr='-ptiles'), ptile_plotname(xvar, iclust), xlist=iclust_plotvals[xvar],
                                    xlabel=xlabel, ylabel=mtitlestr('per-seq', lb_metric), true_inf_str=true_inf_str, iclust=iclust, fnames=wptfns if len(lines)==1 else None, n_per_row=n_per_row)  #iclust_fnames if len(lines)==1 else None)

        if 'among-families' not in lb_metric and not only_csv:
            pt_vals['per-seq']['within-families-mean'] = get_mean_ptile_vals(len(lines), pt_vals['per-seq'], xvar)
            if len(lines) > 1 and not only_distr_plots:
                make_ptile_plot(pt_vals['per-seq']['within-families-mean'], xvar, getplotdir(xvar, extrastr='-ptiles'), ptile_plotname(xvar, None, extra_str='within-cluster-average'),
                                xlabel=xlabel, ylabel=mtitlestr('per-seq', lb_metric), fnames=wptfns, true_inf_str=true_inf_str, n_clusters=len(lines), within_cluster_average=True, xlist=per_seq_plotvals[xvar], n_per_row=n_per_row)

        if 'within-families' not in lb_metric:
            pt_vals['per-seq']['all-clusters'] = get_ptile_vals(lb_metric, per_seq_plotvals, xvar, xlabel, dbgstr='all clusters', debug=debug)  # "averaged" might be a better name than "all", but that's longer
            distr_hists['all-clusters'] = get_distr_hists(per_seq_plotvals, xvar)
            if len(lines) > 1 and not only_csv:
                make_perf_distr_plot(distr_hists['all-clusters'], xvar, tfns=dhfns)
                if not only_distr_plots:
                    make_scatter_plot(per_seq_plotvals, xvar, tfns=nsfns)
                    make_ptile_plot(pt_vals['per-seq']['all-clusters'], xvar, getplotdir(xvar, extrastr='-ptiles'), ptile_plotname(xvar, None), xlist=per_seq_plotvals[xvar],
                                    xlabel=xlabel, ylabel=mtitlestr('per-seq', lb_metric), fnames=aptfns, true_inf_str=true_inf_str, n_clusters=len(lines), n_per_row=n_per_row)

        n_missing = len([hlist for hlist in distr_hists.values() if hlist is None])
        if n_missing > 0:
            print('        %s missing %d / %d h lists for distr hists on %s tree (probably weren\'t any affinity increases in the tree)' % (utils.color('yellow', 'warning'), n_missing, len(distr_hists), 'true' if is_true_line else 'inf'))
        yamlfo = {'percentiles' : pt_vals, 'distr-hists' : {k : {h.title : h.getdict() for h in hlist} for k, hlist in distr_hists.items() if hlist is not None}}
        utils.jsdump('%s/%s.yaml' % (getplotdir(xvar, extrastr='-ptiles'), ptile_plotname(xvar, None)), yamlfo)  # not adding the new correlation keys atm (like in the lb vs affinity fcn)

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
def get_lb_tree_cmd(treestr, outfname, lb_metric, affy_key, subworkdir, metafo=None, tree_style=None, queries_to_include=None, label_all_nodes=False, label_leaf_nodes=False, label_root_node=False, seq_len=None,
                    meta_info_key_to_color=None, meta_info_to_emphasize=None, node_size_key=None, branch_color_key=None, uid_translations=None, node_label_regex=None, min_face_size=None, max_face_size=None):
    treefname = '%s/tree.nwk' % subworkdir
    metafname = '%s/meta.yaml' % subworkdir
    if not os.path.exists(subworkdir):
        os.makedirs(subworkdir)
    with open(treefname, utils.csv_wmode()) as treefile:
        treefile.write(treestr)
    cmdstr = '%s/bin/plot-lb-tree.py --use-node-area --treefname %s' % (utils.get_partis_dir(), treefname)
    if metafo is not None:
        utils.jsdump(metafname, metafo) # had this here, but it was having a crash from a quoting bug (was using " when it should have used '): yaml.dump , Dumper=Dumper)
        cmdstr += ' --metafname %s' % metafname
    if queries_to_include is not None and len(queries_to_include) > 0:
        cmdstr += ' --queries-to-include %s' % ':'.join(queries_to_include)
    if uid_translations is not None and len(uid_translations) > 0:
        tnodes = [n.label for n in treeutils.get_dendro_tree(treestr=treestr)]
        if any(u not in tnodes for u, _ in uid_translations):
            print('  %s %d uid translations missing from tree: %s' % (utils.wrnstr(), len([u for u, _ in uid_translations if u not in tnodes]), ' '.join(u for u, _ in uid_translations if u not in tnodes)))
        cmdstr += ' --uid-translations %s' % ':'.join('%s,%s'%(u, au) for u, au in uid_translations)
    if label_all_nodes:
        cmdstr += ' --label-all-nodes'
    if label_leaf_nodes:
        cmdstr += ' --label-leaf-nodes'
    if label_root_node:
        cmdstr += ' --label-root-node'
    cmdstr += ' --outfname %s' % outfname
    if lb_metric is not None:
        cmdstr += ' --lb-metric %s' % lb_metric
        if 'lbr' in lb_metric:
            cmdstr += ' --log-lbr'
    if affy_key is not None:
        cmdstr += ' --affy-key %s' % utils.reversed_input_metafile_keys[affy_key]
    # cmdstr += ' --lb-tau %f' % lb_tau
    if tree_style is not None:
        cmdstr += ' --tree-style %s' % tree_style
    if meta_info_key_to_color is not None:
        cmdstr += ' --meta-info-key-to-color %s' % meta_info_key_to_color
    if meta_info_to_emphasize is not None:
        cmdstr += ' --meta-info-to-emphasize %s' % ','.join(list(meta_info_to_emphasize.items())[0])
    if node_size_key is not None:
        cmdstr += ' --node-size-key %s' % node_size_key
    if branch_color_key is not None:
        cmdstr += ' --branch-color-key %s' % branch_color_key
    if node_label_regex is not None:
        cmdstr += ' --node-label-regex %s' % node_label_regex
    if min_face_size is not None:
        cmdstr += ' --min-face-size %d' % min_face_size
    if max_face_size is not None:
        cmdstr += ' --max-face-size %d' % max_face_size
    cmdstr = utils.run_ete_script(cmdstr, return_for_cmdfos=True, extra_str='        ')

    return {'cmd_str' : cmdstr, 'workdir' : subworkdir, 'outfname' : outfname, 'workfnames' : [treefname, metafname]}

# ----------------------------------------------------------------------------------------
def plot_lb_trees(args, metric_methods, baseplotdir, lines, base_workdir, is_true_line=False, tree_style=None, fnames=None):
    add_fn(fnames, new_row=True)
    workdir = '%s/ete3-plots' % base_workdir
    plotdir = baseplotdir + '/trees'
    utils.prep_dir(plotdir, wildlings='*.svg')

    if not os.path.exists(workdir):
        os.makedirs(workdir)
    cmdfos = []
    for lb_metric in metric_methods:
        fnames += [['header', '%s (%s)'%(lb_metric, treeutils.default_inference_str if args.tree_inference_method is None else args.tree_inference_method)], []]  # header for html file
        for iclust, line in enumerate(lines):  # note that <min_selection_metric_cluster_size> was already applied in treeutils
            treestr = get_tree_in_line(line, is_true_line, aa='aa-lb' in lb_metric)
            if treestr is None:
                print('    None type tree')
                continue
            qtis = None if args.queries_to_include is None or args.dont_label_queries_to_include else [q for q in args.queries_to_include if q in line['unique_ids']]  # NOTE make sure to *not* modify args.queries_to_include
            altids = [(u, au) for u, au in zip(line['unique_ids'], line['alternate-uids']) if au is not None] if 'alternate-uids' in line else None
            affy_key = 'affinities'  # turning off possibility of using relative affinity for now
            metafo = copy.deepcopy(line['tree-info']['lb'])  # NOTE there's lots of entries in the lb info that aren't observed (i.e. aren't in line['unique_ids'])
            if self.args.infer_trees_with_collapsed_duplicate_seqs:
                print('    %s --infer-trees-with-collapsed-duplicate-seqs was set, but lb tree plotting doesn\'t yet handle multiplicity/collapse info (see how partitionplotter does it' % utils.wrnstr())
            if affy_key in line:  # either 'affinities' or 'relative_affinities'
                metafo[utils.reversed_input_metafile_keys[affy_key]] = {uid : affy for uid, affy in zip(line['unique_ids'], line[affy_key])}
            outfname = '%s/%s-tree-iclust-%d%s.svg' % (plotdir, lb_metric, iclust, '-relative' if 'relative' in affy_key else '')
            cmdfos += [get_lb_tree_cmd(treestr, outfname, lb_metric, affy_key, '%s/sub-%d' % (workdir, len(cmdfos)), metafo=metafo, tree_style=tree_style, queries_to_include=qtis,
                                       label_all_nodes=args.label_tree_nodes, label_leaf_nodes=args.label_leaf_nodes, label_root_node=args.label_root_node, uid_translations=altids, node_label_regex=args.node_label_regex,
                                       seq_len=float(numpy.mean([len(s) for s in line['seqs']])))]
            add_fn(fnames, fn=outfname, n_per_row=4)

    if len(cmdfos) > 0:
        start = time.time()
        utils.run_cmds(cmdfos, clean_on_success=True, shell=True, n_max_procs=utils.auto_n_procs(), proc_limit_str='plot-lb-tree.py')  # I'm not sure what the max number of procs is, but with 21 it's crashing with some of them not able to connect to the X server, and I don't see a big benefit to running them all at once anyways
        print('    made %d ete tree plots (%.1fs)' % (len(cmdfos), time.time() - start))

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
                print('    %s parent \'%s\' not in true affinities, skipping lonr values' % (utils.color('red', 'warning'), lfo['parent']))
                continue
            if lfo['child'] not in true_affinities:
                print('    %s child \'%s\' not in true affinities, skipping lonr values' % (utils.color('red', 'warning'), lfo['child']))
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
