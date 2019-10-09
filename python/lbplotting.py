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

import utils
from hist import Hist
import treeutils
import plotting

# ----------------------------------------------------------------------------------------
# this name is terrible, but it's complicated and I can't think of a better one
lb_metric_axis_stuff = [  # x axis variables against which we plot each lb metric (well, they're the x axis on the scatter plots, not the ptile plots)
    ('lbi', 'affinity', 'affinity'),
    ('lbr', 'n-ancestor', 'N ancestors'),
    ('lbr', 'branch-length', 'branch length')
]

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
        hists, labels = [], []
        for ihist in range(len(numpyhists)):
            nphist = numpyhists[ihist]  # numpy.hist is two arrays: [0] is bin counts, [1] is bin x values (not sure if low, high, or centers)
            obs_time = ihist  #  + 1  # I *think* it's right without the 1 (although I guess it's really a little arbitrary)
            if not plot_this_time(obs_time, numpyhists):
                continue
            if nphist is None:  # time points at which we didn't sample
                hists.append(None)
                labels.append(None)
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
            labels.append('%d (%.0f, %.1f)' % (obs_time, hist.integral(include_overflows=True), hist.get_mean()))

        hists, labels = zip(*[(h, l) for h, l in zip(hists, labels) if h is not None])  # remove None time points
        return hists, labels, xmin, xmax

    # ----------------------------------------------------------------------------------------
    all_hists, all_labels, xmin, xmax = get_hists(histfname)
    if sum(h.integral(include_overflows=True) for h in all_hists) == 0:
        print '  %s no/empty hists in %s' % (utils.color('yellow', 'warning'), histfname)
        return
    jpdata = []
    for hist in all_hists:
        jpdata.append([x for x, y in zip(hist.get_bin_centers(), hist.bin_contents) for _ in range(int(y)) if x > xmin and x < xmax])  # NOTE this is repeating the 'for _ in range()' in the fcn above, but that's because I used to be actually using the Hist()s, and maybe I will again

    pre_fig, pre_ax = plotting.mpl_init()  # not sure to what extent these really get used after joypy is done with things
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # i don't know why it has to warn me that it's clearing the fig/ax I'm passing in, and I don't know how else to stop it
        fig, axes = joypy.joyplot(jpdata, labels=all_labels, fade=True, hist=True, overlap=0.5, ax=pre_ax, x_range=(xmin, xmax), bins=int(xmax - xmin), xlabelsize=15, ylabelsize=15)
    fig.text(0.01, 0.9, 'generation', fontsize=15)
    fig.text(0.01, 0.85, '(N cells, mean)', fontsize=15)
    # NOTE do *not* set your own x ticks/labels in the next line, since they'll be in the wrong place (i.e. not the same as where joypy puts them)
    plotting.mpl_finish(pre_ax, plotdir, plotname, title=title, xlabel=xlabel, ylabel='generation') #, leg_loc=(0.7, 0.45)) #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

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
    from Bio.Seq import Seq

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
def plot_bcr_phylo_simulation(outdir, event, extrastr, metric_for_target_distance):
    utils.prep_dir(outdir + '/plots', wildlings=['*.csv', '*.svg'])

    plot_bcr_phylo_kd_vals(outdir + '/plots', event)
    plot_bcr_phylo_target_attraction(outdir + '/plots', event)

    plot_bcr_phylo_selection_hists('%s/%s_min_aa_target_hdists.p' % (outdir, extrastr), outdir + '/plots', 'min-aa-target-all-cells', title='all cells', xlabel='%s distance to nearest target sequence' % metric_for_target_distance)
    plot_bcr_phylo_selection_hists('%s/%s_sampled_min_aa_target_hdists.p' % (outdir, extrastr), outdir + '/plots', 'min-aa-target-sampled-cells', plot_all=True, title='sampled cells (excluding ancestor sampling)', xlabel='%s distance to nearest target sequence' % metric_for_target_distance)
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
def plot_lb_vs_shm(baseplotdir, lines_to_use, fnames=None, is_true_line=False, add_uids=False, n_per_row=4):  # <is_true_line> is there because we want the true and inferred lines to keep their trees in different places, because the true line just has the one, true, tree, while the inferred line could have a number of them (yes, this means I maybe should have called it the 'true-tree' or something)
    if fnames is None:
        fnames = []
    fnames.append([])
    sorted_lines = sorted([l for l in lines_to_use if get_tree_from_line(l, is_true_line) is not None], key=lambda l: len(l['unique_ids']), reverse=True)

    # note: all clusters together
    subfnames = {lb_metric : [] for lb_metric in treeutils.lb_metrics}
    for lb_metric, lb_label in treeutils.lb_metrics.items():
        plotdir = '%s/%s/%s-vs-shm' % (baseplotdir, lb_metric, lb_metric)
        utils.prep_dir(plotdir, wildlings='*.svg')
        plotvals = {x : {'leaf' : [], 'internal' : []} for x in ['shm', lb_metric]}
        basetitle = '%s %s vs SHM' % ('true' if is_true_line else 'inferred', lb_metric.upper())
        for iclust, line in enumerate(sorted_lines):  # get depth/n_mutations for each node
            iclust_plotvals = {x : {'leaf' : [], 'internal' : []} for x in ['shm', lb_metric, 'uids']}
            dtree = treeutils.get_dendro_tree(treestr=get_tree_from_line(line, is_true_line))
            n_max_mutes = max(line['n_mutations'])  # don't generally have n mutations for internal nodes, so use this to rescale the depth in the tree
            max_depth = max(n.distance_from_root() for n in dtree.leaf_node_iter())
            for node in dtree.preorder_node_iter():
                if lb_metric == 'lbr' and line['tree-info']['lb'][lb_metric][node.taxon.label] == 0:  # lbr equals 0 should really be treated as None/missing
                    continue
                iseq = line['unique_ids'].index(node.taxon.label) if node.taxon.label in line['unique_ids'] else None
                if node.taxon.label not in line['unique_ids']:  # TODO not really sure whether I want to skip them or not (it's nice to have the visual confirmation that I'm not observing any internal nodes, but then again it's also nice to see where the internal nodes fall in the plot)
                    continue
                n_muted = line['n_mutations'][iseq] if node.taxon.label in line['unique_ids'] else node.distance_from_root() * n_max_mutes / float(max_depth)
                tkey = 'leaf' if node.is_leaf() else 'internal'
                iclust_plotvals['shm'][tkey].append(n_muted)
                iclust_plotvals[lb_metric][tkey].append(line['tree-info']['lb'][lb_metric][node.taxon.label])
                affyval = line['affinities'][iseq] if 'affinities' in line and iseq is not None else None
                if add_uids:
                    iclust_plotvals['uids'][tkey].append(node.taxon.label if affyval is not None else None)
            title = '%s (%d observed, %d total)' % (basetitle, len(line['unique_ids']), len(line['tree-info']['lb'][lb_metric]))
            fn = plot_2d_scatter('%s-vs-shm-iclust-%d' % (lb_metric, iclust), plotdir, iclust_plotvals, lb_metric, lb_label, title, xvar='shm', xlabel='N mutations', leg_loc=(0.7, 0.75), log='y' if lb_metric == 'lbr' else '')
            if iclust < n_per_row:  # i.e. only put one row's worth in the html
                subfnames[lb_metric].append(fn)
            for vtype in [vt for vt in plotvals if vt != 'uids']:
                for ltype in plotvals[vtype]:
                    plotvals[vtype][ltype] += iclust_plotvals[vtype][ltype]
        fn = plot_2d_scatter('%s-vs-shm-all-clusters' % lb_metric, plotdir, plotvals, lb_metric, lb_label, '%s (all clusters)' % basetitle, xvar='shm', xlabel='N mutations', leg_loc=(0.7, 0.75), log='y' if lb_metric == 'lbr' else '')
        fnames[-1].append(fn)

    # TODO can't be bothered figuring out how to get these to work with the distributions (below) a.t.m.
    # fnames.append([fn for lbm in treeutils.lb_metrics for fn in subfnames[lbm]])

# ----------------------------------------------------------------------------------------
def plot_lb_distributions(baseplotdir, lines_to_use, fnames=None, plot_str='', n_per_row=4):
    def make_hist(plotvals, n_total, n_skipped, iclust=None):
        if len(plotvals) == 0:
            return
        hist = Hist(30, 0., max(plotvals), value_list=plotvals)
        fig, ax = plotting.mpl_init()
        hist.mpl_plot(ax) #, square_bins=True, errors=False)
        plotname = '%s-%s' % (lb_metric, str(iclust) if iclust is not None else 'all-clusters')
        leafskipstr = ', skipped %d leaves' % n_skipped if n_skipped > 0 else ''  # ok they're not necessarily leaves, but almost all of them are leaves (and not really sure how a non-leaf could get zero, but some of them seem to)
        fn = plotting.mpl_finish(ax, plotdir, plotname, xlabel=lb_label, log='y', ylabel='counts', title='%s %s  (size %d%s)' % (plot_str, lb_metric.upper(), n_total, leafskipstr))
        if iclust is None:
            fnames[-1].append(fn)
        elif iclust < n_per_row:  # i.e. only put one row's worth in the html
            tmpfnames.append(fn)

    sorted_lines = sorted([l for l in lines_to_use if 'tree-info' in l], key=lambda l: len(l['unique_ids']), reverse=True)  # if 'tree-info' is missing, it should be because it's a small cluster we skipped when calculating lb values
    if fnames is None:  # no real effect (except not crashing) since we're not returning it any more
        fnames = []
    if len(fnames) < 1:  # shouldn't really happen
        fnames.append([])
    tmpfnames = []

    for lb_metric, lb_label in treeutils.lb_metrics.items():
        plotvals = []
        n_total_skipped_leaves = 0
        plotdir = '%s/%s/distributions' % (baseplotdir, lb_metric)
        utils.prep_dir(plotdir, wildlings=['*.svg'])
        for iclust, line in enumerate(sorted_lines):
            iclust_plotvals = line['tree-info']['lb'][lb_metric].values()
            cluster_size = len(iclust_plotvals)  # i.e. including leaves
            if lb_metric == 'lbr':
                iclust_plotvals = [v for v in iclust_plotvals if v > 0.]  # don't plot the leaf values, they just make the plot unreadable
            make_hist(iclust_plotvals, cluster_size, cluster_size - len(iclust_plotvals), iclust=iclust)
            plotvals += iclust_plotvals
            n_total_skipped_leaves += cluster_size - len(iclust_plotvals)
        make_hist(plotvals, len(plotvals) + n_total_skipped_leaves, n_total_skipped_leaves)

    # TODO can't be bothered to get this to work with the _vs_shm (above) a.t.m.
    # fnames.append(tmpfnames)

# ----------------------------------------------------------------------------------------
def plot_2d_scatter(plotname, plotdir, plotvals, yvar, ylabel, title, xvar='affinity', xlabel='affinity', log='', leg_loc=None, warn_text=None):
    def getall(k):
        if 'leaf' in plotvals[xvar]:
            return [v for tk in plotvals[k] for v in plotvals[k][tk]]
        else:
            return plotvals[k]

    if len(getall(xvar)) == 0:
        # print '    no %s vs affy info' % yvar
        return '%s/%s.svg' % (plotdir, plotname)
    fig, ax = plotting.mpl_init()
    # cmap, norm = get_normalized_cmap_and_norm()
    # ax.hexbin(plotvals[xvar], plotvals[yvar], gridsize=15, cmap=plt.cm.Blues)
    if 'leaf' not in plotvals[xvar]:  # single plot
        ax.scatter(plotvals[xvar], plotvals[yvar], alpha=0.4)
    else:  # separate plots for leaf/internal nodes
        for tkey, color in zip(plotvals[xvar], (None, 'darkgreen')):
            ax.scatter(plotvals[xvar][tkey], plotvals[yvar][tkey], label=tkey, alpha=0.4, color=color)
    if 'uids' in plotvals:
        for xval, yval, uid in zip(getall(xvar), getall(yvar), getall('uids')):  # note: two ways to signal not to do this: sometimes we have 'uids' in the dict, but don't fill it (so the zip() gives an empty list), but sometimes we populate 'uids' with None values
            if uid is None:
                continue
            ax.plot([xval], [yval], color='red', marker='.', markersize=10)
            ax.text(xval, yval, uid, color='red', fontsize=8)

    if warn_text is not None:
        ax.text(0.6 * ax.get_xlim()[1], 0.75 * ax.get_ylim()[1], warn_text, fontsize=30, fontweight='bold', color='red')
    xmin, xmax = min(getall(xvar)), max(getall(xvar))
    ymin, ymax = min(getall(yvar)), max(getall(yvar))
    xbounds = xmin - 0.02 * (xmax - xmin), xmax + 0.02 * (xmax - xmin)
    if 'y' in log:
        ybounds = 0.75 * ymin, 1.3 * ymax
    else:
        ybounds = ymin - 0.03 * (ymax - ymin), ymax + 0.08 * (ymax - ymin)
    fn = plotting.mpl_finish(ax, plotdir, plotname, title=title, xlabel=xlabel, ylabel=ylabel, xbounds=xbounds, ybounds=ybounds, log=log, leg_loc=leg_loc)
    return fn

# ----------------------------------------------------------------------------------------
def plot_lb_vs_affinity(plot_str, baseplotdir, lines, lb_metric, lb_label, ptile_range_tuple=(50., 100., 1.), is_true_line=False, n_per_row=4, affy_key='affinities', only_csv=False, fnames=None, add_uids=False, debug=False):
    # ----------------------------------------------------------------------------------------
    def get_plotvals(line):
        plotvals = {vt : [] for vt in vtypes + ['uids']}
        # dtree = treeutils.get_dendro_tree(treestr=get_tree_from_line(line, is_true_line))  # keeping this here to remind myself how to get the tree if I need it
        if affy_key not in line:
            return plotvals
        for uid, affy in [(u, a) for u, a in zip(line['unique_ids'], line[affy_key]) if a is not None]:
            plotvals['affinity'].append(affy)
            plotvals[lb_metric].append(line['tree-info']['lb'][lb_metric][uid])
            if add_uids:
                plotvals['uids'].append(uid)
        return plotvals
    # ----------------------------------------------------------------------------------------
    def get_ptile_vals(plotvals):
        ptile_vals = {'lb_ptiles' : [], 'mean_affy_ptiles' : [], 'perfect_vals' : []}  # , 'reshuffled_vals' : []}
        if len(plotvals[lb_metric]) == 0:
            # print '  no affinity values when trying to make lb vs affinity plots'
            return ptile_vals
        if debug:
            print '            %3s         N     mean    mean       |  perfect   perfect' % lb_metric
            print '    ptile  threshold  taken   %saffy    affy ptile |  N taken  mean ptile'  % ('' if affy_key_str == '' else 'r-')
        lbvals = plotvals[lb_metric]  # should really use these shorthands for the previous plot as well
        affyvals = plotvals['affinity']
        sorted_affyvals = sorted(affyvals, reverse=True)
        for percentile in numpy.arange(*ptile_range_tuple):
            lb_ptile_val = numpy.percentile(lbvals, percentile)  # lb value corresponding to <percentile>
            corresponding_affinities = [affy for lb, affy in zip(lbvals, affyvals) if lb > lb_ptile_val]  # affinities corresponding to lb greater than <lb_ptile_val> (i.e. the affinities that you'd get if you took all the lb values greater than that)
            corr_affy_ptiles = [stats.percentileofscore(affyvals, caffy) for caffy in corresponding_affinities]  # affinity percentiles corresponding to each of these affinities  # NOTE this is probably really slow
            if len(corr_affy_ptiles) == 0:
                if debug:
                    print '   %5.0f    no vals' % percentile
                continue
            ptile_vals['lb_ptiles'].append(float(percentile))  # stupid numpy-specific float classes (I only care because I write it to a yaml file below)
            ptile_vals['mean_affy_ptiles'].append(float(numpy.mean(corr_affy_ptiles)))

            # make a "perfect" line using the actual affinities, as opposed to just a straight line (this accounts better for, e.g. the case where the top N affinities are all the same)
            n_to_take = len(corresponding_affinities)  # this used to be (in general) different than the number we took above, hence the weirdness/duplication (could probably clean up at this point)
            corresponding_perfect_affy_vals = sorted_affyvals[:n_to_take]
            corr_perfect_affy_ptiles = [stats.percentileofscore(affyvals, cpaffy) for cpaffy in corresponding_perfect_affy_vals]  # NOTE this is probably really slow
            ptile_vals['perfect_vals'].append(float(numpy.mean(corr_perfect_affy_ptiles)))

            if debug:
                print '   %5.0f   %5.2f     %4d  %8.4f  %5.0f       | %4d    %5.0f' % (percentile, lb_ptile_val, len(corresponding_affinities), numpy.mean(corresponding_affinities), ptile_vals['mean_affy_ptiles'][-1], n_to_take, ptile_vals['perfect_vals'][-1])

            # old way of adding a 'no correlation' line:
            # # add a horizontal line at 50 to show what it'd look like if there was no correlation (this is really wasteful... although it does have a satisfying wiggle to it. Now using a plain flat line [below])
            # shuffled_lb_vals = copy.deepcopy(lbvals)
            # random.shuffle(shuffled_lb_vals)
            # NON_corresponding_affinities = [affy for lb, affy in zip(shuffled_lb_vals, affyvals) if lb > lb_ptile_val]
            # NON_corr_affy_ptiles = [stats.percentileofscore(affyvals, caffy) for caffy in NON_corresponding_affinities]
            # ptile_vals['reshuffled_vals'].append(numpy.mean(NON_corr_affy_ptiles))
        return ptile_vals
    # ----------------------------------------------------------------------------------------
    def getplotdir(extrastr=''):
        return '%s/%s-vs%s-affinity%s' % (baseplotdir, lb_metric, affy_key_str, extrastr)
    # ----------------------------------------------------------------------------------------
    def icstr(iclust):
        return '-all-clusters' if iclust is None else '-iclust-%d' % iclust
    # ----------------------------------------------------------------------------------------
    def make_scatter_plot(plotvals, iclust=None, plotstr=None, vspstuff=None):
        tmplbstr, tmpaffystr, tmpclstr = lb_metric, 'affinity', icstr(iclust)
        xlabel, ylabel = '%s affinity' % affy_key_str.replace('-', ''), lb_label
        title = '%s on %s tree' % (lb_metric.upper(), plot_str)
        if affy_key_str != '':
            tmpaffystr = '%s-%s' % (affy_key_str, tmpaffystr)
        if vspstuff is not None:
            assert iclust is None
            tmplbstr = '%s-%s' % (vspstuff[lb_metric], tmplbstr)
            tmpaffystr = '%s-%s' % (vspstuff['affinity'], tmpaffystr)
            tmpclstr = '-per-cluster'
            title += ' (per family)'
            xlabel = '%s %s' % (vspstuff['affinity'], xlabel)
            ylabel = '%s %s' % (vspstuff[lb_metric], ylabel)
        else:
            if iclust is None:
                title += ' (%d families together)' % len(lines) if iclust is None else ''
        plotname = '%s-vs-%s-%s-tree%s' % (tmplbstr, tmpaffystr, plot_str, tmpclstr)
        warn_text = None
        if len(lines) > 1 and iclust is None and 'relative' in affy_key:  # maybe I should just not make the plot, but then the html would look weird
            warn_text = 'wrong/misleading'
        fn = plot_2d_scatter(plotname, getplotdir(), plotvals, lb_metric, ylabel, title, xlabel=xlabel, warn_text=warn_text)
        if iclust is None: # or iclust < n_per_row:
            fnames[-2 if vspstuff is not None else -1].append(fn)
    # ----------------------------------------------------------------------------------------
    def ptile_plotname(iclust):
        return '%s-vs%s-affinity-%s-tree-ptiles%s' % (lb_metric, affy_key_str, plot_str, icstr(iclust))

    # ----------------------------------------------------------------------------------------
    def make_ptile_plot(ptile_vals, iclust=None):
        fig, ax = plotting.mpl_init()
        ax.plot(ptile_vals['lb_ptiles'], ptile_vals['mean_affy_ptiles'], linewidth=3, alpha=0.7)
        # ax.plot(ax.get_xlim(), [50 + 0.5 * x for x in ax.get_xlim()], linewidth=3, alpha=0.7, color='darkgreen', linestyle='--', label='perfect correlation')  # straight line
        ax.plot(ptile_vals['lb_ptiles'], ptile_vals['perfect_vals'], linewidth=3, alpha=0.7, color='darkgreen', linestyle='--', label='perfect correlation')  # perfect vals
        ax.plot(ax.get_xlim(), (50, 50), linewidth=3, alpha=0.7, color='darkred', linestyle='--', label='no correlation')  # straight line
        # ax.plot(ptile_vals['lb_ptiles'], ptile_vals['reshuffled_vals'], linewidth=3, alpha=0.7, color='darkred', linestyle='--', label='no correlation')  # reshuffled vals
        if len(lines) > 1 and iclust is None:
            ax.text(0.6 * ax.get_xlim()[1], 0.95 * ax.get_ylim()[1], 'choosing among %d families' % len(lines), fontsize=17, fontweight='bold')  # , color='red'
            if 'relative' in affy_key:  # maybe I should just not make the plot, but then the html would look weird
                ax.text(0.6 * ax.get_xlim()[1], 0.75 * ax.get_ylim()[1], 'wrong/misleading', fontsize=30, fontweight='bold', color='red')
        fn = plotting.mpl_finish(ax, getplotdir('-ptiles'), ptile_plotname(iclust), xbounds=(ptile_range_tuple[0], ptile_range_tuple[1]), ybounds=(45, 100), leg_loc=(0.5, 0.2),
                                 title='potential %s thresholds (%s tree)' % (lb_metric.upper(), plot_str),
                                 xlabel='%s threshold (percentile)' % lb_metric.upper(),
                                 ylabel='mean percentile of\nresulting %s' % ' '.join([affy_key_str.replace('-', ''), 'affinities']))
        if iclust is None:
            fnames[-1].append(fn)
    # ----------------------------------------------------------------------------------------
    def getcorr(xvals, yvals):
        return numpy.corrcoef(xvals, yvals)[0, 1]
    # ----------------------------------------------------------------------------------------
    def getcorrkey(xstr, ystr):
        return '-vs-'.join([xstr, ystr])

    # ----------------------------------------------------------------------------------------
    if fnames is None:  # not much point since we're not returning it any more
        fnames = []
    fnames += [[], []]
    affy_key_str = '-relative' if 'relative' in affy_key else ''
    vtypes = ['affinity', lb_metric]  # NOTE this puts relative affinity under the (plain) affinity key, which is kind of bad maybe i think probably
    summary_fcns = {'mean' : numpy.mean, 'max' : max}  # ways in which we summarize the affinity or lb value for all cells in a family
    for estr in ['', '-ptiles']:
        utils.prep_dir(getplotdir(estr), wildlings=['*.svg', '*.yaml'])

    # first plot lb metric vs affinity scatter (all clusters)
    all_plotvals = {vt : [] for vt in vtypes}  # every cell from every cluster plotted together;
    cluster_average_plotvals = {st : {vt : [] for vt in vtypes} for st in summary_fcns}  # each cluster plotted as one point using a summary over its cells (max or mean) for affinity and lb
    per_clust_ptile_vals = {}  # percentile values corresponding to choosing cells only within each cluster (as opposed to among all clusters together)
    all_correlation_vals = {}  # not really analagous to <all_plotvals>, but also not entirely not analagous
    for iclust, line in enumerate(lines):
        iclust_plotvals = get_plotvals(line)
        for vt in vtypes:
            all_plotvals[vt] += iclust_plotvals[vt]
            for st in cluster_average_plotvals:  # store both max and mean for affinity and lb
                cluster_average_plotvals[st][vt].append(summary_fcns[st](iclust_plotvals[vt]))
        iclust_ptile_vals = get_ptile_vals(iclust_plotvals)
        per_clust_ptile_vals['iclust-%d'%iclust] = iclust_ptile_vals
        all_correlation_vals['iclust-%d'%iclust] = {getcorrkey(*vtypes) : getcorr(*[iclust_plotvals[vt] for vt in vtypes])}  # might be cleaner to add this to <per_clust_ptile_vals>, but then I'd have to change the name
        if not only_csv and len(iclust_plotvals['affinity']) > 0:
            make_scatter_plot(iclust_plotvals, iclust=iclust)
            make_ptile_plot(iclust_ptile_vals, iclust=iclust)

    all_correlation_vals[getcorrkey(*vtypes)] = getcorr(*[all_plotvals[vt] for vt in vtypes])
    for st1, st2 in itertools.product(summary_fcns, repeat=2):  # all four combos and orderings of max/mean 
        vspairs = zip(vtypes, (st1, st2))  # assign this (st1, st2) combo to lb and affinity based on their order in <vtypes>
        vspdict = {v : s for v, s in vspairs}  # need to also access it by key
        tmpvals = {vt : cluster_average_plotvals[st][vt] for vt, st in vspairs}  # e.g. 'affinity' : <max affinity value list>, 'lbi' : <mean lbi value list> (this is necessary(ish?) to flip the order/depth of st/vt keys)
        tkey = getcorrkey('%s-affinity%s' % (vspdict['affinity'], affy_key_str), '%s-%s' % (vspdict[lb_metric], lb_metric))
        all_correlation_vals[tkey] = getcorr(tmpvals['affinity'], tmpvals[lb_metric])
        if not only_csv:
            make_scatter_plot(tmpvals, vspstuff=vspdict)  # a little confusing to call them "per-cluster", since it's a slightly different meaning of "per-cluster" than <per_clust_ptile_vals>, but they both make sense
    if not only_csv:
        fnames.append([])
        make_scatter_plot(all_plotvals)

    # then plot potential lb cut thresholds with percentiles
    all_ptile_vals = get_ptile_vals(all_plotvals)  # "averaged" might be a better name than "all", but that's longer
    with open('%s/%s.yaml' % (getplotdir('-ptiles'), ptile_plotname(None)), 'w') as yfile:
        yamlfo = {'percentiles' : {k : v for k, v in all_ptile_vals.items() + per_clust_ptile_vals.items()},
                  'correlations' : all_correlation_vals}
        json.dump(yamlfo, yfile)
    if not only_csv and len(all_plotvals[lb_metric]) > 0:
        make_ptile_plot(all_ptile_vals)

# ----------------------------------------------------------------------------------------
def plot_lb_vs_ancestral_delta_affinity(baseplotdir, true_lines, lb_metric, lb_label, plot_str='true', ptile_range_tuple=(50., 100., 1.), min_affinity_change=1e-6, n_max_steps=15, only_csv=False, fnames=None, n_per_row=4, debug=False):
    # plot lb[ir] vs number of ancestors to nearest affinity decrease (well, decrease as you move upwards in the tree/backwards in time)
    # NOTE now that I've done a huge refactor, this fcn is very similar to plot_lb_vs_affinity(), so they could be eventually combined to clean up quite a bit
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
        if abs(max_diff) / numpy.mean(affinity_changes) > 0.2:
            print'      %s not all affinity increases were the same size (min: %.4f   max: %.4f   abs(diff) / mean: %.4f' % (utils.color('yellow', 'warning'), affinity_changes[0], affinity_changes[-1], abs(max_diff) / numpy.mean(affinity_changes))
    # ----------------------------------------------------------------------------------------
    def get_n_ancestor_vals(node, dtree, line, affinity_changes):
        # find number of steps/ancestors to the nearest ancestor with lower affinity than <node>'s
        #   - also finds the corresponding distance, which is to the lower end of the branch containing the corresponding affinity-increasing mutation
        #   - this is chosen so that <n_steps> and <branch_len> are both 0 for the node at the bottom of a branch on which affinity increases, and are *not* the distance *to* the lower-affinity node
        #   - because it's so common for affinity to get worse from ancestor to descendent, it's important to remember that here we are looking for the first ancestor with lower affinity than the node in question, which is *different* to looking for the first ancestor that has lower affinity than one of its immediate descendents (which we could also plot, but it probably wouldn't be significantly different to the metric performance, since for the metric performance we only really care about the left side of the plot, but this only affects the right side)
        #   - <min_affinity_change> is just to eliminate floating point precision issues (especially since we're deriving affinity by inverting kd) (note that at least for now, and with default settings, the affinity changes should all be pretty similar, and not small)
        this_affinity = utils.per_seq_val(line, 'affinities', node.taxon.label)
        if debug:
            print '     %12s %12s %8s %9.4f' % (node.taxon.label, '', '', this_affinity)
    
        ancestor_node = node
        chosen_ancestor_affinity = None
        n_steps, branch_len  = 0, 0.
        while n_steps < n_max_steps:  # note that if we can't find an ancestor with worse affinity, we don't plot the node
            if ancestor_node is dtree.seed_node:
                break
            ancestor_distance = ancestor_node.edge_length  # distance from current <ancestor_node> to its parent (who in the next line becomes <ancestor_node>)
            ancestor_node = ancestor_node.parent_node  #  move one more step up the tree
            ancestor_uid = ancestor_node.taxon.label
            if ancestor_uid not in line['unique_ids']:
                print '    %s ancestor %s of %s not in true line' % (utils.color('yellow', 'warning'), ancestor_uid, node.taxon.label)
                break
            ancestor_affinity = utils.per_seq_val(line, 'affinities', ancestor_uid)
            if this_affinity - ancestor_affinity > min_affinity_change:  # if we found an ancestor with lower affinity, we're done
                chosen_ancestor_affinity = ancestor_affinity
                affinity_changes.append(this_affinity - ancestor_affinity)
                break
            if debug:
                print '     %12s %12s %8.4f %9.4f%s' % ('', ancestor_uid, branch_len, ancestor_affinity, utils.color('green', ' x') if ancestor_node is dtree.seed_node else '')
            n_steps += 1
            branch_len += ancestor_distance
    
        if chosen_ancestor_affinity is None:  # couldn't find ancestor with lower affinity
            return None, None
        if debug:
            print '     %12s %12s %8.4f %9.4f  %s%-9.4f' % ('', ancestor_uid, branch_len, chosen_ancestor_affinity, utils.color('red', '+'), this_affinity - chosen_ancestor_affinity)
        return n_steps, branch_len
    # ----------------------------------------------------------------------------------------
    def get_plotvals(line, xvar):
        plotvals = {vt : [] for vt in [lb_metric, xvar]}  # , 'uids']}
        dtree = treeutils.get_dendro_tree(treestr=line['tree'])
        affinity_changes = []
        for uid in line['unique_ids']:
            node = dtree.find_node_with_taxon_label(uid)
            if node is dtree.seed_node:  # root doesn't have any ancestors
                continue
            lbval = line['tree-info']['lb'][lb_metric][uid]
            if lb_metric == 'lbr' and lbval == 0:  # lbr equals 0 should really be treated as None/missing
                continue
            n_steps, branch_len = get_n_ancestor_vals(node, dtree, line, affinity_changes)  # also adds to <affinity_changes>
            if n_steps is None:
                continue
            plotvals[xvar].append(n_steps if xvar == 'n-ancestor' else branch_len)
            plotvals[lb_metric].append(lbval)
            # plotvals['uids'].append(uid)
        check_affinity_changes(affinity_changes)
        return plotvals

    # ----------------------------------------------------------------------------------------
    def get_ptile_vals(plotvals, xvar, xlabel, dbgstr=None):
        # NOTE xvar and xlabel refer to the x axis on the scatter plot from which we make this ptile plot. On this ptile plot it's the y axis. (I tried calling it something else, but it was more confusing)
        xkey = 'mean_%s_ptiles' % xvar
        ptile_vals = {'lb_ptiles' : [], xkey : [], 'perfect_vals' : []}  # , 'reshuffled_vals' : []}
        if len(plotvals[xvar]) == 0:
            return ptile_vals
        # plot potential lb cut thresholds with percentiles
        if debug:
            print '    getting ptile vals%s' % ('' if dbgstr is None else (' for %s' % utils.color('blue', dbgstr)))
            print '               %s       N      mean        |  perfect   perfect' % lb_metric
            print '      ptile  threshold  taken  %-13s|  N taken   mean %-13s' % (xlabel, xlabel)
        sorted_xvals = sorted(plotvals[xvar])
        for percentile in numpy.arange(*ptile_range_tuple):
            lb_ptile_val = numpy.percentile(plotvals[lb_metric], percentile)  # lb value corresponding to <percentile>
            corresponding_xvals = [n_anc for lb, n_anc in zip(plotvals[lb_metric], plotvals[xvar]) if lb >= lb_ptile_val]  # n-anc vals corresponding to lb greater than <lb_ptile_val> (i.e. the n-anc vals that you'd get if you took all the lb values greater than that)
            if len(corresponding_xvals) == 0:  # calling it n_anc even though it can now also be branch length
                if debug:
                    print '     %5.0f    no vals' % percentile
                continue
            ptile_vals['lb_ptiles'].append(float(percentile))  # cast to float is so the yaml file written below isn't so damn ugly
            ptile_vals[xkey].append(float(numpy.mean(corresponding_xvals)))

            # values for perfect line
            n_to_take = len(corresponding_xvals)  # this used to be (in general) different than the number we took above, hence the weirdness/duplication (could probably clean up at this point)
            perfect_xvals = sorted_xvals[:n_to_take]
            ptile_vals['perfect_vals'].append(float(numpy.mean(perfect_xvals)))

            if debug:
                print '     %5.0f    %5.2f    %4d   %6.2f        | %4d      %-8s' % (percentile, lb_ptile_val, len(corresponding_xvals), ptile_vals[xkey][-1], n_to_take, ('%8.2f' if xvar == 'n-ancestor' else '%8.6f') % ptile_vals['perfect_vals'][-1])
        return ptile_vals
    # ----------------------------------------------------------------------------------------
    def getplotdir(xvar, extrastr=''):
        return '%s/%s-vs-%s%s' % (baseplotdir, lb_metric, xvar, extrastr)
    # ----------------------------------------------------------------------------------------
    def icstr(iclust):
        return '-all-clusters' if iclust is None else '-iclust-%d' % iclust
    # ----------------------------------------------------------------------------------------
    def make_scatter_plot(plotvals, xvar, iclust=None):
        fn = plot_2d_scatter('%s-vs-%s-%s-tree%s' % (lb_metric, xvar, plot_str, icstr(iclust)), getplotdir(xvar), plotvals, lb_metric, lb_label, '%s (true tree)' % lb_metric.upper(), xvar=xvar, xlabel='%s since affinity increase' % xlabel, log='y' if lb_metric == 'lbr' else '')
        if iclust is None: # or iclust < n_per_row:
            fnames[-1].append(fn)
    # ----------------------------------------------------------------------------------------
    def ptile_plotname(xvar, iclust):
        return '%s-vs-%s-%s-tree-ptiles%s' % (lb_metric, xvar, plot_str, icstr(iclust))
    # ----------------------------------------------------------------------------------------
    def make_ptile_plot(plotvals, ptile_vals, xvar, xlabel, iclust=None):
        fig, ax = plotting.mpl_init()
        xmean = numpy.mean(plotvals[xvar])  # NOTE this is mean of "xvar", which is the x axis on the scatter plot, but here it's the y axis on the ptile plot
        xkey = 'mean_%s_ptiles' % xvar
        ymax = max([xmean] + ptile_vals[xkey] + ptile_vals['perfect_vals'])
        ax.plot(ptile_vals['lb_ptiles'], ptile_vals[xkey], linewidth=3, alpha=0.7)
        ax.plot(ax.get_xlim(), (xmean, xmean), linewidth=3, alpha=0.7, color='darkred', linestyle='--', label='no correlation')
        ax.plot(ptile_vals['lb_ptiles'], ptile_vals['perfect_vals'], linewidth=3, alpha=0.7, color='darkgreen', linestyle='--', label='perfect correlation')
        fn = plotting.mpl_finish(ax, getplotdir(xvar, extrastr='-ptiles'), ptile_plotname(xvar, iclust), xbounds=(ptile_range_tuple[0], ptile_range_tuple[1]), ybounds=(-0.02*ymax, 1.1*ymax), leg_loc=(0.5, 0.6),
                                 title='potential %s thresholds (%s tree)' % (lb_metric.upper(), plot_str),
                                 xlabel='%s threshold (percentile)' % lb_metric.upper(),
                                 ylabel='mean %s\nsince affinity increase' % xlabel)
        if iclust is None:
            fnames[-1].append(fn)

    # ----------------------------------------------------------------------------------------
    if fnames is None:  # no real effect (except not crashing) since we're not returning it any more
        fnames = []
    fnames += [[]]
    xvar_list = collections.OrderedDict([(xvar, xlabel) for metric, xvar, xlabel in lb_metric_axis_stuff if metric == 'lbr'])
    for xvar, estr in itertools.product(xvar_list, ['', '-ptiles']):
        utils.prep_dir(getplotdir(xvar, extrastr=estr), wildlings=['*.svg', '*.yaml'])
    if debug:
        print 'finding ancestors with most recent affinity increases'
    for xvar, xlabel in xvar_list.items():
        all_plotvals = {vt : [] for vt in [lb_metric, xvar]}  # , 'uids']}
        per_clust_ptile_vals = {}
        for iclust, line in enumerate(true_lines):
            if debug:
                if iclust == 0:
                    print ' %s' % utils.color('green', xvar)
                print '  %s' % utils.color('blue', 'iclust %d' % iclust)
                print '         node        ancestors  distance   affinity (%sX: change for chosen ancestor, %s: reached root without finding lower-affinity ancestor)' % (utils.color('red', '+'), utils.color('green', 'x'))
            iclust_plotvals = get_plotvals(line, xvar)
            for vtype in all_plotvals:
                all_plotvals[vtype] += iclust_plotvals[vtype]
            iclust_ptile_vals = get_ptile_vals(iclust_plotvals, xvar, xlabel, dbgstr='iclust %d'%iclust)
            per_clust_ptile_vals['iclust-%d'%iclust] = iclust_ptile_vals
            if not only_csv and len(iclust_plotvals[xvar]) > 0:
                make_scatter_plot(iclust_plotvals, xvar, iclust=iclust)
                make_ptile_plot(iclust_plotvals, iclust_ptile_vals, xvar, xlabel, iclust=iclust)
        if not only_csv:
            make_scatter_plot(all_plotvals, xvar)
        all_ptile_vals = get_ptile_vals(all_plotvals, xvar, xlabel, dbgstr='all clusters')  # "averaged" might be a better name than "all", but that's longer
        with open('%s/%s.yaml' % (getplotdir(xvar, extrastr='-ptiles'), ptile_plotname(xvar, None)), 'w') as yfile:
            json.dump({'percentiles' : {k : v for k, v in all_ptile_vals.items() + per_clust_ptile_vals.items()}}, yfile)  # not adding the new correlation keys atm (like in the lb vs affinity fcn)
        if not only_csv and len(all_plotvals[lb_metric]) > 0:
            make_ptile_plot(all_plotvals, all_ptile_vals, xvar, xlabel)

# ----------------------------------------------------------------------------------------
def plot_true_vs_inferred_lb(plotdir, true_lines, inf_lines, lb_metric, lb_label, debug=False):
    plotvals = {val_type : {uid : l['tree-info']['lb'][lb_metric][uid] for l in lines for uid in l['unique_ids']}
                for val_type, lines in (('true', true_lines), ('inf', inf_lines))}
    common_uids = set(plotvals['true']) & set(plotvals['inf'])  # there should/may be a bunch of internal nodes in the simulation lines but not in the inferred lines, but otherwise they should have the same uids
    plotvals = {val_type : [plotvals[val_type][uid] for uid in common_uids] for val_type in plotvals}
    plotname = '%s-true-vs-inferred' % lb_metric
    fn = plot_2d_scatter(plotname, plotdir, plotvals, 'inf', '%s on inferred tree' % lb_metric.upper(), 'true vs inferred %s' % lb_metric.upper(), xvar='true', xlabel='%s on true tree' % lb_metric.upper())
    return [fn]

# ----------------------------------------------------------------------------------------
def get_lb_tree_cmd(treestr, outfname, lb_metric, affy_key, ete_path, subworkdir, metafo=None, tree_style=None):
    treefname = '%s/tree.nwk' % subworkdir
    metafname = '%s/meta.yaml' % subworkdir
    if not os.path.exists(subworkdir):
        os.makedirs(subworkdir)
    with open(treefname, 'w') as treefile:
        treefile.write(treestr)
    cmdstr = './bin/plot-lb-tree.py --treefname %s' % treefname
    if metafo is not None:
        with open(metafname, 'w') as metafile:
            yaml.dump(metafo, metafile)
        cmdstr += ' --metafname %s' % metafname
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
def plot_lb_trees(baseplotdir, lines, ete_path, base_workdir, is_true_line=False, tree_style=None):
    workdir = '%s/ete3-plots' % base_workdir
    plotdir = baseplotdir + '/trees'
    utils.prep_dir(plotdir, wildlings='*.svg')

    if not os.path.exists(workdir):
        os.makedirs(workdir)
    cmdfos = []
    for lb_metric, lb_label in treeutils.lb_metrics.items():
        for iclust, line in enumerate(lines):  # note that <min_tree_metric_cluster_size> was already applied in treeutils
            treestr = get_tree_from_line(line, is_true_line)
            for affy_key in treeutils.affy_keys[lb_metric]:
                metafo = copy.deepcopy(line['tree-info']['lb'])
                if affy_key in line:  # either 'affinities' or 'relative_affinities'
                    metafo[utils.reversed_input_metafile_keys[affy_key]] = {uid : affy for uid, affy in zip(line['unique_ids'], line[affy_key])}
                outfname = '%s/%s-tree-iclust-%d%s.svg' % (plotdir, lb_metric, iclust, '-relative' if 'relative' in affy_key else '')
                cmdfos += [get_lb_tree_cmd(treestr, outfname, lb_metric, affy_key, ete_path, '%s/sub-%d' % (workdir, len(cmdfos)), metafo=metafo, tree_style=tree_style)]

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
