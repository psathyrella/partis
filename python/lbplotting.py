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

import utils
from hist import Hist
import treeutils
import plotting

# ----------------------------------------------------------------------------------------
lb_metric_axis_stuff = [('lbi', 'affinity', 'affinity'), ('lbr', 'n-ancestor', 'N ancestors'), ('lbr', 'branch-length', 'branch length')]

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
                bin_center = (hi_edge + lo_edge) / 2.
                for _ in range(bin_contents[ibin]):
                    hist.fill(bin_center)
                    xmin = lo_edge if xmin is None else min(xmin, lo_edge)
                    xmax = hi_edge if xmax is None else max(xmax, hi_edge)
            hists.append(hist)
            labels.append('%d (%.1f)' % (obs_time, hist.get_mean()))

        # hists = [Hist(1, xmin, xmax) if h is None else h for h in hists]  # replace the None values with empty hists
        hists, labels = zip(*[(h, l) for h, l in zip(hists, labels) if h is not None])
        return hists, labels, xmin, xmax

    # ----------------------------------------------------------------------------------------
    all_hists, all_labels, xmin, xmax = get_hists(histfname)
    jpdata = []
    for hist in all_hists:
        jpdata.append([x for x, y in zip(hist.get_bin_centers(), hist.bin_contents) for _ in range(int(y)) if x > xmin and x < xmax])  # NOTE this is repeating the 'for _ in range()' in the fcn above, but that's because I used to be actually using the Hist()s, and maybe I will again

    fig, ax = plotting.mpl_init()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')  # i don't know why it has to warn me that it's clearing the fig/ax I'm passing in, and I don't know how else to stop it
        fig, axes = joypy.joyplot(jpdata, labels=all_labels, fade=True, hist=True, overlap=0.5, ax=ax, x_range=(xmin, xmax), bins=int(xmax - xmin))
    # NOTE do *not* set your own x ticks/labels in the next line, since they'll be in the wrong place (i.e. not the same as where joypy puts them)
    plotting.mpl_finish(ax, plotdir, plotname, title=title, xlabel=xlabel, ylabel='generation', leg_loc=(0.7, 0.45)) #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_kd_vals(plotdir, event):
    n_muts, kd_changes = [], []
    dtree = treeutils.get_dendro_tree(treestr=event['tree'])
    for node in dtree.preorder_internal_node_iter():
        for child in node.child_nodes():
            inode = event['unique_ids'].index(node.taxon.label)
            ichild = event['unique_ids'].index(child.taxon.label)
            node_affinity = event['affinities'][inode]
            child_affinity = event['affinities'][ichild]
            n_muts.append(utils.hamming_distance(event['input_seqs'][inode], event['input_seqs'][ichild]))
            kd_changes.append(1./child_affinity - 1./node_affinity)

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
def get_tree_from_line(line, is_simu):
    if is_simu:
        return line['tree']
    if 'tree-info' not in line:  # if 'tree-info' is missing, it should be because it's a small cluster in data that we skipped when calculating lb values
        return None
    return line['tree-info']['lb']['tree']

# ----------------------------------------------------------------------------------------
def plot_lb_vs_shm(baseplotdir, lines_to_use, is_simu=False, n_per_row=4):  # <is_simu> is there because we want the true and inferred lines to keep their trees in different places, because the true line just has the one, true, tree, while the inferred line could have a number of them (yes, this means I maybe should have called it the 'true-tree' or something)
    sorted_lines = sorted([l for l in lines_to_use if get_tree_from_line(l, is_simu) is not None], key=lambda l: len(l['unique_ids']), reverse=True)
    fnames = [[]]

    # note: all clusters together
    subfnames = {lb_metric : [] for lb_metric in treeutils.lb_metrics}
    for lb_metric, lb_label in treeutils.lb_metrics.items():
        plotvals = {x : {'leaf' : [], 'internal' : []} for x in ['shm', lb_metric]}
        basetitle = '%s %s vs SHM' % ('true' if is_simu else 'inferred', lb_metric.upper())
        for iclust, line in enumerate(sorted_lines):  # get depth/n_mutations for each node
            iclust_plotvals = {x : {'leaf' : [], 'internal' : []} for x in ['shm', lb_metric, 'uids']}
            dtree = treeutils.get_dendro_tree(treestr=get_tree_from_line(line, is_simu))
            n_max_mutes = max(line['n_mutations'])  # don't generally have n mutations for internal nodes, so use this to rescale the depth in the tree
            max_depth = max(n.distance_from_root() for n in dtree.leaf_node_iter())
            for node in dtree.preorder_node_iter():
                if lb_metric == 'lbr' and line['tree-info']['lb'][lb_metric][node.taxon.label] == 0:  # lbr equals 0 should really be treated as None/missing
                    continue
                iseq = line['unique_ids'].index(node.taxon.label) if node.taxon.label in line['unique_ids'] else None
                n_muted = line['n_mutations'][iseq] if node.taxon.label in line['unique_ids'] else node.distance_from_root() * n_max_mutes / float(max_depth)
                tkey = 'leaf' if node.is_leaf() else 'internal'
                iclust_plotvals['shm'][tkey].append(n_muted)
                iclust_plotvals[lb_metric][tkey].append(line['tree-info']['lb'][lb_metric][node.taxon.label])
                affyval = line['affinities'][iseq] if 'affinities' in line and iseq is not None else None
                if not is_simu:
                    iclust_plotvals['uids'][tkey].append(node.taxon.label if affyval is not None else None)
            plotname = '%s-vs-shm-iclust-%d' % (lb_metric, iclust)
            title = '%s (%d observed, %d total)' % (basetitle, len(line['unique_ids']), len(line['tree-info']['lb'][lb_metric]))
            fn = plot_2d_scatter(plotname, '%s/%s-vs-shm' % (baseplotdir, lb_metric), iclust_plotvals, lb_metric, lb_label, title, xvar='shm', xlabel='N mutations', leg_loc=(0.7, 0.75), log='y' if lb_metric == 'lbr' else '')
            if iclust < n_per_row:  # i.e. only put one row's worth in the html
                subfnames[lb_metric].append(fn)
            for vtype in [vt for vt in plotvals if vt != 'uids']:
                for ltype in plotvals[vtype]:
                    plotvals[vtype][ltype] += iclust_plotvals[vtype][ltype]
        plotname = '%s-vs-shm' % lb_metric
        plot_2d_scatter(plotname, baseplotdir, plotvals, lb_metric, lb_label, '%s (all clusters)' % basetitle, xvar='shm', xlabel='N mutations', leg_loc=(0.7, 0.75), log='y' if lb_metric == 'lbr' else '')
        fnames[-1].append('%s/%s.svg' % (baseplotdir, plotname))
    fnames += [subfnames[lbm] for lbm in treeutils.lb_metrics]

    return fnames

# ----------------------------------------------------------------------------------------
def plot_lb_distributions(baseplotdir, lines_to_use, n_per_row=4):
    sorted_lines = sorted([l for l in lines_to_use if 'tree-info' in l], key=lambda l: len(l['unique_ids']), reverse=True)  # if 'tree-info' is missing, it should be because it's a small cluster we skipped when calculating lb values
    fnames = []

    for lb_metric, lb_label in treeutils.lb_metrics.items():
        plotdir = baseplotdir + '/' + lb_metric
        utils.prep_dir(plotdir, wildlings=['*.svg'])
        fnames.append([])
        for iclust, line in enumerate(sorted_lines):
            plotvals = line['tree-info']['lb'][lb_metric].values()
            leafskipstr = ''
            if lb_metric == 'lbr':
                plotvals = [v for v in plotvals if v > 0.]  # don't plot the leaf values, they just make the plot unreadable
                leafskipstr = ', skipped %d leaves' % len([v for v in line['tree-info']['lb'][lb_metric].values() if v == 0.])  # ok they're not necessarily leaves, but almost all of them are leaves (and not really sure how a non-leaf could get zero, but some of them seem to)
            hist = Hist(30, 0., max(plotvals), value_list=plotvals)
            fig, ax = plotting.mpl_init()
            hist.mpl_plot(ax) #, square_bins=True, errors=False)
            # ax.text(0.45 * ax.get_xlim()[1], 0.85 * ax.get_ylim()[1], 'size %d' % len(line['unique_ids']), fontsize=17, color='red', fontweight='bold')  # omfg this is impossible to get in the right place
            plotname = '%s-%d' % (lb_metric, iclust)
            plotting.mpl_finish(ax, plotdir, plotname, xlabel=lb_label, log='y' if lb_metric == 'lbr' else '', ylabel='counts', title='%s  (size %d%s)' % (lb_metric.upper(), len(line['tree-info']['lb'][lb_metric]), leafskipstr))
            if iclust < n_per_row:  # i.e. only put one row's worth in the html
                fnames[-1].append('%s/%s.svg' % (plotdir, plotname))
            plotting.make_html(plotdir)

    return fnames

# ----------------------------------------------------------------------------------------
def plot_2d_scatter(plotname, plotdir, plotvals, yvar, ylabel, title, xvar='affinity', xlabel='affinity', log='', leg_loc=None):
    def getall(k):
        if 'leaf' in plotvals[xvar]:
            return [v for tk in plotvals[k] for v in plotvals[k][tk]]
        else:
            return plotvals[k]

    if len(getall(xvar)) == 0:
        # print '    no %s vs affy info' % yvar
        return
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

    xmin, xmax = min(getall(xvar)), max(getall(xvar))
    ymin, ymax = min(getall(yvar)), max(getall(yvar))
    xbounds = xmin - 0.02 * (xmax - xmin), xmax + 0.02 * (xmax - xmin)
    if 'y' in log:
        ybounds = 0.75 * ymin, 1.3 * ymax
    else:
        ybounds = ymin - 0.03 * (ymax - ymin), ymax + 0.08 * (ymax - ymin)
    plotting.mpl_finish(ax, plotdir, plotname, title=title, xlabel=xlabel, ylabel=ylabel, xbounds=xbounds, ybounds=ybounds, log=log, leg_loc=leg_loc)
    return '%s/%s.svg' % (plotdir, plotname)

# ----------------------------------------------------------------------------------------
def plot_lb_vs_affinity(plot_str, plotdir, lines, lb_metric, lb_label, make_all_cluster_plot=False, ptile_range_tuple=(50., 100., 1.), is_simu=False, n_per_row=4, affy_key='affinities', only_csv=False, debug=False):
    # ----------------------------------------------------------------------------------------
    def get_plotvals(iclust, line):
        plotvals = {vt : [] for vt in vtypes + ['uids']}
        # dtree = treeutils.get_dendro_tree(treestr=get_tree_from_line(line, is_simu))  # keeping this here to remind myself how to get the tree if I need it
        if affy_key not in line:
            return plotvals
        for uid, affy in [(u, a) for u, a in zip(line['unique_ids'], line[affy_key]) if a is not None]:
            plotvals['affinity'].append(affy)
            plotvals[lb_metric].append(line['tree-info']['lb'][lb_metric][uid])
            if not is_simu:
                plotvals['uids'].append(uid)
        return plotvals
    # ----------------------------------------------------------------------------------------
    def make_scatter_plot(plotvals, plotname, iclust=None):
        fn = plot_2d_scatter(plotname, '%s/%s-vs-affinity' % (plotdir, lb_metric), plotvals, lb_metric, lb_label, '%s (%s tree)' % (lb_metric.upper(), plot_str), xlabel='%s affinity' % affy_key_str.replace('-', ''))
        if iclust is None or iclust < n_per_row:
            fnames.append(fn)
    # ----------------------------------------------------------------------------------------
    def get_ptile_vals(plotvals):
        ptile_vals = {'lb_ptiles' : [], 'mean_affy_ptiles' : [], 'perfect_vals' : []}  # , 'reshuffled_vals' : []}
        if len(plotvals[lb_metric]) == 0:
            # print '  no affinity values when trying to make lb vs affinity plots'
            return ptile_vals
        if debug:
            print '    ptile   %s     mean affy    mean affy ptile' % lb_metric
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
            if debug:
                print '   %5.0f   %5.2f   %8.4f     %5.0f' % (percentile, lb_ptile_val, numpy.mean(corresponding_affinities), ptile_vals['mean_affy_ptiles'][-1])

            # make a "perfect" line from actual affinities, as opposed to just a straight line (this accounts better for, e.g. the case where the top N affinities are all the same)
            n_to_take = int((1. - percentile / 100) * len(sorted_affyvals))
            corresponding_perfect_affy_vals = sorted_affyvals[:n_to_take]
            corr_perfect_affy_ptiles = [stats.percentileofscore(affyvals, cpaffy) for cpaffy in corresponding_perfect_affy_vals]  # NOTE this is probably really slow
            ptile_vals['perfect_vals'].append(float(numpy.mean(corr_perfect_affy_ptiles)) if len(corr_perfect_affy_ptiles) > 0 else 100)  # not really sure just using 100 is right, but I'm pretty sure it doesn't matter (and it gets rid of a stupid numpy warning)

            # old way of adding a 'no correlation' line:
            # # add a horizontal line at 50 to show what it'd look like if there was no correlation (this is really wasteful... although it does have a satisfying wiggle to it. Now using a plain flat line [below])
            # shuffled_lb_vals = copy.deepcopy(lbvals)
            # random.shuffle(shuffled_lb_vals)
            # NON_corresponding_affinities = [affy for lb, affy in zip(shuffled_lb_vals, affyvals) if lb > lb_ptile_val]
            # NON_corr_affy_ptiles = [stats.percentileofscore(affyvals, caffy) for caffy in NON_corresponding_affinities]
            # ptile_vals['reshuffled_vals'].append(numpy.mean(NON_corr_affy_ptiles))
        return ptile_vals
    # ----------------------------------------------------------------------------------------
    def make_ptile_plot(ptile_vals, plotname):
        fig, ax = plotting.mpl_init()
        ax.plot(ptile_vals['lb_ptiles'], ptile_vals['mean_affy_ptiles'], linewidth=3, alpha=0.7)
        # ax.plot(ax.get_xlim(), [50 + 0.5 * x for x in ax.get_xlim()], linewidth=3, alpha=0.7, color='darkgreen', linestyle='--', label='perfect correlation')  # straight line
        ax.plot(ptile_vals['lb_ptiles'], ptile_vals['perfect_vals'], linewidth=3, alpha=0.7, color='darkgreen', linestyle='--', label='perfect correlation')  # perfect vals
        ax.plot(ax.get_xlim(), (50, 50), linewidth=3, alpha=0.7, color='darkred', linestyle='--', label='no correlation')  # straight line
        # ax.plot(ptile_vals['lb_ptiles'], ptile_vals['reshuffled_vals'], linewidth=3, alpha=0.7, color='darkred', linestyle='--', label='no correlation')  # reshuffled vals
        plotting.mpl_finish(ax, plotdir, plotname, xbounds=(ptile_range_tuple[0], ptile_range_tuple[1]), ybounds=(45, 100), leg_loc=(0.5, 0.2),
                            title='potential %s thresholds (%s tree)' % (lb_metric.upper(), plot_str),
                            xlabel='%s threshold (percentile)' % lb_metric.upper(),
                            ylabel='mean percentile of\nresulting %s' % ' '.join([affy_key_str.replace('-', ''), 'affinities']))
        fnames.append('%s/%s.svg' % (plotdir, plotname))
    # ----------------------------------------------------------------------------------------
    def dump_plot_info_to_yaml(fname, yamlfo):
        with open(fname, 'w') as yfile:
            json.dump(yamlfo, yfile)

    # ----------------------------------------------------------------------------------------
    fnames = []
    affy_key_str = '-relative' if 'relative' in affy_key else ''
    vtypes = [lb_metric, 'affinity']

    # first plot lb metric vs affinity scatter (all clusters)
    all_plotvals = {val_type : [] for val_type in [lb_metric, 'affinity']}  # , 'uids']}  # NOTE this puts relative affinity under the (plain) affinity key, which is kind of bad maybe i think probably
    for iclust, line in enumerate(lines):
        iclust_plotvals = get_plotvals(iclust, line)
        if not make_all_cluster_plot and len(iclust_plotvals['affinity']) > 0 and not only_csv:
            make_scatter_plot(iclust_plotvals, 'iclust-%d%s' % (iclust, affy_key_str), iclust=iclust)
        for vtype in vtypes:
            all_plotvals[vtype] += iclust_plotvals[vtype]
    if make_all_cluster_plot and not only_csv:
        make_scatter_plot(all_plotvals, '%s-vs-affinity-%s-tree%s' % (lb_metric, plot_str, affy_key_str))

    # then plot potential lb cut thresholds with percentiles
    ptile_vals = get_ptile_vals(all_plotvals)
    ptile_plotname = '%s-vs-affinity-%s-tree-ptiles%s' % (lb_metric, plot_str, affy_key_str)
    dump_plot_info_to_yaml('%s/%s.yaml' % (plotdir, ptile_plotname), ptile_vals)
    if not only_csv and len(all_plotvals[lb_metric]) > 0:
        make_ptile_plot(ptile_vals, ptile_plotname)

    return [fnames]

# ----------------------------------------------------------------------------------------
def plot_lb_vs_delta_affinity(plotdir, true_lines, lb_metric, lb_label, debug=False):
    fnames = []

    delta_affinity_vals = {val_type : [] for val_type in [lb_metric, 'delta-affinity']}
    for line in true_lines:
        dtree = treeutils.get_dendro_tree(treestr=line['tree'])
        for uid, affinity in zip(line['unique_ids'], line['affinities']):
            node = dtree.find_node_with_taxon_label(uid)
            if node is dtree.seed_node:  # root won't have a parent
                continue
            parent_uid = node.parent_node.taxon.label
            if parent_uid not in line['unique_ids']:
                print '    %s parent %s of %s not in true line' % (utils.color('yellow', 'warning'), parent_uid, uid)
                continue
            iparent = line['unique_ids'].index(parent_uid)
            parent_affinity = line['affinities'][iparent]
            delta_affinity_vals['delta-affinity'].append(affinity - parent_affinity)
            delta_affinity_vals[lb_metric].append(line['tree-info']['lb'][lb_metric][uid])
    fig, ax = plotting.mpl_init()
    ax.scatter(delta_affinity_vals['delta-affinity'], delta_affinity_vals[lb_metric], alpha=0.4)
    sorted_xvals = sorted(delta_affinity_vals['delta-affinity'])  # not sure why, but ax.scatter() is screwing up the x bounds
    xmin, xmax = sorted_xvals[0], sorted_xvals[-1]
    # ax.hexbin(delta_affinity_vals['delta-affinity'], delta_affinity_vals[lb_metric], gridsize=15, cmap=plt.cm.Blues)
    plotname = '%s-vs-delta-affinity' % lb_metric
    plotting.mpl_finish(ax, plotdir, plotname, title='%s (true tree)' % lb_metric.upper(), xlabel='affinity change (from parent)', ylabel=lb_label, xbounds=(1.05 * xmin, 1.05 * xmax))  # NOTE factor on <xmin> is only right if xmin is negative, but it should always be
    fnames.append('%s/%s.svg' % (plotdir, plotname))

    return [fnames]

# ----------------------------------------------------------------------------------------
def make_ancestral_ptile_plot(plotvals, vstr, vstr_label, lb_metric, ptile_range_tuple, plotdir, plot_str, only_csv=False, debug=False):
    def dump_plot_info_to_yaml(fname, yamlfo):
        with open(fname, 'w') as yfile:
            json.dump(yamlfo, yfile)

    plotname = '%s-vs-%s-%s-tree-ptiles' % (lb_metric, vstr, plot_str)
    if len(plotvals[vstr]) == 0:
        print '    no %s vals' % vstr
        dump_plot_info_to_yaml('%s/%s.yaml' % (plotdir, plotname), {'lb_ptiles' : [], 'mean_%s_ptiles' % vstr : [], 'perfect_vals' : []})
        return plotname

    # plot potential lb cut thresholds with percentiles
    if debug:
        print '    ptile   %s     mean %s' % (lb_metric, vstr_label)
    lb_ptile_vals = {'lb_ptiles' : [], 'mean_%s' % vstr: []}
    lbvals = plotvals[lb_metric]
    n_anc_vals = plotvals[vstr]  # leaving this as <n_anc_vals>, but it can now also be branch length
    for percentile in numpy.arange(*ptile_range_tuple):
        lb_ptile_val = numpy.percentile(lbvals, percentile)  # lb value corresponding to <percentile>
        corresponding_n_anc_vals = [n_anc for lb, n_anc in zip(lbvals, n_anc_vals) if lb >= lb_ptile_val]  # n-anc vals corresponding to lb greater than <lb_ptile_val> (i.e. the n-anc vals that you'd get if you took all the lb values greater than that)
        if len(corresponding_n_anc_vals) == 0:
            if debug:
                print '   %5.0f    no vals' % percentile
            continue
        lb_ptile_vals['lb_ptiles'].append(float(percentile))  # cast to float is so the yaml file written below isn't so damn ugly
        lb_ptile_vals['mean_%s' % vstr].append(float(numpy.mean(corresponding_n_anc_vals)))
        if debug:
            print '   %5.0f   %5.2f   %8.4f' % (percentile, lb_ptile_val, lb_ptile_vals['mean_%s' % vstr][-1])

    # make the "perfect" line
    if debug:
        print '  perfect:'
        print '    ptile  n_taken   %s     mean %s' % ('affy', vstr_label)
    perfect_ptile_vals = {'ptiles' : [], 'mean_%s' % vstr: []}
    sorted_n_anc_vals = sorted(plotvals[vstr])
    for percentile in numpy.arange(*ptile_range_tuple):
        n_to_take = int((1. - percentile / 100) * len(sorted_n_anc_vals))
        corresponding_n_anc_vals = sorted_n_anc_vals[:n_to_take]
        if len(corresponding_n_anc_vals) == 0:
            if debug:
                print '   %5.0f   %5d   no vals' % (percentile, n_to_take)
            continue
        perfect_ptile_vals['ptiles'].append(float(percentile))
        perfect_ptile_vals['mean_%s' % vstr].append(float(numpy.mean(corresponding_n_anc_vals)))
        if debug:
            print '   %5.0f   %5d  %8.4f' % (percentile, n_to_take, perfect_ptile_vals['mean_%s' % vstr][-1])

    if not only_csv:
        fig, ax = plotting.mpl_init()
        mean_n_anc = numpy.mean(n_anc_vals)
        ymax = max([mean_n_anc] + lb_ptile_vals['mean_%s' % vstr] + perfect_ptile_vals['mean_%s' % vstr])
        ax.plot(lb_ptile_vals['lb_ptiles'], lb_ptile_vals['mean_%s' % vstr], linewidth=3, alpha=0.7)
        ax.plot(ax.get_xlim(), (mean_n_anc, mean_n_anc), linewidth=3, alpha=0.7, color='darkred', linestyle='--', label='no correlation')
        ax.plot(perfect_ptile_vals['ptiles'], perfect_ptile_vals['mean_%s' % vstr], linewidth=3, alpha=0.7, color='darkgreen', linestyle='--', label='perfect correlation')
        plotting.mpl_finish(ax, plotdir, plotname, xbounds=(ptile_range_tuple[0], ptile_range_tuple[1]), ybounds=(-0.02*ymax, 1.1*ymax), leg_loc=(0.5, 0.6),
                            title='potential %s thresholds (%s tree)' % (lb_metric.upper(), plot_str),
                            xlabel='%s threshold (percentile)' % lb_metric.upper(),
                            ylabel='mean %s\nsince affinity increase' % vstr_label)

    # it would be nice to rewrite the above so that the stuff going into the yaml dict looked more like the analagous lines in the affy fcn
    dump_plot_info_to_yaml('%s/%s.yaml' % (plotdir, plotname), {'lb_ptiles' : lb_ptile_vals['lb_ptiles'], 'mean_%s_ptiles' % vstr : lb_ptile_vals['mean_%s' % vstr], 'perfect_vals' : perfect_ptile_vals['mean_%s' % vstr]})

    return plotname

# ----------------------------------------------------------------------------------------
def add_n_ancestor_vals(node, dtree, line, plotvals, affinity_changes, lb_metric, min_affinity_change=1e-6, n_max_steps=15, debug=False):
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
        return

    plotvals['n-ancestor'].append(n_steps)
    plotvals['branch-length'].append(branch_len)
    plotvals[lb_metric].append(line['tree-info']['lb'][lb_metric][node.taxon.label])
    # plotvals['uids'].append(node.taxon.label)

    if debug:
        print '     %12s %12s %8.4f %9.4f  %s%-9.4f' % ('', ancestor_uid, branch_len, chosen_ancestor_affinity, utils.color('red', '+'), this_affinity - chosen_ancestor_affinity)

# ----------------------------------------------------------------------------------------
def plot_lb_vs_ancestral_delta_affinity(plotdir, true_lines, lb_metric, lb_label, plot_str='true', ptile_range_tuple=(50., 100., 1.), only_csv=False, debug=False):
    # plot lb[ir] vs number of ancestors to nearest affinity decrease (well, decrease as you move upwards in the tree/backwards in time)

    vstr_list = collections.OrderedDict([(xstr, xlabel) for metric, xstr, xlabel in lb_metric_axis_stuff if metric == 'lbr'])
    # first plot lb metric vs number of ancestors since last affinity increase (all clusters)
    if debug:
        print '  finding N ancestors since last affinity increase'
        print '         node        ancestors  distance   affinity (%sX: change for chosen ancestor, %s: reached root without finding lower-affinity ancestor)' % (utils.color('red', '+'), utils.color('green', 'x'))
    plotvals = {val_type : [] for val_type in [lb_metric] + vstr_list.keys()}  # , 'uids']}
    for line in true_lines:
        dtree = treeutils.get_dendro_tree(treestr=line['tree'])
        affinity_changes = []
        for uid in line['unique_ids']:
            node = dtree.find_node_with_taxon_label(uid)
            if node is dtree.seed_node:  # root doesn't have any ancestors
                continue
            if lb_metric == 'lbr' and line['tree-info']['lb'][lb_metric][uid] == 0:  # lbr equals 0 should really be treated as None/missing
                continue
            add_n_ancestor_vals(node, dtree, line, plotvals, affinity_changes, lb_metric, debug=debug)

        # make sure affinity changes are all roughly the same size
        affinity_changes = sorted(affinity_changes)
        if debug:
            print '    chosen affinity changes: %s' % ' '.join(['%.4f' % a for a in affinity_changes])
        if len(affinity_changes) == 0:
            continue
        if len([a for a in affinity_changes if a < 0.]):
            print '  %s negative affinity changes in %s' % (utils.color('red', 'error'), ' '.join(['%.4f' % a for a in affinity_changes]))
        max_diff = affinity_changes[-1] - affinity_changes[0]
        if abs(max_diff) / numpy.mean(affinity_changes) > 0.2:
            print'      %s not all affinity increases were the same size (min: %.4f   max: %.4f   abs(diff) / mean: %.4f' % (utils.color('yellow', 'warning'), affinity_changes[0], affinity_changes[-1], abs(max_diff) / numpy.mean(affinity_changes))

    fnames = []
    for vstr, vstr_label in vstr_list.items():
        if not only_csv:
            scatter_plotname = '%s-vs-%s-%s-tree' % (lb_metric, vstr, plot_str)  # 'nearest ancestor with lower affinity' would in some ways be a better xlabel, since it clarifies the note at the top of the loop, but it's also less clear in other ways
            plot_2d_scatter(scatter_plotname, plotdir, plotvals, lb_metric, lb_label, '%s (true tree)' % lb_metric.upper(), xvar=vstr, xlabel='%s since affinity increase' % vstr_label, log='y' if lb_metric == 'lbr' else '')
            fnames.append('%s/%s.svg' % (plotdir, scatter_plotname))

        pname = make_ancestral_ptile_plot(plotvals, vstr, vstr_label, lb_metric, ptile_range_tuple, plotdir, plot_str, only_csv=only_csv, debug=debug)
        fnames.append('%s/%s.svg' % (plotdir, pname))

    return [fnames]

# ----------------------------------------------------------------------------------------
def plot_true_vs_inferred_lb(plotdir, true_lines, inf_lines, lb_metric, lb_label, debug=False):
    plotvals = {val_type : {uid : l['tree-info']['lb'][lb_metric][uid] for l in lines for uid in l['unique_ids']}
                for val_type, lines in (('true', true_lines), ('inf', inf_lines))}
    common_uids = set(plotvals['true']) & set(plotvals['inf'])  # there should/may be a bunch of internal nodes in the simulation lines but not in the inferred lines, but otherwise they should have the same uids
    plotvals = {val_type : [plotvals[val_type][uid] for uid in common_uids] for val_type in plotvals}
    plotname = '%s-true-vs-inferred' % lb_metric
    plot_2d_scatter(plotname, plotdir, plotvals, 'inf', '%s on inferred tree' % lb_metric.upper(), 'true vs inferred %s' % lb_metric.upper(), xvar='true', xlabel='%s on true tree' % lb_metric.upper())
    return ['%s/%s.svg' % (plotdir, plotname)]

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
    cmdstr += ' --log-lbr'
    if tree_style is not None:
        cmdstr += ' --tree-style %s' % tree_style
    cmdstr, _ = utils.run_ete_script(cmdstr, ete_path, return_for_cmdfos=True, tmpdir=subworkdir, extra_str='        ')

    return {'cmd_str' : cmdstr, 'workdir' : subworkdir, 'outfname' : outfname, 'workfnames' : [treefname, metafname]}

# ----------------------------------------------------------------------------------------
def plot_lb_trees(plotdir, lines, ete_path, base_workdir, is_simu=False, tree_style=None):
    workdir = '%s/ete3-plots' % base_workdir
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    cmdfos = []
    for lb_metric, lb_label in treeutils.lb_metrics.items():
        for iclust, line in enumerate(lines):  # note that <min_tree_metric_cluster_size> was already applied in treeutils
            treestr = get_tree_from_line(line, is_simu)
            for affy_key in treeutils.affy_keys[lb_metric]:
                metafo = copy.deepcopy(line['tree-info']['lb'])
                if affy_key in line:  # either 'affinities' or 'relative_affinities'
                    metafo[utils.reversed_input_metafile_keys[affy_key]] = {uid : affy for uid, affy in zip(line['unique_ids'], line[affy_key])}
                outfname = '%s/trees/%s-tree-iclust-%d%s.svg' % (plotdir, lb_metric, iclust, '-relative' if 'relative' in affy_key else '')
                cmdfos += [get_lb_tree_cmd(treestr, outfname, lb_metric, affy_key, ete_path, '%s/sub-%d' % (workdir, len(cmdfos)), metafo=metafo, tree_style=tree_style)]

    start = time.time()
    utils.run_cmds(cmdfos, clean_on_success=True, shell=True) #, debug='print')
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
