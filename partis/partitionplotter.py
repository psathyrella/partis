from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import csv
import os
import shutil
import math
import numpy
import itertools
import sys
import time
import collections
import operator
import copy
import glob
import re

from .hist import Hist
from . import hutils
from . import utils
from .clusterpath import ClusterPath
from . import mds
from . import treeutils
from io import open

# ----------------------------------------------------------------------------------------
class PartitionPlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, args, glfo=None):
        self.args = args
        self.glfo = glfo
        from . import plotting
        self.plotting = sys.modules['partis.plotting']
        from . import lbplotting
        self.lbplotting = sys.modules['partis.lbplotting']

        self.n_clusters_per_joy_plot = 50 if self.args.meta_info_key_to_color is None else 30
        self.n_max_joy_plots = 12
        self.n_max_mutations = 65
        self.n_joyplots_in_html = {'shm-vs-size' : self.n_max_joy_plots, 'overview' : 2}  # the rest of them are still in the dir, they're just not displayed in the html (note this is just
        self.min_high_mutation_cluster_size = 1
        self.n_biggest_to_plot = 25  # this functions as the number of plots for trees, mds, sfs, and laplacian spectra NOTE this is overridden by self.args.queries_to_include (i.e. those queries always get plotted), but *not* by self.args.meta_info_to_emphasize
        self.n_plots_per_row = 4

        self.size_vs_shm_min_cluster_size = 3  # don't plot singletons and pairs for really big repertoires
        self.min_pairwise_cluster_size = 3  # note that singletons are meaningless for pairwise diversity0
        self.min_tree_cluster_size = 5  # we also use min_selection_metric_cluster_size
        # self.n_max_tree_plots = 30  # don't put more than this many in the overview html (the rest will still be in trees.html)
        self.n_max_bubbles = 100  # circlify is really slow
        self.mds_max_cluster_size = 50000  # it's way tf too slow NOTE also max_cluster_size in make_mds_plots() (should combine them or at least put them in the same place)
        self.min_mds_cluster_size = 3
        self.laplacian_spectra_min_clusters_size = 4
        self.min_n_clusters_to_apply_size_vs_shm_min_cluster_size = 350  # don't apply the previous thing unless the repertoire's actually pretty large

        self.n_mds_components = 2
        self.indexing = 1  # for mutation plotting labels
        self.n_max_alternatives = None  # 3 # just for testing, this lets you limit the number of alternative/sampled trees that we loop over

        self.sclusts, self.antn_dict, self.treefos, self.mut_info, self.base_plotdir = None, None, None, None, None

    # ----------------------------------------------------------------------------------------
    def init_subd(self, subd, allow_other_files=False):
        plotdir = self.base_plotdir + '/' + subd
        utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg', '*.png'], allow_other_files=allow_other_files)
        return subd, plotdir

    # ----------------------------------------------------------------------------------------
    def get_cdr3_title(self, annotation):
        naive_cdr3_seq, _ = utils.subset_iseq(annotation, 0, restrict_to_region='cdr3')
        title = ''
        if len(naive_cdr3_seq) % 3 != 0:
            # print '  out of frame: adding %s' % ((3 - len(naive_cdr3_seq) % 3) * 'N')
            naive_cdr3_seq += (3 - len(naive_cdr3_seq) % 3) * 'N'
            title += ' (out of frame)'
        title = utils.ltranslate(naive_cdr3_seq) + title
        return title

    # ----------------------------------------------------------------------------------------
    # colorfcn should return the *value* (not color), which then gets mapped to a color, e.g. viridis (note colorfcn is overridden by self.args.meta_info_key_to_color)
    def make_cluster_scatter(self, plotdir, plotname, xfcn, yfcn, clusters_to_use, repertoire_size, min_csize, dxy=0.1, log_cluster_size=False, xlabel='?', xbounds=None, repfrac_ylabel=False, colorfcn=None, leg_title=None,
                             alpha=0.65, title=''):
        # ----------------------------------------------------------------------------------------
        def add_emph_clusters(xvals, yvals):
            for iclust, cluster in enumerate(clusters_to_use):
                tqtis = set()  # queries to emphasize in this cluster
                if self.args.queries_to_include is not None and not self.args.dont_label_queries_to_include:
                    tqtis |= set(cluster) & set(self.args.queries_to_include)
                if self.args.meta_info_to_emphasize is not None and self.any_meta_emph(cluster):
                    key, val = list(self.args.meta_info_to_emphasize.items())[0]
                    tqtis.add(utils.meta_emph_str(key, val, formats=self.args.meta_emph_formats))
                if len(tqtis) == 0:
                    continue
                xval, yval = xvals[iclust], yvals[iclust]
                ax.plot([xval], [yval], color='red', marker='.', markersize=2)
                ax.text(xval + dxy, yval + (0.01 if log_cluster_size else 0.1), ' '.join(tqtis), color='red', fontsize=8)

        # ----------------------------------------------------------------------------------------
        skipped_small_clusters = False
        if len(clusters_to_use) > self.min_n_clusters_to_apply_size_vs_shm_min_cluster_size:  # if repertoire is really big, ignore smaller clusters to keep the plots from being huge
            if any(len(c) < min_csize for c in clusters_to_use):
                clusters_to_use = [cluster for cluster in clusters_to_use if len(cluster) >= min_csize]
                skipped_small_clusters = True
        if len(clusters_to_use) == 0:
            print('  %s no clusters to plot for cluster size scatter' % utils.color('yellow', 'warning'))
            return
        xvals, yvals = zip(*[[xfcn(cluster), yfcn(cluster)] for cluster in clusters_to_use])
        original_yvals = yvals
        if log_cluster_size:
            yvals = [math.log(yv) for yv in yvals]
        colors = None
        if self.args.meta_info_key_to_color is None and colorfcn is not None:  # this is mostly copied from lbplotting.plot_2d_scatter()
            colorvals = [colorfcn(c) for c in clusters_to_use]  # values (numerical/string) that will tell us which color to use
            smap = self.plotting.get_normalized_scalar_map(colorvals, 'viridis')
            smfcn = lambda x: 'grey' if x is None else self.plotting.get_smap_color(smap, None, val=x)
            colors = [smfcn(c) for c in colorvals]
            leg_entries = self.plotting.get_leg_entries(vals=colorvals, colorfcn=smfcn)
            self.plotting.plot_legend_only(leg_entries, plotdir, plotname+'-legend', title=leg_title, n_digits=2)
        fig, ax = self.plotting.mpl_init()
        import matplotlib.pyplot as plt
        if len(clusters_to_use) >= self.min_n_clusters_to_apply_size_vs_shm_min_cluster_size:
            hb = ax.hexbin(xvals, yvals, gridsize=8, cmap=plt.cm.Greys) #, alpha=alpha)
        if self.args.meta_info_key_to_color is None:
            sct_xvals, sct_yvals = xvals, yvals
        else:
            total_bubble_seqs, fake_cluster, bubfos = self.get_cluster_bubfos(clusters_to_use, min_cluster_size=min_csize, dont_sort=True)
            bfos_to_plot, sct_xvals, sct_yvals = [], [], []  # scatter x/y vals
            for bfo in [b for b in bubfos if b['id']!='fake']:
                iclust = int(bfo['id'])
                tclust = clusters_to_use[iclust]
                assert len(tclust) == bfo['radius']  # make sure we have the right cluster
                bkg_frac = None
                if self.args.meta_info_bkg_vals is not None and any(f['label'] in self.args.meta_info_bkg_vals for f in bfo['fracs']):
                    bkg_ffos = [f for f in bfo['fracs'] if f['label'] in self.args.meta_info_bkg_vals]
                    bkg_frac = sum(f['fraction'] for f in bkg_ffos)
                if bkg_frac == 1:
                    sct_xvals.append(xvals[iclust])
                    sct_yvals.append(yvals[iclust])
                else:
                    bfos_to_plot.append(bfo)  # have to plot them *after* making the scatter plot, so they end up on top
        ax.scatter(sct_xvals, sct_yvals, s=(360 * 0.01)**2, marker='.', color='grey', alpha=0.4)  # colors  # eh, turning off colors for now
        if self.args.meta_info_key_to_color is not None:
            for bfo in bfos_to_plot:
                iclust = int(bfo['id'])
                self.plotting.plot_pie_chart_marker(ax, xvals[iclust], yvals[iclust], 0.03, bfo['fracs'], alpha=alpha)

        ymin, ymax = min(yvals), max(yvals)
        if yfcn == len:  # special stuff if y axis is cluster size
            ymin = min(ymin, math.log(1) if log_cluster_size else 1)  # make ymin 1, even if we aren't plotting small clusters, to make it more obvious that we skipped them
            yticks, _ = self.plotting.get_cluster_size_xticks(xmin=min([1] + list(original_yvals)), xmax=max(original_yvals))
            if log_cluster_size:
                yticks = [math.log(y) for y in yticks]
        else:
            nticks = 5
            yticks = [ymax + itick * (ymin - ymax) / float(nticks - 1) for itick in range(nticks)]
        ytlfcn = (lambda yt: utils.get_repfracstr(yt, repertoire_size)) if repfrac_ylabel else (lambda yt: '%.0f'%yt)
        yticklabels = [ytlfcn(math.exp(yt) if log_cluster_size else yt) for yt in yticks]

        # if necessary, add red labels to clusters
        if self.args.queries_to_include is not None or self.args.meta_info_to_emphasize is not None:
            add_emph_clusters(xvals, yvals)

        ylabel = ('family size\n(frac. of %d)' % repertoire_size) if repfrac_ylabel else 'clonal family size'
        if log_cluster_size:
            ylabel = '(log) ' + ylabel
            plotname += '-log'
        if skipped_small_clusters:
            fig.text(0.8, 0.25, 'skipping clusters\nsmaller than %d' % min_csize, color='green', fontsize=8)
        self.plotting.mpl_finish(ax, plotdir, plotname, xlabel=xlabel, ylabel=ylabel, xbounds=xbounds, ybounds=(ymin - 0.03*(ymax-ymin), 1.05 * ymax), yticks=yticks, yticklabels=yticklabels, title=title, leg_title=leg_title)

    # ----------------------------------------------------------------------------------------
    def make_single_hexbin_size_vs_shm_plot(self, repertoire_size, plotdir, plotname, log_cluster_size=False, debug=False):  # NOTE not using <repertoire_size> any more, but don't remember if there was a reason I should leave it
        def mean_mutations(cluster):
            return numpy.mean(self.antn_dict[':'.join(cluster)]['n_mutations'])

        clusters_to_use = [cluster for cluster in self.sclusts if mean_mutations(cluster) < self.n_max_mutations]  # have to do it as a separate line so the zip/* don't crash if no clusters pass the criterion
        self.make_cluster_scatter(plotdir, plotname, mean_mutations, len, clusters_to_use, repertoire_size, self.size_vs_shm_min_cluster_size, log_cluster_size=log_cluster_size, xlabel='mean N mutations', xbounds=(0, self.n_max_mutations))

    # ----------------------------------------------------------------------------------------
    def addfname(self, fnames, fname, force_new_row=False, suffix='.svg'):
        if fname is None:
            return
        fname += suffix
        if force_new_row or len(fnames[-1]) >= self.n_plots_per_row:
            fnames.append([fname])
        else:
            fnames[-1].append(fname)

    # ----------------------------------------------------------------------------------------
    def meta_emph(self, cluster, uid):  # return True if <uid> from <cluster> satisfies criteria in self.args.meta_info_to_emphasize
        antn = self.antn_dict[':'.join(cluster)]
        key, val = list(self.args.meta_info_to_emphasize.items())[0]
        return key in antn and utils.meta_info_equal(key, val, utils.per_seq_val(antn, key, uid), formats=self.args.meta_emph_formats)

    # ----------------------------------------------------------------------------------------
    def any_meta_emph(self, cluster):
        if self.args.meta_info_to_emphasize is None:
            return False
        return any(self.meta_emph(cluster, u) for u in cluster)

    # ----------------------------------------------------------------------------------------
    def plot_this_cluster(self, iclust, plottype=None):
        if len(self.sclusts[iclust]) == 1:
            return False
        if plottype == 'mds' and len(self.sclusts[iclust]) > self.mds_max_cluster_size:
            print('     skipping mds plots for cluster with size %d > %d' % (len(self.sclusts[iclust]), self.mds_max_cluster_size))
            return False
        if plottype == 'trees' and len(self.sclusts[iclust]) < self.args.min_selection_metric_cluster_size:
            return False
        if self.args.cluster_indices is not None and iclust not in self.args.cluster_indices:
            return False
        if iclust < self.n_biggest_to_plot:
            return True
        if self.args.queries_to_include is not None and len(set(self.args.queries_to_include) & set(self.sclusts[iclust])) > 0:  # seed is added to <args.queries_to_include> in bin/partis
            return True
        if plottype == 'trees' and len(self.sclusts[iclust]) >= self.args.min_selection_metric_cluster_size:  # ok, ugh, it sucks to have this twice in both directions
            return True
        return False  # falls through if <iclust> is too big, or if there's no --queries-to-include (which includes the seed)

    # ----------------------------------------------------------------------------------------
    def make_shm_vs_cluster_size_plots(self, no_slug=False, debug=False):
        def get_fname(iclustergroup=None, high_mutation=False, hexbin=False):
            if iclustergroup is not None:  # index of this group of clusters
                return 'size-vs-shm-%d' % iclustergroup
            elif high_mutation:
                return 'size-vs-shm-high-mutation'
            elif hexbin:
                return 'size-vs-shm-hexbin'
            else:
                assert False
        subd, plotdir = self.init_subd('shm-vs-size')

        repertoire_size = sum([len(c) for c in self.sclusts])
        cluster_indices = {':'.join(self.sclusts[i]) : i for i in range(len(self.sclusts))}  # index over all clusters, in the order that the mds plots will appear (compare to the two other indices I need within plotting.make_single_joyplot())

        # size vs shm joy plots
        iclustergroup = 0
        fnd = {'joy' : [], 'hex' : []}
        high_mutation_clusters = []
        sorted_cluster_groups = [self.sclusts[i : i + self.n_clusters_per_joy_plot] for i in range(0, len(self.sclusts), self.n_clusters_per_joy_plot)]
        if debug:
            print('  shm vs size joyplots: divided repertoire of size %d with %d clusters into %d cluster groups' % (repertoire_size, len(self.sclusts), len(sorted_cluster_groups)))
        all_emph_vals, emph_colors = None, None
        if self.args.meta_info_key_to_color is not None:  # have to do this out here before the loop so that the colors are synchronized (and all plots include all possible values)
            all_emph_vals, emph_colors = self.plotting.meta_emph_init(self.args.meta_info_key_to_color, clusters=self.sclusts, antn_dict=self.antn_dict, formats=self.args.meta_emph_formats)
        for subclusters in sorted_cluster_groups:
            if iclustergroup > self.n_max_joy_plots:  # note that when this is activated, the high mutation plot is no longer guaranteed to have every high mutation cluster (but it should have every high mutation cluster that was bigger than the cluster size when we started skipping here)
                continue
            if no_slug:
                continue
            if debug:
                print('    %d: making joyplot with %d clusters' % (iclustergroup, len(subclusters)))
            title = 'per-family SHM (%d / %d)' % (iclustergroup + 1, len(sorted_cluster_groups))  # NOTE it's important that this denominator is still right even when we don't make plots for all the clusters (which it is, now)
            high_mutation_clusters += self.plotting.make_single_joyplot(subclusters, self.antn_dict, repertoire_size, plotdir, get_fname(iclustergroup=iclustergroup), cluster_indices=cluster_indices, title=title, high_x_val=self.n_max_mutations,
                                                                        queries_to_include=self.args.queries_to_include, meta_info_to_emphasize=self.args.meta_info_to_emphasize, meta_info_key_to_color=self.args.meta_info_key_to_color, meta_emph_formats=self.args.meta_emph_formats, all_emph_vals=all_emph_vals, emph_colors=emph_colors,
                                                                        make_legend=self.args.meta_info_key_to_color is not None, dont_label_queries_to_include=self.args.dont_label_queries_to_include, debug=debug)  # add_clone_id=True, add_index=True have to make legend for every plot
            if len(fnd['joy']) < self.n_joyplots_in_html['shm-vs-size']:
                fnd['joy'].append(get_fname(iclustergroup=iclustergroup))
            iclustergroup += 1
        if self.args.meta_info_key_to_color is not None:
            fnd['leg'] = [get_fname(iclustergroup=0)+'-legend']
        if len(high_mutation_clusters) > 0 and len(high_mutation_clusters[0]) > self.min_high_mutation_cluster_size:
            high_mutation_clusters = [cluster for cluster in high_mutation_clusters if len(cluster) > self.min_high_mutation_cluster_size]
            self.plotting.make_single_joyplot(high_mutation_clusters, self.antn_dict, repertoire_size, plotdir, get_fname(high_mutation=True), plot_high_x=True, cluster_indices=cluster_indices, title='families with mean > %d mutations' % self.n_max_mutations,
                                              high_x_val=self.n_max_mutations, queries_to_include=self.args.queries_to_include, meta_info_to_emphasize=self.args.meta_info_to_emphasize, meta_info_key_to_color=self.args.meta_info_key_to_color, meta_emph_formats=self.args.meta_emph_formats, all_emph_vals=all_emph_vals, emph_colors=emph_colors,
                                              make_legend=self.args.meta_info_key_to_color is not None, dont_label_queries_to_include=self.args.dont_label_queries_to_include, debug=debug) # add_clone_id=True, add_index=True
            fnd['high'] = [get_fname(high_mutation=True)]

        # size vs shm hexbin plots
        self.make_single_hexbin_size_vs_shm_plot(repertoire_size, plotdir, get_fname(hexbin=True))
        fnd['hex'].append(get_fname(hexbin=True))
        self.make_single_hexbin_size_vs_shm_plot(repertoire_size, plotdir, get_fname(hexbin=True), log_cluster_size=True)
        fnd['hex'].append(get_fname(hexbin=True) + '-log')

        fnames, rfnames = [[]], [[]]
        if 'leg' in fnd:
            self.addfname(fnames, fnd['leg'][0])
            self.addfname(rfnames, subd + '/' + fnd['leg'][0])
        for ifn, fn in enumerate(fnd['joy']):
            self.addfname(fnames, fn)
            if ifn < self.n_joyplots_in_html['overview']:
                self.addfname(rfnames, subd + '/' + fn)
        if 'high' in fnd:
            self.addfname(fnames, fnd['high'][0])
            self.addfname(rfnames, subd + '/' + fnd['high'][0])
        for ifn, fn in enumerate(fnd['hex']):
            self.addfname(fnames, fn, force_new_row=ifn==0)
            self.addfname(rfnames, subd + '/' + fn, force_new_row=ifn==0)

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames, new_table_each_row=True, extra_links=[('all shm-vs-size plots', subd)])

        return rfnames

    # ----------------------------------------------------------------------------------------
    def get_cluster_bubfos(self, clusters_to_use, n_to_write_size=9999999, n_max_bubbles=None, min_cluster_size=None, dont_sort=False):
        mekey = self.args.meta_info_key_to_color
        fake_cluster, fake_antn = [], {'unique_ids' : [], mekey : []} # make a fake cluster with all sequences from all skipped clusters (circlify is too slow to run on all the smaller clusters)
        bubfos = []
        for iclust, cluster in enumerate(clusters_to_use):
            if n_max_bubbles is not None and iclust < n_max_bubbles or min_cluster_size is not None and len(cluster) >= min_cluster_size:
                bubfos.append({'id' : str(iclust), 'radius' : len(cluster)})
            else:
                fake_cluster += cluster
                fake_antn['unique_ids'] += cluster
                if mekey is not None:
                    antn = self.antn_dict.get(':'.join(cluster))
                    fake_antn[mekey] += antn[mekey] if antn is not None and mekey in antn else [None for _ in cluster]
        if len(fake_cluster) > 0:
            bubfos.append({'id' : 'fake', 'radius' : len(fake_cluster)})

        if mekey is not None:
            all_emph_vals, emph_colors = self.plotting.meta_emph_init(mekey, clusters=clusters_to_use, antn_dict=self.antn_dict, formats=self.args.meta_emph_formats)
            hcolors = {v : c for v, c in emph_colors}
        def getclust(idl): return clusters_to_use[int(idl)] if idl!='fake' else fake_cluster
        if not dont_sort:
            bubfos = sorted(bubfos, key=lambda b: len(getclust(b['id'])), reverse=True)
        for bfo in bubfos:
            cluster = getclust(bfo['id'])
            if bfo['id']=='fake' or int(bfo['id']) < n_to_write_size:
                tstr, fsize, tcol, dx = ('small clusters', 12, 'red', -0.3) if bfo['id']=='fake' else (len(cluster), 5, 'black', -0.04)
                bfo['texts'] = [{'tstr' : tstr, 'fsize' : fsize, 'tcol' : tcol, 'dx' : dx}]
            antn = self.antn_dict.get(':'.join(cluster)) if bfo['id'] != 'fake' else fake_antn
            if antn is None or mekey is None:
                bfo['fracs'] = None
            else:
                emph_fracs = utils.get_meta_emph_fractions(mekey, all_emph_vals, cluster, antn, formats=self.args.meta_emph_formats)
                bfo['fracs'] = [{'label' : v, 'fraction' : f, 'color' : hcolors[v]} for v, f in emph_fracs.items()]
            if bfo['id']!='fake' and self.any_meta_emph(cluster):
                mkey, mval = list(self.args.meta_info_to_emphasize.items())[0]
                bfo['texts'].append({'tstr' : utils.meta_emph_str(mkey, mval, formats=self.args.meta_emph_formats), 'fsize' : 5, 'tcol' : 'red', 'dx' : -0.04, 'dy' : -0.03})
        total_bubble_seqs = sum(len(getclust(b['id'])) for b in bubfos if b['id']!='fake')
        return total_bubble_seqs, fake_cluster, bubfos

    # ----------------------------------------------------------------------------------------
    def make_cluster_bubble_plots(self, alpha=0.6, debug=False):
        import matplotlib.pyplot as plt
        subd, plotdir = self.init_subd('cluster-bubble')
        total_bubble_seqs, fake_cluster, bubfos = self.get_cluster_bubfos(self.sclusts, n_max_bubbles=self.n_max_bubbles)
        nbub = len(bubfos)
        if len(fake_cluster) > 0:
            nbub -= 1
        xtra_text = None
        if len(fake_cluster) > 0:
            xtra_text = {'x' : 0.3, 'y' : 0.85, 'color' : 'green', 'text' : '%d seqs in %d clusters with size %d or smaller' % (len(fake_cluster), len(self.sclusts) - nbub, len(self.sclusts[self.n_max_bubbles - 1]))}
        repertoire_size = sum([len(c) for c in self.sclusts])
        fn = self.plotting.bubble_plot('bubbles', plotdir, [b for b in bubfos if b['id']!='fake'], title='bubbles for %d largest clusters (%d/%d seqs)'%(nbub, total_bubble_seqs, repertoire_size), alpha=alpha)
        fn2 = self.plotting.bubble_plot('small-bubbles', plotdir, [b for b in bubfos if b['id']=='fake'], title='single bubble for %d smaller clusters'%(len(self.sclusts) - self.n_max_bubbles), xtra_text=xtra_text, alpha=alpha)
        fnames = [[]]
        self.addfname(fnames, fn)
        self.addfname(fnames, fn2)
        if self.args.meta_info_key_to_color is not None:
            all_emph_vals, emph_colors = self.plotting.meta_emph_init(self.args.meta_info_key_to_color, clusters=self.sclusts, antn_dict=self.antn_dict, formats=self.args.meta_emph_formats)
            lfn = self.plotting.make_meta_info_legend(plotdir, 'cluster-bubbles', self.args.meta_info_key_to_color, emph_colors, all_emph_vals, meta_emph_formats=self.args.meta_emph_formats, alpha=alpha)
            self.addfname(fnames, lfn)

        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_pairwise_diversity_plots(self, individual_plots=True, n_max_seqs=50):
        def antn(c): return self.antn_dict.get(':'.join(c))
        subd, plotdir = self.init_subd('diversity')

        repertoire_size = sum([len(c) for c in self.sclusts])
        for cluster in self.sclusts:
            utils.add_seqs_aa(antn(cluster))

        mekey = self.args.meta_info_key_to_color
        if mekey is not None:
            all_emph_vals, emph_colors = self.plotting.meta_emph_init(mekey, clusters=self.sclusts, antn_dict=self.antn_dict, formats=self.args.meta_emph_formats, lgroups=self.plotting.legend_groups)
            hcolors = {v : c for v, c in emph_colors}

        fnames = [[]]
        subfns = [[]] if individual_plots else None
        for tstr in ['nuc']: #, 'aa']:
            skey = '' if tstr=='nuc' else '_aa'
            def pwfcn(tclust, slist=None):
                if slist is None:
                    slist = antn(tclust)['seqs%s'%skey]
                return utils.mean_pairwise_hdist(slist, amino_acid=tstr=='aa', n_max_seqs=n_max_seqs)

            # 2d scatter plots with diversity vs mutation/size:
            def mean_mutations(c): return numpy.mean(antn(c)['n_mutations'] if tstr=='nuc' else [utils.shm_aa(antn(c), iseq=i) for i in range(len(antn(c)['unique_ids']))])

            start = time.time()
            fname = 'mean-pairwise-hdist-%s' % tstr
            self.make_cluster_scatter(plotdir, fname, pwfcn, len, self.sclusts, repertoire_size, self.min_pairwise_cluster_size, dxy=0.003, log_cluster_size=True,
                                      xlabel='mean pairwise %s distance'%tstr, colorfcn=mean_mutations, title='pairwise %s diversity'%tstr, leg_title='N %s muts'%tstr) #, xbounds=(0, self.n_max_mutations))
            for fstr in ['legend', 'log']:
                if os.path.exists('%s/%s-%s.svg'%(plotdir, fname, fstr)):
                    self.lbplotting.add_fn(fnames, fn='%s/%s-%s.svg' % (subd, fname, fstr))

            if 'seaborn' not in sys.modules:
                import seaborn
            sns = sys.modules['seaborn']
            thist = self.plotting.make_csize_hist(self.sclusts, xbins=self.args.cluster_size_bins)
            xticklabels = thist.bin_labels
            ytitle = 'mean pairwise %s distance' % tstr
            def hblist(): return [[] for _ in thist.bin_contents]
            if mekey is None:
                plotvals = hblist()
                for iclust, tclust in enumerate(self.sclusts):
                    plotvals[thist.find_bin(len(tclust))].append(pwfcn(tclust))
                ax = sns.boxplot(data=plotvals, color='grey')
                ax.set_xticklabels(xticklabels, rotation='vertical', size=15)
                fn = self.plotting.mpl_finish(ax, plotdir, 'diversity-box-%s'%tstr, xlabel='family size', ylabel=ytitle) #, yticks=xticks  # , xticklabels=xticklabels
                self.lbplotting.add_fn(fnames, fn=fn)
            else:
                plotvals = {v : hblist() for v in all_emph_vals}
                for iclust, tclust in enumerate(self.sclusts):
                    def psfcn(u): return utils.meta_emph_str(mekey, utils.per_seq_val(antn(tclust), mekey, u, use_default=True), formats=self.args.meta_emph_formats, legend_groups=self.plotting.legend_groups)
                    # NOTE that the scatter plots above use the *whole* cluster, incorporating all seqs, whereas here we split the cluster apart by meta info key (so e.g. separate [sub-]clusters for pbmc and skin)
                    for v_emph, vgroup in utils.group_seqs_by_value(tclust, psfcn, return_values=True):
                        vgroup = list(vgroup)  # some fcns i might pass it to might iterate more than once
                        pdist = pwfcn(None, slist=[antn(tclust)['seqs%s'%skey][antn(tclust)['unique_ids'].index(u)] for u in vgroup])
                        plotvals[v_emph][thist.find_bin(len(tclust))].append(pdist)
                median_pvals = {}
                for v_emph in sorted(all_emph_vals):
                    mvals = [None if len(pvals)==0 else numpy.median(pvals) for pvals in plotvals[v_emph]]  # average over clusters in this bin
                    median_pvals[v_emph] = mvals
                    ax = sns.boxplot(data=plotvals[v_emph], color=hcolors[v_emph], boxprops={'alpha' : 0.4}, whiskerprops={'color' : hcolors[v_emph]}, flierprops={'markerfacecolor' : hcolors[v_emph]}, medianprops={'color' : hcolors[v_emph]})
                    # ax = sns.stripplot(data=plotvals[v_emph], color=hcolors[v_emph], size=10, alpha=0.4)
                    ax.set_xticklabels(xticklabels, rotation='vertical', size=15)
                    if individual_plots:
                        fn = self.plotting.mpl_finish(ax, plotdir, 'diversity-box-%s-%s'%(tstr, v_emph.replace('/', '_slash_')), xlabel='family size', ylabel=ytitle, title='%s-only seqs'%v_emph) #, yticks=xticks  # , xticklabels=xticklabels
                        self.lbplotting.add_fn(fnames, fn=fn)
                # # uncomment this to get the swarm/box plots for each color overlaid on each other (as opposed to in separate files)
                # fn = self.plotting.mpl_finish(ax, plotdir, 'diversity-box-%s'%tstr, xlabel='family size', ylabel=ytitle) #, yticks=xticks  # , xticklabels=xticklabels
                # self.lbplotting.add_fn(fnames, fn=fn)
                # fn = self.plotting.stack_meta_hists('diveristy-stacked-%s'%tstr, plotdir, mekey, median_pvals, template_hist=thist, is_bin_contents=True, log='x', figsize=(8, 5), ytitle='%s\n(median over clusters)'%ytitle, rotation='vertical', plottitle='diversity (within seqs with each %s)'%mekey, remove_empty_bins=False)
                # self.lbplotting.add_fn(fnames, fn=fn)
                # if individual_plots:
                #     self.plotting.make_html(plotdir, fnames=subfns, extra_links=[('individual diversity plots', subd)], bgcolor='#FFFFFF')  # , new_table_each_row=True
                fn = self.plotting.make_meta_info_legend(plotdir, 'diversity', mekey, emph_colors, all_emph_vals, meta_emph_formats=self.args.meta_emph_formats, alpha=0.4)
                self.lbplotting.add_fn(fnames, fn='%s/%s.svg'%(subd, fn))

            print('    calculated %s pairwise diversity on %d clusters%s (%.1f sec)' % (tstr, len(self.sclusts), '' if n_max_seqs is None else ', subsampling each cluster to %d seqs'%n_max_seqs, time.time() - start))

        return fnames

    # ----------------------------------------------------------------------------------------
    def make_mds_plots(self, max_cluster_size=10000, reco_info=None, color_rule=None, run_in_parallel=False, aa=False, debug=False):
        debug = True
        # ----------------------------------------------------------------------------------------
        def get_fname(ic):
            return 'icluster-%d' % ic
        # ----------------------------------------------------------------------------------------
        def get_cluster_info(full_cluster, iclust):
            # ----------------------------------------------------------------------------------------
            def addseq(name, seq):
                found = False
                for sfo in seqfos:  # mds barfs if we have duplicate sequences, so if the sequence is already in there with a different name we just rename it (ick)
                    if sfo['seq'] == seq:
                        found = True
                        if sfo['name'] in tqtis:
                            name += '@(%s)' % tqtis[sfo['name']].lstrip('_')
                        sfo['name'] = name
                        break
                if not found:
                    seqfos.append({'name' : name, 'seq' : seq})
                color_scale_vals[name] = 0
                tqtis[name] = name  # i think if we're calling this fcn we always want to label with its actual name
            # ----------------------------------------------------------------------------------------
            def addqtis(tqtis, new_queries):  # NOTE that this can't do anything to prevent two points being on top of each other
                for uid, val in new_queries.items():
                    if uid in tqtis:
                        tqtis[uid] += '@%s' % val  # doesn't nicely condense duplicates like the joy plots do, but they should be much much less likely here (also note it can't use ',' since it gets passed on the command line)
                    else:
                        tqtis[uid] = '_%s' % val
            # ----------------------------------------------------------------------------------------
            full_info = self.antn_dict[':'.join(full_cluster)]
            # title = '%s   (size: %d)' % (self.get_cdr3_title(full_info), len(full_cluster))
            title = '%sMDS (index %d size %d)' % ('AA ' if aa else '', iclust, len(full_cluster))
            leg_title = 'N %s mutations' % ('AA' if aa else 'nuc')
            tstr = '_aa' if aa else ''

            kept_indices = list(range(len(full_cluster)))
            if len(kept_indices) > max_cluster_size:
                uids_to_choose_from = set([full_cluster[i] for i in kept_indices])  # note similarity to code in seqfileopener.post_process()
                if self.args.queries_to_include is not None:
                    uids_to_choose_from -= set(self.args.queries_to_include)
                if self.args.meta_info_to_emphasize is not None:
                    uids_to_choose_from -= set([u for u in uids_to_choose_from if self.meta_emph(full_cluster, u)])
                n_to_remove = len(kept_indices) - max_cluster_size
                if n_to_remove >= len(uids_to_choose_from):  # i.e. if we'd have to start removing queries that are in <self.args.queries_to_include>
                    removed_uids = uids_to_choose_from
                else:
                    removed_uids = numpy.random.choice(list(uids_to_choose_from), n_to_remove, replace=False)  # i think this'll still crash if len(uids_to_choose_from) is zero, but, meh
                kept_indices = sorted(set(kept_indices) - set([full_cluster.index(uid) for uid in removed_uids]))
                title += ' (subset: %d / %d)' % (len(kept_indices), len(full_cluster))

            if aa:
                utils.add_seqs_aa(full_info)
                utils.add_naive_seq_aa(full_info)
            seqfos = [{'name' : full_info['unique_ids'][iseq], 'seq' : full_info['seqs%s'%tstr][iseq]} for iseq in kept_indices]
            color_scale_vals = {full_cluster[iseq] : (utils.shm_aa(full_info, iseq) if aa else full_info['n_mutations'][iseq]) for iseq in kept_indices}

            tqtis = {}
            if self.args.queries_to_include is not None:
                addqtis(tqtis, {u : u for u in self.args.queries_to_include})
            if self.args.meta_info_to_emphasize is not None:
                key, val = list(self.args.meta_info_to_emphasize.items())[0]
                addqtis(tqtis, {full_cluster[i] : utils.meta_emph_str(key, val, formats=self.args.meta_emph_formats) for i in kept_indices if self.meta_emph(full_cluster, full_cluster[i])})  # leading '_' is so dot doesn't cover up label
            addseq('_naive', full_info['naive_seq%s'%tstr])  # note that if any naive sequences that were removed above are in self.args.queries_to_include, they won't be labeled in the plot (but, screw it, who's going to ask to specifically label a sequence that's already specifically labeled?)
            addseq('_consensus', utils.cons_seq_of_line(full_info, aa=aa))  # leading underscore is 'cause the mds will crash if there's another sequence with the same name, and e.g. christian's simulation spits out the naive sequence with name 'naive'. No, this is not a good long term fix

            # remove duplicates, since they crash mds
            new_seqfos, all_seqs = [], {}
            for sfo in seqfos:
                # ----------------------------------------------------------------------------------------
                def add_to_label(nid, oid):  # add label for new id <nid> to any label for old id <oid>
                    if oid in tqtis:
                        tqtis[oid] += '@'
                        tqtis[nid] = tqtis[nid].lstrip('_')
                    else:
                        tqtis[oid] = ''
                    tqtis[oid] += tqtis[nid]
                    del tqtis[nid]
                # ----------------------------------------------------------------------------------------
                if sfo['seq'] in all_seqs:  # if we already added this seq with a different uid, we skip it
                    if sfo['name'] in tqtis:  # *and* if we need this seq to be labeled, then we have to add its label under the previous/other uid
                        add_to_label(sfo['name'], all_seqs[sfo['seq']])
                    continue
                new_seqfos.append(sfo)
                all_seqs[sfo['seq']] = sfo['name']
            seqfos = new_seqfos

            return seqfos, color_scale_vals, tqtis, title, leg_title

        # ----------------------------------------------------------------------------------------
        def get_labels_for_coloring(full_cluster, color_rule):
            full_info = self.antn_dict[':'.join(full_cluster)]
            if color_rule == 'nearest-target':  # color by the index of the nearest cluster index (bcr-phylo simulation only)
                if 'target_seqs' not in reco_info[full_cluster[0]]:
                    return
                labels = {uid : str(reco_info[uid]['nearest_target_indices'][0]) for uid in full_cluster}
                labels['_naive'] = 'foop'
            elif color_rule == 'wtf':
                labels = {uid : uid.split('@')[1] for uid in full_cluster}
                labels['_naive'] = 'foop'
            else:
                assert False

            return labels

        # ----------------------------------------------------------------------------------------
        def prep_cmdfo(iclust, seqfos, tqtis, color_scale_vals, title, leg_title):
            subworkdir = '%s/mds-%d' % (self.args.workdir, iclust)
            utils.prep_dir(subworkdir)
            tmpfname = '%s/seqs.fa' % subworkdir
            with open(tmpfname, 'w') as tmpfile:
                for sfo in seqfos:
                    csval = None
                    if sfo['name'] in color_scale_vals:
                        csval = color_scale_vals[sfo['name']]
                    tmpfile.write('>%s%s\n%s\n' % (sfo['name'], (' %d' % csval) if csval is not None else '' , sfo['seq']))
            cmdstr = '%s/bin/mds-run.py %s --aligned --plotdir %s --plotname %s --workdir %s --seed %d' % (utils.get_partis_dir(), tmpfname, plotdir, get_fname(iclust), subworkdir, self.args.random_seed)
            if tqtis is not None:
                cmdstr += ' --queries-to-include %s' % ':'.join(','.join([u, l]) for u, l in tqtis.items())
            if title is not None:
                cmdstr += ' --title=%s' % title.replace(' ', '@')
            if leg_title is not None:
                cmdstr += ' --leg-title=%s' % leg_title.replace(' ', '@')
            return {'cmd_str' : cmdstr, 'workdir' : subworkdir, 'outfname' : '%s/%s.svg' % (plotdir, get_fname(iclust)), 'workfnames' : [tmpfname]}

        # ----------------------------------------------------------------------------------------
        subd, plotdir = self.init_subd('mds')

        start = time.time()
        if debug:
            if not run_in_parallel:
                print('    making mds plots starting with %d clusters' % len(self.sclusts))
                print('       size (+naive)   mds    plot   total')
        plotted_cluster_lengths, skipped_cluster_lengths = [], []
        fnames = [[]]
        cmdfos = []
        self.addfname(fnames, 'mds-legend')
        for iclust in range(len(self.sclusts)):
            if not self.plot_this_cluster(iclust, plottype='mds'):
                skipped_cluster_lengths.append(len(self.sclusts[iclust]))
                continue
            plotted_cluster_lengths.append(len(self.sclusts[iclust]))

            seqfos, color_scale_vals, tqtis, title, leg_title = get_cluster_info(self.sclusts[iclust], iclust)
            if len(seqfos) < self.min_mds_cluster_size:
                print('    %s cluster index %d too small (%d < %d) after e.g. removing duplicates, skipping' % (utils.wrnstr(), iclust, len(seqfos), self.min_mds_cluster_size))
                continue

            labels = None
            if color_rule is not None:
                labels = get_labels_for_coloring(self.sclusts[iclust], color_rule)
                # print '   %s setting color_scale_vals to None so we can use colors for nearest target seq index' % utils.color('red', 'note')
                color_scale_vals = None  # not sure this is really the best way to do this

            if debug and not run_in_parallel:
                substart = time.time()
                subset_str = '' if len(self.sclusts[iclust]) <= max_cluster_size else utils.color('red', '/%d' % len(self.sclusts[iclust]), width=6, padside='right')  # -1 is for the added naive seq
                tmpfo = self.antn_dict[':'.join(self.sclusts[iclust])]
                # n_naive_in_cluster = len([iseq for iseq in range(len(self.sclusts[iclust])) if tmpfo['n_mutations'][iseq] == 0])  # work out if there was a sequence already in the cluster that was the same as the naive sequence
                # print '      %4d%6s' % (len(seqfos) - 1 + n_naive_in_cluster, subset_str),
                print('      %4d%6s' % (len(seqfos), subset_str), end=' ')

            if run_in_parallel:
                assert labels is None  # would need to implement this (or just switch to non-parallel version if you need to run with labels set)
                cmdfos.append(prep_cmdfo(iclust, seqfos, tqtis, color_scale_vals, title, leg_title))
            else:
                mds.run_bios2mds(self.n_mds_components, None, seqfos, self.args.workdir, self.args.random_seed,
                                 aligned=True, plotdir=plotdir, plotname=get_fname(iclust),
                                 queries_to_include=tqtis, color_scale_vals=color_scale_vals, labels=labels, title=title, leg_title=leg_title)
                if debug:
                    print('  %5.1f' % (time.time() - substart))
            self.addfname(fnames, '%s' % get_fname(iclust))

        if run_in_parallel and len(cmdfos) > 0:
            utils.run_cmds(cmdfos, clean_on_success=True)  #, debug='print')

        if len(skipped_cluster_lengths) > 0:
            print('    mds: skipped %d clusters with lengths: %s' % (len(skipped_cluster_lengths), utils.cluster_size_str(skipped_cluster_lengths, only_passing_lengths=True)))

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames)

        print('      made %d mds plots (%.1fs) with sizes %s' % (sum(len(x) for x in fnames), time.time() - start, utils.cluster_size_str(plotted_cluster_lengths, only_passing_lengths=True)))

        return [[subd + '/' + fn for fn in fnames[i]] for i in range(min(2, len(fnames)))]

    # ----------------------------------------------------------------------------------------
    def set_treefos(self, set_alternatives=False):
        # ----------------------------------------------------------------------------------------
        def find_atn(clust):
            # ----------------------------------------------------------------------------------------
            def obs_nodes(tclst):  # try to hackily ignore iqtree inferred ancestral nodes, since they have the same name in all clusters
                return [u for u in tclst if len(re.findall('Node[0-9][0-9]*', u)) == 0]
            # ----------------------------------------------------------------------------------------
            def dupstr():
                dcluststrs = [' '.join(utils.color('red' if u in clust else None, u) for u in c) for c in ovlp_clusts]
                return 'more than one cluster overlaps with %s:\n        %s' % (' '.join(clust), '\n        '.join(dcluststrs))
            # ----------------------------------------------------------------------------------------
            ovlp_clusts = [l['unique_ids'] for l in antn_list if len(set(l['unique_ids']) & set(clust)) > 0]
            if len(ovlp_clusts) > 1:
                print('    %s %s' % (utils.wrnstr(), dupstr()))
                print('       (trying to remove inferred internal nodes named e.g. NodeXX, if this works it\'s probably fine)')
                ovlp_clusts = [c for c in ovlp_clusts if len(set(obs_nodes(c)) & set(clust)) > 0]
                if len(ovlp_clusts) > 1:
                    raise Exception(dupstr())
            return utils.get_single_entry(ovlp_clusts)
        # ----------------------------------------------------------------------------------------
        if self.treefos is not None:
            return
        antn_list = [self.antn_dict.get(':'.join(c)) for c in self.sclusts]  # NOTE we *need* cluster indices here to match those in all the for loops in this file
        cibak = self.args.cluster_indices  # ick
        self.args.cluster_indices = [i for i in range(len(self.sclusts)) if self.plot_this_cluster(i, plottype='trees')]
        self.treefos = treeutils.get_treefos(self.args, antn_list, glfo=self.glfo)
        self.args.cluster_indices = cibak
        if self.args.tree_inference_method is not None:  # if the inference method infers ancestors, those sequences get added to the annotation during inference (e.g gctree, iqtree, linearham), but here we have to update the antn_dict keys and sclusts for these new seqs
            self.antn_dict = utils.get_annotation_dict(antn_list)
            for iclust, clust in enumerate(self.sclusts):  # NOTE can't just sort the 'unique_ids', since we need the indices to match up with the other plots here
                self.sclusts[iclust] = find_atn(clust)
        if set_alternatives:
            for iclust, clust in enumerate(self.sclusts):
                if self.treefos[iclust] is None:
                    continue
                annotation = self.antn_dict[':'.join(clust)]
                self.treefos[iclust]['alternatives'] = []
                for alt_atn in annotation['alternative-annotations']:
                    self.treefos[iclust]['alternatives'].append(treeutils.get_dendro_tree(treestr=alt_atn['tree-info']['lb']['tree']))
                    if self.n_max_alternatives is not None and len(self.treefos[iclust]['alternatives']) >= self.n_max_alternatives:
                        break

    # ----------------------------------------------------------------------------------------
    def set_mut_infos(self, set_alternatives=False):
        # ----------------------------------------------------------------------------------------
        def get_mutfo(dtree, antn):
            utils.add_seqs_aa(antn)
            aa_mutations, nuc_mutations = {}, {}
            treeutils.get_aa_tree(dtree, antn, nuc_mutations=nuc_mutations, aa_mutations=aa_mutations, quiet=True)
            treeutils.re_index_mut_info(self.indexing, antn, aa_mutations, nuc_mutations, is_fake_paired=antn.get('is_fake_paired', False))
            return {'mcounts' : treeutils.collect_common_mutations(aa_mutations, nuc_mutations, is_fake_paired=antn.get('is_fake_paired', False)),
                    'aa_muts' : aa_mutations,
                    'nuc_muts' : nuc_mutations}
        # ----------------------------------------------------------------------------------------
        if self.mut_info is not None:
            return
        self.mut_info = [None for _ in self.sclusts]
        for iclust, clust in enumerate(self.sclusts):  # NOTE can't just sort the 'unique_ids', since we need the indices to match up with the other plots here
            if self.treefos[iclust] is None:
                continue
            annotation = self.antn_dict[':'.join(clust)]
            self.mut_info[iclust] = get_mutfo(self.treefos[iclust]['tree'], annotation)
            if set_alternatives:
                self.mut_info[iclust]['alternatives'] = []
                for ialt, alt_atn in enumerate(annotation['alternative-annotations']):
                    dtree = self.treefos[iclust]['alternatives'][ialt]
                    self.mut_info[iclust]['alternatives'].append(get_mutfo(dtree, alt_atn))
                    if self.n_max_alternatives is not None and len(self.mut_info[iclust]['alternatives']) >= self.n_max_alternatives:
                        break

    # ----------------------------------------------------------------------------------------
    def get_treestr(self, iclust):  # at some point i'll probably need to handle None type trees, but that day is not today
        return self.treefos[iclust]['tree'].as_string(schema='newick').strip()

    # ----------------------------------------------------------------------------------------
    def make_laplacian_spectra_plots(self, debug=False):  # NOTE it's kind of weird to have this here, but all the other tree-dependent plotting in treeutils, but it's because this is for comparing clusters, whereas the stuff in treeutils is all about lb values, which are mostly useful within clusters
        subd, plotdir = self.init_subd('laplacian-spectra')
        self.set_treefos()
        if self.treefos.count(None) == len(self.treefos):
            return [['x.svg']]

        fnames = [[]]
        for iclust in range(len(self.sclusts)):
            if not self.plot_this_cluster(iclust):
                continue
            annotation = self.antn_dict[':'.join(self.sclusts[iclust])]
            if len(annotation['unique_ids']) < self.laplacian_spectra_min_clusters_size:
                continue
            if len(set(self.sclusts[iclust])) < len(self.sclusts[iclust]):
                repeated_uids = [u for u, count in collections.Counter(self.sclusts[iclust]).items() if count > 1]
                print('  skipping laplacian spectra plotting for cluster with %d duplicate uids (%s)' % (len(repeated_uids), ' '.join(repeated_uids)))
                continue
            treeutils.run_laplacian_spectra(self.get_treestr(iclust), plotdir=plotdir, plotname='icluster-%d' % iclust, title='size %d' % len(annotation['unique_ids']))
            if len(fnames[-1]) < self.n_plots_per_row:
                self.addfname(fnames, 'icluster-%d' % iclust)

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames)

        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_sfs_plots(self, restrict_to_region=None, debug=False):
        def addplot(oindexlist, ofracslist, n_seqs, fname, cdr3titlestr, red_text=None):
            hist = Hist(30, 0., 1.)
            for ofracs in ofracslist:
                hist.fill(ofracs)
            fig, ax = self.plotting.mpl_init()
            hist.mpl_plot(ax, remove_empty_bins=True)
            ax.text(0.65, 0.8 * ax.get_ylim()[1], 'h: %+.1f' % utils.fay_wu_h(line=None, restrict_to_region=restrict_to_region, occurence_indices=oindexlist, n_seqs=n_seqs), fontsize=17, fontweight='bold')
            if red_text is not None:
                ax.text(0.65, 0.7 * ax.get_ylim()[1], red_text, fontsize=17, color='red', fontweight='bold')
            titlestr = '%s (size: %d)' % (cdr3titlestr, n_seqs)

            regionstr = restrict_to_region + ' ' if restrict_to_region is not None else ''
            self.plotting.mpl_finish(ax, plotdir, fname, title=titlestr, xlabel=regionstr + 'mutation frequency', ylabel=regionstr + 'density of mutations', xticks=[0, 1], log='')  # xticks=[min(occurence_fractions), max(occurence_fractions)], 
            self.addfname(fnames, fname)

        subd, plotdir = self.init_subd('sfs')

        fnames = [[]]
        for iclust in range(len(self.sclusts)):
            if not self.plot_this_cluster(iclust):
                continue
            annotation = self.antn_dict[':'.join(self.sclusts[iclust])]
            occurence_indices, occurence_fractions = utils.get_sfs_occurence_info(annotation, restrict_to_region=restrict_to_region)
            red_text = None
            assert self.args.meta_info_to_emphasize is None  # would need to be implemented
            if self.args.queries_to_include is not None and len(set(self.args.queries_to_include) & set(self.sclusts[iclust])) > 0:
                red_text = '%s' % ' '.join(set(self.args.queries_to_include) & set(self.sclusts[iclust]))
            addplot(occurence_indices, occurence_fractions, len(self.sclusts[iclust]), 'icluster-%d' % iclust, self.get_cdr3_title(annotation), red_text=red_text)

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames)

        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_cluster_size_distribution(self):
        # ----------------------------------------------------------------------------------------
        def plot_me_color_hists():  # plot mean fraction of cluster that's X for each cluster size
            mekey = self.args.meta_info_key_to_color
            all_emph_vals, emph_colors = self.plotting.meta_emph_init(mekey, clusters=self.sclusts, antn_dict=self.antn_dict, formats=self.args.meta_emph_formats)
            hcolors = {v : c for v, c in emph_colors}
            plotvals = {v : [] for v in all_emph_vals}  # for each possible value, a list of (cluster size, fraction of seqs in cluster with that val) for clusters that contain seqs with that value
            for csize, cluster_group in itertools.groupby(self.sclusts, key=len):
                for tclust in cluster_group:
                    antn = self.antn_dict.get(':'.join(tclust))
                    if antn is None:
                        continue
                    emph_fracs = utils.get_meta_emph_fractions(mekey, all_emph_vals, tclust, antn, formats=self.args.meta_emph_formats)
                    for v, frac in emph_fracs.items():
                        plotvals[v].append((csize, frac, frac * len(tclust)))
            csize_hists.update({v : Hist(template_hist=csize_hists['best']) for v in all_emph_vals})  # for each possible value, a list of (cluster size, fraction of seqs in cluster with that val) for clusters that contain seqs with that value
            del csize_hists['best']
            ctot_hists = copy.deepcopy(csize_hists)
            for e_val, cvals in plotvals.items():
                ehist = csize_hists[e_val] #utils.meta_emph_str(mekey, e_val, formats=self.args.meta_emph_formats)]
                for ibin in ehist.ibiniter(include_overflows=True):
                    ib_vals = [(f, c) for s, f, c in cvals if ehist.find_bin(s)==ibin]  # fracs (+total seqs) whose cluster sizes fall in this bin (should all be quite similar in size if our bins are sensible, so shouldn't need to do an average weighted for cluster size)
                    if len(ib_vals) == 0:
                        continue
                    ib_fracs, ib_counts = zip(*ib_vals)
                    mval = numpy.mean(ib_fracs)
                    err = mval / math.sqrt(2) if len(ib_fracs) == 1 else numpy.std(ib_fracs, ddof=1) / math.sqrt(len(ib_fracs))  # that isn't really right for len 1, but whatever
                    ehist.set_ibin(ibin, mval, err)
                    ctot_hists[e_val].set_ibin(ibin, sum(ib_counts), math.sqrt(mval))
            ytitle = 'mean fraction of each cluster'
            bfn = 'cluster-size-fractions'
            for bstr, hlist in zip(['', '-total'], [csize_hists, ctot_hists]):
                for hname, thist in hlist.items():
                    thist.write('%s/%s%s.csv' % (plotdir, bfn+bstr, '' if hname=='best' else '-'+hname))
            self.plotting.plot_cluster_size_hists(plotdir, bfn, csize_hists, hcolors=hcolors, ytitle=ytitle, log='x', no_legend=True)
            self.plotting.plot_cluster_size_hists(plotdir, bfn+'-tot', ctot_hists, hcolors=hcolors, ytitle='total N seqs', log='', stacked_bars=True, no_legend=True)
            lfn = self.plotting.make_meta_info_legend(plotdir, bfn, self.args.meta_info_key_to_color, emph_colors, all_emph_vals, meta_emph_formats=self.args.meta_emph_formats, alpha=0.6)
            fnlist.extend([lfn, bfn, bfn+'-tot'])
        # ----------------------------------------------------------------------------------------
        subd, plotdir = self.init_subd('sizes', allow_other_files=True)  # it's too hard to rm the cluster-size-fractions subdirs

        csize_hists = {'best' : self.plotting.make_csize_hist(self.sclusts, xbins=self.args.cluster_size_bins)}
        bfn = 'cluster-sizes'
        csize_hists['best'].write('%s/%s.csv' % (plotdir, bfn))
        self.plotting.plot_cluster_size_hists(plotdir, bfn, csize_hists)
        fnlist = [bfn]

        if True: #weight_by_n_seqs:
            shist = self.plotting.make_csize_hist(self.sclusts, xbins=self.args.cluster_size_bins, weight_by_n_seqs=True)
            sfn = 'seqs-in-cluster-sizes'
            self.plotting.plot_cluster_size_hists(plotdir, sfn, {'best' : shist}, ytitle='total N seqs')
            fnlist.append(sfn)

        if self.args.meta_info_key_to_color is not None:  # would need to update this function to weight by N seqs
            plot_me_color_hists()

        return [['%s/%s.svg' % (subd, f) for f in fnlist]]

    # ----------------------------------------------------------------------------------------
    def make_tree_plots(self):
        # ----------------------------------------------------------------------------------------
        def add_mut_labels(annotation, metafo, iclust):
            # NOTE that the label/str machinations in here interact in complicated ways (i.e. it's hackey) with those in plot-lb-tree.py
            # ----------------------------------------------------------------------------------------
            def add_aa_mutstr(uid):
                if annotation.get('is_fake_paired', False):
                    chstrs = []
                    for tch in 'hl':
                        chmuts = [m for m in mutfo['aa'][uid] if m['str'].find(tch+':')==0]
                        if len(chmuts) > 0:
                             chstrs.append('%s: %s' % (tch, ', '.join(m['str'].replace(tch+':', '') for m in chmuts)))
                    metafo['labels'][uid] += '\n' + '\n'.join(chstrs)
                else:
                    metafo['labels'][uid] += '\n' + ', '.join(m['str'] for m in mutfo['aa'][uid])
            # ----------------------------------------------------------------------------------------
            mutfo = {tstr : self.mut_info[iclust]['%s_muts'%tstr] for tstr in ['aa', 'nuc']}
            metafo['labels'] = {}
            is_short = 'short' in self.args.mutation_label_cfg
            if 'leaf' in self.args.mutation_label_cfg:
                dtree = self.treefos[iclust]['tree']
            for uid in annotation['unique_ids']:
                if uid not in mutfo['nuc']:
                    continue
                has_aa_mutfo = 'mut-strs' in self.args.mutation_label_cfg and uid in mutfo['aa'] and len(mutfo['aa'][uid]) > 0
                divstr = ', ' if has_aa_mutfo else '\n'
                if 'leaf' in self.args.mutation_label_cfg:
                    node = dtree.find_node_with_taxon_label(uid)
                    if node is None or not node.is_leaf():
                        continue
                    metafo['labels'][uid] = '%d nuc%s%d aa' % (utils.per_seq_val(annotation, 'n_mutations', uid), divstr, utils.shm_aa(annotation, uid=uid))
                elif 'all' in self.args.mutation_label_cfg:
                    if len(mutfo['nuc'][uid]) == 0:
                        metafo['labels'][uid] = '' if is_short else '0'
                    # elif annotation.get('is_fake_paired', False):
                    #     metafo['labels'][uid] = '%d nuc, %d aa' % (len(mutfo['nuc'][uid]), len(mutfo['aa'][uid])) #''
                    else:
                        n_nuc, n_aa = len(mutfo['nuc'][uid]), len(mutfo['aa'][uid])
                        metafo['labels'][uid] = '%d (%d)' % (n_nuc, n_aa) if is_short else '%d nuc%s%d aa' % (n_nuc, divstr, n_aa)
                else:
                    raise Exception('expected either \'leaf\' or \'all\' in --mutation-label-cfg but got: %s' % self.args.mutation_label_cfg)
                if has_aa_mutfo:
                    add_aa_mutstr(uid)
        # ----------------------------------------------------------------------------------------
        def get_metafo(annotation, iclust):
            # ----------------------------------------------------------------------------------------
            def vmuts(vclass_muts, mutfo, uid):
                if uid not in mutfo:
                    return None
                obs_muts = [m['str'] for m in mutfo[uid]]
                obs_vclass_muts = [m for m in obs_muts if m in vclass_muts]
                return 0 if len(obs_muts) == 0 else len(obs_vclass_muts) #/ float(len(obs_muts))
            # ----------------------------------------------------------------------------------------
            malist = [self.args.meta_info_key_to_color, self.args.meta_info_to_emphasize, self.args.node_size_key, self.args.branch_color_key, self.args.mutation_label_cfg]
            if all(a is None for a in malist):
                return None, None
            metafo, cdr3fo = {}, {}
            tklist = [self.args.meta_info_key_to_color, self.args.node_size_key, self.args.branch_color_key]
            if self.args.meta_info_to_emphasize is not None:
                tklist.append(list(self.args.meta_info_to_emphasize.items())[0][0])
            if 'vrc01-muts' in tklist:  # have to set the values, and don't want to incorporate into utils.py since it's complicated and requires the tree
                import partis.vrc01 as vrc01
                vc_muts = vrc01.vrc01_class_mutation_set()  # somehow this is way faster if i don't put it in the fcn call in the next line???
                annotation['vrc01-muts'] = [vmuts(vc_muts, self.mut_info[iclust]['aa_muts'], u) for u in annotation['unique_ids']]
            for tk in [k for k in tklist if k is not None and k in annotation]:
                metafo[tk] = {u : f for u, f in zip(annotation['unique_ids'], annotation[tk])}
            if 'multiplicities' in tklist and self.args.infer_trees_with_collapsed_duplicate_seqs:
                if 'multiplicities' in metafo:  # if it was in the actual annotation
                   assert False  # needs implementing (probably can add together any existing multiplicity to the new values)
                else:
                    collapsed_seqfos = utils.collapse_seqfos_with_identical_seqs(utils.seqfos_from_line(annotation, extra_keys=[self.args.meta_info_key_to_color]), keys_not_to_collapse=[self.args.meta_info_key_to_color], debug=True)
                    metafo['multiplicities'] = {s['name'] : s['multiplicity'] for s in collapsed_seqfos}
                    metafo['duplicates'] = {s['name'] : s['duplicates'] for s in collapsed_seqfos}
                assert 'labels' not in metafo  # would need to add the new labels to the old ones, or something
                metafo['labels'] = {u : str(n) for u, n in metafo['multiplicities'].items() if n > 1}
            if self.args.mutation_label_cfg is not None:
                add_mut_labels(annotation, metafo, iclust)
                if annotation.get('is_fake_paired', False):
                    for tch in 'hl':
                        offset = annotation['%s_offset'%tch]
                        cbounds = [(b-offset)//3 + self.indexing for b in annotation['%s_cdr3_bounds'%tch]]
                        cdr3fo[tch] = '[%d - %d]' % tuple(cbounds)
            return metafo, cdr3fo
        # ----------------------------------------------------------------------------------------
        mekey = self.args.meta_info_key_to_color
        if len(self.sclusts) == 0:
            print('  %s no clusters to plot' % utils.wrnstr())
            return [['x.svg']]
        from . import lbplotting  # this is really slow because of the scipy stats import
        subd, plotdir = self.init_subd('trees')
        self.set_treefos()
        self.set_mut_infos()
        if self.treefos.count(None) == len(self.treefos):
            return [['x.svg']]
        workdir = '%s/ete3-plots' % self.args.workdir
        fnames = [['header', '%s trees'%utils.non_none([self.args.tree_inference_method, treeutils.default_inference_str])], []]  # header for html file
        bkgfns = [[]]  # for extra html file
        cmdfos = []
        n_skipped = {'not-to-plot' : 0, 'too-small' : 0, 'bkg-val' : 0}
        for iclust in range(len(self.sclusts)):
            if not self.plot_this_cluster(iclust, plottype='trees'):
                n_skipped['not-to-plot'] += 1  # this fcn has a ton of clauses, so can't really summarize it well
                continue
            annotation = self.antn_dict[':'.join(self.sclusts[iclust])]
            if len(annotation['unique_ids']) < self.min_tree_cluster_size:
                n_skipped['too-small'] += 1
                continue
            fnlist = fnames
            # find clusters that are entirely made up of any bkg val (used to skip them, now putting in different html file)
            if self.args.meta_info_bkg_vals is not None and self.args.meta_info_key_to_color in annotation:
                if all(any(utils.meta_info_equal(mekey, bval, mval) for bval in self.args.meta_info_bkg_vals) for mval in annotation[mekey]):
                    fnlist = bkgfns
                    # n_skipped['bkg-val'] += 1
                    # continue
            plotname = 'tree-iclust-%d-%s' % (iclust, utils.get_clone_id(annotation['unique_ids']))
            qtis = None if self.args.queries_to_include is None else [q for q in self.args.queries_to_include if q in annotation['unique_ids']]  # NOTE make sure to *not* modify args.queries_to_include
            altids = [(u, au) for u, au in zip(annotation['unique_ids'], annotation['alternate-uids']) if au is not None] if 'alternate-uids' in annotation else None
            mfo, cdr3fo = get_metafo(annotation, iclust)
            cfo = lbplotting.get_lb_tree_cmd(self.get_treestr(iclust), '%s/%s.svg'%(plotdir, plotname), None, None, '%s/sub-%d'%(workdir, len(cmdfos)), metafo=mfo,
                                             queries_to_include=None if self.args.dont_label_queries_to_include else qtis, meta_info_key_to_color=mekey, meta_info_to_emphasize=self.args.meta_info_to_emphasize, uid_translations=altids,
                                             label_all_nodes=self.args.label_tree_nodes, label_leaf_nodes=self.args.label_leaf_nodes, label_root_node=self.args.label_root_node, node_size_key=self.args.node_size_key, branch_color_key=self.args.branch_color_key, node_label_regex=self.args.node_label_regex, min_face_size=self.args.min_face_size, max_face_size=self.args.max_face_size)
            cmdfos.append(cfo)
            self.addfname(fnlist, plotname)
            if mekey is not None:
                self.addfname(fnlist, '%s-legend'%plotname)

            # cdr3 position info plots
            if annotation.get('is_fake_paired', False) and cdr3fo is not None and len(cdr3fo) > 0:
                self.plotting.plot_legend_only(collections.OrderedDict([('%s %s'%(c, cfo), {'color' : 'blue' if c=='h' else 'green'}) for c, cfo in cdr3fo.items()]), plotdir, '%s-cdr3'%plotname, title='CDR3')
                self.addfname(fnlist, '%s-cdr3'%plotname)
        if sum(n_skipped.values()) > 0:
            skip_dbg_strs = {'bkg-val' : 'entirely --meta-info-bkg-vals %s' % self.args.meta_info_bkg_vals}
            skstrs = ',  '.join('%d (%s)'%(n_skipped[k], skip_dbg_strs.get(k, k)) for k in [k for k in sorted(n_skipped) if n_skipped[k]>0])
            print('      skipped %d / %d trees (ran %d): %s' % (sum(n_skipped.values()), len(self.sclusts), len(cmdfos), skstrs))

        clean_up = True
        if len(cmdfos) > 0:
# TODO keep all i/o files, save command line
            # sys.exit()
            start = time.time()
            utils.run_cmds(cmdfos, clean_on_success=clean_up, shell=True, n_max_procs=utils.auto_n_procs(), proc_limit_str='plot-lb-tree.py')  # I'm not sure what the max number of procs is, but with 21 it's crashing with some of them not able to connect to the X server, and I don't see a big benefit to running them all at once anyways
            print('    made %d ete tree plots (%.1fs)' % (len(cmdfos), time.time() - start))
        if clean_up and os.path.exists(workdir):
            os.rmdir(workdir)

        if self.args.tree_inference_method == 'linearham' and self.args.outfname is not None:
            fnames.append([])
            for iclust in range(len(self.sclusts)):
                lin_plot_dir = '%s/%s/iclust-%d/lineage-plots' % (utils.fpath(utils.getprefix(self.args.outfname)), self.args.tree_inference_method, iclust)  # this duplicates code in treeutils.get_treefos() and treeutils.get_trees_for_annotations()
                for fn in glob.glob('%s/*.png'%lin_plot_dir):
                    utils.makelink(plotdir, fn, '%s/%s'%(plotdir, os.path.basename(fn)))  # ok this is two levels of links, which kind of sucks but oh well
                    self.addfname(fnames, utils.getprefix(os.path.basename(fn)), suffix='.png')

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames, extra_links=[('all tree plots', subd)], bgcolor='#FFFFFF')  # , new_table_each_row=True
            if any(len(fnl)>0 for fnl in bkgfns):
                self.plotting.make_html(plotdir, fnames=bkgfns, extra_links=[('all tree plots', subd)], bgcolor='#FFFFFF', htmlfname='%s/meta-bkg.html'%os.path.dirname(plotdir), subdname='trees')  # , new_table_each_row=True

        assert fnames[0][0] == 'header'
        # # shenanigans to put fewer trees in the main html file, maybe not worth the trouble, but also don't want to delete quite yet:
        # shorter_fnames, n_tot = [], 0  # UGH
        # for fnl in fnames:
        #     if 'header' in fnl:
        #         shorter_fnames.append(fnl)
        #     else:
        #         new_fnl = []
        #         for fn in fnl:
        #             new_fnl.append(subd+'/'+fn)
        #             if 'legend' not in fn:
        #                 n_tot += 1
        #             if n_tot > self.n_max_tree_plots:
        #                 break
        #         shorter_fnames.append(new_fnl)
        #         if n_tot > self.n_max_tree_plots:
        #             break
        # fnames[0][1] = '%d / %d %s' % (self.n_max_tree_plots, len(cmdfos), fnames[0][1])
        # fnames = [fnl if 'header' in fnl else [subd + '/' + fn for fn in fnl] for il, fnl in enumerate(fnames) if il*self.n_plots_per_row < self.n_max_tree_plots]
        return [fnl if 'header' in fnl else [subd + '/' + fn for fn in fnl] for fnl in fnames]

    # ----------------------------------------------------------------------------------------
    def make_mut_bubble_plots(self, min_n_muts=3, debug=False):
        # ----------------------------------------------------------------------------------------
        def process_single_tree(tkey, mcounts, dtree, antn, tdbg=False):
            dtree.resolve_node_ages()
            meta_vals = {u : utils.meta_emph_str(self.args.meta_info_key_to_color, v, formats=self.args.meta_emph_formats) for u, v in zip(antn['unique_ids'], antn[self.args.meta_info_key_to_color])}
            mutfos, ncounts = {}, {}
            for mstr, mct in sorted(list(mcounts[tkey].items()), key=lambda x: x[1]['count'], reverse=True):  # loop over observed mutations/positions
                # if len(mct['nodes']) < min_n_muts:  # ignore mutations that we saw fewer than <min_n_muts> times (can't do this here if we're averaging over multiple trees since we'd ignore mutations in some trees but not others
                #     continue
                if tdbg:
                    print('  %-6s   %3d occurences' % (str(mstr), len(mct['nodes'])))
                assert len(mct['nodes']) == len(set(mct['nodes']))  # shouldn't be possible for a node to be in the list twice (and wouldn't make sense)
                nodevals = collections.Counter()  # meta values for all nodes under any instance of this mutation
                for nlabel in mct['nodes']:
                    mval_counts = {}
                    tnode = dtree.find_node_with_taxon_label(nlabel)
                    if tnode is None:
                        raise Exception('couldn\'t find node with label \'%s\' in tree: %s' % (nlabel, dtree.as_string(schema='newick')))
                    treeutils.get_clade_purity(meta_vals, tnode, meta_vals[nlabel], mval_counts=mval_counts, exclude_vals=['None'])
                    nodevals.update(mval_counts)
                    if tdbg:
                        print('        %s       %s' % ('  '.join('%s: %d'%(k, v) for k, v in mval_counts.items()), utils.antnval(antn, 'alternate-uids', antn['unique_ids'].index(nlabel), default_val=nlabel, use_default=True)))
                ntot = sum(nodevals.values())
                nodevals = {k : v / float(ntot) for k, v in nodevals.items()}
                if tdbg:
                    srtd_vals = sorted(list(nodevals.items()), key=operator.itemgetter(1))
                    print('    %d total nodes:  %s  (%s)' % (ntot, '  '.join('%s %d'%(k, v*ntot) for k, v in srtd_vals), '  '.join('%s %.3f'%(k, v) for k, v in srtd_vals)))
                mutfos[mstr] = nodevals
                ncounts[mstr] = {'n-obs' : len(mct['nodes']), 'n-total-nodes' : ntot}
            return mutfos, ncounts
        # ----------------------------------------------------------------------------------------
        use_alts = self.args.tree_inference_method == 'linearham'
        if len(self.sclusts) == 0:
            print('  %s no clusters to plot' % utils.wrnstr())
            return [['x.svg']]
        subd, plotdir = self.init_subd('mut-bubble')
        self.set_treefos(set_alternatives=use_alts)
        self.set_mut_infos(set_alternatives=use_alts)
        if self.treefos.count(None) == len(self.treefos):
            return [['x.svg']]
        fnames = []
        for tkey, tlabel in [['str', 'mutation'], ['pos', 'position']]:
            fnames.append([])
            for iclust in range(len(self.sclusts)):
                if not self.plot_this_cluster(iclust, plottype='trees'):
                    continue
                annotation = self.antn_dict[':'.join(self.sclusts[iclust])]
                if len(annotation['unique_ids']) < self.min_tree_cluster_size:
                    continue

                if self.args.meta_info_key_to_color is not None:
                    all_emph_vals, emph_colors = self.plotting.meta_emph_init(self.args.meta_info_key_to_color, clusters=[self.sclusts[iclust]], antn_dict=self.antn_dict, formats=self.args.meta_emph_formats)
                    hcolors = {v : c for v, c in emph_colors}

                if use_alts:
                    talist, mclist, dtrees = annotation['alternative-annotations'], [d['mcounts'] for d in self.mut_info[iclust]['alternatives']], self.treefos[iclust]['alternatives']
                else:
                    talist, mclist, dtrees = [annotation], [self.mut_info[iclust]['mcounts']], [self.treefos[iclust]['tree']]
                all_mfos, all_ncts = [], [] # list of <mutfos> for each tree, where each <mutfos> is dict from <mstr> to {dict of <meta val> : <fraction of nodes with that meta val>}
                for tatn, mcounts, dtree in zip(talist, mclist, dtrees):
                    mutfos, ncounts = process_single_tree(tkey, mcounts, dtree, tatn)
                    all_mfos.append(mutfos)
                    all_ncts.append(ncounts)
                all_mstrs = set(k for mfs in all_mfos for k in mfs)  # all mutation/position strs
                if debug:
                    all_mvals = sorted(set([k for mfs in all_mfos for nvals in mfs.values() for k in nvals]))
                    nlen, flen = 4 * len(all_mfos), 4 * len(all_mfos)
                    def fstr(v): return '1' if v==1 else ('-' if v is None else '%.2f'%v)
                    print('  %s  %s  %s   %s' % (tlabel, utils.wfmt('N obs', 5+nlen), utils.wfmt('N total nodes', 5+nlen), '  '.join(utils.wfmt(v, 10+flen) for v in all_mvals)))
                bubfos = []
                for mstr in all_mstrs:
                    mstr_mvals = set([k for mfs in all_mfos for k in mfs.get(mstr, {})])
                    klists = {k : [mfs[mstr].get(k, 0) if mstr in mfs else None for mfs in all_mfos] for k in mstr_mvals}
                    avg_nodevals = {k : numpy.mean([v for v in klists[k] if v is not None]) for k in klists}  # {} gives 0 if mstr wasn't in that tree, 0 gives 0 if that meta value wasn't in the clade below the mutation in that tree
                    n_obs_list = [ncts.get(mstr, {'n-obs' : 0})['n-obs'] for ncts in all_ncts]
                    n_tn_list = [ncts.get(mstr, {'n-total-nodes' : 0})['n-total-nodes'] for ncts in all_ncts]
                    n_obs_mean, n_tn_mean = [numpy.mean(l) for l in [n_obs_list, n_tn_list]]
                    n_obs_std, n_tn_std = [numpy.std(l, ddof=1) / math.sqrt(len(l)) for l in [n_obs_list, n_tn_list]]
                    def mestr(m, e): return '%s+/-%s' % (('%.0f' if int(m)==m else '%.1f')%m, '0' if e==0 else '%.1f'%e)
                    if debug:
                        def klstr(k): return (' '*(10+flen)) if k not in klists else '[' + ', '.join(fstr(v) for v in klists[k]) + ']'
                        print('    %-6s  %10s %s   %10s %s      %s' % (mstr, mestr(n_obs_mean, n_obs_std), utils.wfmt(n_obs_list, nlen, jfmt='-'), mestr(n_tn_mean, n_tn_std), utils.wfmt(n_tn_list, nlen, jfmt='-'),
                                                                       '  '.join('%s %s'%(fstr(avg_nodevals[k]) if k in avg_nodevals else '', utils.wfmt(klstr(k), flen, jfmt='-')) for k in all_mvals)))

                    if n_obs_mean < min_n_muts:  # ignore mutations that we saw fewer than <min_n_muts> times (it sucks we have to get all the info for all mutations for every tree, but otherwise we'd miss some that were below the threshold in some but not others)
                        continue
                    bfo = {'id' : str(mstr), 'radius' : n_obs_mean} #int(n_obs_mean)}  # truncate to integer since a) we don't care exactly how big the bubbles are and b) i think something further along requires an integer
                    cntstr = '%s\n%s'%(mestr(n_obs_mean, n_obs_std), mestr(n_tn_mean, n_tn_std)) if use_alts else '%d\n%d'%(n_obs_mean, n_tn_mean)
                    bfo['texts'] = [{'tstr' : str(mstr), 'fsize' : 6, 'tcol' : 'black', 'dx' : -0.07 if 'is_fake_paired' in annotation else -0.05, 'dy' : 0.035},
                                    {'tstr' : cntstr, 'fsize' : 6, 'tcol' : 'black', 'dx' : -0.05, 'dy' : -0.08 if use_alts else -0.05}]
                    bfo['fracs'] = [{'label' : v, 'fraction' : f, 'color' : hcolors[utils.meta_emph_str(self.args.meta_info_key_to_color, v, formats=self.args.meta_emph_formats)]} for v, f in avg_nodevals.items()]
                    bubfos.append(bfo)

                plotname = '%s-bubbles-iclust-%d' % (tlabel, iclust)
                title = 'bubbles for %d/%d %ss (at least %d observations)' % (len(bubfos), len(all_mstrs), tlabel, min_n_muts)
                xtra_text = {'x' : 0.1, 'y' : 0.75, 'color' : 'black', 'text' : '<%s>\n<N obs>\n<N total nodes below obs>%s\n(size: N obs)' % (tlabel, ('\nmean+/-std err over %d trees'%len(all_mfos)) if use_alts else '')}
                fn = self.plotting.bubble_plot(plotname, plotdir, bubfos, title=title, xtra_text=xtra_text) #, alpha=alpha)
                self.addfname(fnames, plotname) #os.path.basename(utils.getprefix(fn)))
                lfn = self.plotting.make_meta_info_legend(plotdir, plotname, self.args.meta_info_key_to_color, emph_colors, all_emph_vals, meta_emph_formats=self.args.meta_emph_formats, alpha=0.6)
                self.addfname(fnames, lfn)

        return [[subd + '/' + fn for fn in fnl] for fnl in fnames]

    # ----------------------------------------------------------------------------------------
    def make_subtree_purity_plots(self):
        if len(self.sclusts) == 0:
            print('  %s no clusters to plot' % utils.wrnstr())
            return [['x.svg']]
        from . import lbplotting  # this is really slow because of the scipy stats import
        subd, plotdir = self.init_subd('subtree-purity')
        self.set_treefos()
        if self.treefos.count(None) == len(self.treefos):
            return [['x.svg']]
        fnames = [[]]
        for iclust in range(len(self.sclusts)):
            if not self.plot_this_cluster(iclust, plottype='trees'):
                continue
            annotation = self.antn_dict[':'.join(self.sclusts[iclust])]
            if len(annotation['unique_ids']) < self.min_tree_cluster_size:
                continue

            fnlists = lbplotting.plot_subtree_purity(plotdir, 'subtree-purity-iclust-%d' % iclust, self.treefos[iclust]['tree'], annotation, self.args.meta_info_key_to_color, meta_emph_formats=self.args.meta_emph_formats, only_csv=self.args.only_csv_plots, max_size=21)
            for fnl in fnlists:
                fnames.append(fnl)

        return [[subd + '/' + fn + '.svg' for fn in fnl] for fnl in fnames]

    # ----------------------------------------------------------------------------------------
    def make_vrc01_class_mut_plots(self):
        if len(self.sclusts) == 0:
            print('  %s no clusters to plot' % utils.wrnstr())
            return [['x.svg']]
        subd, plotdir = self.init_subd('vrc01-muts')
        fnames = [[]]
        for iclust in range(len(self.sclusts)):
            if not self.plot_this_cluster(iclust):
                continue
            annotation = self.antn_dict[':'.join(self.sclusts[iclust])]
            fnlists = self.plotting.plot_vrc01_class_muts(plotdir, 'vrc01-muts-iclust-%d' % iclust, annotation, mekey=self.args.meta_info_key_to_color, formats=self.args.meta_emph_formats, only_csv=self.args.only_csv_plots)
            for fnl in fnlists:
                fnames.append(fnl)

        return [[subd + '/' + fn + '.svg' for fn in fnl] for fnl in fnames]

    # ----------------------------------------------------------------------------------------
    def remove_failed_clusters(self, partition, annotations):
        # remove clusters with failed annotations
        failed_clusters = []
        for cluster in partition:
            if ':'.join(cluster) not in annotations:
                print('    %s cluster %s not in annotations' % (utils.color('red', 'warning'), ':'.join(cluster)))
                # # this prints overlap/extra/missing:
                # for clstr in annotations:
                #     x = clstr.split(':')
                #     if len(set(x) & set(cluster)) > 0:
                #         print '    common: %s' % ' '.join(set(x) & set(cluster))
                #         print '    not in cluster in partition: %s' % ' '.join(set(x) - set(cluster))
                #         print '    not in annotation: %s' % ' '.join(set(cluster) - set(x))
                failed_clusters.append(cluster)
        for fclust in failed_clusters:
            partition.remove(fclust)

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, partition, annotations, reco_info=None, args=None):
        if self.args.only_csv_plots:
            print('  --only-csv-plots not implemented for partition plots, so returning without plotting')
            return
        if args is not None and args.sub_plotdir is not None:
            plotdir += '/' + args.sub_plotdir
        print('    plotting partition with %d cluster%s to %s' % (len(partition), utils.plural(len(partition)), plotdir))
        sys.stdout.flush()
        start = time.time()
        if args is None or args.partition_plot_cfg is None:
            plot_cfg = utils.default_ptn_plot_cfg
        else:
            plot_cfg = args.partition_plot_cfg
        subdirs = copy.deepcopy(plot_cfg)

        fnames = []
        self.remove_failed_clusters(partition, annotations)
        self.sclusts, self.antn_dict, self.base_plotdir = sorted(partition, key=lambda c: len(c), reverse=True), annotations, plotdir
        if len(self.sclusts) == 0:
            print('    no partitions to plot')
            return

        if 'shm-vs-size' in plot_cfg:
            fnames += self.make_shm_vs_cluster_size_plots(no_slug='no-slug' in plot_cfg)
        if 'diversity' in plot_cfg:
            fnames += self.make_pairwise_diversity_plots()
        if 'cluster-bubble' in plot_cfg:
            fnames += self.make_cluster_bubble_plots()
        if 'mut-bubble' in plot_cfg:  # NOTE this needs to be before 'trees' so that we get the alternatives stuff if we're going to need it for 'mut-bubble'
            fnames += self.make_mut_bubble_plots()
        if 'trees' in plot_cfg:
            fnames += self.make_tree_plots()
        if 'subtree-purity' in plot_cfg:
            assert args is not None and args.meta_info_key_to_color is not None
            fnames += self.make_subtree_purity_plots()
        if 'mds' in plot_cfg:
            if utils.check_cmd('R', options=['--slave', '--version'], return_bool=True):
                fnames += self.make_mds_plots(reco_info=reco_info, run_in_parallel=True, aa=True) #, color_rule='wtf')
            else:
                print('  note: R does not seem to be installed, so skipping mds partition plots')
        if 'sizes' in plot_cfg:
            csfns = self.make_cluster_size_distribution()
            if 'shm-vs-size' in plot_cfg:
                if len(fnames[0]) < 4 and len(csfns[0]) < 2:
                    fnames[0] += csfns[0]
                else:
                    fnames.insert(1, csfns[0])
            else:
                fnames.append(csfns[0])

        # these probably needs testing/updating
        if 'laplacian-spectra' in plot_cfg:
            fnames += self.make_laplacian_spectra_plots()
        if 'sfs' in plot_cfg:
            fnames += self.make_sfs_plots()
        if 'vrc01-muts' in plot_cfg and self.args.locus == 'igh':
            fnames += self.make_vrc01_class_mut_plots()

        if not self.args.only_csv_plots and os.path.exists(self.base_plotdir):
            self.plotting.make_html(self.base_plotdir, fnames=fnames, new_table_each_row=True, htmlfname=self.base_plotdir + '/overview.html', extra_links=[(subd, '%s.html'%subd if os.path.exists('%s/%s.html'%(self.base_plotdir,subd)) else subd) for subd in subdirs], bgcolor='#FFFFFF')

        print('    partition plotting time: %.1f sec' % (time.time()-start))
