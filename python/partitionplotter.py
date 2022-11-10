import csv
import os
import math
import numpy
import itertools
import sys
import time
import collections
import operator
import copy

from hist import Hist
import hutils
import utils
from clusterpath import ClusterPath
import mds
import treeutils

# ----------------------------------------------------------------------------------------
class PartitionPlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, args, glfo=None):
        self.args = args
        self.glfo = glfo
        import plotting
        self.plotting = sys.modules['plotting']

        self.n_clusters_per_joy_plot = 50 if self.args.meta_info_key_to_color is None else 30
        self.n_max_joy_plots = 12
        self.n_max_mutations = 65
        self.n_joyplots_in_html = {'shm-vs-size' : self.n_max_joy_plots, 'overview' : 2}  # the rest of them are still in the dir, they're just not displayed in the html (note this is just
        self.min_high_mutation_cluster_size = 1
        self.n_biggest_to_plot = 24  # this functions as the number of plots for mds, sfs, and laplacian spectra NOTE this is overridden by self.args.queries_to_include (i.e. those queries always get plotted), but *not* by self.args.meta_info_to_emphasize
        self.n_plots_per_row = 4

        self.size_vs_shm_min_cluster_size = 3  # don't plot singletons and pairs for really big repertoires
        self.min_pairwise_cluster_size = 7
        self.min_tree_cluster_size = 5
        self.n_max_bubbles = 100  # circlify is really slow
        self.mds_max_cluster_size = 50000  # it's way tf too slow NOTE also max_cluster_size in make_mds_plots() (should combine them or at least put them in the same place)
        self.laplacian_spectra_min_clusters_size = 4
        self.min_n_clusters_to_apply_size_vs_shm_min_cluster_size = 500  # don't apply the previous thing unless the repertoire's actually pretty large

        self.n_mds_components = 2

        self.sclusts, self.antn_dict, self.treefos, self.base_plotdir = None, None, None, None

    # ----------------------------------------------------------------------------------------
    def init_subd(self, subd):
        plotdir = self.base_plotdir + '/' + subd
        utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'])
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
    def make_cluster_scatter(self, plotdir, plotname, xfcn, yfcn, clusters_to_use, repertoire_size, min_csize, dxy=0.1, log_cluster_size=False, xlabel='?', xbounds=None, repfrac_ylabel=True, colorfcn=None, leg_title=None,
                             alpha=0.65, title=''):
        # ----------------------------------------------------------------------------------------
        def add_emph_clusters():
            for cluster in clusters_to_use:
                tqtis = set()  # queries to emphasize in this cluster
                if self.args.queries_to_include is not None:
                    tqtis |= set(cluster) & set(self.args.queries_to_include)
                if self.args.meta_info_to_emphasize is not None and self.any_meta_emph(cluster):
                    key, val = self.args.meta_info_to_emphasize.items()[0]
                    tqtis.add(utils.meta_emph_str(key, val, formats=self.args.meta_emph_formats))
                if len(tqtis) == 0:
                    continue
                xval, yval = xfcn(cluster), yfcn(cluster)  # it kind of sucks to recalulate the x and y vals here, but by default (i.e. no emphasis) we don't keep track of which cluster goes with which x val, and it's nice to keep that simple
                if log_cluster_size:
                    yval = math.log(yval)
                ax.plot([xval], [yval], color='red', marker='.', markersize=2)
                ax.text(xval + dxy, yval + (0.01 if log_cluster_size else 0.1), ' '.join(tqtis), color='red', fontsize=8)

        # ----------------------------------------------------------------------------------------

        skipped_small_clusters = False
        if len(clusters_to_use) > self.min_n_clusters_to_apply_size_vs_shm_min_cluster_size:  # if repertoire is really big, ignore smaller clusters to keep the plots from being huge
            clusters_to_use = [cluster for cluster in clusters_to_use if len(cluster) >= min_csize]
            skipped_small_clusters = True
        if len(clusters_to_use) == 0:
            print '  %s no clusters to plot for cluster size scatter' % utils.color('yellow', 'warning')
            return
        xvals, yvals = zip(*[[xfcn(cluster), yfcn(cluster)] for cluster in clusters_to_use])
        if log_cluster_size:
            yvals = [math.log(yv) for yv in yvals]
        colors = None
        if colorfcn is not None:  # this is mostly copied from lbplotting.plot_2d_scatter()
            colorvals = [colorfcn(c) for c in clusters_to_use]
            smap = self.plotting.get_normalized_scalar_map(colorvals, 'viridis')
            smfcn = lambda x: 'grey' if x is None else self.plotting.get_smap_color(smap, None, val=x)
            colors = [smfcn(c) for c in colorvals]
            leg_entries = self.plotting.get_leg_entries(vals=colorvals, colorfcn=smfcn)
            self.plotting.plot_legend_only(leg_entries, plotdir, plotname+'-legend', title=leg_title, n_digits=2)
        fig, ax = self.plotting.mpl_init()
        hb = ax.scatter(xvals, yvals, c=colors, alpha=alpha)

        nticks = 5
        ymin, ymax = yvals[-1], yvals[0]
        ymin = 1  # make it 1, even if we aren't plotting small clusters, to make it more obvious that we skipped them
        yticks = [ymax + itick * (ymin - ymax) / float(nticks - 1) for itick in range(nticks)]
        if repfrac_ylabel:
            ytlfcn = lambda yt: utils.get_repfracstr(yt, repertoire_size)
        else:
            ytlfcn = lambda yt: ('%.0f' % yt) if yt > 5 else ('%.1f' % yt)
        yticklabels = [ytlfcn(math.exp(yt) if log_cluster_size else yt) for yt in yticks]

        # if necessary, add red labels to clusters
        if self.args.queries_to_include is not None or self.args.meta_info_to_emphasize is not None:
            add_emph_clusters()

        ylabel = ('family size\n(frac. of %d)' % repertoire_size) if repfrac_ylabel else 'clonal family size'
        if log_cluster_size:
            ylabel = '(log) ' + ylabel
            plotname += '-log'
        if skipped_small_clusters:
            fig.text(0.8, 0.25, 'skipping clusters\nsmaller than %d' % min_csize, color='green', fontsize=8)
        self.plotting.mpl_finish(ax, plotdir, plotname, xlabel=xlabel, ylabel=ylabel, xbounds=xbounds, ybounds=(ymin, 1.05 * ymax), yticks=yticks, yticklabels=yticklabels, title=title, leg_title=leg_title)

    # ----------------------------------------------------------------------------------------
    def make_single_hexbin_size_vs_shm_plot(self, repertoire_size, plotdir, plotname, log_cluster_size=False, debug=False):  # NOTE not using <repertoire_size> any more, but don't remember if there was a reason I should leave it
        def mean_mutations(cluster):
            return numpy.mean(self.antn_dict[':'.join(cluster)]['n_mutations'])

        clusters_to_use = [cluster for cluster in self.sclusts if mean_mutations(cluster) < self.n_max_mutations]  # have to do it as a separate line so the zip/* don't crash if no clusters pass the criterion
        self.make_cluster_scatter(plotdir, plotname, mean_mutations, len, clusters_to_use, repertoire_size, self.size_vs_shm_min_cluster_size, log_cluster_size=log_cluster_size, xlabel='mean N mutations', xbounds=(0, self.n_max_mutations))

    # ----------------------------------------------------------------------------------------
    def addfname(self, fnames, fname, force_new_row=False):
        fname += '.svg'
        if force_new_row or len(fnames[-1]) >= self.n_plots_per_row:
            fnames.append([fname])
        else:
            fnames[-1].append(fname)

    # ----------------------------------------------------------------------------------------
    def meta_emph(self, cluster, uid):  # return True if <uid> from <cluster> satisfies criteria in self.args.meta_info_to_emphasize
        antn = self.antn_dict[':'.join(cluster)]
        key, val = self.args.meta_info_to_emphasize.items()[0]
        return key in antn and utils.meta_info_equal(key, val, utils.per_seq_val(antn, key, uid), formats=self.args.meta_emph_formats)

    # ----------------------------------------------------------------------------------------
    def any_meta_emph(self, cluster):
        return any(self.meta_emph(cluster, u) for u in cluster)

    # ----------------------------------------------------------------------------------------
    def plot_this_cluster(self, iclust, plottype=None):
        if len(self.sclusts[iclust]) == 1:
            return False
        if plottype == 'mds' and len(self.sclusts[iclust]) > self.mds_max_cluster_size:
            print '     skipping mds plots for cluster with size %d > %d' % (len(self.sclusts[iclust]), self.mds_max_cluster_size)
            return False
        if plottype == 'trees' and len(self.sclusts[iclust]) < self.args.min_selection_metric_cluster_size:
            return False
        if self.args.cluster_indices is not None and iclust not in self.args.cluster_indices:
            return False
        if iclust < self.n_biggest_to_plot:
            return True
        if self.args.queries_to_include is not None and len(set(self.args.queries_to_include) & set(self.sclusts[iclust])) > 0:  # seed is added to <args.queries_to_include> in bin/partis
            return True
        return False  # falls through if <iclust> is too big, or if there's no --queries-to-include (which includes the seed)

    # ----------------------------------------------------------------------------------------
    def make_shm_vs_cluster_size_plots(self, debug=False):
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
            print '  shm vs size joyplots: divided repertoire of size %d with %d clusters into %d cluster groups' % (repertoire_size, len(self.sclusts), len(sorted_cluster_groups))
        all_emph_vals, emph_colors = None, None
        if self.args.meta_info_key_to_color is not None:  # have to do this out here before the loop so that the colors are synchronized (and all plots include all possible values)
            all_emph_vals, emph_colors = self.plotting.meta_emph_init(self.args.meta_info_key_to_color, clusters=self.sclusts, antn_dict=self.antn_dict, formats=self.args.meta_emph_formats)
        for subclusters in sorted_cluster_groups:
            if iclustergroup > self.n_max_joy_plots:  # note that when this is activated, the high mutation plot is no longer guaranteed to have every high mutation cluster (but it should have every high mutation cluster that was bigger than the cluster size when we started skipping here)
                continue
            if debug:
                print '    %d: making joyplot with %d clusters' % (iclustergroup, len(subclusters))
            title = 'per-family SHM (%d / %d)' % (iclustergroup + 1, len(sorted_cluster_groups))  # NOTE it's important that this denominator is still right even when we don't make plots for all the clusters (which it is, now)
            high_mutation_clusters += self.plotting.make_single_joyplot(subclusters, self.antn_dict, repertoire_size, plotdir, get_fname(iclustergroup=iclustergroup), cluster_indices=cluster_indices, title=title, high_x_val=self.n_max_mutations,
                                                                        queries_to_include=self.args.queries_to_include, meta_info_to_emphasize=self.args.meta_info_to_emphasize, meta_info_key_to_color=self.args.meta_info_key_to_color, meta_emph_formats=self.args.meta_emph_formats, all_emph_vals=all_emph_vals, emph_colors=emph_colors,
                                                                        make_legend=self.args.meta_info_key_to_color is not None, debug=debug)  # have to make legend for every plot
            if len(fnd['joy']) < self.n_joyplots_in_html['shm-vs-size']:
                fnd['joy'].append(get_fname(iclustergroup=iclustergroup))
            iclustergroup += 1
        if self.args.meta_info_key_to_color is not None:
            fnd['leg'] = [get_fname(iclustergroup=0)+'-legend']
        if len(high_mutation_clusters) > 0 and len(high_mutation_clusters[0]) > self.min_high_mutation_cluster_size:
            high_mutation_clusters = [cluster for cluster in high_mutation_clusters if len(cluster) > self.min_high_mutation_cluster_size]
            self.plotting.make_single_joyplot(high_mutation_clusters, self.antn_dict, repertoire_size, plotdir, get_fname(high_mutation=True), plot_high_x=True, cluster_indices=cluster_indices, title='families with mean > %d mutations' % self.n_max_mutations,
                                              high_x_val=self.n_max_mutations, queries_to_include=self.args.queries_to_include, meta_info_to_emphasize=self.args.meta_info_to_emphasize, meta_info_key_to_color=self.args.meta_info_key_to_color, meta_emph_formats=self.args.meta_emph_formats, all_emph_vals=all_emph_vals, emph_colors=emph_colors,
                                              make_legend=self.args.meta_info_key_to_color is not None, debug=debug)
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
    def make_bubble_plots(self, alpha=0.4, n_to_write_size=10, debug=False):
        import matplotlib.pyplot as plt
        subd, plotdir = self.init_subd('bubble')
        rfn = '%s/csize-radii.csv' % self.base_plotdir
        bpfn = '%s/bubble-positions.csv' % self.base_plotdir
        workfnames = [rfn, bpfn]
        mekey = self.args.meta_info_key_to_color
        fake_cluster, fake_antn = [], {'unique_ids' : [], mekey : []} # make a fake cluster with all sequences from all skipped clusters (circlify is too slow to run on all the smaller clusters)
        with open(rfn, 'w') as rfile:
            writer = csv.DictWriter(rfile, ['id', 'radius'])
            writer.writeheader()
            for iclust, cluster in enumerate(self.sclusts):
                if iclust < self.n_max_bubbles:
                    writer.writerow({'id' : iclust, 'radius' : len(cluster)})
                else:
                    fake_cluster += cluster
                    fake_antn['unique_ids'] += cluster
                    if mekey is not None:
                        antn = self.antn_dict.get(':'.join(cluster))
                        fake_antn[mekey] += antn[mekey] if antn is not None and mekey in antn else [None for _ in cluster]
            if len(fake_cluster) > 0:
                writer.writerow({'id' : 'fake', 'radius' : len(fake_cluster)})

        cmd = '%s/bin/circle-plots.py %s %s' % (utils.get_partis_dir(), rfn, bpfn)
        utils.simplerun(cmd)
        bubble_positions = []
        with open(bpfn) as bpfile:
            def cfn(k, v): return v if k=='id' else float(v)
            reader = csv.DictReader(bpfile)
            for line in reader:
                bubble_positions.append({k : cfn(k, v) for k, v in line.items()})
        fig, ax = self.plotting.mpl_init()
        if len(bubble_positions) == 0:
            print '  %s no bubble positions, returning' % utils.wrnstr()
            return 'not-plotted.svg'
        lim = max(max(abs(bfo['x']) + bfo['r'], abs(bfo['y']) + bfo['r']) for bfo in bubble_positions)
        plt.xlim(-lim, lim)
        plt.ylim(-lim, lim)
        ax.axis('off')
        plt.gca().set_aspect('equal')
        if self.args.meta_info_key_to_color is not None:
            all_emph_vals, emph_colors = self.plotting.meta_emph_init(mekey, clusters=self.sclusts, antn_dict=self.antn_dict, formats=self.args.meta_emph_formats)
            hcolors = {v : c for v, c in emph_colors}
        def getclust(idl): return self.sclusts[int(idl)] if idl!='fake' else fake_cluster
        for bfo in sorted(bubble_positions, key=lambda b: len(getclust(b['id'])), reverse=True):
            cluster = getclust(bfo['id'])
            if bfo['id']=='fake' or int(bfo['id']) < n_to_write_size:
                tstr, fsize, tcol, dx = ('small clusters', 12, 'red', -0.3) if bfo['id']=='fake' else (len(cluster), 8, 'black', -0.04)
                ax.text(bfo['x']+dx, bfo['y'], tstr, fontsize=fsize, alpha=0.4, color=tcol)
            antn = self.antn_dict.get(':'.join(cluster)) if bfo['id'] != 'fake' else fake_antn
            if antn is None or mekey is None:
                ax.add_patch(plt.Circle((bfo['x'], bfo['y']), bfo['r'], alpha=alpha, linewidth=2, fill=True))  # plain circle
            else:
                emph_fracs = utils.get_meta_emph_fractions(mekey, all_emph_vals, cluster, antn, formats=self.args.meta_emph_formats)
                self.plotting.plot_pie_chart_marker(ax, bfo['x'], bfo['y'], bfo['r'], [{'label' : v, 'fraction' : f, 'color' : hcolors[v]} for v, f in emph_fracs.items()], alpha=alpha)
        fname = 'bubbles'
        total_bubble_seqs = sum(len(getclust(b['id'])) for b in bubble_positions if b['id']!='fake')
        repertoire_size = sum([len(c) for c in self.sclusts])
        nbub = len(bubble_positions)
        if len(fake_cluster) > 0:
            nbub -= 1
        title = 'bubbles for %d/%d clusters (%d/%d seqs)' % (nbub, len(self.sclusts), total_bubble_seqs, repertoire_size)
        if len(fake_cluster) > 0:
            fig.text(0.2, 0.85, 'small clusters: %d seqs in %d clusters smaller than %d' % (len(fake_cluster), len(self.sclusts) - len(bubble_positions), len(self.sclusts[self.n_max_bubbles - 1])), fontsize=8, color='green')
        self.plotting.mpl_finish(ax, plotdir, fname, title=title)
        fnames = [[]]
        self.addfname(fnames, fname)
        if mekey is not None:
            lfn = self.plotting.make_meta_info_legend(plotdir, fname, mekey, emph_colors, all_emph_vals, meta_emph_formats=self.args.meta_emph_formats, alpha=alpha)
            self.addfname(fnames, lfn)

        for wfn in workfnames:
            os.remove(wfn)
        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_pairwise_diversity_plots(self):
        def antn(c): return self.antn_dict[':'.join(c)]
        subd, plotdir = self.init_subd('diversity')

        repertoire_size = sum([len(c) for c in self.sclusts])
        for cluster in self.sclusts:
            utils.add_seqs_aa(antn(cluster))

        fnlist = []
        for tstr in ['nuc', 'aa']:
            skey = '' if tstr=='nuc' else '_aa'
            def pwfcn(c): return utils.mean_pairwise_hdist(antn(c)['seqs%s'%skey], amino_acid=tstr=='aa')
            def mean_mutations(c): return numpy.mean(antn(c)['n_mutations'] if tstr=='nuc' else [utils.shm_aa(antn(c), iseq=i) for i in range(len(antn(c)['unique_ids']))])
            fname = 'mean-pairwise-hdist-%s' % tstr
            self.make_cluster_scatter(plotdir, fname, pwfcn, len, self.sclusts, repertoire_size, self.min_pairwise_cluster_size, dxy=0.003, log_cluster_size=True,
                                      xlabel='mean pairwise %s distance'%tstr, colorfcn=mean_mutations, title='pairwise %s diversity'%tstr, leg_title='N %s muts'%tstr) #, xbounds=(0, self.n_max_mutations))
            fnlist += ['%s/%s-%s.svg' % (subd, fname, fstr) for fstr in ['legend', 'log']]
        return [fnlist]

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
                key, val = self.args.meta_info_to_emphasize.items()[0]
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
                print '    making mds plots starting with %d clusters' % len(self.sclusts)
                print '       size (+naive)   mds    plot   total'
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
                print '      %4d%6s' % (len(seqfos), subset_str),

            if run_in_parallel:
                assert labels is None  # would need to implement this (or just switch to non-parallel version if you need to run with labels set)
                cmdfos.append(prep_cmdfo(iclust, seqfos, tqtis, color_scale_vals, title, leg_title))
            else:
                mds.run_bios2mds(self.n_mds_components, None, seqfos, self.args.workdir, self.args.random_seed,
                                 aligned=True, plotdir=plotdir, plotname=get_fname(iclust),
                                 queries_to_include=tqtis, color_scale_vals=color_scale_vals, labels=labels, title=title, leg_title=leg_title)
                if debug:
                    print '  %5.1f' % (time.time() - substart)
            self.addfname(fnames, '%s' % get_fname(iclust))

        if run_in_parallel and len(cmdfos) > 0:
            utils.run_cmds(cmdfos, clean_on_success=True)  #, debug='print')

        if len(skipped_cluster_lengths) > 0:
            print '    mds: skipped %d clusters with lengths: %s' % (len(skipped_cluster_lengths), utils.cluster_size_str(skipped_cluster_lengths, only_passing_lengths=True))

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames)

        print '      made %d mds plots (%.1fs) with sizes %s' % (sum(len(x) for x in fnames), time.time() - start, utils.cluster_size_str(plotted_cluster_lengths, only_passing_lengths=True))

        return [[subd + '/' + fn for fn in fnames[i]] for i in range(min(2, len(fnames)))]

    # ----------------------------------------------------------------------------------------
    def set_treefos(self):
        if self.treefos is not None:
            return
        antn_list = [self.antn_dict.get(':'.join(c)) for c in self.sclusts]  # NOTE we *need* cluster indices here to match those in all the for loops in this file
        self.treefos = treeutils.get_treefos(self.args, antn_list, cpath=self.cpath, glfo=self.glfo)

    # ----------------------------------------------------------------------------------------
    def get_treestr(self, iclust):  # at some point i'll probably need to handle None type trees, but that day is not today
        return self.treefos[iclust]['tree'].as_string(schema='newick').strip()

    # ----------------------------------------------------------------------------------------
    def make_laplacian_spectra_plots(self, cpath=None, debug=False):  # NOTE it's kind of weird to have this here, but all the other tree-dependent plotting in treeutils, but it's because this is for comparing clusters, whereas the stuff in treeutils is all about lb values, which are mostly useful within clusters
        subd, plotdir = self.init_subd('laplacian-spectra')
        self.set_treefos()

        fnames = [[]]
        for iclust in range(len(self.sclusts)):
            if not self.plot_this_cluster(iclust):
                continue
            annotation = self.antn_dict[':'.join(self.sclusts[iclust])]
            if len(annotation['unique_ids']) < self.laplacian_spectra_min_clusters_size:
                continue
            if len(set(self.sclusts[iclust])) < len(self.sclusts[iclust]):
                repeated_uids = [u for u, count in collections.Counter(self.sclusts[iclust]).items() if count > 1]
                print '  skipping laplacian spectra plotting for cluster with %d duplicate uids (%s)' % (len(repeated_uids), ' '.join(repeated_uids))
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
        subd, plotdir = self.init_subd('sizes')
        fname = 'cluster-sizes'
        hcolors = None
        cslist = [len(c) for c in self.sclusts]
        # is_log_x = len(cslist) > 100 # and len([s for s in cslist if s>50]) > 30
        csize_hists = {'best' : hutils.make_hist_from_list_of_values(cslist, 'int', fname, is_log_x=True)}  # seems kind of wasteful to make a bin for every integer (as here), but it's not going to be *that* many, and we want to be able to sample from them, and it's always a hassle getting the bins you want UPDATE ok now sometimes using log bins, for aesthetic plotting reasons, but maybe also ok for sampling?
        self.plotting.plot_cluster_size_hists(plotdir, fname, csize_hists)
        csize_hists['best'].write('%s/%s.csv' % (plotdir, fname))
        fnlist = [subd + '/' + fname + '.svg']
        ytitle = None
        if self.args.meta_info_key_to_color is not None:  # plot mean fraction of cluster that's X for each cluster size
            fname = 'cluster-size-fractions'
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
                        plotvals[v].append((csize, frac))
            bhist = csize_hists['best']
            csize_hists.update({v : Hist(template_hist=bhist) for v in all_emph_vals})  # for each possible value, a list of (cluster size, fraction of seqs in cluster with that val) for clusters that contain seqs with that value
            del csize_hists['best']
            for e_val, cvals in plotvals.items():
                ehist = csize_hists[e_val] #utils.meta_emph_str(mekey, e_val, formats=self.args.meta_emph_formats)]
                for ibin in ehist.ibiniter(include_overflows=True):
                    ib_vals = [f for s, f in cvals if ehist.find_bin(s)==ibin]  # fracs whose cluster sizes fall in this bin (should all be quite similar in size if our bins are sensible, so shouldn't need to do an average weighted for cluster size)
                    if len(ib_vals) == 0:
                        continue
                    mval = numpy.mean(ib_vals)
                    err = mval / math.sqrt(2) if len(ib_vals) == 1 else numpy.std(ib_vals, ddof=1) / math.sqrt(len(ib_vals))  # that isn't really right for len 1, but whatever
                    ehist.set_ibin(ibin, mval, err)
            ytitle = 'mean fraction of each cluster'

        self.plotting.plot_cluster_size_hists(plotdir, fname, csize_hists, hcolors=hcolors, ytitle=ytitle)
        for hname, thist in csize_hists.items():
            thist.write('%s/%s%s.csv' % (plotdir, fname, '' if hname=='best' else '-'+hname))
        fnlist.append(subd + '/' + fname + '.svg')
        return [fnlist]

    # ----------------------------------------------------------------------------------------
    def make_tree_plots(self, cpath=None):
        # ----------------------------------------------------------------------------------------
        def get_metafo(annotation):
            if self.args.meta_info_key_to_color is None and self.args.node_size_key is None and not self.args.label_mutations:
                return None
            metafo = {}
            for tk in [k for k in [self.args.meta_info_key_to_color, self.args.node_size_key] if k is not None and k in annotation]:
                metafo[tk] = {u : f for u, f in zip(annotation['unique_ids'], annotation[tk])}
            if self.args.label_mutations:
                utils.add_seqs_aa(annotation)
                aa_mutations, nuc_mutations = {}, {}
                treeutils.get_aa_tree(treeutils.get_dendro_tree(treestr=treestr), annotation, nuc_mutations=nuc_mutations, aa_mutations=aa_mutations)
                tkeys = {'pos' : 'position', 'str' : 'mutation'}
                mcounts = {k : {} for k in tkeys}
                for mfos in aa_mutations.values():
                    for mfo in mfos:
                        for tkey in mcounts:
                            if mfo[tkey] not in mcounts[tkey]:
                                mcounts[tkey][mfo[tkey]] = 0
                            mcounts[tkey][mfo[tkey]] += 1
                n_to_print = 10
                for tkey in mcounts:
                    print '      count   %s' % tkeys[tkey]
                    for mstr, mct in sorted(mcounts[tkey].items(), key=operator.itemgetter(1), reverse=True)[:n_to_print]:
                        print '       %3d     %s' % (mct, mstr)
                metafo['labels'] = {}
                for uid in annotation['unique_ids']:
                    if uid not in nuc_mutations:
                        continue
                    if len(nuc_mutations[uid]) == 0:
                        metafo['labels'][uid] = '0'
                    else:
                        metafo['labels'][uid] = '%d nuc, %d aa' % (len(nuc_mutations[uid]), len(aa_mutations[uid]))
                    if uid in aa_mutations and len(aa_mutations[uid]) > 0:
                        metafo['labels'][uid] += '\n' + ', '.join(mfo['str'] for mfo in aa_mutations[uid])
            return metafo
        # ----------------------------------------------------------------------------------------
        if len(self.sclusts) == 0:
            print '  %s no clusters to plot' % utils.wrnstr()
            return [['x.svg']]
        import lbplotting  # this is really slow because of the scipy stats import
        subd, plotdir = self.init_subd('trees')
        self.set_treefos()
        workdir = '%s/ete3-plots' % self.args.workdir
        fnames = [[]]
        cmdfos = []
        for iclust in range(len(self.sclusts)):
            if not self.plot_this_cluster(iclust, plottype='trees'):
                continue
            annotation = self.antn_dict[':'.join(self.sclusts[iclust])]
            if len(annotation['unique_ids']) < self.min_tree_cluster_size:
                continue
            plotname = 'tree-iclust-%d' % iclust
            cmdfos += [lbplotting.get_lb_tree_cmd(self.get_treestr(iclust), '%s/%s.svg'%(plotdir, plotname), None, None, self.args.ete_path, '%s/sub-%d'%(workdir, len(cmdfos)), metafo=get_metafo(annotation), queries_to_include=self.args.queries_to_include, meta_info_key_to_color=self.args.meta_info_key_to_color, label_all_nodes=self.args.label_tree_nodes, label_root_node=self.args.label_root_node, node_size_key=self.args.node_size_key)]
            self.addfname(fnames, plotname)
        if len(cmdfos) > 0:
            start = time.time()
            utils.run_cmds(cmdfos, clean_on_success=True, shell=True, n_max_procs=utils.auto_n_procs(), proc_limit_str='plot-lb-tree.py')  # I'm not sure what the max number of procs is, but with 21 it's crashing with some of them not able to connect to the X server, and I don't see a big benefit to running them all at once anyways
            print '    made %d ete tree plots (%.1fs)' % (len(cmdfos), time.time() - start)
        if os.path.exists(workdir):
            os.rmdir(workdir)

        if self.args.meta_info_key_to_color is not None:
            mekey = self.args.meta_info_key_to_color
            all_emph_vals, emph_colors = self.plotting.meta_emph_init(mekey, clusters=self.sclusts, antn_dict=self.antn_dict, formats=self.args.meta_emph_formats)
            self.plotting.make_meta_info_legend(plotdir, 'tree', mekey, emph_colors, all_emph_vals, meta_emph_formats=self.args.meta_emph_formats, alpha=0.6)
            self.addfname(fnames, 'tree-legend')

        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_subtree_purity_plots(self, cpath=None):

        if len(self.sclusts) == 0:
            print '  %s no clusters to plot' % utils.wrnstr()
            return [['x.svg']]
        import lbplotting  # this is really slow because of the scipy stats import
        subd, plotdir = self.init_subd('subtree-purity')
        if self.treefos is None:
            self.set_treefos()
        fnames = []
        for iclust in range(len(self.sclusts)):
            if not self.plot_this_cluster(iclust, plottype='trees'):
                continue
            annotation = self.antn_dict[':'.join(self.sclusts[iclust])]
            if len(annotation['unique_ids']) < self.min_tree_cluster_size:
                continue
            ifns = lbplotting.plot_subtree_purity(plotdir, 'subtree-purity-iclust-%d' % iclust, self.treefos[iclust]['tree'], annotation, self.args.meta_info_key_to_color, meta_emph_formats=self.args.meta_emph_formats, only_csv=self.args.only_csv_plots)
            fnames += ifns
            # for fn in ifns:
            #     self.addfname(fnames, fn)

        return [[subd + '/' + fn + '.svg' for fn in fnl] for fnl in fnames]

    # ----------------------------------------------------------------------------------------
    def remove_failed_clusters(self, partition, annotations):
        # remove clusters with failed annotations
        failed_clusters = []
        for cluster in partition:
            if ':'.join(cluster) not in annotations:
                print '    %s cluster %s not in annotations' % (utils.color('red', 'warning'), ':'.join(cluster))
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
    def plot(self, plotdir, partition, annotations, reco_info=None, cpath=None, args=None):
        if self.args.only_csv_plots:
            print '  --only-csv-plots not implemented for partition plots, so returning without plotting'
            return
        print '  plotting partitions to %s' % plotdir
        sys.stdout.flush()
        start = time.time()
        if args is None or args.partition_plot_cfg is None:
            plot_cfg = utils.default_ptn_plot_cfg
        else:
            plot_cfg = args.partition_plot_cfg
        subdirs = copy.deepcopy(plot_cfg)

        fnames = []
        self.remove_failed_clusters(partition, annotations)
        self.sclusts, self.antn_dict, self.cpath, self.base_plotdir = sorted(partition, key=lambda c: len(c), reverse=True), annotations, cpath, plotdir

        if 'shm-vs-size' in plot_cfg:
            fnames += self.make_shm_vs_cluster_size_plots()
        if 'diversity' in plot_cfg:
            fnames += self.make_pairwise_diversity_plots()
        if 'bubble' in plot_cfg:
            fnames += self.make_bubble_plots()
        if 'trees' in plot_cfg: # and self.args.meta_info_key_to_color is not None:
            fnames += self.make_tree_plots()
            if args is not None and args.meta_info_key_to_color is not None:
                fnames += self.make_subtree_purity_plots()
        if 'mds' in plot_cfg:
            if utils.check_cmd('R', options=['--slave', '--version'], return_bool=True):
                fnames += self.make_mds_plots(reco_info=reco_info, run_in_parallel=True, aa=True) #, color_rule='wtf')
            else:
                print '  note: R does not seem to be installed, so skipping mds partition plots'
        if 'sizes' in plot_cfg:
            csfns = self.make_cluster_size_distribution()
            if len(fnames) == 0: fnames.append([])
            fnames[0] += csfns[0]

        # these probably needs testing/updating
        if 'laplacian-spectra' in plot_cfg:
            fnames += self.make_laplacian_spectra_plots()
        if 'sfs' in plot_cfg:
            fnames += self.make_sfs_plots()

        if not self.args.only_csv_plots:
            self.plotting.make_html(self.base_plotdir, fnames=fnames, new_table_each_row=True, htmlfname=self.base_plotdir + '/overview.html', extra_links=[(subd, '%s.html'%subd if os.path.exists('%s/%s.html'%(self.base_plotdir,subd)) else subd) for subd in subdirs], bgcolor='#FFFFFF')

        print '    partition plotting time: %.1f sec' % (time.time()-start)
