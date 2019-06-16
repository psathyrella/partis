import os
import math
import numpy
import itertools
import sys
import time
import collections

from hist import Hist
import utils
from clusterpath import ClusterPath
import mds
import treeutils

# ----------------------------------------------------------------------------------------
class PartitionPlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, args):
        self.args = args
        import plotting
        self.plotting = sys.modules['plotting']
        from Bio.Seq import Seq
        self.Seq = sys.modules['Bio.Seq']

        self.n_clusters_per_joy_plot = 50
        self.n_max_joy_plots = 30
        self.n_max_mutations = 65
        self.n_joyplots_in_html = 4  # the rest of them are still in the dir, they're just not displayed in the html
        self.min_high_mutation_cluster_size = 1
        self.n_biggest_to_plot = 24
        self.n_plots_per_row = 4

        self.size_vs_shm_min_cluster_size = 3  # don't plot singletons and pairs for really big repertoires
        self.laplacian_spectra_min_clusters_size = 4
        self.max_clusters_to_apply_size_vs_shm_min_cluster_size = 500  # don't apply the previous thing unless the repertoire's actually pretty large

        self.n_mds_components = 2

    # ----------------------------------------------------------------------------------------
    def init_subd(self, subd, base_plotdir):
        plotdir = base_plotdir + '/' + subd
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
        title = self.Seq.Seq(naive_cdr3_seq).translate() + title
        return title

    # # ----------------------------------------------------------------------------------------
    # def make_single_hexbin_shm_vs_identity_plot(self, cluster, annotation, plotdir, plotname, debug=False):
    #     """ shm (identity to naive sequence) vs identity to some reference sequence """
    #     import matplotlib.pyplot as plt

    #     fig, ax = self.plotting.mpl_init()

    #     if self.args.seed_unique_id is not None and self.args.seed_unique_id in cluster:
    #         seed_index = cluster.index(self.args.seed_unique_id)
    #         ref_seq = annotation['seqs'][seed_index]
    #         ref_label = 'seed seq'
    #         xref = annotation['n_mutations'][seed_index]
    #     else:
    #         ref_seq = utils.cons_seq(0.1, aligned_seqfos=[{'name' : cluster[iseq], 'seq' : annotation['seqs'][iseq]} for iseq in range(len(cluster))])
    #         ref_label = 'consensus seq'
    #         xref = utils.hamming_distance(annotation['naive_seq'], ref_seq)

    #     xvals, yvals = zip(*[[annotation['n_mutations'][iseq], utils.hamming_distance(ref_seq, annotation['seqs'][iseq])] for iseq in range(len(cluster))])
    #     hb = ax.hexbin(xvals, yvals, gridsize=40, cmap=plt.cm.Blues, bins='log')

    #     # add a red mark for the reference sequence
    #     yref = 0
    #     ax.plot([xref], [yref], color='red', marker='.', markersize=10)
    #     ax.text(xref, yref, ref_label, color='red', fontsize=8)

    #     if self.args.queries_to_include is not None:  # note similarity to code in make_single_hexbin_size_vs_shm_plot()
    #         queries_to_include_in_this_cluster = set(cluster) & set(self.args.queries_to_include)  # TODO merge with similar code in make_single_hexbin_shm_vs_identity_plot
    #         for uid in queries_to_include_in_this_cluster:
    #             iseq = cluster.index(uid)
    #             xval = annotation['n_mutations'][iseq]
    #             yval = utils.hamming_distance(ref_seq, annotation['seqs'][iseq])
    #             ax.plot([xval], [yval], color='red', marker='.', markersize=10)
    #             ax.text(xval, yval, uid, color='red', fontsize=8)

    #     ylabel = 'identity to %s' % ref_label
    #     self.plotting.mpl_finish(ax, plotdir, plotname, xlabel='N mutations', ylabel=ylabel, title='%d sequences'  % len(cluster))

    # ----------------------------------------------------------------------------------------
    def make_single_hexbin_size_vs_shm_plot(self, sorted_clusters, annotations, repertoire_size, plotdir, plotname, log_cluster_size=False, debug=False):  # NOTE not using <repertoire_size> any more, but don't remember if there was a reason I should leave it
        import matplotlib.pyplot as plt
        def getnmutelist(cluster):
            return annotations[':'.join(cluster)]['n_mutations']

        fig, ax = self.plotting.mpl_init()

        clusters_to_use = [cluster for cluster in sorted_clusters if numpy.mean(getnmutelist(cluster)) < self.n_max_mutations]  # have to do it as a separate line so the zip/* don't crash if no clusters pass the criterion
        skipped_small_clusters = False
        if len(clusters_to_use) > self.max_clusters_to_apply_size_vs_shm_min_cluster_size:  # if repertoire is really big, ignore smaller clusters to keep the plots from being huge
            clusters_to_use = [cluster for cluster in clusters_to_use if len(cluster) >= self.size_vs_shm_min_cluster_size]
            skipped_small_clusters = True
        if len(clusters_to_use) == 0:
            print '  %s no clusters to plot' % utils.color('yellow', 'warning')
            return
        xvals, yvals = zip(*[[numpy.mean(getnmutelist(cluster)), len(cluster)] for cluster in clusters_to_use])
        if log_cluster_size:
            yvals = [math.log(yv) for yv in yvals]
        # hb = ax.hexbin(xvals, yvals, gridsize=self.n_max_mutations, cmap=plt.cm.Blues, bins='log')
        hb = ax.scatter(xvals, yvals, alpha=0.4)

        nticks = 5
        ymin, ymax = yvals[-1], yvals[0]
        ymin = 1  # make it 1, even if we aren't plotting small clusters, to make it more obvious that we skipped them
        yticks = [ymax + itick * (ymin - ymax) / float(nticks - 1) for itick in range(nticks)]
        if log_cluster_size:
            yticklabels = [math.exp(yt) for yt in yticks]
            yticklabels = [('%.0f' % yt) if yt > 5 else ('%.1f' % yt) for yt in yticklabels]
        else:
            yticklabels = [int(yt) for yt in yticks]

        if self.args.queries_to_include is not None:  # TODO merge with similar code in make_single_hexbin_shm_vs_identity_plot
            for cluster in sorted_clusters:  # NOTE just added <clusters_to_use>, and I'm not sure if I should use <sorted_clusters> or <clusters_to_use> here, but I think it's ok how it is
                queries_to_include_in_this_cluster = set(cluster) & set(self.args.queries_to_include)
                if len(queries_to_include_in_this_cluster) == 0:
                    continue
                xval = numpy.mean(getnmutelist(cluster))
                yval = len(cluster)
                if log_cluster_size:
                    yval = math.log(yval)
                ax.plot([xval], [yval], color='red', marker='.', markersize=10)
                ax.text(xval, yval, ' '.join(queries_to_include_in_this_cluster), color='red', fontsize=8)

        ylabel = 'clonal family size'
        if log_cluster_size:
            ylabel += ' (log)'
            plotname += '-log'
        if skipped_small_clusters:
            ax.text(0.7 * self.n_max_mutations, math.log(ymin + 5) if log_cluster_size else ymin + 1, 'skipping clusters smaller than %d' % self.size_vs_shm_min_cluster_size, color='green', fontsize=8)
        self.plotting.mpl_finish(ax, plotdir, plotname, xlabel='mean N mutations', ylabel=ylabel, xbounds=(0, self.n_max_mutations), ybounds=(ymin, 1.05 * ymax), yticks=yticks, yticklabels=yticklabels)

    # ----------------------------------------------------------------------------------------
    def get_repfracstr(self, csize, repertoire_size):
        repfrac = float(csize) / repertoire_size
        denom = int(1. / repfrac)
        estimate = 1. / denom
        frac_error = (estimate - repfrac) / repfrac
        if frac_error > 0.10:  # if it's more than 10% off just use the float str
            # print 'pretty far off: (1/denom - repfrac) / repfrac = (1./%d - %f) / %f = %f' % (denom, repfrac, repfrac, frac_error)
            repfracstr = '%.2f' % repfrac
        elif denom > 1000:
            repfracstr = '%.0e' % repfrac
        else:
            repfracstr = '1/%d' % denom
        return repfracstr

    # ----------------------------------------------------------------------------------------
    def make_single_joyplot(self, sorted_clusters, annotations, repertoire_size, plotdir, plotname, plot_high_mutation=False, cluster_indices=None, title=None, debug=False):
        def gety(minval, maxval, xmax, x):
            slope = (maxval - minval) / xmax
            return slope * x + minval
        def getnmutelist(cluster):
            return annotations[':'.join(cluster)]['n_mutations']

        colors = ['#006600', '#3399ff', '#ffa500']
        # goldenrod '#daa520'
        # red '#cc0000',
        # dark red '#990012'
        # purple '#a821c7'
        # grey '#808080'

        dpi = 80
        xpixels = 450
        ypixels = max(400, 10 * len(sorted_clusters))
        fig, ax = self.plotting.mpl_init(figsize=(xpixels / dpi, ypixels / dpi))

        min_linewidth = 0.3
        max_linewidth = 12
        # min_alpha = 0.1
        # max_alpha = 1.
        # linewidth = 7
        alpha = 0.55

        ymin, ymax = 9999, 0
        iclust_global = 0  # index within this plot
        yticks, yticklabels = [], []

        high_mutation_clusters = []
        biggest_n_mutations = None

        if debug:
            print '  %s   %d x %d   %s' % (plotname, xpixels, ypixels, utils.color('red', 'high mutation') if plot_high_mutation else '')
            print '      size   frac      yval    median   mean'

        for csize, cluster_group in itertools.groupby(sorted_clusters, key=lambda c: len(c)):
            cluster_group = sorted(list(cluster_group), key=lambda c: numpy.median(getnmutelist(c)))
            n_clusters = len(cluster_group)
            repfracstr = self.get_repfracstr(csize, repertoire_size)
            for iclust in range(len(cluster_group)):  # index within the clusters of this size
                cluster = cluster_group[iclust]
                nmutelist = sorted(getnmutelist(cluster))
                nmedian = numpy.median(nmutelist)
                nmean = numpy.mean(nmutelist)  # maybe should use this instead of the median?
                if biggest_n_mutations is None or nmutelist[-1] > biggest_n_mutations:
                    biggest_n_mutations = nmutelist[-1]

                if nmedian > self.n_max_mutations and not plot_high_mutation:
                    high_mutation_clusters.append(cluster)
                    continue

                yval = len(sorted_clusters) - iclust_global
                if yval < ymin:
                    ymin = yval
                if yval > ymax:
                    ymax = yval
                yticks.append(yval)
                # yticklabels.append('%d' % csize)
                yticklabels.append(repfracstr)

                base_color = colors[iclust_global % len(colors)]
                qti_n_muted = {}
                if self.args.queries_to_include is not None:
                    queries_to_include_in_this_cluster = set(cluster) & set(self.args.queries_to_include)
                    if len(queries_to_include_in_this_cluster) > 0:
                        unsorted_nmutelist = getnmutelist(cluster)
                        qti_n_muted = {uid : unsorted_nmutelist[cluster.index(uid)] for uid in queries_to_include_in_this_cluster}  # add a red line for each of 'em (i.e. color that hist bin red)
                        if plot_high_mutation:
                            xtext = 1.1
                        elif float(nmedian) / self.n_max_mutations < 0.5:
                            xtext = 0.75
                        else:
                            xtext = 0.1
                        ax.text(xtext * self.n_max_mutations, yval, ' '.join(sorted(queries_to_include_in_this_cluster, key=lambda q: qti_n_muted[q])), color='red', fontsize=8)

                if debug:
                    print '     %5s  %-10s  %4.1f  %6.1f  %6.1f' % ('%d' % csize if iclust == 0 else '', repfracstr if iclust == 0 else '', yval, nmedian, nmean)

                nbins = nmutelist[-1] - nmutelist[0] + 1
                hist = Hist(nbins, nmutelist[0] - 0.5, nmutelist[-1] + 0.5)
                for nm in nmutelist:
                    hist.fill(nm)
                assert hist.overflow_contents() == 0.  # includes underflows
                xmax = max(hist.bin_contents)  # NOTE no relation to <ymax> above
                for ibin in range(1, hist.n_bins + 1):
                    linewidth = gety(min_linewidth, max_linewidth, xmax, hist.bin_contents[ibin])
                    color = base_color
                    # alpha = gety(min_alpha, max_alpha, xmax, hist.bin_contents[ibin])
                    for nmuted in qti_n_muted.values():
                        if hist.find_bin(nmuted) == ibin:
                            color = 'red'
                    if hist.bin_contents[ibin] == 0.:
                        color = 'grey'
                        linewidth = min_linewidth
                        alpha = 0.4
                    ax.plot([hist.low_edges[ibin], hist.low_edges[ibin+1]], [yval, yval], color=color, linewidth=linewidth, alpha=alpha, solid_capstyle='butt')

                if cluster_indices is not None:
                    xtext = nmutelist[-1] if plot_high_mutation else self.n_max_mutations  # NOTE reuse of <xtext> (arg)
                    xwidth = ax.get_xlim()[1] - ax.get_xlim()[0] if plot_high_mutation else self.n_max_mutations
                    ax.text(0.05 * xwidth + xtext, yval, str(cluster_indices[':'.join(cluster)]), color=base_color, fontsize=6, alpha=alpha, fontdict={'weight' : 'bold'})
                    ax.text(0.12 * xwidth + xtext, yval, str(csize), color=base_color, fontsize=6, alpha=alpha, fontdict={'weight' : 'bold'})

                iclust_global += 1

        xbounds = [-0.2, self.n_max_mutations] if not plot_high_mutation else [self.n_max_mutations, biggest_n_mutations]
        ybounds = [0.95 * ymin, 1.05 * ymax]
        n_ticks = 5
        if len(yticks) > n_ticks:
            yticks = [yticks[i] for i in range(0, len(yticks), int(len(yticks) / float(n_ticks - 1)))]
            yticklabels = [yticklabels[i] for i in range(0, len(yticklabels), int(len(yticklabels) / float(n_ticks - 1)))]
        self.plotting.mpl_finish(ax, plotdir, plotname, xlabel='N mutations', ylabel='fraction of repertoire', title=title,  # ylabel = 'clonal family size'
                                 xbounds=xbounds, ybounds=ybounds, yticks=yticks, yticklabels=yticklabels, adjust={'left' : 0.25})

        return high_mutation_clusters

    # ----------------------------------------------------------------------------------------
    def addfname(self, fnames, fname, force_new_row=False):
        fname += '.svg'
        if force_new_row or len(fnames[-1]) >= self.n_plots_per_row:
            fnames.append([fname])
        else:
            fnames[-1].append(fname)

    # ----------------------------------------------------------------------------------------
    def plot_this_cluster(self, sorted_clusters, iclust):
        if len(sorted_clusters[iclust]) == 1:
            return False
        if iclust < self.n_biggest_to_plot:
            return True
        if self.args.queries_to_include is not None and len(set(self.args.queries_to_include) & set(sorted_clusters[iclust])) > 0:  # seed is added to <args.queries_to_include> in bin/partis
            return True
        return False  # falls through if <iclust> is too big, or if there's no --queries-to-include (which includes the seed)

    # ----------------------------------------------------------------------------------------
    def make_shm_vs_cluster_size_plots(self, sorted_clusters, annotations, base_plotdir, debug=False):
        def get_fname(iclustergroup=None, high_mutation=False, hexbin=False):
            if iclustergroup is not None:  # index of this group of clusters
                return 'size-vs-shm-%d' % iclustergroup
            elif high_mutation:
                return 'size-vs-shm-high-mutation'
            elif hexbin:
                return 'size-vs-shm-hexbin'
            else:
                assert False
        subd, plotdir = self.init_subd('shm-vs-size', base_plotdir)

        repertoire_size = sum([len(c) for c in sorted_clusters])
        cluster_indices = {':'.join(sorted_clusters[i]) : i for i in range(len(sorted_clusters))}  # index over all clusters, in the order that the mds plots will appear (compare to the two other indices I need within make_single_joyplot())

        # size vs shm joy plots
        iclustergroup = 0
        fnames = [[]]
        high_mutation_clusters = []
        sorted_cluster_groups = [sorted_clusters[i : i + self.n_clusters_per_joy_plot] for i in range(0, len(sorted_clusters), self.n_clusters_per_joy_plot)]
        if debug:
            print 'divided repertoire of size %d with %d clusters into %d cluster groups' % (repertoire_size, len(sorted_clusters), len(sorted_cluster_groups))
        for subclusters in sorted_cluster_groups:
            if iclustergroup > self.n_max_joy_plots:  # note that when this is activated, the high mutation plot is no longer guaranteed to have every high mutation cluster (but it should have every high mutation cluster that was bigger than the cluster size when we started skipping here)
                continue
            title = 'per-family SHM (%d / %d)' % (iclustergroup + 1, len(sorted_cluster_groups))  # NOTE it's important that this denominator is still right even when we don't make plots for all the clusters (which it is, now)
            high_mutation_clusters += self.make_single_joyplot(subclusters, annotations, repertoire_size, plotdir, get_fname(iclustergroup=iclustergroup), cluster_indices=cluster_indices, title=title, debug=debug)
            if len(fnames[-1]) < self.n_joyplots_in_html:
                self.addfname(fnames, get_fname(iclustergroup=iclustergroup))
            iclustergroup += 1
        if len(high_mutation_clusters) > self.n_clusters_per_joy_plot and len(high_mutation_clusters[0]) > self.min_high_mutation_cluster_size:
            high_mutation_clusters = [cluster for cluster in high_mutation_clusters if len(cluster) > self.min_high_mutation_cluster_size]
        self.make_single_joyplot(high_mutation_clusters, annotations, repertoire_size, plotdir, get_fname(high_mutation=True), plot_high_mutation=True, cluster_indices=cluster_indices, title='families with mean > %d mutations' % self.n_max_mutations, debug=debug)
        self.addfname(fnames, get_fname(high_mutation=True))

        # size vs shm hexbin plots
        self.make_single_hexbin_size_vs_shm_plot(sorted_clusters, annotations, repertoire_size, plotdir, get_fname(hexbin=True))
        self.addfname(fnames, get_fname(hexbin=True), force_new_row=True)
        self.make_single_hexbin_size_vs_shm_plot(sorted_clusters, annotations, repertoire_size, plotdir, get_fname(hexbin=True), log_cluster_size=True)
        self.addfname(fnames, get_fname(hexbin=True) + '-log')

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames, new_table_each_row=True)

        return [[subd + '/' + fn for fn in [fnames[0][0]] + fnames[-1]]]  # take the first joy plot, and the two hexbin plots

    # # ----------------------------------------------------------------------------------------
    # def make_shm_vs_inverse_identity_plots(self, sorted_clusters, annotations, base_plotdir, debug=False):
    #     def get_fname(cluster_rank):
    #         return 'shm-vs-identity-icluster-%d' % cluster_rank
    #     plotdir = base_plotdir + '/shm-vs-identity'
    #     utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'])

    #     skipped_cluster_lengths = []
    #     fnames = [[]]
    #     for iclust in range(len(sorted_clusters)):
    #         if not self.plot_this_cluster(sorted_clusters, iclust):
    #             skipped_cluster_lengths.append(len(sorted_clusters[iclust]))
    #             continue
    #         self.make_single_hexbin_shm_vs_identity_plot(sorted_clusters[iclust], annotations[':'.join(sorted_clusters[iclust])], plotdir, get_fname(cluster_rank=iclust))
    #         self.addfname(fnames, get_fname(cluster_rank=iclust))

    #     print '    skipped %d clusters with lengths: %s' % (len(skipped_cluster_lengths), ' '.join(['%d' % l for l in skipped_cluster_lengths]))

    #     if not self.args.only_csv_plots:
    #         self.plotting.make_html(plotdir, fnames=fnames)

    #     return fnames

    # ----------------------------------------------------------------------------------------
    def make_mds_plots(self, sorted_clusters, annotations, base_plotdir, max_cluster_size=10000, reco_info=None, color_rule=None, run_in_parallel=False, debug=False):
        debug = True
        # ----------------------------------------------------------------------------------------
        def get_fname(ic):
            return 'icluster-%d' % ic
        # ----------------------------------------------------------------------------------------
        def get_cluster_info(full_cluster):
            full_info = annotations[':'.join(full_cluster)]
            title = '%s   (size: %d)' % (self.get_cdr3_title(full_info), len(full_cluster))

            all_seqs = set()
            kept_indices = []
            for iseq in range(len(full_cluster)):
                if full_info['seqs'][iseq] in all_seqs:  # duplicates are from shm indels (I think I did this on purpose in sw)
                    continue
                if full_info['n_mutations'][iseq] == 0:  # remove unmutated sequences since a) they'll crash mds after we add the naive seq below and b) they'd show up in the same spot anyway (note that the only way there can be more than one is if there's Ns either within the sequences or on either end)
                    continue
                kept_indices.append(iseq)
                all_seqs.add(full_info['seqs'][iseq])

            if len(kept_indices) > max_cluster_size:
                uids_to_choose_from = set([full_cluster[i] for i in kept_indices])  # note similarity to code in seqfileopener.post_process()
                if self.args.queries_to_include is not None:
                    uids_to_choose_from -= set(self.args.queries_to_include)
                n_to_remove = len(kept_indices) - max_cluster_size
                if n_to_remove >= len(uids_to_choose_from):  # i.e. if we'd have to start removing queries that are in <queries_to_include>
                    removed_uids = uids_to_choose_from
                else:
                    removed_uids = numpy.random.choice(list(uids_to_choose_from), n_to_remove, replace=False)  # i think this'll still crash if len(uids_to_choose_from) is zero, but, meh
                kept_indices = sorted(set(kept_indices) - set([full_cluster.index(uid) for uid in removed_uids]))
                title += ' (subset: %d / %d)' % (len(kept_indices), len(full_cluster))

            seqfos = [{'name' : full_info['unique_ids'][iseq], 'seq' : full_info['seqs'][iseq]} for iseq in kept_indices]
            color_scale_vals = {full_cluster[iseq] : full_info['n_mutations'][iseq] for iseq in kept_indices}

            seqfos.append({'name' : '_naive', 'seq' : full_info['naive_seq']})  # note that if any naive sequences that were removed above are in self.args.queries_to_include, they won't be labeled in the plot (but, screw it, who's going to ask to specifically label a sequence that's already specifically labeled?)
            color_scale_vals['_naive'] = 0  # leading underscore is 'cause the mds will crash if there's another sequence with the same name, and e.g. christian's simulation spits out the naive sequence with name 'naive'. No, this is not a good long term fix
            queries_to_include = ['_naive']
            if self.args.queries_to_include is not None:
                queries_to_include += self.args.queries_to_include

            return seqfos, color_scale_vals, queries_to_include, title

        # ----------------------------------------------------------------------------------------
        def get_labels_for_coloring(full_cluster, color_rule):
            full_info = annotations[':'.join(full_cluster)]
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
        def prep_cmdfo(iclust, seqfos, queries_to_include, color_scale_vals, title):
            subworkdir = '%s/mds-%d' % (self.args.workdir, iclust)
            utils.prep_dir(subworkdir)
            tmpfname = '%s/seqs.fa' % subworkdir
            with open(tmpfname, 'w') as tmpfile:
                for sfo in seqfos:
                    csval = None
                    if sfo['name'] in color_scale_vals:
                        csval = color_scale_vals[sfo['name']]
                    tmpfile.write('>%s%s\n%s\n' % (sfo['name'], (' %d' % csval) if csval is not None else '' , sfo['seq']))
            cmdstr = './bin/mds-run.py %s --aligned --plotdir %s --plotname %s --workdir %s --seed %d' % (tmpfname, plotdir, get_fname(iclust), subworkdir, self.args.seed)
            if queries_to_include is not None:
                cmdstr += ' --queries-to-include %s' % ':'.join(queries_to_include)
            if title is not None:
                cmdstr += ' --title=%s' % title.replace(' ', '@')
            return {'cmd_str' : cmdstr, 'workdir' : subworkdir, 'outfname' : '%s/%s.svg' % (plotdir, get_fname(iclust)), 'workfnames' : [tmpfname]}

        # ----------------------------------------------------------------------------------------
        subd, plotdir = self.init_subd('mds', base_plotdir)

        start = time.time()
        if debug:
            if not run_in_parallel:
                print '    making mds plots starting with %d clusters' % len(sorted_clusters)
                print '       size (+naive)   mds    plot   total'
        skipped_cluster_lengths = []
        fnames = [[]]
        cmdfos = []
        for iclust in range(len(sorted_clusters)):
            if not self.plot_this_cluster(sorted_clusters, iclust):
                skipped_cluster_lengths.append(len(sorted_clusters[iclust]))
                continue

            seqfos, color_scale_vals, queries_to_include, title = get_cluster_info(sorted_clusters[iclust])

            labels = None
            if color_rule is not None:
                labels = get_labels_for_coloring(sorted_clusters[iclust], color_rule)
                # print '   %s setting color_scale_vals to None so we can use colors for nearest target seq index' % utils.color('red', 'note')
                color_scale_vals = None  # not sure this is really the best way to do this

            if debug and not run_in_parallel:
                substart = time.time()
                subset_str = '' if len(sorted_clusters[iclust]) <= max_cluster_size else utils.color('red', '/%d' % len(sorted_clusters[iclust]), width=6, padside='right')  # -1 is for the added naive seq
                tmpfo = annotations[':'.join(sorted_clusters[iclust])]
                # n_naive_in_cluster = len([iseq for iseq in range(len(sorted_clusters[iclust])) if tmpfo['n_mutations'][iseq] == 0])  # work out if there was a sequence already in the cluster that was the same as the naive sequence
                # print '      %4d%6s' % (len(seqfos) - 1 + n_naive_in_cluster, subset_str),
                print '      %4d%6s' % (len(seqfos), subset_str),

            if run_in_parallel:
                assert labels is None  # would need to implement this (or just switch to non-parallel version if you need to run with labels set)
                cmdfos.append(prep_cmdfo(iclust, seqfos, queries_to_include, color_scale_vals, title))
            else:
                mds.run_bios2mds(self.n_mds_components, None, seqfos, self.args.workdir, self.args.seed,
                                 aligned=True, plotdir=plotdir, plotname=get_fname(iclust),
                                 queries_to_include=queries_to_include, color_scale_vals=color_scale_vals, labels=labels, title=title)
                if debug:
                    print '  %5.1f' % (time.time() - substart)
            self.addfname(fnames, '%s' % get_fname(iclust))

        if run_in_parallel:
            utils.run_cmds(cmdfos, clean_on_success=True) #, debug='print')

        if debug and len(skipped_cluster_lengths) > 0:
            print '    skipped %d clusters with lengths: %s (+%d singletons)' % (len(skipped_cluster_lengths), ' '.join(['%d' % l for l in skipped_cluster_lengths if l > 1]), skipped_cluster_lengths.count(1))

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames)

        print '    made %d mds plots (%.1fs)' % (sum(len(x) for x in fnames), time.time() - start)

        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_laplacian_spectra_plots(self, sorted_clusters, annotations, plotdir, cpath=None, debug=False):  # NOTE it's kind of weird to have this here, but all the other tree-dependent plotting in treeutils, but it's because this is for comparing clusters, whereas the stuff in treeutils is all about lb values, which are mostly useful within clusters
        subd, plotdir = self.init_subd('laplacian-spectra', plotdir)

        fnames = [[]]
        for iclust in range(len(sorted_clusters)):
            if not self.plot_this_cluster(sorted_clusters, iclust):
                continue
            annotation = annotations[':'.join(sorted_clusters[iclust])]
            if len(annotation['unique_ids']) < self.laplacian_spectra_min_clusters_size:
                continue
            if len(set(sorted_clusters[iclust])) < len(sorted_clusters[iclust]):
                repeated_uids = [u for u, count in collections.Counter(sorted_clusters[iclust]).items() if count > 1]
                print '  skipping laplacian spectra plotting for cluster with %d duplicate uids (%s)' % (len(repeated_uids), ' '.join(repeated_uids))
                continue
            if 'tree-info' in annotation and 'lb' in annotation['tree-info']:
                treestr = annotation['tree-info']['lb']['tree']
            else:  # if this is simulation, and calculate_tree_metrics() was called with use_true_clusters=True, then we probably have to get our own trees here for the actual clusters in the best partition
                treefo = treeutils.get_tree_for_line(annotation, cpath=cpath, annotations=annotations, debug=debug)
                print '  %s no tree in annotation, so getting new tree from/with \'%s\'' % (utils.color('yellow', 'warning'), treefo['origin'])
                treestr = treefo['tree'].as_string(schema='newick').strip()
            treeutils.run_laplacian_spectra(treestr, plotdir=plotdir, plotname='icluster-%d' % iclust, title='size %d' % len(annotation['unique_ids']))
            if len(fnames[-1]) < self.n_plots_per_row:
                self.addfname(fnames, 'icluster-%d' % iclust)

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames)

        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_sfs_plots(self, sorted_clusters, annotations, base_plotdir, restrict_to_region=None, debug=False):
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

        subd, plotdir = self.init_subd('sfs', base_plotdir)

        fnames = [[]]
        for iclust in range(len(sorted_clusters)):
            if not self.plot_this_cluster(sorted_clusters, iclust):
                continue
            annotation = annotations[':'.join(sorted_clusters[iclust])]
            occurence_indices, occurence_fractions = utils.get_sfs_occurence_info(annotation, restrict_to_region=restrict_to_region)
            red_text = None
            if self.args.queries_to_include is not None and len(set(self.args.queries_to_include) & set(sorted_clusters[iclust])) > 0:
                red_text = '%s' % ' '.join(set(self.args.queries_to_include) & set(sorted_clusters[iclust]))
            addplot(occurence_indices, occurence_fractions, len(sorted_clusters[iclust]), 'icluster-%d' % iclust, self.get_cdr3_title(annotation), red_text=red_text)

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames)

        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_cluster_size_distribution(self, base_plotdir, partition=None, infiles=None):
        subd, plotdir = self.init_subd('sizes', base_plotdir)

        if partition is not None:  # one partition
            csize_hists = {'best' : self.plotting.get_cluster_size_hist(partition)}
        elif infiles is not None:  # plot the mean of a partition from each file
            subset_hists = []
            for fname in infiles:
                cp = ClusterPath(fname=fname)
                subset_hists.append(self.plotting.get_cluster_size_hist(cp.partitions[cp.i_best]))
            csize_hists = {'best' : self.plotting.make_mean_hist(subset_hists)}
            for ih in range(len(subset_hists)):
                subset_hists[ih].write(plotdir + ('/subset-%d-cluster-sizes.csv' % ih))
        else:
            assert False

        self.plotting.plot_cluster_size_hists(plotdir + '/cluster-sizes.svg', csize_hists, title='', log='x')
        return [[subd + '/cluster-sizes.svg']]

    # ----------------------------------------------------------------------------------------
    def remove_failed_clusters(self, partition, annotations):
        # remove clusters with failed annotations
        failed_clusters = []
        for cluster in partition:
            if ':'.join(cluster) not in annotations:
                print '    %s cluster %s not in annotations' % (utils.color('red', 'warning'), ':'.join(cluster))
                failed_clusters.append(cluster)
        for fclust in failed_clusters:
            partition.remove(fclust)

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, partition=None, infiles=None, annotations=None, reco_info=None, cpath=None):
        if self.args.only_csv_plots:
            print '  --only-csv-plots not implemented for partition plots, so skipping'
            return
        assert (partition is None and annotations is not None) or infiles is None
        print '  plotting partitions'
        sys.stdout.flush()
        start = time.time()

        fnames = []
        if partition is not None:  # only one partition
            self.remove_failed_clusters(partition, annotations)
            sorted_clusters = sorted(partition, key=lambda c: len(c), reverse=True)
            fnames += self.make_shm_vs_cluster_size_plots(sorted_clusters, annotations, plotdir)
            fnames += self.make_mds_plots(sorted_clusters, annotations, plotdir, reco_info=reco_info, run_in_parallel=True) #, color_rule='wtf')
            # fnames += self.make_laplacian_spectra_plots(sorted_clusters, annotations, plotdir, cpath=cpath)
            # fnames += self.make_sfs_plots(sorted_clusters, annotations, plotdir)
        self.make_cluster_size_distribution(plotdir, partition=partition, infiles=infiles)

        if not self.args.only_csv_plots:
            self.plotting.make_html(plotdir, fnames=fnames, new_table_each_row=True, htmlfname=plotdir + '/overview.html', extra_links=[(subd, '%s/%s.html' % (plotdir, subd)) for subd in ['shm-vs-size', 'mds', 'laplacian-spectra']])  # , 'sfs

        print '    partition plotting time: %.1f sec' % (time.time()-start)
