import math
import numpy
import itertools
import sys
import time

from hist import Hist
import utils
from clusterpath import ClusterPath
import mds

# ----------------------------------------------------------------------------------------
class PartitionPlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, args):
        self.args = args
        import plotting
        self.plotting = sys.modules['plotting']

        self.n_clusters_per_joy_plot = 50
        self.n_max_mutations = 65
        self.n_joyplots_in_html = 4  # the rest of them are still in the dir, they're just not displayed in the html
        self.min_high_mutation_cluster_size = 1
        self.n_biggest_to_plot = 12
        self.n_plots_per_row = 4

        self.n_mds_components = 2

    # self.n_bins = 30
    # # ----------------------------------------------------------------------------------------
    # def get_cdr3_length_classes(self, partition, annotations):
    #     # NOTE this should be replaced by the fcn from utils
    #     classes = {}
    #     for cluster in partition:
    #         if ':'.join(cluster) not in annotations:  # i.e. if this cluster's annotation failed in the hmm when getting annotations for the final partition (which is weird, I'm surprised there's messed-up-enough sequences at this stage that that would happen)
    #             continue
    #         info = annotations[':'.join(cluster)]
    #         if info['cdr3_length'] not in classes:
    #             classes[info['cdr3_length']] = []
    #         classes[info['cdr3_length']].append(cluster)
    #     return classes

    # # ----------------------------------------------------------------------------------------
    # def plot_each_within_vs_between_hist(self, distances, plotdir, plotname, plottitle):
    #     xmax = 1.2 * max([d for dtype in distances for d in distances[dtype]])
    #     hists = {}
    #     for dtype in distances:
    #         hists[dtype] = Hist(self.n_bins, 0., xmax, title=dtype)
    #         for mut_freq in distances[dtype]:
    #             hists[dtype].fill(mut_freq)
    #     self.plotting.draw_no_root(hists['within'], plotname=plotname, plotdir=plotdir, more_hists=[hists['between']], plottitle=plottitle, xtitle='hamming distance', errors=True)

    # # ----------------------------------------------------------------------------------------
    # def plot_within_vs_between_hists(self, partition, annotations, base_plotdir):
    #     classes = self.get_cdr3_length_classes(partition, annotations)

    #     overall_distances = {'within' : [mut_freq for info in annotations.values() for mut_freq in info['mut_freqs']],
    #                          'between' : []}
    #     sub_distances = {}
    #     def nseq(cl):
    #         return annotations[':'.join(cl)]['naive_seq']
    #     for cdr3_length, clusters in classes.items():  # for each cdr3 length, loop over each pair of clusters that have that cdr3 length
    #         # NOTE/TODO I'm extremely unhappy that I have to put the naive seq length check here. But we pad cdr3 length subclasses to the same length during smith waterman, and by the time we get to here, in very rare cases the the cdr3 length has changed.
    #         hfracs = [utils.hamming_fraction(nseq(cl_a), nseq(cl_b)) for cl_a, cl_b in itertools.combinations(clusters, 2) if len(nseq(cl_a)) == len(nseq(cl_b))]  # hamming fractions for each pair of clusters with this cdr3 length
    #         sub_distances[cdr3_length] = {'within' : [mut_freq for cluster in clusters for mut_freq in annotations[':'.join(cluster)]['mut_freqs']],
    #                                       'between' : hfracs}
    #         overall_distances['between'] += hfracs

    #     self.plot_each_within_vs_between_hist(overall_distances, base_plotdir + '/overall', 'within-vs-between', '')
    #     for cdr3_length, subd in sub_distances.items():
    #         self.plot_each_within_vs_between_hist(subd, base_plotdir + '/within-vs-between', 'cdr3-length-%d' % cdr3_length, 'CDR3 %d' % cdr3_length)

    # ----------------------------------------------------------------------------------------
    def make_single_hexbin_shm_vs_identity_plot(self, cluster, annotation, plotdir, plotname, debug=False):
        """ shm (identity to naive sequence) vs identity to some reference sequence """
        import matplotlib.pyplot as plt

        fig, ax = self.plotting.mpl_init()

        if self.args.seed_unique_id is not None and self.args.seed_unique_id in cluster:
            seed_index = cluster.index(self.args.seed_unique_id)
            ref_seq = annotation['seqs'][seed_index]
            ref_label = 'seed seq'
            xref = annotation['n_mutations'][seed_index]
        else:
            ref_seq = utils.cons_seq(0.1, aligned_seqfos=[{'name' : cluster[iseq], 'seq' : annotation['seqs'][iseq]} for iseq in range(len(cluster))])
            ref_label = 'consensus seq'
            xref = utils.hamming_distance(annotation['naive_seq'], ref_seq)

        xvals, yvals = zip(*[[annotation['n_mutations'][iseq], utils.hamming_distance(ref_seq, annotation['seqs'][iseq])] for iseq in range(len(cluster))])
        hb = ax.hexbin(xvals, yvals, gridsize=40, cmap=plt.cm.Blues, bins='log')

        # add a red mark for the reference sequence
        yref = 0
        ax.plot([xref], [yref], color='red', marker='.', markersize=10)
        ax.text(xref, yref, ref_label, color='red', fontsize=8)

        if self.args.queries_to_include is not None:  # note similarity to code in make_single_hexbin_size_vs_shm_plot()
            queries_to_include_in_this_cluster = set(cluster) & set(self.args.queries_to_include)  # TODO merge with similar code in make_single_hexbin_shm_vs_identity_plot
            for uid in queries_to_include_in_this_cluster:
                iseq = cluster.index(uid)
                xval = annotation['n_mutations'][iseq]
                yval = utils.hamming_distance(ref_seq, annotation['seqs'][iseq])
                ax.plot([xval], [yval], color='red', marker='.', markersize=10)
                ax.text(xval, yval, uid, color='red', fontsize=8)

        ylabel = 'identity to %s' % ref_label
        self.plotting.mpl_finish(ax, plotdir, plotname, xlabel='N mutations', ylabel=ylabel, title='%d sequences'  % len(cluster))

    # ----------------------------------------------------------------------------------------
    def make_single_hexbin_size_vs_shm_plot(self, sorted_clusters, annotations, repertoire_size, plotdir, plotname, log_cluster_size=False, debug=False):  # NOTE not using <repertoire_size> any more, but don't remember if there was a reason I should leave it
        import matplotlib.pyplot as plt
        def getnmutelist(cluster):
            return annotations[':'.join(cluster)]['n_mutations']

        fig, ax = self.plotting.mpl_init()

        clusters_to_use = [cluster for cluster in sorted_clusters if numpy.mean(getnmutelist(cluster)) < self.n_max_mutations]  # have to do it as a separate line so the zip/* don't crash if no clusters pass the criterion
        if len(clusters_to_use) == 0:
            print '  %s no clusters to plot' % utils.color('yellow', 'warning')
            return
        xvals, yvals = zip(*[[numpy.mean(getnmutelist(cluster)), len(cluster)] for cluster in clusters_to_use])
        if log_cluster_size:
            yvals = [math.log(yv) for yv in yvals]
        hb = ax.hexbin(xvals, yvals, gridsize=self.n_max_mutations, cmap=plt.cm.Blues, bins='log')

        nticks = 5
        yticks = [yvals[0] + itick * (yvals[-1] - yvals[0]) / float(nticks - 1) for itick in range(nticks)]
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
        self.plotting.mpl_finish(ax, plotdir, plotname, xlabel='mean N mutations', ylabel=ylabel, xbounds=[0, self.n_max_mutations], yticks=yticks, yticklabels=yticklabels)

    # ----------------------------------------------------------------------------------------
    def get_repfracstr(self, csize, repertoire_size):
        repfrac = float(csize) / repertoire_size
        denom = int(1. / repfrac)
        estimate = 1. / denom
        frac_error = (estimate - repfrac) / repfrac
        if frac_error > 0.10:  # if it's more than 10% off just use the float str
            repfracstr = '%f' % repfrac
        elif denom > 10000:
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
                yticklabels.append('%d' % csize)
                # yticklabels.append(repfracstr)

                base_color = colors[iclust_global % len(colors)]
                if self.args.queries_to_include is not None:
                    queries_to_include_in_this_cluster = set(cluster) & set(self.args.queries_to_include)
                    if len(queries_to_include_in_this_cluster) > 0:
                        base_color = 'red'
                        if plot_high_mutation:
                            xtext = 1.1
                        elif float(nmedian) / self.n_max_mutations < 0.5:
                            xtext = 0.75
                        else:
                            xtext = 0.1
                        ax.text(xtext * self.n_max_mutations, yval, ' '.join(queries_to_include_in_this_cluster), color='red', fontsize=8)

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
                    if hist.bin_contents[ibin] == 0.:
                        color = 'grey'
                        linewidth = min_linewidth
                        alpha = 0.4
                    ax.plot([hist.low_edges[ibin], hist.low_edges[ibin+1]], [yval, yval], color=color, linewidth=linewidth, alpha=alpha, solid_capstyle='butt')

                if cluster_indices is not None:
                    xtext = nmutelist[-1] if plot_high_mutation else self.n_max_mutations  # NOTE reuse of <xtext> (arg)
                    xwidth = ax.get_xlim()[1] - ax.get_xlim()[0] if plot_high_mutation else self.n_max_mutations
                    ax.text(0.05 * xwidth + xtext, yval, str(cluster_indices[':'.join(cluster)]), color=base_color, fontsize=6, alpha=alpha, fontdict={'weight' : 'bold'})
                    ax.text(0.10 * xwidth + xtext, yval, str(csize), color=base_color, fontsize=6, alpha=alpha, fontdict={'weight' : 'bold'})

                iclust_global += 1

        xbounds = [-0.2, self.n_max_mutations] if not plot_high_mutation else [self.n_max_mutations, biggest_n_mutations]
        ybounds = [0.95 * ymin, 1.05 * ymax]
        n_ticks = 5
        if len(yticks) > n_ticks:
            yticks = [yticks[i] for i in range(0, len(yticks), int(len(yticks) / float(n_ticks - 1)))]
            yticklabels = [yticklabels[i] for i in range(0, len(yticklabels), int(len(yticklabels) / float(n_ticks - 1)))]
        self.plotting.mpl_finish(ax, plotdir, plotname, xlabel='N mutations', ylabel='clonal family size', title=title,
                            xbounds=xbounds, ybounds=ybounds, yticks=yticks, yticklabels=yticklabels, adjust={'left' : 0.18})

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
        plotting = sys.modules['plotting']
        subd = 'shm-vs-size'
        plotdir = base_plotdir + '/' + subd
        utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'])

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
            title = 'per-family SHM (%d / %d)' % (iclustergroup + 1, len(sorted_cluster_groups))
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

        self.plotting.make_html(plotdir, fnames=fnames, new_table_each_row=True)

        return [[subd + '/' + fn for fn in [fnames[0][0]] + fnames[-1]]]  # take the first joy plot, and the two hexbin plots

    # # ----------------------------------------------------------------------------------------
    # def make_shm_vs_inverse_identity_plots(self, sorted_clusters, annotations, base_plotdir, debug=False):
    #     def get_fname(cluster_rank):
    #         return 'shm-vs-identity-icluster-%d' % cluster_rank
    #     plotting = sys.modules['plotting']
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

    #     self.plotting.make_html(plotdir, fnames=fnames)

    #     return fnames

    # ----------------------------------------------------------------------------------------
    def make_mds_plots(self, sorted_clusters, annotations, base_plotdir, debug=False):
        def get_fname(ic):
            return 'icluster-%d' % ic
        plotting = sys.modules['plotting']
        subd = 'mds'
        plotdir = base_plotdir + '/' + subd
        utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'])

        skipped_cluster_lengths = []
        fnames = [[]]
        for iclust in range(len(sorted_clusters)):
            cluster = sorted_clusters[iclust]
            if not self.plot_this_cluster(sorted_clusters, iclust):
                skipped_cluster_lengths.append(len(sorted_clusters[iclust]))
                continue
            info = annotations[':'.join(cluster)]
            seqfos = [{'name' : info['unique_ids'][iseq], 'seq' : info['seqs'][iseq]} for iseq in range(len(cluster))]
            seqfos.append({'name' : 'naive', 'seq' : info['naive_seq']})
            color_vals = {cluster[iseq] : info['n_mutations'][iseq] for iseq in range(len(cluster))}
            color_vals['naive'] = 0
            queries_to_include = ['naive']
            if self.args.queries_to_include is not None:
                queries_to_include += self.args.queries_to_include
            mds.bios2mds_kmeans_cluster(self.n_mds_components, None, seqfos, self.args.workdir, self.args.seed, plotdir=plotdir, plotname=get_fname(iclust), queries_to_include=queries_to_include, color_vals=color_vals)
            self.addfname(fnames, '%s' % get_fname(iclust))

        print '    skipped %d clusters with lengths: %s' % (len(skipped_cluster_lengths), ' '.join(['%d' % l for l in skipped_cluster_lengths]))

        self.plotting.make_html(plotdir, fnames=fnames)

        return [[subd + '/' + fn for fn in fnames[0]]]

    # ----------------------------------------------------------------------------------------
    def make_cluster_size_distribution(self, base_plotdir, partition=None, infiles=None):
        plotting = sys.modules['plotting']
        subd = 'sizes'
        plotdir = base_plotdir + '/' + subd
        utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'])

        if partition is not None:  # one partition
            csize_hists = {'best' : self.plotting.get_cluster_size_hist(partition)}
        elif infiles is not None:  # plot the mean of a partition from each file
            subset_hists = []
            for fname in infiles:
                cp = ClusterPath()
                cp.readfile(fname)
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
    def plot(self, plotdir, partition=None, infiles=None, annotations=None, only_csv=None):
        assert (partition is None and annotations is not None) or infiles is None
        print '  plotting partitions'
        sys.stdout.flush()
        start = time.time()

        fnames = []
        if partition is not None:  # only one partition
            self.remove_failed_clusters(partition, annotations)
            sorted_clusters = sorted(partition, key=lambda c: len(c), reverse=True)
            fnames += self.make_shm_vs_cluster_size_plots(sorted_clusters, annotations, plotdir)
            fnames += self.make_mds_plots(sorted_clusters, annotations, plotdir)
            # fnames += self.make_shm_vs_inverse_identity_plots(sorted_clusters, annotations, plotdir)
        fnames += self.make_cluster_size_distribution(plotdir, partition=partition, infiles=infiles)

        assert not only_csv  # meh, too much trouble to deal with this now
        self.plotting.make_html(plotdir, fnames=fnames, new_table_each_row=True, htmlfname=plotdir + '/overview.html')

        print '(%.1f sec)' % (time.time()-start)
