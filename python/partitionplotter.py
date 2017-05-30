import numpy
import itertools
import sys
import time

from hist import Hist
import utils
from clusterpath import ClusterPath

# ----------------------------------------------------------------------------------------
class PartitionPlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, args):
        # between vs within stuff:
        self.n_bins = 30
        self.subplotdirs = ['overall', 'within-vs-between']
        self.args = args

    # ----------------------------------------------------------------------------------------
    def get_cdr3_length_classes(self, partition, annotations):
        classes = {}
        for cluster in partition:
            if ':'.join(cluster) not in annotations:  # i.e. if this cluster's annotation failed in the hmm when getting annotations for the final partition (which is weird, I'm surprised there's messed-up-enough sequences at this stage that that would happen)
                continue
            info = annotations[':'.join(cluster)]
            if info['cdr3_length'] not in classes:
                classes[info['cdr3_length']] = []
            classes[info['cdr3_length']].append(cluster)
        return classes

    # ----------------------------------------------------------------------------------------
    def plot_each_within_vs_between_hist(self, distances, plotdir, plotname, plottitle):
        import plotting
        xmax = 1.2 * max([d for dtype in distances for d in distances[dtype]])
        hists = {}
        for dtype in distances:
            hists[dtype] = Hist(self.n_bins, 0., xmax, title=dtype)
            for mut_freq in distances[dtype]:
                hists[dtype].fill(mut_freq)
        plotting.draw_no_root(hists['within'], plotname=plotname, plotdir=plotdir, more_hists=[hists['between']], plottitle=plottitle, xtitle='hamming distance', errors=True)

    # ----------------------------------------------------------------------------------------
    def plot_within_vs_between_hists(self, partition, annotations, base_plotdir):
        classes = self.get_cdr3_length_classes(partition, annotations)

        overall_distances = {'within' : [mut_freq for info in annotations.values() for mut_freq in info['mut_freqs']],
                             'between' : []}
        sub_distances = {}
        def nseq(cl):
            return annotations[':'.join(cl)]['naive_seq']
        for cdr3_length, clusters in classes.items():  # for each cdr3 length, loop over each pair of clusters that have that cdr3 length
            # NOTE/TODO I'm extremely unhappy that I have to put the naive seq length check here. But we pad cdr3 length subclasses to the same length during smith waterman, and by the time we get to here, in very rare cases the the cdr3 length has changed.
            hfracs = [utils.hamming_fraction(nseq(cl_a), nseq(cl_b)) for cl_a, cl_b in itertools.combinations(clusters, 2) if len(nseq(cl_a)) == len(nseq(cl_b))]  # hamming fractions for each pair of clusters with this cdr3 length
            sub_distances[cdr3_length] = {'within' : [mut_freq for cluster in clusters for mut_freq in annotations[':'.join(cluster)]['mut_freqs']],
                                          'between' : hfracs}
            overall_distances['between'] += hfracs

        self.plot_each_within_vs_between_hist(overall_distances, base_plotdir + '/overall', 'within-vs-between', '')
        for cdr3_length, subd in sub_distances.items():
            self.plot_each_within_vs_between_hist(subd, base_plotdir + '/within-vs-between', 'cdr3-length-%d' % cdr3_length, 'CDR3 %d' % cdr3_length)

    # ----------------------------------------------------------------------------------------
    def make_single_size_vs_shm_plot(self, sorted_clusters, annotations, base_plotdir, plotname, plot_high_mutation=False, debug=False):
        import plotting
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
        fig, ax = plotting.mpl_init(figsize=(xpixels / dpi, ypixels / dpi))

        n_max_mutations = 50

        min_linewidth = 0.3
        max_linewidth = 12
        # min_alpha = 0.1
        # max_alpha = 1.
        # linewidth = 7
        alpha = 0.55

        ymin, ymax = 9999, 0
        iclust_global = 0
        yticks, yticklabels = [], []

        high_mutation_clusters = []
        biggest_n_mutations = None

        if debug:
            print '  %s' % plotname
            print '    size    yval    median'

        for csize, cluster_group in itertools.groupby(sorted_clusters, key=lambda c: len(c)):
            cluster_group = sorted(list(cluster_group), key=lambda c: numpy.median(getnmutelist(c)))
            n_clusters = len(cluster_group)
            for iclust in range(len(cluster_group)):
                cluster = cluster_group[iclust]
                base_color = colors[iclust_global % len(colors)]
                if self.args.plotting_seed_ids is not None and len(set(cluster) & set(self.args.plotting_seed_ids)) > 0:
                    base_color = 'red'
                nmutelist = sorted(getnmutelist(cluster))
                if biggest_n_mutations is None or nmutelist[-1] > biggest_n_mutations:
                    biggest_n_mutations = nmutelist[-1]

                yval = len(sorted_clusters) - iclust_global
                if yval < ymin:
                    ymin = yval
                if yval > ymax:
                    ymax = yval
                yticks.append(yval)
                yticklabels.append('%d' % csize)

                if debug:
                    print '     %3s    %4.1f  %6.1f' % ('%d' % csize if iclust == 0 else '', yval, numpy.median(nmutelist)),

                if numpy.median(nmutelist) > n_max_mutations and not plot_high_mutation:
                    if debug:
                        print '%s' % utils.color('red', 'high mutation')
                    high_mutation_clusters.append(cluster)
                    continue

                if debug:
                    print ''

                nbins = nmutelist[-1] - nmutelist[0] + 1
                hist = Hist(nbins, nmutelist[0] - 0.5, nmutelist[-1] + 0.5)
                for nm in nmutelist:
                    hist.fill(nm)
                assert hist.overflow_contents() == 0.  # includes underflows
                xmax = max(hist.bin_contents)  # float(csize)
                for ibin in range(1, hist.n_bins + 1):
                    linewidth = gety(min_linewidth, max_linewidth, xmax, hist.bin_contents[ibin])
                    color = base_color
                    # alpha = gety(min_alpha, max_alpha, xmax, hist.bin_contents[ibin])
                    if hist.bin_contents[ibin] == 0.:
                        color = 'grey'
                        linewidth = min_linewidth
                        alpha = 0.4
                    ax.plot([hist.low_edges[ibin], hist.low_edges[ibin+1]], [yval, yval], color=color, linewidth=linewidth, alpha=alpha, solid_capstyle='butt')

                iclust_global += 1

        xbounds = [-0.2, n_max_mutations] if not plot_high_mutation else [n_max_mutations, biggest_n_mutations]
        ybounds = [0.95 * ymin, 1.05 * ymax]
        n_ticks = 5
        if len(yticks) > n_ticks:
            yticks = [yticks[i] for i in range(0, len(yticks), int(len(yticks) / float(n_ticks - 1)))]
            yticklabels = [yticklabels[i] for i in range(0, len(yticklabels), int(len(yticklabels) / float(n_ticks - 1)))]
        plotting.mpl_finish(ax, base_plotdir + '/overall', plotname, xlabel='N mutations', ylabel='clonal family size',
                            xbounds=xbounds, ybounds=ybounds, yticks=yticks, yticklabels=yticklabels, adjust={'left' : 0.25})

        return high_mutation_clusters

    # ----------------------------------------------------------------------------------------
    def plot_size_vs_shm(self, partition, annotations, base_plotdir, debug=False):
        fnames = []

        sorted_clusters = sorted(partition, key=lambda c: len(c), reverse=True)
        # sorted_clusters = [clust for clust in sorted_clusters if len(clust) > self.args.plotting_min_cluster_size]

        high_mutation_clusters = []

        n_clusters_per_plot = 75
        iclust = 0
        for subclusters in [sorted_clusters[i : i + n_clusters_per_plot] for i in range(0, len(sorted_clusters), n_clusters_per_plot)]:
            fnames.append('size-vs-shm-%d' % iclust)
            high_mutation_clusters += self.make_single_size_vs_shm_plot(subclusters, annotations, base_plotdir, fnames[-1], debug=debug)
            iclust += 1
        fnames.append('size-vs-shm-high-mutation')
        self.make_single_size_vs_shm_plot(high_mutation_clusters, annotations, base_plotdir, fnames[-1], plot_high_mutation=True, debug=debug)

        return [fn + '.svg' for fn in fnames]

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, partition=None, infiles=None, annotations=None, only_csv=None):
        import plotting
        print '  plotting partitions'
        sys.stdout.flush()
        start = time.time()
        for subdir in self.subplotdirs:
            utils.prep_dir(plotdir + '/' + subdir, wildlings=['*.csv', '*.svg'])

        fnames = []

        if partition is not None:  # one partition
            assert infiles is None
            assert annotations is not None
            csize_hists = {'best' : plotting.get_cluster_size_hist(partition)}
            # self.plot_within_vs_between_hists(partition, annotations, plotdir)
            fnames.append(self.plot_size_vs_shm(partition, annotations, plotdir))
        elif infiles is not None:  # plot the mean of a partition from each file
            subset_hists = []
            for fname in infiles:
                cp = ClusterPath()
                cp.readfile(fname)
                subset_hists.append(plotting.get_cluster_size_hist(cp.partitions[cp.i_best]))
            csize_hists = {'best' : plotting.make_mean_hist(subset_hists)}
            for ih in range(len(subset_hists)):
                subset_hists[ih].write(plotdir + ('/subset-%d-cluster-sizes.csv' % ih))
        else:
            assert False

        plotting.plot_cluster_size_hists(plotdir + '/overall/cluster-sizes.svg', csize_hists, title='', log='x')
        fnames.append(['cluster-sizes.svg'])

        if not only_csv:
            for subdir in self.subplotdirs:
                plotting.make_html(plotdir + '/' + subdir, fnames=fnames)

        print '(%.1f sec)' % (time.time()-start)
