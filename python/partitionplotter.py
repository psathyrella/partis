import itertools
import sys
import time

from hist import Hist
import utils
from clusterpath import ClusterPath

# ----------------------------------------------------------------------------------------
class PartitionPlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self):
        # between vs within stuff:
        self.n_bins = 30
        self.subplotdirs = ['overall', 'within-vs-between']

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
    def plot_size_vs_mfreq(self, partition, annotations, base_plotdir):
        import plotting

        binwidth = 0.0003
        min_alpha = 0.1
        max_alpha = 0.8
        min_linewidth = 1
        max_linewidth = 12

        def gety(minval, maxval, xmax, x):
            slope = (maxval - minval) / xmax
            return slope * x + minval

        fig, ax = plotting.mpl_init()
        colors = ['#006600', '#990012', '#cc0000', '#3399ff', '#a821c7', '#808080']  # plotting.default_colors
        last_cluster_size = None
        last_color = 0
        y_adjust = 0
        sorted_clusters = sorted(partition, key=lambda c: len(c), reverse=True)
        for cluster in sorted_clusters:
            mfreqs = sorted(annotations[':'.join(cluster)]['mut_freqs'])

            if last_cluster_size is None or last_cluster_size != len(cluster):
                print ' %3d' % len(cluster)
                color = colors[0]
                y_adjust = 0.
            else:
                color = colors[(colors.index(last_color) + 1) % len(colors)]
                y_adjust = (-1)**(colors.index(color)) * colors.index(color) * 0.15

            print '      %8s   %5.2f   %s' % (color, y_adjust, ' '.join(['%.3f' % mf for mf in mfreqs]))

            nbins = int((mfreqs[-1] - mfreqs[0]) / binwidth) + 2
            hist = Hist(nbins, mfreqs[0] - binwidth, mfreqs[-1] + binwidth)
            for mf in mfreqs:
                hist.fill(mf)
            assert hist.overflow_contents() == 0.
            max_bin_contents = max(hist.bin_contents)
            xmax = float(len(cluster))  # max_bin_contents
            yval = len(cluster) + y_adjust
            for ibin in range(1, hist.n_bins + 1):
                alpha = gety(min_alpha, max_alpha, xmax, hist.bin_contents[ibin])
                linewidth = gety(min_linewidth, max_linewidth, xmax, hist.bin_contents[ibin])
                ax.plot([hist.low_edges[ibin], hist.low_edges[ibin+1]], [yval, yval], color=color, linewidth=linewidth, alpha=alpha, solid_capstyle='butt')

            last_cluster_size = len(cluster)
            last_color = color

        ybounds = [0.9 * len(sorted_clusters[-1]), 1.05 * len(sorted_clusters[0])]
        plotnames = {'size-vs-mfreq' : [-0.001, 0.25], 'size-vs-mfreq-zoomed' : [-0.001, 0.05]}
        for name, xbounds in plotnames.items():
            plotting.mpl_finish(ax, base_plotdir + '/overall', name, xlabel='mut freq', ylabel='clonal family size', xbounds=xbounds, ybounds=ybounds)

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, partition=None, infiles=None, annotations=None, only_csv=None):
        import plotting
        print '  plotting partitions'
        sys.stdout.flush()
        start = time.time()
        for subdir in self.subplotdirs:
            utils.prep_dir(plotdir + '/' + subdir, wildlings=['*.csv', '*.svg'])

        if partition is not None:  # one partition
            assert infiles is None
            assert annotations is not None
            csize_hists = {'best' : plotting.get_cluster_size_hist(partition)}
            # self.plot_within_vs_between_hists(partition, annotations, plotdir)
            self.plot_size_vs_mfreq(partition, annotations, plotdir)
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

        if not only_csv:
            for subdir in self.subplotdirs:
                plotting.make_html(plotdir + '/' + subdir)

        print '(%.1f sec)' % (time.time()-start)
