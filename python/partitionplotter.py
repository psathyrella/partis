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

        colors = ['#006600', '#3399ff', '#ffa500']
        # goldenrod '#daa520'
        # red '#cc0000',
        # dark red '#990012'
        # purple '#a821c7'
        # grey '#808080'

        def gety(minval, maxval, xmax, x):
            slope = (maxval - minval) / xmax
            return slope * x + minval
        def getnmutelist(cluster):
            return annotations[':'.join(cluster)]['n_mutations']

        fig, ax = plotting.mpl_init()

        sorted_clusters = sorted(partition, key=lambda c: len(c), reverse=True)
        biggest_cluster_size = len(sorted_clusters[0])
        min_linewidth = 0.3  # * int(float(biggest_cluster_size) / 30.)
        max_linewidth = 12 * int(float(biggest_cluster_size) / 30.)
        for csize, cluster_group in itertools.groupby(sorted_clusters, key=lambda c: len(c)):
            # if csize < 30:
            #     break
            cluster_group = sorted(list(cluster_group), key=lambda c: numpy.median(getnmutelist(c)))
            n_clusters = len(cluster_group)
            print ' %3d' % csize
            for iclust in range(len(cluster_group)):
                cluster = cluster_group[iclust]
                nmutelist = sorted(getnmutelist(cluster))

                nbins = nmutelist[-1] - nmutelist[0] + 1
                hist = Hist(nbins, nmutelist[0] - 0.5, nmutelist[-1] + 0.5)
                for nm in nmutelist:
                    hist.fill(nm)
                assert hist.overflow_contents() == 0.  # includes underflows
                xmax = float(csize)
                yval = csize
                if iclust > 0:
                    yval += (-1)**(iclust) * 0.15

                print '      %8.2f %s' % (numpy.median(nmutelist), ' '.join(['%3d' % nm for nm in nmutelist]))

                for ibin in range(1, hist.n_bins + 1):
                    linewidth = gety(min_linewidth, max_linewidth, xmax, hist.bin_contents[ibin])
                    color = colors[iclust % len(colors)]
                    alpha = 0.55
                    if hist.bin_contents[ibin] == 0.:
                        color = 'grey'
                        alpha = 0.4
                    ax.plot([hist.low_edges[ibin], hist.low_edges[ibin+1]], [yval, yval], color=color, linewidth=linewidth, alpha=alpha, solid_capstyle='butt')

        ybounds = [0.9 * len(sorted_clusters[-1]), 1.05 * len(sorted_clusters[0])]
        plotnames = {'size-vs-shm' : [-0.2, 50], 'size-vs-shm-zoomed' : [-0.2, 60]}
        for name, xbounds in plotnames.items():
            plotting.mpl_finish(ax, base_plotdir + '/overall', name, xlabel='N mutations', ylabel='clonal family size', xbounds=xbounds, ybounds=ybounds)

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
