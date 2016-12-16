import itertools
import sys
import time

from hist import Hist
import utils
import plotting
from clusterpath import ClusterPath

# ----------------------------------------------------------------------------------------
class PartitionPlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self):
        # between vs within stuff:
        self.n_bins = 30

    # ----------------------------------------------------------------------------------------
    def get_cdr3_length_classes(self, partition, annotations):
        classes = {}
        for cluster in partition:
            info = annotations[':'.join(cluster)]
            if info['cdr3_length'] not in classes:
                classes[info['cdr3_length']] = []
            classes[info['cdr3_length']].append(cluster)
        return classes

    # ----------------------------------------------------------------------------------------
    def get_within_vs_between_hists(self, partition, annotations):
        classes = self.get_cdr3_length_classes(partition, annotations)

        distances = {'within' : [mut_freq for info in annotations.values() for mut_freq in info['mut_freqs']],
                     'between' : []}
        # for each cdr3 length, loop over each pair of clusters that have that cdr3 length
        for cdr3_length, clusters in classes.items():
            for cluster_a, cluster_b in itertools.combinations(clusters, 2):
                nseq_a = annotations[':'.join(cluster_a)]['naive_seq']
                nseq_b = annotations[':'.join(cluster_b)]['naive_seq']
                distances['between'].append(utils.hamming_fraction(nseq_a, nseq_b))

        xmax = 1.2 * max(distances['between'] + distances['within'])
        hists = {}
        for dtype in distances:
            hists[dtype] = Hist(self.n_bins, 0., xmax, title=dtype)
            for mut_freq in distances[dtype]:
                hists[dtype].fill(mut_freq)

        return hists['within'], hists['between']

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, partition=None, infiles=None, annotations=None, only_csv=None):
        print '  plotting partitions'
        sys.stdout.flush()
        start = time.time()
        utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'])

        if partition is not None:  # one partition
            assert infiles is None
            assert annotations is not None
            csize_hists = {'best' : plotting.get_cluster_size_hist(partition)}
            hwithin, hbetween = self.get_within_vs_between_hists(partition, annotations)
            plotting.draw_no_root(hwithin, plotname='within-vs-between', plotdir=plotdir, more_hists=[hbetween], xtitle='hamming distance', errors=True)
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

        plotting.plot_cluster_size_hists(plotdir + '/cluster-sizes.svg', csize_hists, title='')

        if not only_csv:
            plotting.make_html(plotdir)

        print '(%.1f sec)' % (time.time()-start)
