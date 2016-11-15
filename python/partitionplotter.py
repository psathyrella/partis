import sys
import time

import utils
import plotting
from clusterpath import ClusterPath

# ----------------------------------------------------------------------------------------
class PartitionPlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, cpath=None, infiles=None):
        self.cpaths = []
        if cpath is not None:
            self.cpaths.append(cpath)
        elif infiles is not None:  # specify one clusterpath per input file
            for fname in infiles:
                self.cpaths.append(ClusterPath())
                self.cpaths[-1].readfile(fname)

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, only_csv):
        if len(self.cpaths) == 0:
            print '  no partitions to plot'
            return

        print '  plotting partitions'
        sys.stdout.flush()
        start = time.time()
        
        utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'])

        if len(self.cpaths) == 1:
            csize_hists = {'best' : plotting.get_cluster_size_hist(self.cpaths[0].partitions[self.cpaths[0].i_best])}
        else:
            subset_hists = [plotting.get_cluster_size_hist(cpath.partitions[cpath.i_best]) for cpath in self.cpaths]
            csize_hists = {'best' : plotting.make_mean_hist(subset_hists)}
            for ih in range(len(subset_hists)):
                htmp.write(plotdir + ('/subset-%d-cluster-sizes.csv' % ih))

        plotting.plot_cluster_size_hists(plotdir + '/cluster-sizes.svg', csize_hists, title='')

        if not only_csv:
            plotting.make_html(plotdir)

        print '(%.1f sec)' % (time.time()-start)
