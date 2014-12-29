import csv
import math
from opener import opener

# ----------------------------------------------------------------------------------------
class Hist(object):
    """ a simple histogram """
    def __init__(self, n_bins, xmin, xmax):
        self.n_bins = n_bins
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.low_edges = []  # lower edge of each bin
        self.centers = []  # center of each bin
        self.bin_contents = []
        dx = (self.xmax - self.xmin) / self.n_bins
        for ib in range(self.n_bins + 2):  # ROOT conventions: zero is underflow and last bin is overflow
            self.low_edges.append(self.xmin + (ib-1)*dx)  # subtract one from ib so underflow bin has upper edge xmin
            self.centers.append(self.low_edges[-1] + 0.5*dx)
            self.bin_contents.append(0.0)

    # ----------------------------------------------------------------------------------------
    def fill(self, value, weight=1.0):
        if value < self.low_edges[0]:  # underflow
            self.bin_contents[0] += weight
        elif value >= self.low_edges[self.n_bins + 1]:  # overflow
            self.bin_contents[self.n_bins + 1] += weight
        else:
            for ib in range(self.n_bins + 2):  # loop over the rest of the bins
                if value >= self.low_edges[ib] and value < self.low_edges[ib+1]:
                    self.bin_contents[ib] += weight

    # ----------------------------------------------------------------------------------------
    def normalize(self):
        sum_value = 0.0
        for ib in range(1, self.n_bins + 1):  # don't include under/overflows in sum_value
            sum_value += self.bin_contents[ib]
        if sum_value == 0.0:
            print 'WARNING sum zero in Hist::normalize, returning without doing anything'
            return
        # make sure there's not too much stuff in the under/overflows
        if self.bin_contents[0]/sum_value > 1e-10 or self.bin_contents[self.n_bins+1]/sum_value > 1e-10:
            print 'WARNING under/overflows'
        for ib in range(1, self.n_bins + 1):
            self.bin_contents[ib] /= sum_value
        check_sum = 0.0
        for ib in range(1, self.n_bins + 1):  # check it
            print 
            check_sum += self.bin_contents[ib]
        assert math.fabs(check_sum - 1.0) < 1e-10

    # ----------------------------------------------------------------------------------------
    def write(self, outfname):
        with opener('w')(outfname) as outfile:
            writer = csv.DictWriter(outfile, ('bin_low_edge', 'contents'))
            writer.writeheader()
            for ib in range(self.n_bins + 2):
                writer.writerow({'bin_low_edge':self.low_edges[ib], 'contents':self.bin_contents[ib]})
