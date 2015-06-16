import csv
import math
from opener import opener
from utils import is_normed

# ----------------------------------------------------------------------------------------
class Hist(object):
    """ a simple histogram """
    def __init__(self, n_bins=None, xmin=None, xmax=None, sumw2=False, xbins=None, fname=None):  # if <sumw2>, keep track of errors with <sum_weights_squared>
        self.low_edges, self.bin_contents, self.bin_labels = [], [], []
        self.xtitle, self.ytitle = '', ''

        if fname is None:
            self.scratch_init(n_bins, xmin, xmax, sumw2=sumw2, xbins=xbins)
        else:
            self.file_init(fname)

    # ----------------------------------------------------------------------------------------
    def scratch_init(self, n_bins, xmin, xmax, sumw2=None, xbins=None):
        self.n_bins = int(n_bins)
        self.xmin, self.xmax = float(xmin), float(xmax)
        self.errors = None if sumw2 else []
        self.sum_weights_squared = [] if sumw2 else None

        if xbins is not None:  # check validity of handmade bins
            assert len(xbins) == self.n_bins + 1
            assert self.xmin == xbins[0]
            assert self.xmax == xbins[-1]

        dx = 0.0 if self.n_bins == 0 else (self.xmax - self.xmin) / self.n_bins
        for ib in range(self.n_bins + 2):  # using ROOT conventions: zero is underflow and last bin is overflow
            self.bin_labels.append('')
            if xbins is None:  # uniform binning
                self.low_edges.append(self.xmin + (ib-1)*dx)  # subtract one from ib so underflow bin has upper edge xmin. NOTE this also means that <low_edges[-1]> is the lower edge of the overflow
            else:  # handmade bins
                if ib == 0:
                    self.low_edges.append(xbins[0] - dx)  # low edge of underflow needs to be less than xmin, but is otherwise arbitrary, so just choose something that kinda makes sense
                else:
                    self.low_edges.append(xbins[ib-1])
            self.bin_contents.append(0.0)
            if sumw2:
                self.sum_weights_squared.append(0.)
            else:
                self.errors.append(0.)  # don't set the error values until we <write> (that is unless you explicitly set them with <set_ibin()>

    # ----------------------------------------------------------------------------------------
    def file_init(self, fname):
        self.errors, self.sum_weights_squared = [], []  # kill the unused one after reading file
        with opener('r')(fname) as infile:
            reader = csv.DictReader(infile)
            for line in reader:
                self.low_edges.append(float(line['bin_low_edge']))
                self.bin_contents.append(float(line['contents']))
                if 'sum-weights-squared' in line:
                    self.sum_weights_squared.append(float(line['sum-weights-squared']))
                if 'error' in line or 'binerror' in line:  # in theory I should go find all the code that writes these files and make 'em use the same header for this
                    assert 'sum-weights-squared' not in line
                    tmp_error = float(line['error']) if 'error' in line else float(line['binerror'])
                    self.errors.append(tmp_error)
                if 'binlabel' in line:
                    self.bin_labels.append(line['binlabel'])
                else:
                    self.bin_labels.append('')
                if 'xtitle' in line:  # should be the same for every line in the file... but this avoids complicating the file format
                    self.xtitle = line['xtitle']

        self.n_bins = len(self.low_edges) - 2  # file should have a line for the under- and overflow bins
        self.xmin, self.xmax = self.low_edges[1], self.low_edges[-1]  # *upper* bound of underflow, *lower* bound of overflow

        assert sorted(self.low_edges) == self.low_edges
        assert len(self.bin_contents) == len(self.low_edges)
        assert len(self.low_edges) == len(self.bin_labels)
        if len(self.errors) == 0:  # (re)set to None if the file didn't have errors listed
            self.errors = None
            assert len(self.sum_weights_squared) == len(self.low_edges)
        if len(self.sum_weights_squared) == 0:
            self.sum_weights_squared = None
            assert len(self.errors) == len(self.low_edges)

    # ----------------------------------------------------------------------------------------
    def set_ibin(self, ibin, value, error=None, label=''):
        """ set <ibin>th bin to <value> """
        self.bin_contents[ibin] = value
        if error is not None:
            assert self.errors is not None  # you shouldn't have set sumw2 to True in the constructor
            self.errors[ibin] = error
        self.bin_labels[ibin] = label

    # ----------------------------------------------------------------------------------------
    def fill_ibin(self, ibin, weight=1.0):
        """ fill <ibin>th bin with <weight> """
        self.bin_contents[ibin] += weight
        if self.sum_weights_squared is not None:
            self.sum_weights_squared[ibin] += weight*weight
        if self.errors is not None:
            if weight != 1.0:
                print 'WARNING using errors instead of sumw2 with weight != 1.0 in Hist::fill_ibin()'
            self.errors[ibin] = math.sqrt(self.bin_contents[ibin])

    # ----------------------------------------------------------------------------------------
    def find_bin(self, value):
        """ find <ibin> corresponding to <value>. NOTE boundary is owned by the upper bin. """
        if value < self.low_edges[0]:  # is it below the low edge of the underflow?
            return 0
        elif value >= self.low_edges[self.n_bins + 1]:  # or above the low edge of the overflow?
            return self.n_bins + 1
        else:
            for ib in range(self.n_bins + 2):  # loop over all the bins (including under/overflow)
                if value >= self.low_edges[ib] and value < self.low_edges[ib+1]:  # NOTE <ib> never gets to <n_bins> + 1 because we already get all the overflows above (which is good 'cause this'd fail with an IndexError)
                    return ib

    # ----------------------------------------------------------------------------------------
    def fill(self, value, weight=1.0):
        """ fill bin corresponding to <value> with <weight> """
        self.fill_ibin(self.find_bin(value), weight)

    # ----------------------------------------------------------------------------------------
    def normalize(self, overflow_warn=True):  # since when you normalize hists you have to make the arbitrary decision whether you're going to include the under/overflow bins (we don't include them here), in general we prefer to avoid having under/overflow entries
        """ NOTE does not multiply/divide by bin widths """
        sum_value = 0.0
        for ib in range(1, self.n_bins + 1):  # don't include under/overflows
            sum_value += self.bin_contents[ib]
        if sum_value == 0.0:
            print 'WARNING sum zero in Hist::normalize(), returning without doing anything'
            return
        # make sure there's not too much stuff in the under/overflows
        if overflow_warn and (self.bin_contents[0]/sum_value > 1e-10 or self.bin_contents[self.n_bins+1]/sum_value > 1e-10):
            print 'WARNING under/overflows in Hist::normalize()'
        for ib in range(1, self.n_bins + 1):
            self.bin_contents[ib] /= sum_value
            if self.sum_weights_squared is not None:
                self.sum_weights_squared[ib] /= sum_value*sum_value
            if self.errors is not None:
                self.errors[ib] /= sum_value
        check_sum = 0.0
        for ib in range(1, self.n_bins + 1):  # check it
            check_sum += self.bin_contents[ib]
        assert is_normed(check_sum, this_eps=1e-10)

    # ----------------------------------------------------------------------------------------
    def divide_by(self, denom_hist, debug=False):
        """ NOTE doesn't check bin edges are the same, only that they've got the same number of bins """
        if self.n_bins != denom_hist.n_bins or self.xmin != denom_hist.xmin or self.xmax != denom_hist.xmax:
            raise Exception('ERROR bad limits in Hist::divide_by')
        for ib in range(0, self.n_bins + 2):
            if debug:
                print ib, self.bin_contents[ib], float(denom_hist.bin_contents[ib])
            if denom_hist.bin_contents[ib] == 0.0:
                self.bin_contents[ib] = 0.0
            else:
                self.bin_contents[ib] /= float(denom_hist.bin_contents[ib])

    # ----------------------------------------------------------------------------------------
    def add(self, h2, debug=False):
        """ NOTE doesn't check bin edges are the same, only that they've got the same number of bins """
        if self.n_bins != h2.n_bins or self.xmin != h2.xmin or self.xmax != h2.xmax:
            raise Exception('ERROR bad limits in Hist::add')
        for ib in range(0, self.n_bins + 2):
            if debug:
                print ib, self.bin_contents[ib], float(h2.bin_contents[ib])
            self.bin_contents[ib] += h2.bin_contents[ib]

    # ----------------------------------------------------------------------------------------
    def write(self, outfname):
        with opener('w')(outfname) as outfile:
            header = [ 'bin_low_edge', 'contents', 'binlabel' ]
            if self.errors is not None:
                header.append('error')
            else:
                header.append('sum-weights-squared')
            writer = csv.DictWriter(outfile, header)
            writer.writeheader()
            for ib in range(self.n_bins + 2):
                row = {'bin_low_edge':self.low_edges[ib], 'contents':self.bin_contents[ib], 'binlabel':self.bin_labels[ib] }
                if self.errors is not None:
                    row['error'] = self.errors[ib]
                else:
                    row['sum-weights-squared'] = self.sum_weights_squared[ib]
                writer.writerow(row)
