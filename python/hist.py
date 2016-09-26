import csv
import math
import os
from opener import opener
from utils import is_normed

# ----------------------------------------------------------------------------------------
class Hist(object):
    """ a simple histogram """
    def __init__(self, n_bins=None, xmin=None, xmax=None, sumw2=False, xbins=None, fname=None, xtitle='', ytitle='', title=''):  # if <sumw2>, keep track of errors with <sum_weights_squared>
        self.low_edges, self.bin_contents, self.bin_labels = [], [], []
        self.xtitle, self.ytitle, self.title = xtitle, ytitle, title

        self.all_data = None

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
        for ib in range(self.n_bins + 2):  # using root conventions: zero is underflow and last bin is overflow
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
    def get_maximum(self, xbounds=None):
        # NOTE includes under/overflows by default
        ibin_start = 0 if xbounds is None else self.find_bin(xbounds[0])
        ibin_end = self.n_bins + 2 if xbounds is None else self.find_bin(xbounds[1])

        if ibin_start == ibin_end:
            return self.bin_contents[ibin_start]

        ymax = None
        for ibin in range(ibin_start, ibin_end):
            if ymax is None or self.bin_contents[ibin] > ymax:
                ymax = self.bin_contents[ibin]

        assert ymax is not None
        return ymax

    # ----------------------------------------------------------------------------------------
    def get_bounds(self, include_overflows=False):
        if include_overflows:
            imin, imax = 0, self.n_bins + 2
        else:
            imin, imax = 1, self.n_bins + 1
        return imin, imax

    # ----------------------------------------------------------------------------------------
    def integral(self, include_overflows):
        """ NOTE does not multiply/divide by bin widths """
        imin, imax = self.get_bounds(include_overflows)
        sum_value = 0.0
        for ib in range(imin, imax):
            sum_value += self.bin_contents[ib]
        return sum_value

    # ----------------------------------------------------------------------------------------
    def normalize(self, include_overflows=True, expect_empty=False, expect_overflows=False, overflow_eps_to_ignore=1e-15):
        sum_value = self.integral(include_overflows)
        imin, imax = self.get_bounds(include_overflows)
        if sum_value == 0.0:
            return
        if sum_value == 0.0:
            if not expect_empty:
                print 'WARNING sum zero in Hist::normalize()'
            return
        if not expect_overflows and not include_overflows and (self.bin_contents[0]/sum_value > overflow_eps_to_ignore or self.bin_contents[self.n_bins+1]/sum_value > overflow_eps_to_ignore):
            print 'WARNING under/overflows in Hist::normalize()'
        for ib in range(imin, imax):
            self.bin_contents[ib] /= sum_value
            if self.sum_weights_squared is not None:
                self.sum_weights_squared[ib] /= sum_value*sum_value
            if self.errors is not None:
                self.errors[ib] /= sum_value
        check_sum = 0.0
        for ib in range(imin, imax):  # check it
            check_sum += self.bin_contents[ib]
        if not is_normed(check_sum, this_eps=1e-10):
            raise Exception('not normalized: %f' % check_sum)

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
        if not os.path.exists(os.path.dirname(outfname)):
            os.makedirs(os.path.dirname(outfname))
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

    # ----------------------------------------------------------------------------------------
    def get_bin_centers(self, ignore_overflows=False):
        bin_centers = []
        for ibin in range(len(self.low_edges)):
            low_edge = self.low_edges[ibin]
            if ibin < len(self.low_edges) - 1:
                high_edge = self.low_edges[ibin + 1]
            else:
                high_edge = low_edge + (self.low_edges[ibin] - self.low_edges[ibin - 1])  # overflow bin has undefined upper limit, so just use the next-to-last bin width
            bin_centers.append(0.5 * (low_edge + high_edge))
        if ignore_overflows:
            return bin_centers[1:-1]
        else:
            return bin_centers

    # ----------------------------------------------------------------------------------------
    def bin_contents_no_zeros(self, value):
        """ replace any zeros with <value> """
        bin_centers_no_zeros = list(self.bin_contents)
        for ibin in range(len(bin_centers_no_zeros)):
            if bin_centers_no_zeros[ibin] == 0.:
                bin_centers_no_zeros[ibin] = value
        return bin_centers_no_zeros

    # ----------------------------------------------------------------------------------------
    def get_mean(self, ignore_overflows=False):
        if ignore_overflows:
            imin, imax = 1, self.n_bins + 1
        else:
            imin, imax = 0, self.n_bins + 2
        centers = self.get_bin_centers()
        total, integral = 0.0, 0.0
        for ib in range(imin, imax):
            # print '    ', centers[ib], self.bin_contents[ib]
            total += self.bin_contents[ib] * centers[ib]
            integral += self.bin_contents[ib]
        # print total, integral
        if integral > 0.:
            return total / integral
        else:
            return 0.

    # ----------------------------------------------------------------------------------------
    def rebin(self, factor):
        print 'TODO implement Hist::rebin()'

    # ----------------------------------------------------------------------------------------
    def __str__(self):
        str_list = []
        for ib in range(len(self.low_edges)):
            str_list += ['    %7.4f  %12.3f'  % (self.low_edges[ib], self.bin_contents[ib]), ]
            if ib == 0:
                str_list += ['   under']
            if ib == len(self.low_edges) - 1:
                str_list += ['   over']
            str_list += ['\n']
        return ''.join(str_list)

    # ----------------------------------------------------------------------------------------
    def mpl_plot(self, ax, ignore_overflows=False, label=None, color=None, alpha=None, linewidth=None, linestyle=None, markersize=None, errors=True):
        if self.integral(include_overflows=(not ignore_overflows)) == 0.0:
            # print '   integral is zero in hist::mpl_plot'
            return None
        if ignore_overflows:
            xvals = self.get_bin_centers()[1:-1]
            yvals = self.bin_contents[1:-1]
            yerrs = self.errors[1:-1]
        else:
            xvals = self.get_bin_centers()
            yvals = self.bin_contents
            yerrs = self.errors

        defaults = {'color' : 'black',
                    'alpha' : 1.,
                    'linewidth' : 2,
                    'linestyle' : '-',
                    'marker' : '.',
                    'markersize' : 13}
        kwargs = {}
        argvars = locals()
        for arg in defaults:
            if arg in argvars and argvars[arg] is not None:
                kwargs[arg] = argvars[arg]
            else:
                kwargs[arg] = defaults[arg]
        if label is not None:
            kwargs['label'] = label
        elif self.title != '':
            kwargs['label'] = self.title
        if errors:
            kwargs['yerr'] = yerrs
            return ax.errorbar(xvals, yvals, **kwargs)  #, fmt='-o')
        else:
            return ax.plot(xvals, yvals, **kwargs)  #, fmt='-o')
