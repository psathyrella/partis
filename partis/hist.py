from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import csv
import math
import os
import numpy
import sys

from . import utils
from io import open

# ----------------------------------------------------------------------------------------
class Hist(object):
    """ a simple histogram """
    # NOTE see autobins() or set_bins() in hutils.py for log/str binning (don't need them for linear/uniform float bins)
    def __init__(self, n_bins=None, xmin=None, xmax=None, sumw2=False, xbins=None, template_hist=None, fname=None, value_list=None, weight_list=None, init_int_bins=False, xtitle='', ytitle='counts', title=''):  # if <sumw2>, keep track of errors with <sum_weights_squared>
        # <xbins>: low edge of all non-under/overflow bins, plus low edge of overflow bin (i.e. low edge of every bin except underflow). Weird, but it's root conventions and it's fine
        self.low_edges, self.bin_contents, self.bin_labels = [], [], []
        self.xtitle, self.ytitle, self.title = xtitle, ytitle, title

        if fname is None:
            if init_int_bins:  # use value_list to initialize integer bins NOTE adding this very late, and i tink there's lots of places where it could be used
                assert value_list is not None
                xmin, xmax = min(value_list) - 0.5, max(value_list) + 0.5
                n_bins = xmax - xmin
            if template_hist is not None:
                assert n_bins is None and xmin is None and xmax is None and xbins is None  # specify *only* the template hist
                n_bins = template_hist.n_bins
                xmin, xmax = template_hist.xmin, template_hist.xmax
                xbins = template_hist.low_edges[1:]
            assert n_bins is not None
            assert xmin is not None and xmax is not None
            self.scratch_init(n_bins, xmin, xmax, sumw2=sumw2, xbins=xbins)
            if template_hist is not None and template_hist.bin_labels.count('') != len(template_hist.bin_labels):
                self.bin_labels = [l for l in template_hist.bin_labels]
        else:
            self.file_init(fname)

        if value_list is not None:
            if any(math.isnan(v) for v in value_list):
                raise Exception('nan value in value_list: %s' % value_list)
            if any(v < xmin or v >= xmax for v in value_list):  # maybe because you forgot that xmax is low edge of overflow bin, so it's included in that
                # NOTE it would be nice to integrate this with hutils.make_hist_from_list_of_values() and hutils.make_hist_from_dict_of_counts(), but there's just too much bin/bound infrastructure in there to make it worthwhile
                obvals = sorted([v for v in value_list if v < xmin or v >= xmax])
                print('        %s %d values outside bounds [%s, %s] in hist list fill: %s' % (utils.color('yellow', 'warning'), len(obvals), xmin, xmax, ' '.join('%.2f'%v for v in obvals)))
            self.list_fill(value_list, weight_list=weight_list)

    # ----------------------------------------------------------------------------------------
    def scratch_init(self, n_bins, xmin, xmax, sumw2=False, xbins=None):
        self.n_bins = int(n_bins)
        self.xmin, self.xmax = float(xmin), float(xmax)
        self.errors = None if sumw2 else []
        self.sum_weights_squared = [] if sumw2 else None

        if xbins is not None:  # check validity of handmade bins
            if len(xbins) != self.n_bins + 1:
                raise Exception('misspecified xbins: should be n_bins + 1 (%d, i.e. the low edges of each non-under/overflow bin plus the low edge of the overflow bin) but got %d' % (self.n_bins + 1, len(xbins)))
            assert self.xmin == xbins[0]
            assert self.xmax == xbins[-1]
            if len(set(xbins)) != len(xbins):
                raise Exception('xbins has duplicate entries: %s' % xbins)

        dx = 0.0 if self.n_bins == 0 else (self.xmax - self.xmin) / float(self.n_bins)
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
        with open(fname, 'r') as infile:
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
    def getdict(self):  # get a dict suitable for writing to json/yaml file (ick! but i don't always want the hists to be in their own file) NOTE code reversing this is in test/cf-tree-metrics.py
        return {'n_bins' : self.n_bins, 'xmin' : self.xmin, 'xmax' : self.xmax, 'bin_contents' : self.bin_contents}

    # ----------------------------------------------------------------------------------------
    def is_overflow(self, ibin):  # return true if <ibin> is either the under or over flow bin
        return ibin in [0, self.n_bins + 1]

    # ----------------------------------------------------------------------------------------
    def overflow_contents(self):
        return self.bin_contents[0] + self.bin_contents[-1]

    # ----------------------------------------------------------------------------------------
    def ibiniter(self, include_overflows, reverse=False):  # return iterator over ibins (adding this late, so could probably be used in a lot of places that it isn't)
        if include_overflows:
            istart, istop = 0, self.n_bins + 2
        else:
            istart, istop = 1, self.n_bins + 1
        step = 1
        if reverse:
            itmp = istart
            istart = istop - 1
            istop = itmp - 1
            step = -1
        return list(range(istart, istop, step))

    # ----------------------------------------------------------------------------------------
    def set_ibin(self, ibin, value, error, label=None):
        """ set <ibin>th bin to <value> """
        self.bin_contents[ibin] = value
        if error is not None:
            if self.errors is None:
                raise Exception('attempted to set ibin error with none type <self.errors>')
            else:
                self.errors[ibin] = error
        if label is not None:
            self.bin_labels[ibin] = label

    # ----------------------------------------------------------------------------------------
    def fill_ibin(self, ibin, weight=1.0):
        """ fill <ibin>th bin with <weight> """
        self.bin_contents[ibin] += weight
        if self.sum_weights_squared is not None:
            self.sum_weights_squared[ibin] += weight*weight
        if self.errors is not None:
            if weight != 1.0:
                print('WARNING using errors instead of sumw2 with weight != 1.0 in Hist::fill_ibin()')
            self.errors[ibin] = math.sqrt(self.bin_contents[ibin])

    # ----------------------------------------------------------------------------------------
    def find_bin(self, value, label=None):  # if <label> is set, find ibin corresponding to bin label <label>
        """ find <ibin> corresponding to <value>. NOTE boundary is owned by the upper bin. """
        if label is not None:
            if label not in self.bin_labels:
                raise Exception('asked for label \'%s\' that isn\'t among bin labels: %s' % (label, self.bin_labels))
            return utils.get_single_entry([i for i, l in enumerate(self.bin_labels) if l==label])
        if value < self.low_edges[0]:  # is it below the low edge of the underflow?
            return 0
        elif value >= self.low_edges[self.n_bins + 1]:  # or above the low edge of the overflow?
            return self.n_bins + 1
        else:
            for ib in range(self.n_bins + 2):  # loop over all the bins (including under/overflow)
                if value >= self.low_edges[ib] and value < self.low_edges[ib+1]:  # NOTE <ib> never gets to <n_bins> + 1 because we already get all the overflows above (which is good 'cause this'd fail with an IndexError)
                    return ib
        print(self)
        raise Exception('couldn\'t find bin for value %f (see lines above)' % value)

    # ----------------------------------------------------------------------------------------
    def fill(self, value, weight=1.0):
        """ fill bin corresponding to <value> with <weight> """
        self.fill_ibin(self.find_bin(value), weight)

    # ----------------------------------------------------------------------------------------
    def list_fill(self, value_list, weight_list=None):
        if weight_list is None:
            for value in value_list:
                self.fill(value)
        else:
            for value, weight in zip(value_list, weight_list):
                self.fill(value, weight=weight)

    # ----------------------------------------------------------------------------------------
    def get_extremum(self, mtype, xbounds=None, exclude_empty=False):  # NOTE includes under/overflows by default for max, but *not* for min
        if xbounds is None:
            if mtype == 'min':
                ibin_start, ibin_end = 1, self.n_bins + 1
            if mtype == 'max':
                ibin_start, ibin_end = 0, self.n_bins + 2
        else:
            ibin_start, ibin_end = [self.find_bin(b) for b in xbounds]

        if ibin_start == ibin_end:
            return self.bin_contents[ibin_start]

        ymin, ymax = None, None
        for ibin in range(ibin_start, ibin_end):
            if exclude_empty and self.bin_contents[ibin] == 0:
                continue
            if ymin is None or self.bin_contents[ibin] < ymin:
                ymin = self.bin_contents[ibin]
        for ibin in range(ibin_start, ibin_end):
            if ymax is None or self.bin_contents[ibin] > ymax:
                ymax = self.bin_contents[ibin]
        if ymin is None and exclude_empty:
            print('  %s couldn\'t find ymin for hist, maybe because <exclude_empty> was set (setting ymin arbitrarily to -99999)' % utils.wrnstr())
            ymin = -99999
        assert ymin is not None and ymax is not None

        if mtype == 'min':
            return ymin
        elif mtype == 'max':
            return ymax
        else:
             assert False

    # ----------------------------------------------------------------------------------------
    def get_maximum(self, xbounds=None):  # NOTE includes under/overflows by default
        return self.get_extremum('max', xbounds=xbounds)

    # ----------------------------------------------------------------------------------------
    def get_minimum(self, xbounds=None, exclude_empty=False):  # NOTE does *not* include under/overflows by default (unlike previous fcn, since we expect under/overflows to be zero)
        return self.get_extremum('min', xbounds=xbounds, exclude_empty=exclude_empty)

    # ----------------------------------------------------------------------------------------
    def get_filled_ibins(self):  # return indices of bins with nonzero contents
        return [i for i, c in enumerate(self.bin_contents) if c > 0.]

    # ----------------------------------------------------------------------------------------
    def get_filled_bin_xbounds(self, extra_pads=0):  # low edge of lowest filled bin, high edge of highest filled bin (for search: "ignores empty bins")
        fbins = self.get_filled_ibins()
        if len(fbins) > 0:
            imin, imax = fbins[0], fbins[-1]
        else:
            imin, imax = 1, self.n_bins + 1  # not sure this is really the best default, but i'm adding it long after writing the rest of the fcn and it seems ok? (note after afterwards: it might be better to return None, but this is better than what i originally had which included under/overflows)
        if extra_pads > 0:  # give a little extra space on either side
            imin = max(0, imin - extra_pads)
            imax = min(len(self.low_edges) - 1, imax + extra_pads)
        return self.low_edges[imin], self.low_edges[imax + 1] if imax + 1 < len(self.low_edges) else self.xmax  # if it's not already the overflow bin, we want the next low edge, otherwise <self.xmax>

    # ----------------------------------------------------------------------------------------
    def get_bounds(self, include_overflows=False):
        if include_overflows:
            imin, imax = 0, self.n_bins + 2
        else:
            imin, imax = 1, self.n_bins + 1
        return imin, imax

    # ----------------------------------------------------------------------------------------
    def binwidth(self, ibin):
        if ibin == 0:  # use width of first bin for underflow
            ibin += 1
        elif ibin == self.n_bins + 1:  # and last bin for overflow
            ibin -= 1
        return self.low_edges[ibin+1] - self.low_edges[ibin]

    # ----------------------------------------------------------------------------------------
    def integral(self, include_overflows, ibounds=None, multiply_by_bin_width=False, multiply_by_bin_center=False):
        """ NOTE by default does not multiply by bin widths """
        if ibounds is None:
            imin, imax = self.get_bounds(include_overflows)
        else:
            imin, imax = ibounds
        sum_value = 0.0
        for ib in range(imin, imax):
            sum_value += self.bin_contents[ib] * (self.binwidth(ib) if multiply_by_bin_width else 1) * (self.get_bin_centers()[ib] if multiply_by_bin_center else 1)
        return sum_value

    # ----------------------------------------------------------------------------------------
    def normalize(self, include_overflows=True, expect_overflows=False, overflow_eps_to_ignore=1e-15, multiply_by_bin_width=False):
        sum_value = self.integral(include_overflows, multiply_by_bin_width=multiply_by_bin_width)
        if multiply_by_bin_width and any(abs(self.binwidth(i)-self.binwidth(1)) > utils.eps for i in self.ibiniter(False)):
            print('  %s normalizing with multiply_by_bin_width set, but bins aren\'t all the same width, which may not work' % utils.wrnstr())  # it would be easy to add but i don't want to test it now
        imin, imax = self.get_bounds(include_overflows)
        if sum_value == 0.0:
            return
        if not expect_overflows and not include_overflows and (self.bin_contents[0]/sum_value > overflow_eps_to_ignore or self.bin_contents[self.n_bins+1]/sum_value > overflow_eps_to_ignore):
            print('WARNING under/overflows in Hist::normalize()')
        for ib in range(imin, imax):
            self.bin_contents[ib] /= sum_value
            if self.sum_weights_squared is not None:
                self.sum_weights_squared[ib] /= sum_value*sum_value
            if self.errors is not None:
                self.errors[ib] /= sum_value
        check_sum = 0.0
        for ib in range(imin, imax):  # check it
            check_sum += self.bin_contents[ib] * (self.binwidth(ib) if multiply_by_bin_width else 1)
        if not utils.is_normed(check_sum, this_eps=1e-10):
            raise Exception('not normalized: %f' % check_sum)
        self.ytitle = 'fraction of %.0f' % sum_value

    # ----------------------------------------------------------------------------------------
    def sample(self, n_vals, include_overflows=False, debug_plot=False):  # draw <n_vals> random numbers from the x axis, according to the probabilities given by the bin contents NOTE similarity to recombinator.choose_vdj_combo()
        assert not include_overflows  # probably doesn't really make sense (since contents of overflows could've been from anywhere below/above, but we'd only return bin center), this is just a way to remind that it doesn't make sense
        self.normalize(include_overflows=include_overflows)  # if this is going to get called a lot with n_vals of 1, this would be slow, but otoh we *really* want to make sure things are normalized with include_overflows the same as it is here
        centers = self.get_bin_centers()
        pvals = numpy.random.uniform(0, 1, size=n_vals)
        return_vals = [None for _ in pvals]
        sum_prob, last_sum_prob = 0., 0.
        for ibin in self.ibiniter(include_overflows):
            sum_prob += self.bin_contents[ibin]
            for iprob, pval in enumerate(pvals):
                if pval < sum_prob and pval >= last_sum_prob:
                    return_vals[iprob] = centers[ibin]
            last_sum_prob = sum_prob
        assert return_vals.count(None) == 0

        if debug_plot:
            from . import plotting
            fig, ax = plotting.mpl_init()
            self.mpl_plot(ax, label='original')
            shist = Hist(value_list=return_vals, init_int_bins=True)
            shist.normalize(include_overflows=False)
            shist.mpl_plot(ax, label='sampled', color='red')
            plotting.mpl_finish(ax, '', 'tmp')

        return return_vals

    # ----------------------------------------------------------------------------------------
    def logify(self, factor):
        for ib in self.ibiniter(include_overflows=True):
            if self.bin_contents[ib] > 0:
                if self.bin_contents[ib] <= factor:
                    raise Exception('factor %f passed to hist.logify() must be less than all non-zero bin entries, but found a bin with %f' % (factor, self.bin_contents[ib]))
                self.bin_contents[ib] = math.log(self.bin_contents[ib] / float(factor))
            if self.errors[ib] > 0:  # I'm not actually sure this makes sense
                self.errors[ib] = math.log(self.errors[ib] / float(factor))

    # ----------------------------------------------------------------------------------------
    def divide_by(self, denom_hist, debug=False):
        """ NOTE doesn't check bin edges are the same, only that they've got the same number of bins """
        if self.n_bins != denom_hist.n_bins or self.xmin != denom_hist.xmin or self.xmax != denom_hist.xmax:
            raise Exception('ERROR bad limits in Hist::divide_by')
        for ib in range(0, self.n_bins + 2):
            if debug:
                print(ib, self.bin_contents[ib], float(denom_hist.bin_contents[ib]))
            if denom_hist.bin_contents[ib] == 0.0:
                self.bin_contents[ib] = 0.0
            else:
                self.bin_contents[ib] /= float(denom_hist.bin_contents[ib])

    # ----------------------------------------------------------------------------------------
    # NOTE if you're here, you may be looking for plotting.make_mean_hist()
    def add(self, h2, debug=False):
        """ NOTE doesn't check bin edges are the same, only that they've got the same number of bins """
        if self.n_bins != h2.n_bins or self.xmin != h2.xmin or self.xmax != h2.xmax:
            raise Exception('ERROR bad limits in Hist::add')
        for ib in range(0, self.n_bins + 2):
            if debug:
                print(ib, self.bin_contents[ib], float(h2.bin_contents[ib]))
            self.bin_contents[ib] += h2.bin_contents[ib]

    # ----------------------------------------------------------------------------------------
    def write(self, outfname):
        if not os.path.exists(os.path.dirname(outfname)):
            os.makedirs(os.path.dirname(outfname))
        with open(outfname, utils.csv_wmode()) as outfile:
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
        bin_contents_no_zeros = list(self.bin_contents)
        for ibin in range(len(bin_contents_no_zeros)):
            if bin_contents_no_zeros[ibin] == 0.:
                bin_contents_no_zeros[ibin] = value
        return bin_contents_no_zeros

    # ----------------------------------------------------------------------------------------
    def get_mean(self, ignore_overflows=False, absval=False, ibounds=None):
        if ibounds is not None:
            imin, imax = ibounds
            if self.integral(False, ibounds=(0, imin)) > 0 or self.integral(False, ibounds=(imax, self.n_bins + 2)) > 0:
                print('  %s called hist.get_mean() with ibounds %s that exclude bins with nonzero entries:  below %.3f   above %.3f' % (utils.color('yellow', 'warning'), ibounds, self.integral(False, ibounds=(0, imin)), self.integral(False, ibounds=(imax, self.n_bins + 2))))
            if imin < 0:
                print('  %s increasing specified imin %d to 0' % (utils.wrnstr(), imin))
                imin = 0
            if imax > self.n_bins + 2:
                print('  %s decreasing specified imax %d to %d' % (utils.wrnstr(), imax, self.n_bins + 2))
                imax = self.n_bins + 2
        elif ignore_overflows:
            imin, imax = 1, self.n_bins + 1
        else:
            imin, imax = 0, self.n_bins + 2
        centers = self.get_bin_centers()
        total, integral = 0.0, 0.0
        for ib in range(imin, imax):
            total += self.bin_contents[ib] * (abs(centers[ib]) if absval else centers[ib])
            integral += self.bin_contents[ib]
        if integral > 0.:
            return total / integral
        else:
            return 0.

    # ----------------------------------------------------------------------------------------
    def rebin(self, factor):
        print('TODO implement Hist::rebin()')

    # ----------------------------------------------------------------------------------------
    def horizontal_print(self, bin_centers=False, bin_decimals=4, contents_decimals=3):
        bin_format_str = '%7.' + str(bin_decimals) + 'f'
        contents_format_str = '%7.' + str(contents_decimals) + 'f'
        binlist = self.get_bin_centers() if bin_centers else self.low_edges
        binline = ''.join([bin_format_str % b for b in binlist])
        contentsline = ''.join([contents_format_str % c for c in self.bin_contents])
        return [binline, contentsline]

    # ----------------------------------------------------------------------------------------
    def __str__(self, print_ibin=False):
        str_list = ['   %s %10s%12s%s'  % ('ibin ' if print_ibin else '', 'low edge', 'contents', '' if self.errors is None else '     err'), '\n', ]
        for ib in range(len(self.low_edges)):
            str_list += ['   %s %10.4f%12.3f'  % ('%4d'%ib if print_ibin else '', self.low_edges[ib], self.bin_contents[ib]), ]
            if self.errors is not None:
                str_list += ['%9.2f' % self.errors[ib]]
            if self.bin_labels.count('') != len(self.bin_labels):
                str_list += ['%12s' % self.bin_labels[ib]]
            if ib == 0:
                str_list += ['   (under)']
            if ib == len(self.low_edges) - 1:
                str_list += ['   (over)']
            str_list += ['\n']
        return ''.join(str_list)

    # ----------------------------------------------------------------------------------------
    # NOTE remove_empty_bins can be a bool (remove/not all empty bins) or a list of length two (remove empty bins outside range)
    def mpl_plot(self, ax, ignore_overflows=False, label=None, color=None, alpha=None, linewidth=None, linestyle=None, markersize=None, errors=True, remove_empty_bins=False,
                 square_bins=False, no_vertical_bin_lines=False):
        # ----------------------------------------------------------------------------------------
        def keep_bin(xvlist, yvlist, ib):
            if isinstance(remove_empty_bins, list):
                xmin, xmax = remove_empty_bins
                if xvlist[ib] > xmin and xvlist[ib] < xmax:
                    return True  # always keep within range
            return yvlist[ib] != 0.
        # ----------------------------------------------------------------------------------------
        def sqbplot(kwargs):
            kwargs['markersize'] = 0
            for ibin in self.ibiniter(include_overflows=False):
                if not keep_bin(self.get_bin_centers(), self.bin_contents, ibin):
                    continue
                tplt = ax.plot([self.low_edges[ibin], self.low_edges[ibin+1]], [self.bin_contents[ibin], self.bin_contents[ibin]], **kwargs)  # horizontal line for this bin
                kwargs['label'] = None  # make sure there's only one legend entry for each hist
                if not no_vertical_bin_lines:
                    ax.plot([self.low_edges[ibin], self.low_edges[ibin]], [self.bin_contents[ibin-1], self.bin_contents[ibin]], **kwargs)  # vertical line from last bin contents
                if errors:
                    bcenter = self.get_bin_centers()[ibin]
                    tplt = ax.plot([bcenter, bcenter], [self.bin_contents[ibin] - self.errors[ibin], self.bin_contents[ibin] + self.errors[ibin]], **kwargs)  # horizontal line for this bin
            if not no_vertical_bin_lines:
                tplt = ax.plot([self.low_edges[ibin+1], self.low_edges[ibin+1]], [self.bin_contents[ibin], 0], **kwargs)  # vertical line for right side of last bin
            return tplt  # not sure if this gets used anywhere?
        # ----------------------------------------------------------------------------------------
        if linewidth is not None:  # i'd really rather this wasn't here, but the error message mpl kicks is spectacularly uninformative so you have to catch it beforehand (when writing the svg file, it throws TypeError: Cannot cast array data from dtype('<U1') to dtype('float64') according to the rule 'safe')
            if not isinstance(linewidth, int):
                raise Exception('have to pass linewidth as int, not %s' % type(linewidth))
        # note: bin labels are/have to be handled elsewhere
        if self.integral(include_overflows=(not ignore_overflows)) == 0.0:
            # print '   integral is zero in hist::mpl_plot'
            return None
        if ignore_overflows:
            xvals = self.get_bin_centers()[1:-1]
            yvals = self.bin_contents[1:-1]
            yerrs = self.errors[1:-1] if self.errors is not None else [math.sqrt(w) for w in self.sum_weights_squared[1:-1]]
        else:
            xvals = self.get_bin_centers()
            yvals = self.bin_contents
            yerrs = self.errors if self.errors is not None else [math.sqrt(w) for w in self.sum_weights_squared]

        defaults = {'color' : 'black',
                    'alpha' : 0.6,
                    'linewidth' : 3,
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
        if remove_empty_bins is not False:  # NOTE can be bool, but can also be list of length two (remove bins only outside those bounds)
            xvals, yvals, yerrs = zip(*[(xvals[iv], yvals[iv], yerrs[iv]) for iv in range(len(xvals)) if keep_bin(xvals, yvals, iv)])
        if errors and not square_bins:
            kwargs['yerr'] = yerrs
            return ax.errorbar(xvals, yvals, **kwargs)  #, fmt='-o')
        else:
            if square_bins:
                return sqbplot(kwargs)
            else:
                return ax.plot(xvals, yvals, **kwargs)  #, fmt='-o')

    # ----------------------------------------------------------------------------------------
    def fullplot(self, plotdir, plotname, pargs={}, fargs={}, texts=None, only_csv=False): #**kwargs):  # i.e. full plotting process, not just the ax.plot type stuff above
        self.write('%s/%s.csv'%(plotdir, plotname))
        if only_csv:
            return
        from . import plotting
        fig, ax = plotting.mpl_init()  # this'll need to be updated when i want to use a kwarg for this fcn
        self.mpl_plot(ax, **pargs)
        if texts is not None:
            for xv, yv, tx in texts:
                fig.text(xv, yv, tx, fontsize=15)
        if 'xticks' not in fargs and any(l != '' for l in self.bin_labels):
            fargs['xticks'] = self.get_bin_centers()
            fargs['xticklabels'] = self.bin_labels
        return plotting.mpl_finish(ax, plotdir, plotname, **fargs)
