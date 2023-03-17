import math
import sys

from hist import Hist
import utils

# ----------------------------------------------------------------------------------------
def get_expanded_bounds(values, dxmin, dxmax=None, only_down=False):  # NOTE see also plotting.expand_bounds()
    def dfcn(): return max(1e-5, abs(values[0]))  # hackey, but effective
    values = sorted(values)  # should already be sorted, but this is the only way to enforce it
    if dxmax is None:
        dxmax = dxmin
    if dxmin == 0:
        dxmin = dfcn()
    if dxmax == 0:
        dxmax = dfcn()
    xmin = values[0] - 0.1*dxmin  # expand a little to avoid underflows
    if only_down:
        xmax = values[-1]
    else:
        xmax = values[-1] + 0.1*dxmax  # and overflows
    return [xmin, xmax]

# ----------------------------------------------------------------------------------------
def set_bins(values, n_bins, is_log_x, xbins, var_type='float'):  # NOTE this fcn/signature is weird because it mimics an old root fcn (for search: log_bins log bins) UPDATE see autobins() fcn below
    assert len(values) > 0
    values = sorted(values)
    if is_log_x:  # NOTE similarity to get_auto_y_ticks() and get_cluster_size_xticks()
        log_xmin = math.log(pow(10.0, math.floor(math.log10(values[0]))))  # round down to the nearest power of 10 from the first entry
        log_xmax = math.log(pow(10.0, math.ceil(math.log10(values[-1]))))  # round up to the nearest power of 10 from the last entry
        log_dx = (log_xmax - log_xmin) / float(n_bins)
        log_xmin -= 0.1*log_dx  # expand a little to avoid overflows
        log_xmax += 0.1*log_dx
        log_dx = (log_xmax - log_xmin) / float(n_bins)
        for ib in range(n_bins+1):
            low_edge = math.exp(log_xmin + ib*log_dx)
            xbins[ib] = low_edge
    else:
        xmin, xmax = 0, 0
        if var_type == 'string':
            xmin = 0.5
            xmax = n_bins + 0.5
        else:
            xmin, xmax = get_expanded_bounds(values, float(values[-1] - values[0]) / n_bins)
        dx = float(xmax - xmin) / n_bins  # then recalculate dx
        if dx < 1e-10:
            print '  %s very small difference %f between xmin %f and xmax %f in hutils.set_bins() (values: %s)' % (utils.wrnstr(), dx, xmin, xmax, values)
        for ib in range(n_bins+1):
            xbins[ib] = xmin + ib*dx

# ----------------------------------------------------------------------------------------
def autobins(values, n_bins, is_log=False, var_type='float'):  # just making a better interface for set_bins()
    xbins = [0. for _ in range(n_bins+1)]  # the +1 is for the lower edge of the overflow bin (hist.scratch_init() adds low edge of underflow)
    set_bins(values, n_bins, is_log, xbins, var_type=var_type)
    return xbins

# ----------------------------------------------------------------------------------------
def binprint(xbins, values):
    print '    chose %d bins (%d bin low edges, including low edge of overflow) using %d values with min/max %f %f' % (len(xbins) - 1, len(xbins), len(values), min(values), max(values))
    print '      %s' % ' '.join('%8d'%(i+1) for i in range(len(xbins)))  # start at 1, to match hist.py/root convention that 0 is underflow (whose low edge is added in hist scratch init fcn, since its value has no effect)
    print '      %s' % ' '.join('%8.4f'%x for x in xbins)
    print '      %s' % ' '.join('%8d'%(len([v for v in values if v >= xbins[i] and v < xbins[i+1]])) for i in range(len(xbins) - 1))

# ----------------------------------------------------------------------------------------
def auto_bin_expand(values, xbins, int_bins=False, debug=False):
    padval = 0.5 if int_bins else 0.01 * abs(max(values) - min(values))
    if xbins[0] > values[0] - padval:  # make sure the hi edge of the underflow bin (lo edge of first bin) is below the smallest value
        print '      expanding lower edge of first bin (hi edge of underflow): %.3f --> %.3f' % (xbins[0], values[0] - padval)
        xbins[0] = values[0] - padval
    if xbins[-1] < values[-1] + padval:  # make sure the lo edge of the overflow bin is above the largest value
        print '      expanding upper edge of last bin (lo edge of overflow): %.3f --> %.3f' % (xbins[-1], values[-1] + padval)
        xbins[-1] = values[-1] + padval
    if debug:
        binprint(xbins, values)

# ----------------------------------------------------------------------------------------
# attempts to make <n_bins> bins with equal numbers of entries per bin (may return fewer)
def auto_volume_bins(values, n_bins, int_bins=False, min_xdist=None, debug=False):
    n_per_bin = int(len(values) / float(n_bins))
    values = sorted([v for v in values if v is not None])
    if debug > 1:
        print '  trying to divide %d values%s from %s to %s into %d bins (%d per bin) %s' % (len(values), ' (with int bins)' if int_bins else '', values[0], values[-1], n_bins, n_per_bin, values if len(values) < 150 else '')
    if int_bins:
        if len(set(values)) < n_bins:
            print '    %s SHOULD PROBABLY MAYBE reduce n_bins %d to the number of different values %d' % (utils.wrnstr(), n_bins, len(set(values)))
        xbins = [values[0] - 0.5]
        n_this_bin = 0
        for iv, val in enumerate(values):  # similar to the ibins = below, but we have to deal with ties here
            n_this_bin += 1
            potential_hi_edge = val + 0.5
            if n_this_bin < n_per_bin:
                if debug > 1: print '  %3d %5.1f  too few %d'  % (iv, potential_hi_edge, n_this_bin)
                continue
            if min_xdist is not None and potential_hi_edge - xbins[-1] < min_xdist:
                if debug > 1: print '  %3d %5.1f  too close %.1f'  % (iv, potential_hi_edge, potential_hi_edge - xbins[-1])
                continue
            xbins.append(potential_hi_edge)
            n_this_bin = 0
            if debug > 1: print '  %3d %5.1f  adding hi edge for bin with %d entries'  % (iv, potential_hi_edge, n_this_bin)
        if xbins[-1] != values[-1] + 0.5:  # add the lo edge of the overflow bin
            xbins.append(values[-1] + 0.5)
            if debug > 1: print '  %3d %5.1f  adding hi edge for bin with %d entries'  % (len(values) - 1, values[-1] + 0.5, n_this_bin)
        n_bins = len(xbins) - 1
    else:
        if len(set(values)) > 2 * n_bins:  # if there's many more values than bins, do actual auto volume bins (~same number of entries per bin)
            ibins = [min(i * n_per_bin, len(values) - 1) for i in range(n_bins + 1)]  # indices (in sorted values) of bin boundaries that have ~equal entries per bin
            xbins = [values[i] for i in ibins]
        else:
            print '    %s number of distinct values %d less than or equal to requested number of bins %d, so using one bin (n_bins = 2) in/instead of auto volume bins' % (utils.wrnstr(), len(set(values)), n_bins)
            xbins = [values[0], values[-1]]
            n_bins = 1
            debug = True  # turn debug on just to make more clear what's going on
        dxmin, dxmax = [abs(float(xbins[ist+1] - xbins[ist])) for ist in (0, len(xbins) - 2)]  # width of (first, last) bin [not under/overflows]
        xbins[0], xbins[-1] = get_expanded_bounds(values, dxmin, dxmax=dxmax)
    if len(set(xbins)) != len(xbins):
        print '    %s duplicate xbins in auto volume bins, so removing duplicates (and reducing n_bins)' % utils.wrnstr()
        xbins = sorted(set(xbins))
        n_bins = len(xbins) - 1
    if debug:
        binprint(xbins, values)
    assert len(xbins) == n_bins + 1  # will cause the n bin reduction to crash
    return xbins, n_bins

# ----------------------------------------------------------------------------------------
def make_hist_from_list_of_values(vlist, var_type, hist_label, is_log_x=False, xmin_force=0.0, xmax_force=0.0, sort_by_counts=False, arg_bins=None):
    vdict = {v : vlist.count(v) for v in set(vlist)}
    return make_hist_from_dict_of_counts(vdict, var_type, hist_label, is_log_x=is_log_x, xmin_force=xmin_force, xmax_force=xmax_force, sort_by_counts=sort_by_counts, arg_bins=arg_bins)

# ----------------------------------------------------------------------------------------
# <values> is of form {<bin 1>:<counts 1>, <bin 2>:<counts 2>, ...}
def make_hist_from_dict_of_counts(values, var_type, hist_label, is_log_x=False, xmin_force=0.0, xmax_force=0.0, sort_by_counts=False, no_sort=False, default_n_bins=30, arg_bins=None):  # default_n_bins is only used if is_log_x set we're doing auto log bins
    """ Fill a histogram with values from a dictionary (each key will correspond to one bin) """
    assert var_type == 'int' or var_type == 'string'  # floats should be handled by Hist class in hist.py

    if len(values) == 0:
        print 'WARNING no values for %s in make_hist' % hist_label
        return Hist(1, 0, 1)

    if not no_sort:
        bin_labels = sorted(values)  # by default sort by keys in dict (i.e. these aren't usually actually string "labels")
    else:
        bin_labels = values.keys()
    if sort_by_counts:  # instead sort by counts
        bin_labels = sorted(values, key=values.get, reverse=True)

    if arg_bins is not None:
        n_bins = len(arg_bins) - 1
    elif var_type == 'string':
        n_bins = len(values)
    else:
        n_bins = bin_labels[-1] - bin_labels[0] + 1 if not is_log_x else default_n_bins

    hist = None
    xbins = [0. for _ in range(n_bins+1)]  # NOTE the +1 is 'cause you need the lower edge of the overflow bin (hist.scratch_init() adds low edge of underflow)
    if xmin_force == xmax_force:  # if boundaries aren't set explicitly, work out what they should be
        if var_type == 'string':
            set_bins(bin_labels, n_bins, is_log_x, xbins, var_type)
            hist = Hist(n_bins, xbins[0], xbins[-1], xbins=xbins)
        else:
            if arg_bins is not None:
                hist = Hist(n_bins, arg_bins[0], arg_bins[-1], xbins=arg_bins)
            elif is_log_x:  # get automatic log-spaced bins
                set_bins(bin_labels, n_bins, is_log_x, xbins, var_type)
                hist = Hist(n_bins, xbins[0], xbins[-1], xbins=xbins)
            else:
                hist = Hist(n_bins, bin_labels[0] - 0.5, bin_labels[-1] + 0.5)  # for integers, just go from the first to the last bin label (they're sorted)
    else:
      hist = Hist(n_bins, xmin_force, xmax_force)

    for ival in range(len(values)):
        if var_type == 'string':
            label = bin_labels[ival]
            ibin = ival + 1
        else:
            label = ''
            ibin = hist.find_bin(bin_labels[ival])
        hist.set_ibin(ibin, values[bin_labels[ival]], error=math.sqrt(values[bin_labels[ival]]), label=label)

    # make sure there's no overflows
    if hist.bin_contents[0] != 0.0 or hist.bin_contents[-1] != 0.0:
        for ibin in range(hist.n_bins + 2):
            print '%d %f %f' % (ibin, hist.low_edges[ibin], hist.bin_contents[ibin])
        raise Exception('overflows in ' + hist_label)

    return hist
