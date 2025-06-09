from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import math
import sys
import collections

from .hist import Hist
from . import utils
from . import plotting

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
            print('  %s very small difference %f between xmin %f and xmax %f in hutils.set_bins() (values: %s)' % (utils.wrnstr(), dx, xmin, xmax, values))
        for ib in range(n_bins+1):
            xbins[ib] = xmin + ib*dx

# ----------------------------------------------------------------------------------------
def autobins(values, n_bins, is_log=False, var_type='float'):  # just making a better interface for set_bins()
    xbins = [0. for _ in range(n_bins+1)]  # the +1 is for the lower edge of the overflow bin (hist.scratch_init() adds low edge of underflow)
    set_bins(values, n_bins, is_log, xbins, var_type=var_type)
    return xbins

# ----------------------------------------------------------------------------------------
def binprint(xbins, values):
    print('    chose %d bins (%d bin low edges, including low edge of overflow) using %d values with min/max %f %f' % (len(xbins) - 1, len(xbins), len(values), min(values), max(values)))
    print('      %s' % ' '.join('%8d'%(i+1) for i in range(len(xbins))))  # start at 1, to match hist.py/root convention that 0 is underflow (whose low edge is added in hist scratch init fcn, since its value has no effect)
    print('      %s' % ' '.join('%8.4f'%x for x in xbins))
    print('      %s' % ' '.join('%8d'%(len([v for v in values if v >= xbins[i] and v < xbins[i+1]])) for i in range(len(xbins) - 1)))

# ----------------------------------------------------------------------------------------
def auto_bin_expand(values, xbins, int_bins=False, debug=False):
    padval = 0.5 if int_bins else 0.01 * abs(max(values) - min(values))
    if xbins[0] > values[0] - padval:  # make sure the hi edge of the underflow bin (lo edge of first bin) is below the smallest value
        print('      expanding lower edge of first bin (hi edge of underflow): %.3f --> %.3f' % (xbins[0], values[0] - padval))
        xbins[0] = values[0] - padval
    if xbins[-1] < values[-1] + padval:  # make sure the lo edge of the overflow bin is above the largest value
        print('      expanding upper edge of last bin (lo edge of overflow): %.3f --> %.3f' % (xbins[-1], values[-1] + padval))
        xbins[-1] = values[-1] + padval
    if debug:
        binprint(xbins, values)

# ----------------------------------------------------------------------------------------
# attempts to make <n_bins> bins with equal numbers of entries per bin (may return fewer)
def auto_volume_bins(values, n_bins, int_bins=False, min_xdist=None, debug=False):
    if len(values) == 0:
        raise Exception('zero length values passed to auto volume bins')
    n_per_bin = int(len(values) / float(n_bins))
    values = sorted([v for v in values if v is not None])
    if debug > 1:
        print('  trying to divide %d values%s from %s to %s into %d bins (%d per bin) %s' % (len(values), ' (with int bins)' if int_bins else '', values[0], values[-1], n_bins, n_per_bin, values if len(values) < 150 else ''))
    if int_bins:
        if len(set(values)) < n_bins:
            print('        %s reducing n_bins from %d to the number of different values %d in auto volume bins' % (utils.wrnstr(), n_bins, len(set(values))))
            n_bins = len(set(values))
        xbins = [values[0] - 0.5]
        n_this_bin = 0
        for iv, val in enumerate(values):  # similar to the ibins = below, but we have to deal with ties here
            n_this_bin += 1
            potential_hi_edge = val + 0.5
            if n_this_bin < n_per_bin:
                if debug > 1: print('  %3d %5.1f  too few %d'  % (iv, potential_hi_edge, n_this_bin))
                continue
            if min_xdist is not None and potential_hi_edge - xbins[-1] < min_xdist:
                if debug > 1: print('  %3d %5.1f  too close %.1f'  % (iv, potential_hi_edge, potential_hi_edge - xbins[-1]))
                continue
            xbins.append(potential_hi_edge)
            n_this_bin = 0
            if debug > 1: print('  %3d %5.1f  adding hi edge for bin with %d entries'  % (iv, potential_hi_edge, n_this_bin))
        if xbins[-1] != values[-1] + 0.5:  # add the lo edge of the overflow bin
            xbins.append(values[-1] + 0.5)
            if debug > 1: print('  %3d %5.1f  adding hi edge for bin with %d entries'  % (len(values) - 1, values[-1] + 0.5, n_this_bin))
        n_bins = len(xbins) - 1
    else:
        if len(set(values)) > 2 * n_bins:  # if there's many more values than bins, do actual auto volume bins (~same number of entries per bin)
            ibins = [min(i * n_per_bin, len(values) - 1) for i in range(n_bins + 1)]  # indices (in sorted values) of bin boundaries that have ~equal entries per bin
            xbins = [values[i] for i in ibins]
        else:
            print('        %s number of distinct values %d less than or equal to requested number of bins %d, so using one bin (n_bins = 2) in/instead of auto volume bins' % (utils.wrnstr(), len(set(values)), n_bins))
            xbins = [values[0], values[-1]]
            n_bins = 1
            debug = True  # turn debug on just to make more clear what's going on
        dxmin, dxmax = [abs(float(xbins[ist+1] - xbins[ist])) for ist in (0, len(xbins) - 2)]  # width of (first, last) bin [not under/overflows]
        xbins[0], xbins[-1] = get_expanded_bounds(values, dxmin, dxmax=dxmax)
    if len(set(xbins)) != len(xbins):
        print('        %s %d duplicate xbins in auto volume bins, so removing duplicates (and reducing n_bins): %s' % (utils.wrnstr(), len(xbins) - len(set(xbins)), ' '.join('%.1f'%b for b in xbins)))
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
def make_hist_from_dict_of_counts(valdict, var_type, hist_label, is_log_x=False, xmin_force=0.0, xmax_force=0.0, sort_by_counts=False, no_sort=False, default_n_bins=30, arg_bins=None):  # default_n_bins is only used if is_log_x set we're doing auto log bins
    """ Fill a histogram with values from a dictionary (each key will correspond to one bin) """
    assert var_type == 'int' or var_type == 'string'  # floats should be handled by Hist class in hist.py
# TODO rename values to value_dict or something

    if len(valdict) == 0:
        print('WARNING no values for %s in make_hist' % hist_label)
        return Hist(1, 0, 1)

    # get <bin_labels>, which are the labels (if 'string') or values (otherwise) for each bin
    if not no_sort:
        bin_labels = sorted(valdict)  # by default sort by keys in dict (i.e. these aren't usually actually string "labels")
    else:
        bin_labels = list(valdict.keys())
    if sort_by_counts:  # instead sort by counts
        bin_labels = sorted(valdict, key=valdict.get, reverse=True)

    if arg_bins is not None:
        n_bins = len(arg_bins) - 1
    elif var_type == 'string':
        n_bins = len(valdict)
    else:
        n_bins = int(bin_labels[-1]) - int(bin_labels[0]) + 1 if not is_log_x else default_n_bins

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

    filled_bins = set()  # kind of weird to keep track of this, but it's nice to use set_ibin() for the (vast majority of) cases where we only set/fill each bin once
    for ival, lbl in enumerate(bin_labels):  # <lbl> is only really a label for 'string', otherwise it's e.g. the value of the bin's low edge
        if var_type == 'string':
            label = lbl
            ibin = ival + 1
        else:
            label = ''
            ibin = hist.find_bin(lbl)
        fill_val = valdict[lbl]
        if ibin in filled_bins:  # if we've already filled this bin (probably uncommon, e.g. if <rebin> is set in calling fcn)
            if fill_val == int(fill_val):  # if it's an integer, it's better to fill it N times
                for _ in range(int(fill_val)):
                    hist.fill_ibin(ibin)
            else:  # this may kick a warning about errors, which you should probably do something about, but really, don't call this fcn if you have non-integer bin values
                hist.fill_ibin(ibin, weight=fill_val)
        else:
            hist.set_ibin(ibin, fill_val, error=math.sqrt(fill_val), label=label)
        filled_bins.add(ibin)

    # make sure there's no overflows
    if hist.bin_contents[0] != 0.0 or hist.bin_contents[-1] != 0.0:
        for ibin in range(hist.n_bins + 2):
            print('%d %f %f' % (ibin, hist.low_edges[ibin], hist.bin_contents[ibin]))
        raise Exception('overflows in ' + hist_label)

    return hist

# ----------------------------------------------------------------------------------------
var_types = {
    'v_gene' : 'string',
    'd_gene' : 'string',
    'j_gene' : 'string',
    'cdr3_length' : 'int',
    'abundance' : 'int',
    'max-abundance' : 'int',
    'cluster-size' : 'cluster-size',
    'mean-muts' : 'float',
    'diversity' : 'int',
}

# ----------------------------------------------------------------------------------------
# call fcn above repeatedly for each list of values in <value_lists> dict, but first finding the OR of all the bins required, so all returned hists can have the same binning
# atm this is only really used in datascripts/meta/anton-hsv/plot-partition-overlaps.py
def make_hists_from_lists_of_values(value_lists, var_name, var_type=None, is_log_x=False, n_bins=10, rebin=None):
    if var_type is None:
        var_type = var_types[var_name]  # if this crashes, you need to either pass in var_type, or add var_name to var_types above
    assert var_type in ['string', 'int', 'float', 'cluster-size']
    xbins = set()
    for label, vlist in value_lists.items():
        if var_type == 'cluster-size':
            htmp = plotting.make_csize_hist(vlist, n_bins=n_bins)
        elif var_type in ['string', 'int']:
            htmp = make_hist_from_list_of_values(vlist, var_type, label, is_log_x=is_log_x)

        if var_type == 'string':
            xbins |= set(vlist)
        elif var_type in ['int', 'cluster-size']:
            xbins |= set(htmp.low_edges[1:])  # don't include the underflow's low edge
        elif var_type == 'float':
            xbins = [mfn(vlist+([xbins[i]] if len(xbins)>0 else [])) for i, mfn in enumerate([min, max])]  # xbins is just (min, max) for 'float'
        else:
            assert False

    xbins = sorted(xbins)
    if rebin is not None:  # keep first and last edge, otherwise keep only every 1/<rebin>th edge
        xbins = [l for i, l in enumerate(xbins) if i==0 or i==len(xbins)-1 or i%rebin==0]
    if var_type == 'float':
        xbins = plotting.expand_bounds(xbins)
    rhists = collections.OrderedDict()
    for label, values in value_lists.items():
        hfcn = make_hist_from_list_of_values
        if var_type == 'string':
            hfcn = make_hist_from_dict_of_counts
            values = {v : values.count(v) for v in set(values)}
            values.update({v : 0 for v in set(xbins) - set(values)})  # add any that are missing from this hist (with zero counts)
        if var_type == 'cluster-size':
            rhists[label] = plotting.make_csize_hist(values, xbins=xbins)
        elif var_type == 'float':
            rhists[label] = Hist(n_bins=n_bins, xmin=xbins[0], xmax=xbins[1], value_list=values)
        else:
            rhists[label] = hfcn(values, var_type, label, is_log_x=is_log_x, arg_bins=None if var_type=='string' else xbins)
    return rhists

# ----------------------------------------------------------------------------------------
def multi_hist_filled_bin_xbounds(hists):
    xmin, xmax = None, None
    for htmp in hists:
        txb = htmp.get_filled_bin_xbounds()
        if xmin is None or txb[0] < xmin:
            xmin = txb[0]
        if xmax is None or txb[1] > xmax:
            xmax = txb[1]
    return (xmin, xmax)
