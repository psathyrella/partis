from __future__ import unicode_literals

import copy
import pickle
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt

import math
import random
from scipy.interpolate import interp1d
from scipy import stats
import os
import glob
import sys
import csv
import numpy
import subprocess
import operator

import utils
import plotconfig
from hist import Hist
import treeutils

default_colors = ['#006600', '#990012', '#2b65ec', '#cc0000', '#3399ff', '#a821c7', '#808080']
default_linewidths = ['5', '3', '2', '2', '2']

plot_ratios = {
    'v' : (30, 3),
    'd' : (8, 4),
    'j' : (8, 3)
}

# ----------------------------------------------------------------------------------------
def get_normalized_cmap_and_norm(vals, cmap=None):
    if cmap is None:
        cmap = plt.cm.Blues
    sorted_vals = sorted(vals)
    vmin = sorted_vals[0] - 0.2 * (sorted_vals[-1] - sorted_vals[0])  # don't want anybody to be white, so set <vmin> to a bit less than the actual min value (i.e. so white corresponds to a value that's a bit less than any of our values)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=sorted_vals[-1])
    return cmap, norm

# ----------------------------------------------------------------------------------------
def set_bins(values, n_bins, is_log_x, xbins, var_type='float'):
    """ NOTE <values> should be sorted """
    assert len(values) > 0
    values = sorted(values)
    if is_log_x:
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
        dx, xmin, xmax = 0, 0, 0
        if var_type == 'string':
            xmin = 0.5
            xmax = n_bins + 0.5
        else:
            dx = float(values[-1] - values[0]) / n_bins
            if dx == 0:  # hackey, but effective
                dx = max(1e-5, values[0])
            xmin = values[0] - 0.1*dx  # expand a little to avoid underflows
            xmax = values[-1] + 0.1*dx  # and overflows
        dx = float(xmax - xmin) / n_bins  # then recalculate dx
        for ib in range(n_bins+1):
            xbins[ib] = xmin + ib*dx

# ----------------------------------------------------------------------------------------
def write_hist_to_file(fname, hist):
    """ see the make_hist_from* functions to reverse this operation """
    with open(fname, 'w') as histfile:
        writer = csv.DictWriter(histfile, ('bin_low_edge', 'contents', 'binerror', 'xtitle', 'binlabel'))  # this is a really crummy way of writing style information, but root files *suck*, so this is what I do for now
        writer.writeheader()
        for ibin in range(hist.GetNbinsX() + 2):
            writer.writerow({
                'bin_low_edge' : hist.GetXaxis().GetBinLowEdge(ibin),
                'contents' : hist.GetBinContent(ibin),
                'binerror' : hist.GetBinError(ibin),
                'xtitle' : hist.GetXaxis().GetTitle(),
                'binlabel' : hist.GetXaxis().GetBinLabel(ibin)
            })

# ----------------------------------------------------------------------------------------
def make_bool_hist(n_true, n_false, hist_label):
    """ fill a two-bin histogram with the fraction false in the first bin and the fraction true in the second """
    if 'fraction_uncertainty' not in sys.modules:
        import fraction_uncertainty

    hist = Hist(2, -0.5, 1.5, ytitle='freq')

    def set_bin(numer, denom, ibin, label):
        frac = float(numer) / denom
        bounds = sys.modules['fraction_uncertainty'].err(numer, denom)
        err = max(abs(frac - bounds[0]), abs(frac - bounds[1]))
        hist.set_ibin(ibin, frac, error=err, label=label)

    set_bin(n_true, n_true + n_false, 1, 'right')
    set_bin(n_false, n_true + n_false, 2, 'wrong')

    return hist

# ----------------------------------------------------------------------------------------
# <values> is of form {<bin 1>:<counts 1>, <bin 2>:<counts 2>, ...}
def make_hist_from_dict_of_counts(values, var_type, hist_label, log='', xmin_force=0.0, xmax_force=0.0, normalize=False, sort=False):
    """ Fill a histogram with values from a dictionary (each key will correspond to one bin) """
    assert var_type == 'int' or var_type == 'string'  # floats should be handled by Hist class in hist.py

    if len(values) == 0:
        print 'WARNING no values for %s in make_hist' % hist_label
        return Hist(1, 0, 1)

    bin_labels = sorted(values)
    if not sort and var_type == 'string':  # for strings, sort so most common value is to left side
        bin_labels = sorted(values, key=values.get, reverse=True)

    if var_type == 'string':
        n_bins = len(values)
    else:
        n_bins = bin_labels[-1] - bin_labels[0] + 1

    hist = None
    xbins = [0. for _ in range(n_bins+1)]  # NOTE the +1 is 'cause you need the lower edge of the overflow bin
    if xmin_force == xmax_force:  # if boundaries aren't set explicitly, work out what they should be
        if var_type == 'string':
            set_bins(bin_labels, n_bins, 'x' in log, xbins, var_type)
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

    if normalize:
        hist.normalize()
        hist.ytitle = 'freq'
    else:
        hist.ytitle = 'counts'
    
    return hist

# ----------------------------------------------------------------------------------------
def add_bin_labels_not_in_all_hists(hists):
    """ find the OR of all bin labels present in <hists>, and remake each hist in <hists> to have zero bins for any that weren't there already """
    # first convert each hist to a map from bin label to entries
    all_labels = []
    histmaps = []
    for hist in hists:
        histmaps.append({})
        for ibin in range(1, hist.n_bins + 1):  # ignore under/over flows, they're kinda useless for bin-labelled hists
            label = hist.bin_labels[ibin]
            histmaps[-1][label] = (hist.bin_contents[ibin], hist.errors[ibin])  # 2-tuple with (content, error)
            if label not in all_labels:
                all_labels.append(label)

    all_labels = sorted(all_labels)

    # then go through and make new histograms for everybody
    finalhists = []
    for ih in range(len(histmaps)):
        original_hist = hists[ih]
        hmap = histmaps[ih]
        finalhists.append(Hist(len(all_labels), 0.5, len(all_labels) + 0.5, title=original_hist.title))
        for ilabel in range(len(all_labels)):
            label = all_labels[ilabel]
            ibin = ilabel + 1  # root conventions
            finalhists[-1].bin_labels[ibin] = label
            if label in hmap:
                finalhists[-1].bin_contents[ibin] = hmap[label][0]
                finalhists[-1].errors[ibin] = hmap[label][1]
            else:
                finalhists[-1].bin_contents[ibin] = 0.0
                finalhists[-1].errors[ibin] = 0.0

    return finalhists

# ----------------------------------------------------------------------------------------
def shift_hist_overflows(hists, xmin, xmax):
    for htmp in hists:
        if htmp is None:
            continue
        underflows, overflows = 0., 0.
        under_err2, over_err2 = 0., 0.  # sum of squared errors
        first_shown_bin, last_shown_bin = -1, -1
        bin_centers = htmp.get_bin_centers(ignore_overflows=False)
        for ib in range(0, htmp.n_bins + 2):
            if bin_centers[ib] <= xmin:
                underflows += htmp.bin_contents[ib]
                under_err2 += htmp.errors[ib]**2
                htmp.set_ibin(ib, 0., error=0.)
            elif first_shown_bin == -1:
                first_shown_bin = ib
            else:
                break
        for ib in reversed(range(0, htmp.n_bins + 2)):
            if bin_centers[ib] >= xmax:
                overflows += htmp.bin_contents[ib]
                over_err2 += htmp.errors[ib]**2
                htmp.set_ibin(ib, 0., error=0.)
            elif last_shown_bin == -1:
                last_shown_bin = ib
            else:
                break

        htmp.set_ibin(first_shown_bin,
                      underflows + htmp.bin_contents[first_shown_bin],
                      error=math.sqrt(under_err2 + htmp.errors[first_shown_bin]**2))
        htmp.set_ibin(last_shown_bin,
                      overflows + htmp.bin_contents[last_shown_bin],
                      error=math.sqrt(over_err2 + htmp.errors[last_shown_bin]**2))

# ----------------------------------------------------------------------------------------
def draw_no_root(hist, log='', plotdir=None, plotname='foop', more_hists=None, scale_errors=None, normalize=False, bounds=None,
                 figsize=None, shift_overflows=False, colors=None, errors=False, write_csv=False, xline=None, yline=None, xyline=None, linestyles=None,
                 linewidths=None, plottitle=None, csv_fname=None, stats='', translegend=(0., 0.), rebin=None,
                 xtitle=None, ytitle=None, markersizes=None, no_labels=False, only_csv=False, alphas=None, remove_empty_bins=False,
                 square_bins=False):
    assert os.path.exists(plotdir)

    hists = [hist,]
    if more_hists is not None:
        hists = hists + more_hists

    xmin, xmax, ymax = None, None, None
    for htmp in hists:
        if htmp.title == 'null':  # empty hists
            continue
        if scale_errors is not None:
            factor = float(scale_errors[0]) if len(scale_errors) == 1 else float(scale_errors[hists.index(htmp)])
            for ibin in range(htmp.n_bins + 2):
                htmp.errors[ibin] *= factor
        if normalize:  # NOTE removed <normalization_bounds> option, hopefully I'm not using it any more
            htmp.normalize()
        if ymax is None or htmp.get_maximum(xbounds=bounds) > ymax:
            ymax = htmp.get_maximum(xbounds=bounds)
        if xmin is None or htmp.xmin < xmin:  # overridden by <bounds> below
            xmin = htmp.xmin
        if xmax is None or htmp.xmax > xmax:
            xmax = htmp.xmax

    if bounds is not None:
        xmin, xmax = bounds

    if shift_overflows:
        if '_vs_per_gene_support' in plotname or '_fraction_correct_vs_mute_freq' in plotname or plotname in [r + '_gene' for r in utils.regions]:
            print '%s overriding overflow shifting for %s' % (utils.color('yellow', 'warning'), plotname)
        else:
            shift_hist_overflows(hists, xmin, xmax)
        # assert '_vs_per_gene_support' not in plotname and '_fraction_correct_vs_mute_freq' not in plotname and plotname.find('_gene') != 1  # really, really, really don't want to shift overflows for these

    if write_csv:
        assert more_hists is None  # can't write a superposition on multiple hists to a single csv
        if csv_fname is None:
            hist.write(plotdir + '/' + plotname + '.csv')
        else:
            hist.write(csv_fname)

    if only_csv:
        return

    # this is the slow part of plotting (well, writing the svg is also slow)
    fig, ax = mpl_init(figsize=figsize)
    mpl.rcParams.update({'legend.fontsize' : 15})

    tmpcolors = copy.deepcopy(colors)  # don't want to modify the arguments
    if tmpcolors is None:  # fiddle here http://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
        tmpcolors = ['royalblue', 'darkred', 'green', 'darkorange']
    n_different_colors = len(tmpcolors)
    while len(tmpcolors) < len(hists):
        tmpcolors += tmpcolors

    tmplinestyles = [] if linestyles is None or len(linestyles) < len(hists) else copy.deepcopy(linestyles)
    itmp = 0
    availstyles = ['-', '--', '-.', ':']
    while len(tmplinestyles) < len(hists):
        tmplinestyles += [availstyles[itmp % len(availstyles)] for _ in range(n_different_colors)]
        itmp += 1

    for ih in range(len(hists)):
        htmp = hists[ih]
        if stats == 'mean':
            htmp.title += ' (mean %.2f)' % htmp.get_mean()
        elif stats == 'absmean':
            htmp.title += ' (abs av %.2f)' % htmp.get_mean(absval=True)
        elif stats == '0-bin':
            htmp.title += ' (right %.2f)' % htmp.bin_contents[1]
        elif stats is not None and stats != '':  # damnit, I ended up with both of the damn things as possible defaults
            raise Exception('unexpected stats str %s' % stats)
        markersize = None
        if markersizes is not None:
            imark = ih if len(markersizes) > 1 else 0
            markersize = markersizes[imark]
        linewidth = None
        if linewidths is None:
            if ih < 6 and len(hists) > 1:
                linewidth = 6-ih
        else:
            ilw = ih if len(linewidths) > 1 and ih < len(linewidths) else 0
            linewidth = linewidths[ilw]
        if rebin is not None:
            htmp.rebin(rebin)
        alpha = 1.
        if alphas is not None:
            alpha = alphas[ih]
        htmp.mpl_plot(ax, color=tmpcolors[ih], linewidth=linewidth, linestyle=tmplinestyles[ih], ignore_overflows=True, errors=errors, alpha=alpha, markersize=markersize, remove_empty_bins=remove_empty_bins, square_bins=square_bins)

    # NOTE it would be nice to combine xline, yline, and xyline (I don't want to go find everwhere that calls this right now)
    if xline is not None:
        ax.plot([xline, xline], [-0.1*ymax, 0.5*ymax], color='black', linestyle='--', linewidth=3)
    if yline is not None:
        print '%s fix y line' % utils.color('red', 'error')
    if xyline is not None:
        assert len(xyline) == 2
        assert len(xyline[0]) == 2 and len(xyline[1]) == 2
        ax.plot([xyline[0][0], xyline[1][0]], [xyline[0][1], xyline[1][1]], color='black', linestyle='--', linewidth=3)
    # if yline is not None:
    #     # if yline < hframe.GetYaxis().GetXmin() or xline > hframe.GetYaxis().GetXmax():  # make sure we got valid a x position for the line
    #     #     print 'WARNING plotting y line at %f out of bounds (%f, %f)' % (float(ymin), hframe.GetYaxis().GetXmin(), hframe.GetYaxis().GetXmax())
    #     yl = TLine(hframe.GetXaxis().GetXmin(), yline, hframe.GetXaxis().GetXmax(), yline)
    #     yl.Draw()

    xticks, xticklabels = None, None
    if not no_labels and hist.bin_labels.count('') != len(hist.bin_labels):
        xticks = hist.get_bin_centers()
        xticklabels = hist.bin_labels

    if plottitle is not None:
        tmptitle = plottitle
    elif plotname in plotconfig.plot_titles:
        tmptitle = plotconfig.plot_titles[plotname]
    else:
        tmptitle = hist.title  # hm, maybe shouldn't be hist.title? I think that's usually supposed to be the legend

    if xtitle is not None:
        tmpxtitle = xtitle
    elif plotname in plotconfig.xtitles:
        tmpxtitle = plotconfig.xtitles[plotname]
    else:
        tmpxtitle = hist.xtitle  # hm, maybe shouldn't be hist.title? I think that's usually supposed to be the legend

    mpl_finish(ax, plotdir, plotname,
               title=tmptitle,
               xlabel=tmpxtitle,
               ylabel=hist.ytitle if ytitle is None else ytitle,
               xbounds=[xmin, xmax],
               ybounds=[-0.03*ymax, 1.15*ymax],
               leg_loc=(0.72 + translegend[0], 0.7 + translegend[1]),
               log=log, xticks=xticks, xticklabels=xticklabels,
               no_legend=(len(hists) <= 1))
    plt.close()

# ----------------------------------------------------------------------------------------
def get_unified_bin_hist(hists):
    """ 
    Unify bins in <hists>.
    Starts from the bins from <hists[0]>, then loops over the rest of 'em adding bins as it goes (with width from <hists[0]>) so we won't have any under/overflows.
    NOTE totally ignores under/overflows in the original hists. That's on purpose, but like everying else in this foolish thing we call life may in fact turn out to be dumb later on.
    """
    assert len(hists) > 0
    dx = hists[0].GetXaxis().GetBinLowEdge(2) - hists[0].GetXaxis().GetBinLowEdge(1)  # always have at least one bin, in which case this'd be the low edge of the overflow bin minus low edge of the first bin
    # print 'dx:', dx
    low_edges = []
    for ib in range(1, hists[0].GetNbinsX()+1):
        low_edges.append(hists[0].GetXaxis().GetBinLowEdge(ib))

    # for d in [ low_edges[i] - low_edges[i-1] for i in range(1, len(low_edges)) ]:
    #     print ' ', d

    for hist in hists[1:]:
        for ib in range(1, hist.GetNbinsX()+1):
            bincenter = hist.GetXaxis().GetBinCenter(ib)
            while bincenter <= low_edges[0]:  # as long as <bincenter> is outside of the current bounds, keep adding bins on the left...
                low_edges.insert(0, low_edges[0] - dx)
            while bincenter >= low_edges[-1] + dx:  # ...and same thing on the right
                low_edges.insert(len(low_edges), low_edges[-1] + dx)

    return Hist(len(low_edges), low_edges[0], low_edges[-1] + dx)

# ----------------------------------------------------------------------------------------
def get_cluster_size_hist(partition, rebin=None):
    sizes = [len(c) for c in partition]
    nbins = max(sizes)
    # if nbins > 30:
    #     rebin = 2
    if rebin is not None:
        nbins = int(float(nbins) / rebin)
    hist = Hist(nbins, 0.5, max(sizes) + 0.5)
    for sz in sizes:
        hist.fill(sz)
    return hist

# ----------------------------------------------------------------------------------------
def make_mean_hist(hists):
    """ return the hist with bin contents the mean over <hists> of each bin """
    binvals = {}
    for hist in hists:  # I could probably do this with list comprehensions or something, but this way handles different bin bounds
        for ib in range(0, hist.n_bins + 2):
            low_edge = hist.low_edges[ib]
            if low_edge not in binvals:
                binvals[low_edge] = []
            binvals[low_edge].append(hist.bin_contents[ib])
    binlist = sorted(binvals.keys())
    meanhist = Hist(len(binlist) - 2, binlist[1], binlist[-1], xbins=binlist[1 :])
    for ib in range(len(binlist)):
        vlist = binvals[binlist[ib]]
        meanhist.set_ibin(ib, numpy.mean(vlist), error=(numpy.std(vlist, ddof=1) / math.sqrt(len(vlist))))
    # meanhist.normalize()
    return meanhist

# ----------------------------------------------------------------------------------------
def interpolate_values(xvals, yvals):
    """ Replace any instances of None in <yvals> which have non-Non values on both sides with a linear interpolation """
    xvals_no_none, yvals_no_none = [], []
    for ip in range(len(yvals)):
        if yvals[ip] is not None:
            xvals_no_none.append(xvals[ip])
            yvals_no_none.append(yvals[ip])

    fcn = interp1d(xvals_no_none, yvals_no_none)
    for ip in range(len(yvals)):
        if yvals[ip] is None:
            try:
                yvals[ip] = int(fcn([xvals[ip], ])[0])
            except ValueError:
                pass

# ----------------------------------------------------------------------------------------
# NOTE annotation stuff is in plotconfig.py
#
legends = {'vollmers-0.9' : 'VJ CDR3 0.9',
           # 'partition partis' : 'full partis',
           'partition' : 'full partis',
           # 'naive-hamming-partition partis' : 'point partis',
           'naive-hamming-partition' : 'point partis',
           # 'vsearch-partition partis' : 'vsearch partis',
           'vsearch-partition' : 'vsearch partis',
           'seed-partition' : 'full partis (seed)',
           'seed-naive-hamming-partition' : 'point partis (seed)',
           'changeo' : 'IMGT + Change-O',
           # '0.1-true-singletons' : '10% random singletons',
           # '0.1-true-reassign' : '10% random reassign',
           'misassign-0.60-singletons' : 'synth. 60%\nsingleton',
           'misassign-0.10-reassign' : 'synth. 10%\nreassign',
           'misassign-distance-0.03' : 'synth.\nneighbor 0.03',
           'mixcr' : 'MiXCR',
           'adj_mi' : 'adj MI',
           'ccf_under' : 'precision',
           'ccf_over' : 'sensitivity',
           'ccf_product' : 'F1 score'
           }

colors = {'true' : '#006600',
          'partition' : '#cc0000',  # 8c001a',
          'vsearch-partition' : '#990012',  #c04000',
          'naive-hamming-partition' : '#990012',
          'seed-partition' : '#990012',
          'seed-naive-hamming-partition' : '#990012',
          'vollmers-0.5' : '#3333ff',
          'vollmers-0.9' : '#3399ff',
          'changeo' :  '#2b65ec',
          'mixcr' : '#2b65ec',
          'misassign-0.60-singletons' : '#808080',
          'misassign-0.10-reassign' : '#808080',
          'misassign-distance-0.03' : '#808080'
}

linewidths = {'true' : 15,
              'vsearch-partition' : 3,
              'naive-hamming-partition' : 3,
              'seed-partition' : 2,
              'seed-naive-hamming-partition' : 2,
              'partition' : 6,
              'vollmers-0.5' : 4,
              'vollmers-0.9' : 6,
              'changeo' : 3,
              'mixcr' : 6,
              'misassign-0.60-singletons' : 4,
              'misassign-0.10-reassign' : 3,
              'misassign-distance-0.03' : 2
}

linestyles = {'naive-hamming-partition' : 'dashed',
              'vsearch-partition' : 'dotted',
              'changeo' : 'dashed',
              'mixcr' : 'dotted',
              'misassign-distance-0.03' : 'dashed'
}

alphas = {'true' : 0.6,
          'vollmers-0.9' : 0.6,
          'misassign-0.60-singletons' : 0.5,
          'misassign-distance-0.03' : 0.8
}

def label_bullshit_transform(label):
    return '-'.join([hex(int(l)) for l in label.split('-')]).replace('0x', '')

# linewidths['v-true'] = 10
# linewidths['cdr3-true'] = 10
# colors['v-true'] = '#006600'
# colors['cdr3-true'] = '#006600'
# colors['v-indels'] = '#cc0000'
# colors['cdr3-indels'] = '#cc0000'

# ----------------------------------------------------------------------------------------
def plot_cluster_size_hists(outfname, hists, title, xmax=None, log='x', normalize=False):
    if 'seaborn' not in sys.modules:
        import seaborn  # really #$*$$*!ing slow to import, but only importing part of it doesn't seem to help
    sys.modules['seaborn'].set_style('ticks')

    fsize = 20
    mpl.rcParams.update({
        # 'font.size': fsize,
        'legend.fontsize': fsize,
        'axes.titlesize': fsize,
        # 'axes.labelsize': fsize,
        'xtick.labelsize': fsize,
        'ytick.labelsize': fsize,
        'axes.labelsize': fsize
    })

    # base_xvals = hists['true'].get_bin_centers()
    # data = {}
    # for name, hist in hists.items():
    #     hist.normalize()

    #     contents = hist.bin_contents
    #     centers = hist.get_bin_centers()
    #     if len(centers) > len(base_xvals):
    #         print centers, base_xvals
    #         assert False
    #     for ic in range(len(base_xvals)):
    #         if ic < len(centers):
    #             if centers[ic] != base_xvals[ic]:
    #                 print centers
    #                 print base_xvals
    #                 assert False
    #         else:
    #             contents.append(0)

    #     data[name] = contents
        
    fig, ax = plt.subplots()
    fig.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.95)
    # dark red '#A52A2A',
    plots = {}
    scplots = {}
    # tmpmax, n_queries = None, None
    for name, hist in hists.items():
        if 'misassign' in name:
            continue
        if 'vollmers' in name:
            if '0.7' in name or '0.8' in name or '0.95' in name or '0.5' in name:
                continue

        # other (old) possibility:
        # ax.bar and ax.scatter also suck
        # plots[name] = ax.plot(base_xvals, data[name], linewidth=linewidth, label=name, color=colors.get(name, 'grey'), linestyle=linestyle, alpha=alpha)
        # scplots[name] = ax.scatter(hists[name].get_bin_centers(), hists[name].bin_contents, linewidth=linewidths.get(name, 4), label=legends.get(name, name), color=colors.get(name, 'grey'), linestyle=linestyle, alpha=alpha)

        kwargs = {'linewidth' : linewidths.get(name, 4),
                  'label' : legends.get(name, name),
                  'color' : colors.get(name, 'grey'),
                  'linestyle' : linestyles.get(name, 'solid'),
                  'alpha' : alphas.get(name, 1.)}

        # was using this:
        if normalize:
            hist.normalize()
        plots[name] = ax.errorbar(hist.get_bin_centers(), hist.bin_contents_no_zeros(1e-8), yerr=hist.errors, **kwargs)

    legend = ax.legend()
    sys.modules['seaborn'].despine()  #trim=True, bottom=True)

    # xmax = tmpmax
    if xmax is None:
        xmax = plt.gca().get_xlim()[1]
    else:
        ax.set_xlim(0.9, xmax)

    if normalize:
        if 'vollmers' in title:
            ymin = 5e-4
        else:
            ymin = 5e-5
        plt.ylim(ymin, 1)
    plt.title(title)
    plt.xlabel('cluster size')
    plt.ylabel('fraction of clusters' if normalize else 'number of clusters')
    plt.subplots_adjust(bottom=0.14, left=0.14)
    if 'x' in log:
        ax.set_xscale('log')
    if 'y' in log:
        ax.set_yscale('log')
        # if n_queries is not None:
        #     plt.ylim(1./n_queries, 1)
    potential_xticks = [1, 2, 3, 9, 30, 75, 200, 500]
    xticks = [xt for xt in potential_xticks if xt < xmax]
    plt.xticks(xticks, [str(xt) for xt in xticks])
    plotdir = os.path.dirname(outfname)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    plt.savefig(outfname)
    plt.close()
    # subprocess.check_call(['chmod', '664', outfname])

# ----------------------------------------------------------------------------------------
def plot_metrics_vs_thresholds(meth, thresholds, info, plotdir, plotfname, title):
    fig, ax = mpl_init()
    if 'adj_mi' in info and meth in info['adj_mi']:
        ax.plot(thresholds, info['adj_mi'][meth], label='adj mi', linewidth=4)
    ax.plot(thresholds, info['ccf_under'][meth], label='clonal fraction', color='#cc0000', linewidth=4)
    ax.plot(thresholds, info['ccf_over'][meth], label='fraction present', color='#cc0000', linestyle='--', linewidth=4)
    ccf_products = [info['ccf_under'][meth][iv] * info['ccf_over'][meth][iv] for iv in range(len(thresholds))]
    ax.plot(thresholds, ccf_products, label='product of fractions', color='#006600', linewidth=4)
    log, xlabel = '', ''
    if meth == 'partition':
        xlabel = 'log prob ratio'
        ymin = 0.3  #0.69
        xticks = [b for b in range(int(thresholds[0]), int(thresholds[-1]), 5)]
        if int(thresholds[-1]) not in xticks:
            xticks.append(int(thresholds[-1]))
    elif meth == 'naive-hamming-partition' or meth == 'vsearch-partition':
        xlabel = 'naive hamming fraction'
        ymin = 0.3
        xticks = [th for th in thresholds if th < 0.12 and th != 0.025]
        log = 'x'
    mpl_finish(ax, plotdir, plotfname, log=log, xticks=xticks, xticklabels=xticks, leg_loc=(0.1, 0.2), xbounds=(xticks[0], xticks[-1]), ybounds=(ymin, 1.01), title=title, xlabel=xlabel, ylabel='metric value')

# ----------------------------------------------------------------------------------------
def plot_adj_mi_and_co(plotname, plotvals, mut_mult, plotdir, valname, xvar, title=''):
    if 'seaborn' not in sys.modules:
        import seaborn  # really #$*$$*!ing slow to import, but only importing part of it doesn't seem to help
    sys.modules['seaborn'].set_style('ticks')

    # ----------------------------------------------------------------------------------------
    def remove_some_duplicates(xyvals):
        hmap = {}
        newvals = []
        for x, y in xyvals:
            if x in hmap and x != 100000 and x != 500000 and x != 1000000:
                continue
            newvals.append((x, y))
            hmap.add(x)
        return newvals

    fig, ax = mpl_init()
    mpl.rcParams.update({
        'legend.fontsize': 15,})
    plots = {}
    for meth, xyvals in plotvals.items():

        # print sorted([xy[0] for xy in xyvals])
        # xyvals = remove_some_duplicates(xyvals)

        xyvals = sorted(xyvals, key=operator.itemgetter(0))
        xvals = [xy[0] for xy in xyvals]  # xyvals.keys()
        yvals = [ve[1][0] for ve in xyvals]
        yerrs = [ve[1][1] for ve in xyvals]
        kwargs = {'linewidth' : linewidths.get(meth, 4),
                  'label' : legends.get(meth, meth),
                  'color' : colors.get(meth, 'grey'),
                  'linestyle' : linestyles.get(meth, 'solid'),
                  'alpha' : alphas.get(meth, 1.),
                  }

        if meth == 'seed-partition':
            kwargs['linewidth'] = 0
            kwargs['alpha'] = 0.5

        if xvar == 'n_leaves':
            kwargs['fmt'] = '-o'
            plots[meth] = ax.errorbar(xvals, yvals, yerr=yerrs, **kwargs)
        else:  # darn it, the order in the legend gets messed up if I do some as .plot and some as .errorbar
            kwargs['marker'] = '.'
            kwargs['markersize'] = 20
            plots[meth] = ax.plot(xvals, yvals, **kwargs)
    
    lx, ly = 1.6, 0.7
    if len(plotvals) != 1:
        legend = ax.legend(bbox_to_anchor=(lx, ly))
    # legend.get_frame().set_facecolor('white')
    ymin = -0.01
    ax.set_ylim(ymin, 1.03)
    sys.modules['seaborn'].despine()  #trim=True, bottom=True)
    plt.title(title)
    xtitle = 'mean N leaves' if xvar == 'n_leaves' else 'sample size'
    plt.xlabel(xtitle)
    plt.ylabel(legends[valname])
    plt.gcf().subplots_adjust(bottom=0.14, left=0.12, right=0.67, top=0.95)

    xticks = xvals

    # put an 'n/a' in the n_leaves=1 column for adj_mi
    # if 1. not in xticks:
    #     xticks = [1, ] + xticks
    #     x1 = 1.
    #     ax.text(x1 - 0.25, 0.5, 'n/a', color='green', fontsize=25)
    #     ax.plot([x1, x1], [0.05, 0.4], color='green', linewidth=3)
    #     ax.plot([x1, x1], [0.6, 0.97], color='green', linewidth=3)
    # ax.set_xlim(xticks[0] - 0.4, xticks[-1])

    xticks = list(xvals)

    if xvar == 'n_leaves':
        ax.set_xscale('log')
        if 100 in xticks and 200 in xticks:
            xticks.remove(100)
        ax.set_xlim(0.95 * xvals[0], 1.05 * xvals[-1])
    elif xvar == 'nseqs':
        # xticks = [xticks[i] for i in range(0, len(xticks), 2)]
        # if 750 in xticks:
        #     xticks.remove(750)
        # xticks += xvals[-1:]
        # xticks = [100, 5000, 10000, 15000]
        xticks = [1000, 10000, 100000, 1000000]
        ax.set_xscale('log')
        ax.set_xlim(0.9 * xvals[0], 1.15 * xvals[-1])

    xticklabels = xticks if xvar == 'n_leaves' else ['%.0e' % xt for xt in xticks]
    plt.xticks(xticks, xticklabels)
    # ax.plot([xticks[0], xticks[-1]], [1., 1.], linewidth=1, color='grey')
    ax.grid(True)

    yticks = [yt for yt in [0., .2, .4, .6, .8, 1.] if yt >= ymin]
    yticklabels = [str(yt) for yt in yticks]
    plt.yticks(yticks, yticklabels)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    plt.savefig(plotdir + '/' + plotname)
    plt.close()

# ----------------------------------------------------------------------------------------
def mpl_init(figsize=None, fontsize=20):
    if 'seaborn' not in sys.modules:
        import seaborn  # really #$*$$*!ing slow to import, but only importing part of it doesn't seem to help
    sys.modules['seaborn'].set_style('ticks')
    fsize = fontsize
    mpl.rcParams.update({
        # 'legend.fontweight': 900,
        'legend.fontsize': fsize,
        'axes.titlesize': fsize,
        # 'axes.labelsize': fsize,
        'xtick.labelsize': fsize,
        'ytick.labelsize': fsize,
        'axes.labelsize': fsize
    })
    fig, ax = plt.subplots(figsize=figsize)
    fig.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.92)

    return fig, ax

# ----------------------------------------------------------------------------------------
def mpl_finish(ax, plotdir, plotname, title='', xlabel='', ylabel='', xbounds=None, ybounds=None, leg_loc=(0.04, 0.6), leg_prop=None, log='',
               xticks=None, xticklabels=None, yticks=None, yticklabels=None, no_legend=False, adjust=None, suffix='svg', leg_title=None):
    if 'seaborn' not in sys.modules:
        import seaborn  # really #$*$$*!ing slow to import, but only importing part of it doesn't seem to help
    if not no_legend:
        handles, labels = ax.get_legend_handles_labels()
        if len(handles) > 0:
            legend = ax.legend(handles, labels, loc=leg_loc, prop=leg_prop, title=leg_title)
    if adjust is None:
        plt.gcf().subplots_adjust(bottom=0.20, left=0.18, right=0.95, top=0.92)
    else:
        plt.gcf().subplots_adjust(**adjust)
    sys.modules['seaborn'].despine()  #trim=True, bottom=True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if 'x' in log:
        ax.set_xscale('symlog')  # 'log' used to work, but now it screws up the x axis labels
    if 'y' in log:
        ax.set_yscale('log')
    if xbounds is not None and xbounds[0] != xbounds[1]:
        plt.xlim(xbounds[0], xbounds[1])
    if ybounds is not None and ybounds[0] != ybounds[1]:
        plt.ylim(ybounds[0], ybounds[1])
    if xticks is not None:
        plt.xticks(xticks)
    if yticks is not None:
        plt.yticks(yticks)
    if xticklabels is not None:
        # mean_length = float(sum([len(xl) for xl in xticklabels])) / len(xticklabels)
        median_length = numpy.median([len(xl) for xl in xticklabels])
        if median_length > 4:
            ax.set_xticklabels(xticklabels, rotation='vertical', size=8)
        else:
            ax.set_xticklabels(xticklabels)
    if yticklabels is not None:
        ax.set_yticklabels(yticklabels)
    plt.title(title, fontweight='bold')
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    fullname = plotdir + '/' + plotname + '.' + suffix
    plt.savefig(fullname)
    plt.close()
    # subprocess.check_call(['chmod', '664', fullname])

# ----------------------------------------------------------------------------------------
def plot_cluster_similarity_matrix(plotdir, plotname, meth1, partition1, meth2, partition2, n_biggest_clusters, title='', debug=False):
    if debug:
        print ''
        print '%s    %s' % (meth1, meth2)
    # partition1 = [['4'], ['7', '8'], ['6', '5'], ['99', '3', '1']]
    # # partition2 = [['1', '2', '3'], ['4'], ['5', '6'], ['7', '8']]
    # partition2 = [['3'], ['5'], ['6'], ['7'], ['8'], ['99', '3', '4']]
    a_cluster_lengths, b_cluster_lengths, smatrix = utils.partition_similarity_matrix(meth1, meth2, partition1, partition2, n_biggest_clusters=n_biggest_clusters, debug=debug)
    if debug:
        print 'a_clusters: ', ' '.join([str(l) for l in a_cluster_lengths])
        print 'b_clusters: ', ' '.join([str(l) for l in b_cluster_lengths])

    fig, ax = plt.subplots()
    plt.gcf().subplots_adjust(bottom=0.14, left=0.18, right=0.95, top=0.92)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    data = numpy.array(smatrix)
    cmap = plt.cm.Blues  #cm.get_cmap('jet')
    cmap.set_under('w')
    heatmap = ax.pcolor(data, cmap=cmap, vmin=0., vmax=1.)
    cbar = plt.colorbar(heatmap)
    
    modulo = 2
    if n_biggest_clusters > 20:
        modulo = 3
    ticks = [n - 0.5 for n in range(1, n_biggest_clusters + 1, modulo)]
    xticklabels = [b_cluster_lengths[it] for it in range(0, len(b_cluster_lengths), modulo)]
    yticklabels = [a_cluster_lengths[it] for it in range(0, len(a_cluster_lengths), modulo)]
    plt.xticks(ticks, xticklabels)
    plt.yticks(ticks, yticklabels)

    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.unicode'] = True

    def boldify(textstr):
        textstr = textstr.replace('%', '\%')
        textstr = textstr.replace('\n', '')
        return r'\textbf{' + textstr + '}'

    plt.xlabel(boldify(legends.get(meth2, meth2)) + ' cluster size')  # I don't know why it's reversed, it just is
    plt.ylabel(boldify(legends.get(meth1, meth1)) + ' cluster size')
    ax.set_xlim(0, n_biggest_clusters)
    ax.set_ylim(0, n_biggest_clusters)

    plt.title(title)
    
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    plt.savefig(plotdir + '/' + plotname + '.svg')
    plt.close()

# ----------------------------------------------------------------------------------------
def make_html(plotdir, n_columns=3, extension='svg', fnames=None, title='foop', new_table_each_row=False, htmlfname=None, extra_links=None):
    if plotdir[-1] == '/':  # remove trailings slash, if present
        plotdir = plotdir[:-1]
    if not os.path.exists(plotdir):
        raise Exception('plotdir %s d.n.e.' % plotdir)
    dirname = os.path.basename(plotdir)
    extra_link_str = ''
    if extra_links is not None:
        extra_link_str = ' '.join(['<a href=%s>%s</a>' % (url, name) for name, url in extra_links])
    lines = ['<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 3.2//EN>',
             '<html>',
             '<head><title>' + title + '</title></head>',
             '<body bgcolor="000000">', 
             '<h3 style="text-align:left; color:DD6600;">' + title + '</h3>',
             extra_link_str,
             '<table>',
             '<tr>']

    def add_newline(lines):
        if new_table_each_row:
            newlines = ['</tr>', '</table>', '<table>', '<tr>']
        else:
            newlines = ['</tr>', '<tr>']
        lines += newlines
    def add_fname(lines, fullfname):  # NOTE <fullname> may, or may not, be a base name (i.e. it might have a subdir tacked on the left side)
        fname = fullfname.replace(plotdir, '').lstrip('/')
        if htmlfname is None:  # dirname screws it up if we're specifying htmfname explicitly, since then the files are in a variety of different subdirs
            fname = dirname + '/' + fname
        line = '<td><a target="_blank" href="' + fname + '"><img src="' + fname + '" alt="' + fname + '" width="100%"></a></td>'
        lines.append(line)

    # if <fnames> wasn't used to tell us how to group them into rows, try to guess based on the file base names
    if fnames is None:
        fnamelist = [os.path.basename(fn) for fn in sorted(glob.glob(plotdir + '/*.' + extension))]
        fnames = []

        # arrange the ones that have '[vdj]_' into group of three
        for v_fn in [fn for fn in fnamelist if os.path.basename(fn).find('v_') == 0]:  # get the ones that start with 'v_', so we can use them as templates for the others
            fstem = v_fn.replace('v_', '')
            fnames.append([rstr + fstem for rstr in plotconfig.rstrings if rstr + fstem in fnamelist])
            for fn in fnames[-1]:
                fnamelist.remove(fn)

        # and group insertion lengths together
        found_bound_fnames = []
        for bound in utils.all_boundaries:
            for fn in fnamelist:
                if bound + '_insertion' in fn:
                    found_bound_fnames.append(fn)
                    break
        if len(found_bound_fnames) == len(utils.all_boundaries):
            fnames.append(found_bound_fnames)
            for fn in found_bound_fnames:
                fnamelist.remove(fn)

        # then do the rest in groups of <n_columns>
        while len(fnamelist) > 0:
            fnames.append(fnamelist[:n_columns])
            fnamelist = fnamelist[n_columns:]

    # write the meat of the html
    for rowlist in fnames:
        for fn in rowlist:
            add_fname(lines, fn)
        add_newline(lines)

    lines += ['</tr>',
              '</table>',
              '</body>',
              '</html>']

    if htmlfname is None:
        htmlfname = os.path.dirname(plotdir) + '/' + dirname + '.html'  # more verbose than necessary
    with open(htmlfname, 'w') as htmlfile:
        htmlfile.write('\n'.join(lines))
    # subprocess.check_call(['chmod', '664', htmlfname])

# ----------------------------------------------------------------------------------------
def make_allele_finding_plot(plotdir, gene, position, values, xmax, fitfos=None, new_gene=None):
    xmin, xmax = -0.3, xmax
    fig, ax = mpl_init()

    ax.errorbar(values['n_mutelist'], values['freqs'], yerr=values['errs'], markersize=15, linewidth=2, marker='.')  #, title='position ' + str(position))

    if fitfos is not None:  # fitted lines
        colors = {'prefo' : 'red', 'postfo' : 'red', 'onefo' : 'green'}
        for ftype in colors:
            if fitfos[ftype]['xvals'] is None:  # not really sure why this happens... probably zero-point fits?
                continue
            linevals = [fitfos[ftype]['slope']*x + fitfos[ftype]['y_icpt'] for x in fitfos[ftype]['xvals']]
            ax.plot(fitfos[ftype]['xvals'], linevals, color=colors[ftype])

    ax.plot([xmin, xmax], [0, 0], linestyle='dashed', alpha=0.5, color='black')
    ymax = max(values['freqs']) + max(values['errs'])
    title = 'position %d in %s' % (position, gene)
    if new_gene is not None:
        ax.text(0.3 * (xmax - xmin), 0.95 * (ymax - 0), 'inferred: %s' % new_gene, color='red', fontsize=15)
    mpl_finish(ax, plotdir, str(position), xlabel='mutations in %s segment' % utils.get_region(gene), ylabel='position\'s mut freq', xbounds=(xmin, xmax), ybounds=(-0.01, ymax), leg_loc=(0.95, 0.1), adjust={'right' : 0.85}, title=title)

# ----------------------------------------------------------------------------------------
def make_fraction_plot(hright, hwrong, plotdir, plotname, xlabel, ylabel, xbounds, only_csv=False, write_csv=False):
    if 'fraction_uncertainty' not in sys.modules:
        import fraction_uncertainty

    # NOTE should really merge this with draw_no_root()
    xvals = hright.get_bin_centers() #ignore_overflows=True)
    right = hright.bin_contents
    wrong = hwrong.bin_contents
    yvals = [float(r) / (r + w) if r + w > 0. else 0. for r, w in zip(right, wrong)]

    # remove values corresponding to bins with no entries
    while yvals.count(0.) > 0:
        iv = yvals.index(0.)
        xvals.pop(iv)
        right.pop(iv)
        wrong.pop(iv)
        yvals.pop(iv)

    tmphilos = [sys.modules['fraction_uncertainty'].err(r, r + w) for r, w in zip(right, wrong)]
    yerrs = [err[1] - err[0] for err in tmphilos]
    # print '%s' % region
    # for iv in range(len(xvals)):
    #     print '   %5.2f     %5.0f / %5.0f  =  %5.2f   +/-  %.3f' % (xvals[iv], right[iv], right[iv] + wrong[iv], yvals[iv], yerrs[iv])

    if write_csv:
        hist_for_csv = Hist(hright.n_bins, hright.xmin, hright.xmax)
        bincenters = hright.get_bin_centers()
        for ibin in range(hright.n_bins):
            bcenter = bincenters[ibin]
            if bcenter in xvals:  # if we didn't remove it
                iy = xvals.index(bcenter)
                hist_for_csv.set_ibin(ibin, yvals[iy], error=yerrs[iy])

        hist_for_csv.write(plotdir + '/' + plotname + '.csv')

    if not only_csv:
        fig, ax = mpl_init()
        ax.errorbar(xvals, yvals, yerr=yerrs, markersize=10, linewidth=1, marker='.')
        if xlabel == 'support':
            ax.plot((0, 1), (0, 1), color='black', linestyle='--', linewidth=3)  # line with slope 1 and intercept 0
        mpl_finish(ax, plotdir, plotname, xlabel=xlabel, ylabel=ylabel, title=plotconfig.plot_titles.get(plotname, plotname), xbounds=xbounds, ybounds=(-0.1, 1.1))

    plt.close()

# ----------------------------------------------------------------------------------------
def plot_gl_inference_fractions(plotdir, plotname, plotvals, labels, xlabel='', ylabel='', leg_title=None, title=None):
    if 'fraction_uncertainty' not in sys.modules:
        import fraction_uncertainty
    fraction_uncertainty = sys.modules['fraction_uncertainty']

    def get_single_vals(pv):
        yvals = [float(c) / t for c, t in zip(pv['ycounts'], pv['ytotals'])]  # total shouldn't be able to be zero
        tmphilos = [fraction_uncertainty.err(c, t) for c, t in zip(pv['ycounts'], pv['ytotals'])]
        yerrs = [err[1] - err[0] for err in tmphilos]
        print '  %s                    %s' % (xlabel, ylabel)
        for iv in range(len(pv['xvals'])):
            print '   %8.0f     %5.0f / %-5.0f  =  %5.2f   +/-  %.3f' % (pv['xvals'][iv], pv['ycounts'][iv], pv['ytotals'][iv], yvals[iv], yerrs[iv])
        return pv['xvals'], yvals, yerrs

    fig, ax = mpl_init()
    mpl.rcParams.update({'legend.fontsize' : 15})

    xmin, xmax, xticks = None, None, None
    for ii in range(len(labels)):
        print labels[ii]
        xvals, yvals, yerrs = get_single_vals(plotvals[ii])
        if xmin is None:
            xmin = xvals[0]
            xmax = xvals[-1]
            xticks = xvals
        kwargs = {
            'markersize' : 13,
            'linewidth' : default_linewidths[min(ii, len(default_linewidths) - 1)],
            'color' : default_colors[min(ii, len(default_colors) - 1)],
            'alpha' : 0.6,
        }
        ax.errorbar(xvals, yvals, yerr=yerrs, label=labels[ii], **kwargs)

    minfrac, maxfrac = 0.95, 1.05
    ax.plot((minfrac * xmin, maxfrac * xmax), (0, 0), color='black', linestyle='--', linewidth=3)  # line at y=0
    ax.plot((minfrac * xmin, maxfrac * xmax), (1, 1), color='black', linestyle='--', linewidth=3)  # line at y=1
    mpl_finish(ax, plotdir, plotname, xlabel=xlabel, ylabel=ylabel, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)
    plt.close()


# ----------------------------------------------------------------------------------------
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_selection_hists(histfname, plotdir, plotname, plot_all=False, n_plots=7, title='', xlabel=''):
    import joypy
    # ----------------------------------------------------------------------------------------
    def plot_this_time(otime, numpyhists):
        if plot_all:
            return True
        if otime == 0:
            return False
        if otime in (len(numpyhists),):
            return True
        if otime % int(len(numpyhists) / float(n_plots)) == 0:
            return True
        return False
    # ----------------------------------------------------------------------------------------
    def get_hists(hfname):
        with open(hfname) as runstatfile:
            numpyhists = pickle.load(runstatfile)
        xmin, xmax = None, None
        hists, labels = [], []
        for ihist in range(len(numpyhists)):
            nphist = numpyhists[ihist]  # numpy.hist is two arrays: [0] is bin counts, [1] is bin x values (not sure if low, high, or centers)
            obs_time = ihist  #  + 1  # I *think* it's right without the 1 (although I guess it's really a little arbitrary)
            if not plot_this_time(obs_time, numpyhists):
                continue
            if nphist is None:  # time points at which we didn't sample
                hists.append(None)
                labels.append(None)
                continue
            bin_contents, bin_edges = nphist
            assert len(bin_contents) == len(bin_edges) - 1
            # print ' '.join('%5.1f' % c for c in bin_contents)
            # print ' '.join('%5.1f' % c for c in bin_edges)
            hist = Hist(len(bin_edges) - 1, bin_edges[0], bin_edges[-1])
            for ibin in range(len(bin_edges) - 1):  # nphist indexing, not Hist indexing
                lo_edge = bin_edges[ibin]
                hi_edge = bin_edges[ibin + 1]
                bin_center = (hi_edge + lo_edge) / 2.
                for _ in range(bin_contents[ibin]):
                    hist.fill(bin_center)
                    xmin = lo_edge if xmin is None else min(xmin, lo_edge)
                    xmax = hi_edge if xmax is None else max(xmax, hi_edge)
            hists.append(hist)
            labels.append('%d (%.1f)' % (obs_time, hist.get_mean()))

        # hists = [Hist(1, xmin, xmax) if h is None else h for h in hists]  # replace the None values with empty hists
        hists, labels = zip(*[(h, l) for h, l in zip(hists, labels) if h is not None])
        return hists, labels, xmin, xmax

    # ----------------------------------------------------------------------------------------
    all_hists, all_labels, xmin, xmax = get_hists(histfname)
    jpdata = []
    for hist in all_hists:
        jpdata.append([x for x, y in zip(hist.get_bin_centers(), hist.bin_contents) for _ in range(int(y)) if x > xmin and x < xmax])  # NOTE this is repeating the 'for _ in range()' in the fcn above, but that's because I used to be actually using the Hist()s, and maybe I will again

    fig, ax = mpl_init()
    fig, axes = joypy.joyplot(jpdata, labels=all_labels, fade=True, hist=True, overlap=0.5, ax=ax, x_range=(xmin, xmax), bins=int(xmax - xmin))
    # NOTE do *not* set your own x ticks/labels in the next line, since they'll be in the wrong place (i.e. not the same as where joypy puts them)
    mpl_finish(ax, plotdir, plotname, title=title, xlabel=xlabel, ylabel='generation', leg_loc=(0.7, 0.45)) #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_kd_vals(plotdir, event):
    n_muts, kd_changes = [], []
    dtree = treeutils.get_dendro_tree(treestr=event['tree'])
    for node in dtree.preorder_internal_node_iter():
        for child in node.child_nodes():
            inode = event['unique_ids'].index(node.taxon.label)
            ichild = event['unique_ids'].index(child.taxon.label)
            node_affinity = event['affinities'][inode]
            child_affinity = event['affinities'][ichild]
            n_muts.append(utils.hamming_distance(event['input_seqs'][inode], event['input_seqs'][ichild]))
            kd_changes.append(1./child_affinity - 1./node_affinity)

    hist = Hist(30, min(kd_changes), max(kd_changes))
    for val in kd_changes:
        hist.fill(val)
    fig, ax = mpl_init()
    hist.mpl_plot(ax, square_bins=True, errors=False)  #remove_empty_bins=True)
    plotname = 'kd-changes'
    mpl_finish(ax, plotdir,  plotname, xlabel='parent-child kd change', ylabel='branches', log='y') #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

    plotvals = {'shm' : [], 'kd_vals' : []}
    for iseq, uid in enumerate(event['unique_ids']):
        plotvals['shm'].append(event['n_mutations'][iseq])
        plotvals['kd_vals'].append(1. / event['affinities'][iseq])
    # new_cmap = truncate_colormap(plt.cm.Blues, 0, 1)
    # ax.hexbin(kd_changes, shms, gridsize=25, cmap=plt.cm.Blues) #, info['ccf_under'][meth], label='clonal fraction', color='#cc0000', linewidth=4)
    fig, ax = mpl_init()
    ax.scatter(plotvals['kd_vals'], plotvals['shm'], alpha=0.4)
    plotname = 'kd-vs-shm'
    mpl_finish(ax, plotdir, plotname, xlabel='Kd', ylabel='N mutations') #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_target_attraction(plotdir, event):  # plots of which sequences are going toward which targets
    from Bio.Seq import Seq

    fig, ax = mpl_init()

    # affinity vs stuff:
    # xvals = [1. / af for line in mutated_events for af in line['affinities']]
    # yvals = [nm for line in mutated_events for nm in line['n_mutations']]

    # # min distance to target:
    # yvals = [hd for line in mutated_events for hd in get_min_target_hdists(line['input_seqs'], line['target_seqs'])]
    # ax.scatter(xvals, yvals, alpha=0.65)

    hist = Hist(len(event['target_seqs']), -0.5, len(event['target_seqs']) - 0.5, value_list=event['nearest_target_indices'])
    hist.mpl_plot(ax, alpha=0.7, ignore_overflows=True)

    plotname = 'nearest-target-identities'
    mpl_finish(ax, plotdir, plotname, xlabel='index (identity) of nearest target sequence', ylabel='counts') #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

# ----------------------------------------------------------------------------------------
def plot_bcr_phylo_simulation(outdir, event, extrastr):
    utils.prep_dir(outdir + '/plots', wildlings=['*.csv', '*.svg'])

    plot_bcr_phylo_kd_vals(outdir + '/plots', event)
    plot_bcr_phylo_target_attraction(outdir + '/plots', event)

    plot_bcr_phylo_selection_hists('%s/%s_min_aa_target_hdists.p' % (outdir, extrastr), outdir + '/plots', 'min-aa-target-all-cells', title='all cells', xlabel='AA distance to nearest target sequence')
    plot_bcr_phylo_selection_hists('%s/%s_sampled_min_aa_target_hdists.p' % (outdir, extrastr), outdir + '/plots', 'min-aa-target-sampled-cells', plot_all=True, title='sampled cells (excluding ancestor sampling)', xlabel='AA distance to nearest target sequence')
    plot_bcr_phylo_selection_hists('%s/%s_n_mutated_nuc_hdists.p' % (outdir, extrastr), outdir + '/plots', 'n-mutated-nuc-all-cells', title='SHM all cells', xlabel='N nucleotide mutations to naive')

    make_html(outdir + '/plots')

# ----------------------------------------------------------------------------------------
def get_tree_from_line(line, is_simu):
    if is_simu:
        return line['tree']
    if 'tree-info' not in line:  # if 'tree-info' is missing, it should be because it's a small cluster in data that we skipped when calculating lb values
        return None
    return line['tree-info']['lb']['tree']

# ----------------------------------------------------------------------------------------
def plot_lb_vs_shm(baseplotdir, lines_to_use, is_simu=False, n_per_row=4):  # <is_simu> is there because we want the true and inferred lines to keep their trees in different places, because the true line just has the one, true, tree, while the inferred line could have a number of them (yes, this means I maybe should have called it the 'true-tree' or something)
    sorted_lines = sorted([l for l in lines_to_use if get_tree_from_line(l, is_simu) is not None], key=lambda l: len(l['unique_ids']), reverse=True)
    fnames = [[]]

    # note: all clusters together
    subfnames = {lb_metric : [] for lb_metric in treeutils.lb_metrics}
    for lb_metric, lb_label in treeutils.lb_metrics.items():
        plotvals = {x : {'leaf' : [], 'internal' : []} for x in ['shm', lb_metric]}
        basetitle = '%s %s vs SHM' % ('true' if is_simu else 'inferred', lb_metric.upper())
        for iclust, line in enumerate(sorted_lines):  # get depth/n_mutations for each node
            iclust_plotvals = {x : {'leaf' : [], 'internal' : []} for x in ['shm', lb_metric, 'uids']}
            dtree = treeutils.get_dendro_tree(treestr=get_tree_from_line(line, is_simu))
            n_max_mutes = max(line['n_mutations'])  # don't generally have n mutations for internal nodes, so use this to rescale the depth in the tree
            max_depth = max(n.distance_from_root() for n in dtree.leaf_node_iter())
            for node in dtree.preorder_node_iter():
                if lb_metric == 'lbr' and line['tree-info']['lb'][lb_metric][node.taxon.label] == 0:  # lbr equals 0 should really be treated as None/missing
                    continue
                iseq = line['unique_ids'].index(node.taxon.label) if node.taxon.label in line['unique_ids'] else None
                n_muted = line['n_mutations'][iseq] if node.taxon.label in line['unique_ids'] else node.distance_from_root() * n_max_mutes / float(max_depth)
                tkey = 'leaf' if node.is_leaf() else 'internal'
                iclust_plotvals['shm'][tkey].append(n_muted)
                iclust_plotvals[lb_metric][tkey].append(line['tree-info']['lb'][lb_metric][node.taxon.label])
                affyval = line['affinities'][iseq] if 'affinities' in line and iseq is not None else None
                if not is_simu:
                    iclust_plotvals['uids'][tkey].append(node.taxon.label if affyval is not None else None)
            plotname = '%s-vs-shm-iclust-%d' % (lb_metric, iclust)
            title = '%s (%d observed, %d total)' % (basetitle, len(line['unique_ids']), len(line['tree-info']['lb'][lb_metric]))
            fn = plot_2d_scatter(plotname, '%s/%s-vs-shm' % (baseplotdir, lb_metric), iclust_plotvals, lb_metric, lb_label, title, xvar='shm', xlabel='N mutations', leg_loc=(0.7, 0.75), log='y' if lb_metric == 'lbr' else '')
            if iclust < n_per_row:  # i.e. only put one row's worth in the html
                subfnames[lb_metric].append(fn)
            for vtype in [vt for vt in plotvals if vt != 'uids']:
                for ltype in plotvals[vtype]:
                    plotvals[vtype][ltype] += iclust_plotvals[vtype][ltype]
        plotname = '%s-vs-shm' % lb_metric
        plot_2d_scatter(plotname, baseplotdir, plotvals, lb_metric, lb_label, '%s (all clusters)' % basetitle, xvar='shm', xlabel='N mutations', leg_loc=(0.7, 0.75), log='y' if lb_metric == 'lbr' else '')
        fnames[-1].append('%s/%s.svg' % (baseplotdir, plotname))
    fnames += [subfnames[lbm] for lbm in treeutils.lb_metrics]

    return fnames

# ----------------------------------------------------------------------------------------
def plot_lb_distributions(baseplotdir, lines_to_use, n_per_row=4):
    sorted_lines = sorted([l for l in lines_to_use if 'tree-info' in l], key=lambda l: len(l['unique_ids']), reverse=True)  # if 'tree-info' is missing, it should be because it's a small cluster we skipped when calculating lb values
    fnames = []

    for lb_metric, lb_label in treeutils.lb_metrics.items():
        plotdir = baseplotdir + '/' + lb_metric
        utils.prep_dir(plotdir, wildlings=['*.svg'])
        fnames.append([])
        for iclust, line in enumerate(sorted_lines):
            plotvals = line['tree-info']['lb'][lb_metric].values()
            leafskipstr = ''
            if lb_metric == 'lbr':
                plotvals = [v for v in plotvals if v > 0.]  # don't plot the leaf values, they just make the plot unreadable
                leafskipstr = ', skipped %d leaves' % len([v for v in line['tree-info']['lb'][lb_metric].values() if v == 0.])  # ok they're not necessarily leaves, but almost all of them are leaves (and not really sure how a non-leaf could get zero, but some of them seem to)
            hist = Hist(30, 0., max(plotvals), value_list=plotvals)
            fig, ax = mpl_init()
            hist.mpl_plot(ax) #, square_bins=True, errors=False)
            # ax.text(0.45 * ax.get_xlim()[1], 0.85 * ax.get_ylim()[1], 'size %d' % len(line['unique_ids']), fontsize=17, color='red', fontweight='bold')  # omfg this is impossible to get in the right place
            plotname = '%s-%d' % (lb_metric, iclust)
            mpl_finish(ax, plotdir, plotname, xlabel=lb_label, log='y' if lb_metric == 'lbr' else '', ylabel='counts', title='%s  (size %d%s)' % (lb_metric.upper(), len(line['tree-info']['lb'][lb_metric]), leafskipstr))
            if iclust < n_per_row:  # i.e. only put one row's worth in the html
                fnames[-1].append('%s/%s.svg' % (plotdir, plotname))
            make_html(plotdir)

    return fnames

# ----------------------------------------------------------------------------------------
def plot_2d_scatter(plotname, plotdir, plotvals, yvar, ylabel, title, xvar='affinity', xlabel='affinity', log='', leg_loc=None):
    def getall(k):
        if 'leaf' in plotvals[xvar]:
            return [v for tk in plotvals[k] for v in plotvals[k][tk]]
        else:
            return plotvals[k]

    if len(plotvals[xvar]) == 0:
        # print '    no %s vs affy info' % yvar
        return
    fig, ax = mpl_init()
    # cmap, norm = get_normalized_cmap_and_norm()
    # ax.hexbin(plotvals[xvar], plotvals[yvar], gridsize=15, cmap=plt.cm.Blues)
    if 'leaf' not in plotvals[xvar]:  # single plot
        ax.scatter(plotvals[xvar], plotvals[yvar], alpha=0.4)
    else:  # separate plots for leaf/internal nodes
        for tkey, color in zip(plotvals[xvar], (None, 'darkgreen')):
            ax.scatter(plotvals[xvar][tkey], plotvals[yvar][tkey], label=tkey, alpha=0.4, color=color)
    if 'uids' in plotvals:
        for xval, yval, uid in zip(getall(xvar), getall(yvar), getall('uids')):  # note: two ways to signal not to do this: sometimes we have 'uids' in the dict, but don't fill it (so the zip() gives an empty list), but sometimes we populate 'uids' with None values
            if uid is None:
                continue
            ax.plot([xval], [yval], color='red', marker='.', markersize=10)
            ax.text(xval, yval, uid, color='red', fontsize=8)

    xmin, xmax = min(getall(xvar)), max(getall(xvar))
    ymin, ymax = min(getall(yvar)), max(getall(yvar))
    xbounds = xmin - 0.01 * (xmax - xmin), 1.05 * xmax
    if 'y' in log:
        ybounds = 0.95 * ymin, 1.05 * ymax
    else:
        ybounds = ymin - 0.01 * (ymax - ymin), 1.05 * ymax
    mpl_finish(ax, plotdir, plotname, title=title, xlabel=xlabel, ylabel=ylabel, xbounds=xbounds, ybounds=ybounds, log=log, leg_loc=leg_loc)
    return '%s/%s.svg' % (plotdir, plotname)

# ----------------------------------------------------------------------------------------
def plot_lb_vs_affinity(plot_str, plotdir, lines, lb_metric, lb_label, all_clusters_together=False, ptile_range_tuple=(50., 100., 1.), is_simu=False, n_per_row=4, debug=False):
    fnames = []

    # first plot lb metric vs affinity scatter (all clusters)
    lb_vs_affinity_vals = {val_type : [] for val_type in [lb_metric, 'affinity']}  # , 'uids']}
    for iclust, line in enumerate(lines):
        iclust_lb_vs_affinity_vals = {val_type : [] for val_type in [lb_metric, 'affinity' , 'uids']}
        # dtree = treeutils.get_dendro_tree(treestr=get_tree_from_line(line, is_simu))
        if 'affinities' not in line:
            continue
        for uid, affy in [(u, a) for u, a in zip(line['unique_ids'], line['affinities']) if a is not None]:
            # node = dtree.find_node_with_taxon_label(uid)
            iclust_lb_vs_affinity_vals['affinity'].append(affy)
            iclust_lb_vs_affinity_vals[lb_metric].append(line['tree-info']['lb'][lb_metric][uid])
            if not is_simu:
                iclust_lb_vs_affinity_vals['uids'].append(uid)
        if not all_clusters_together and len(iclust_lb_vs_affinity_vals['affinity']) > 0:
            fn = plot_2d_scatter('iclust-%d' % iclust, '%s/%s-vs-affinity' % (plotdir, lb_metric), iclust_lb_vs_affinity_vals, lb_metric, lb_label, '%s (%s tree)' % (lb_metric.upper(), plot_str))
            if iclust < n_per_row:
                fnames.append(fn)
        for vtype in [vt for vt in lb_vs_affinity_vals if vt != 'uids']:
            lb_vs_affinity_vals[vtype] += iclust_lb_vs_affinity_vals[vtype]
    if all_clusters_together:
        plotname = '%s-vs-affinity-%s-tree' % (lb_metric, plot_str)
        plot_2d_scatter(plotname, plotdir, lb_vs_affinity_vals, lb_metric, lb_label, '%s (%s tree)' % (lb_metric.upper(), plot_str))
        fnames.append('%s/%s.svg' % (plotdir, plotname))

    if len(lb_vs_affinity_vals[lb_metric]) == 0:
        print '  no affinity values when trying to make lb vs affinity plots'
        return [fnames]

    # then plot potential lb cut thresholds with percentiles
    if debug:
        print '    ptile   %s     mean affy    mean affy ptile' % lb_metric
    ptile_vals = {'lb_ptiles' : [], 'mean_affy_ptiles' : [], 'perfect_vals' : []}  # , 'reshuffled_vals' : []}
    lbvals = lb_vs_affinity_vals[lb_metric]  # should really use these shorthands for the previous plot as well
    affyvals = lb_vs_affinity_vals['affinity']
    sorted_affyvals = sorted(affyvals, reverse=True)
    for percentile in numpy.arange(*ptile_range_tuple):
        lb_ptile_val = numpy.percentile(lbvals, percentile)  # lb value corresponding to <percentile>
        corresponding_affinities = [affy for lb, affy in zip(lbvals, affyvals) if lb > lb_ptile_val]  # affinities corresponding to lb greater than <lb_ptile_val> (i.e. the affinities that you'd get if you took all the lb values greater than that)
        corr_affy_ptiles = [stats.percentileofscore(affyvals, caffy) for caffy in corresponding_affinities]  # affinity percentiles corresponding to each of these affinities  # NOTE this is probably really slow
        if len(corr_affy_ptiles) == 0:
            if debug:
                print '   %5.0f    no vals' % percentile
            continue
        ptile_vals['lb_ptiles'].append(percentile)
        ptile_vals['mean_affy_ptiles'].append(numpy.mean(corr_affy_ptiles))
        if debug:
            print '   %5.0f   %5.2f   %8.4f     %5.0f' % (percentile, lb_ptile_val, numpy.mean(corresponding_affinities), ptile_vals['mean_affy_ptiles'][-1])

        # make a "perfect" line from actual affinities, as opposed to just a straight line (this accounts better for, e.g. the case where the top N affinities are all the same)
        n_to_take = int((1. - percentile / 100) * len(sorted_affyvals))
        corresponding_perfect_affy_vals = sorted_affyvals[:n_to_take]
        corr_perfect_affy_ptiles = [stats.percentileofscore(affyvals, cpaffy) for cpaffy in corresponding_perfect_affy_vals]  # NOTE this is probably really slow
        ptile_vals['perfect_vals'].append(numpy.mean(corr_perfect_affy_ptiles))

        # # add a horizontal line at 50 to show what it'd look like if there was no correlation (this is really wasteful... although it does have a satisfying wiggle to it. Now using a plain flat line [below])
        # shuffled_lb_vals = copy.deepcopy(lbvals)
        # random.shuffle(shuffled_lb_vals)
        # NON_corresponding_affinities = [affy for lb, affy in zip(shuffled_lb_vals, affyvals) if lb > lb_ptile_val]
        # NON_corr_affy_ptiles = [stats.percentileofscore(affyvals, caffy) for caffy in NON_corresponding_affinities]
        # ptile_vals['reshuffled_vals'].append(numpy.mean(NON_corr_affy_ptiles))

    fig, ax = mpl_init()
    ax.plot(ptile_vals['lb_ptiles'], ptile_vals['mean_affy_ptiles'], linewidth=3, alpha=0.7)
    # ax.plot(ax.get_xlim(), [50 + 0.5 * x for x in ax.get_xlim()], linewidth=3, alpha=0.7, color='darkgreen', linestyle='--', label='perfect correlation')  # straight line
    ax.plot(ptile_vals['lb_ptiles'], ptile_vals['perfect_vals'], linewidth=3, alpha=0.7, color='darkgreen', linestyle='--', label='perfect correlation')  # perfect vals
    ax.plot(ax.get_xlim(), (50, 50), linewidth=3, alpha=0.7, color='darkred', linestyle='--', label='no correlation')  # straight line
    # ax.plot(ptile_vals['lb_ptiles'], ptile_vals['reshuffled_vals'], linewidth=3, alpha=0.7, color='darkred', linestyle='--', label='no correlation')  # reshuffled vals
    plotname = '%s-vs-affinity-%s-tree-ptiles' % (lb_metric, plot_str)
    mpl_finish(ax, plotdir, plotname, xbounds=(ptile_range_tuple[0], ptile_range_tuple[1]), ybounds=(45, 100), leg_loc=(0.035, 0.75), title='potential %s thresholds (%s tree)' % (lb_metric.upper(), plot_str), xlabel='%s threshold (percentile)' % lb_metric.upper(), ylabel='mean percentile of resulting affinities')
    fnames.append('%s/%s.svg' % (plotdir, plotname))

    return [fnames]

# ----------------------------------------------------------------------------------------
def plot_lb_vs_delta_affinity(plotdir, true_lines, lb_metric, lb_label, debug=False):
    fnames = []

    delta_affinity_vals = {val_type : [] for val_type in [lb_metric, 'delta-affinity']}
    for line in true_lines:
        dtree = treeutils.get_dendro_tree(treestr=line['tree'])
        for uid, affinity in zip(line['unique_ids'], line['affinities']):
            node = dtree.find_node_with_taxon_label(uid)
            if node is dtree.seed_node:  # root won't have a parent
                continue
            parent_uid = node.parent_node.taxon.label
            if parent_uid not in line['unique_ids']:
                print '    %s parent %s of %s not in true line' % (utils.color('yellow', 'warning'), parent_uid, uid)
                continue
            iparent = line['unique_ids'].index(parent_uid)
            parent_affinity = line['affinities'][iparent]
            delta_affinity_vals['delta-affinity'].append(affinity - parent_affinity)
            delta_affinity_vals[lb_metric].append(line['tree-info']['lb'][lb_metric][uid])
    fig, ax = mpl_init()
    ax.scatter(delta_affinity_vals['delta-affinity'], delta_affinity_vals[lb_metric], alpha=0.4)
    sorted_xvals = sorted(delta_affinity_vals['delta-affinity'])  # not sure why, but ax.scatter() is screwing up the x bounds
    xmin, xmax = sorted_xvals[0], sorted_xvals[-1]
    # ax.hexbin(delta_affinity_vals['delta-affinity'], delta_affinity_vals[lb_metric], gridsize=15, cmap=plt.cm.Blues)
    plotname = '%s-vs-delta-affinity' % lb_metric
    mpl_finish(ax, plotdir, plotname, title='%s (true tree)' % lb_metric.upper(), xlabel='affinity change (from parent)', ylabel=lb_label, xbounds=(1.05 * xmin, 1.05 * xmax))  # NOTE factor on <xmin> is only right if xmin is negative, but it should always be
    fnames.append('%s/%s.svg' % (plotdir, plotname))

    return [fnames]

# ----------------------------------------------------------------------------------------
def plot_lb_vs_ancestral_delta_affinity(plotdir, true_lines, lb_metric, lb_label, plot_str='true', ptile_range_tuple=(50., 100., 1.), debug=False):
    # plot lb[ir] vs number of ancestors to nearest affinity decrease (well, decrease as you move upwards in the tree/backwards in time)
    # NOTE because it's so common for affinity to get worse from ancestor to descendent, it's important to remember that here we are looking for the first ancestor with lower affinity than the node in question, which is *different* to looking for the first ancestor that has lower affinity than one of its immediate descendents (which we could also plot, but it probably wouldn't be significantly different to the metric performance, since for the metric performance we only really care about the left side of the plot, but this only affects the right side)
    fnames = []

    # first plot lb metric vs number of ancestors since last affinity increase (all clusters)
    if debug:
        print '  finding N ancestors since last affinity increase'
        print '         node        ancestors   affinity (%sX: change for chosen ancestor, %s: reached root without finding lower-affinity ancestor)' % (utils.color('red', '+'), utils.color('green', 'x'))
    n_ancestor_vals = {val_type : [] for val_type in [lb_metric, 'n-ancestors']}  # , 'uids']}
    # TMP_plotvals = {x : {'leaf' : [], 'internal' : []} for x in [lb_metric, 'n-ancestors']}  # TODO all this commented stuff is for splitting apart leaf and internal plots. If I end up wanting to do that again, I should combine it with plot_lb_vs_shm().
    for line in true_lines:
        dtree = treeutils.get_dendro_tree(treestr=line['tree'])
        affinity_changes = []
        for this_uid, this_affinity in zip(line['unique_ids'], line['affinities']):
            node = dtree.find_node_with_taxon_label(this_uid)
            if node is dtree.seed_node:  # root doesn't have any ancestors
                continue
            if lb_metric == 'lbr' and line['tree-info']['lb'][lb_metric][this_uid] == 0:  # lbr equals 0 should really be treated as None/missing
                continue

            if debug:
                print '     %12s %12s %9.4f' % (this_uid, '', this_affinity)
            min_affinity_change = 1e-6  # just to eliminate floating point precision issues (especially since we're deriving affinity by inverting kd) (note that at least for now, the affinity changes should all be pretty similar, and not small)
            n_max_steps = 15
            ancestor_node = node
            ancestor_affinity = None
            n_steps = 0
            while True:  # NOTE if we can't find an ancestor with worse affinity, we don't plot the node
                if ancestor_node is dtree.seed_node:
                    break
                ancestor_node = ancestor_node.parent_node  #  move one more step up the tree
                ancestor_uid = ancestor_node.taxon.label
                if ancestor_uid not in line['unique_ids']:
                    print '    %s ancestor %s of %s not in true line' % (utils.color('yellow', 'warning'), ancestor_uid, this_uid)
                    break
                iancestor = line['unique_ids'].index(ancestor_uid)
                if this_affinity - line['affinities'][iancestor] > min_affinity_change:  # if we found an ancestor with lower affinity, we're done
                    ancestor_affinity = line['affinities'][iancestor]
                    affinity_changes.append(this_affinity - line['affinities'][iancestor])
                    break
                if debug:
                    print '     %12s %12s %9.4f%s' % ('', ancestor_uid, line['affinities'][iancestor], utils.color('green', ' x') if ancestor_node is dtree.seed_node else '')
                n_steps += 1
                if n_steps >= n_max_steps:
                    break
            if ancestor_affinity is None:
                # print '    couldn\'t find ancestor with lower affinity for %s within %d steps' % (this_uid, n_max_steps)
                continue

            if debug:
                print '     %12s %12s %9.4f  %s%-9.4f' % ('', ancestor_uid, ancestor_affinity, utils.color('red', '+'), this_affinity - ancestor_affinity)

            n_ancestor_vals['n-ancestors'].append(n_steps)
            n_ancestor_vals[lb_metric].append(line['tree-info']['lb'][lb_metric][this_uid])
            # n_ancestor_vals['uids'].append(this_uid)
            # tkey = 'leaf' if node.is_leaf() else 'internal'
            # TMP_plotvals['n-ancestors'][tkey].append(n_steps)
            # TMP_plotvals[lb_metric][tkey].append(line['tree-info']['lb'][lb_metric][this_uid])

        # make sure affinity changes are all roughly the same size
        affinity_changes = sorted(affinity_changes)
        if debug:
            print '    chosen affinity changes: %s' % ' '.join(['%.4f' % a for a in affinity_changes])
        if len([a for a in affinity_changes if a < 0.]):
            print '  %s negative affinity changes in %s' % (utils.color('red', 'error'), ' '.join(['%.4f' % a for a in affinity_changes]))
        max_diff = affinity_changes[-1] - affinity_changes[0]
        if abs(max_diff) / numpy.mean(affinity_changes) > 0.2:
            print'  %s not all affinity increases were the same size (min: %.4f   max: %.4f   abs(diff) / mean: %.4f' % (utils.color('yellow', 'warning'), affinity_changes[0], affinity_changes[-1], abs(max_diff) / numpy.mean(affinity_changes))

    plotname = '%s-vs-n-ancestors-%s-tree' % (lb_metric, plot_str)  # 'nearest ancestor with lower affinity' would in some ways be a better xlabel, since it clarifies the note at the top of the loop, but it's also less clear in other ways
    # fig, ax = mpl_init()
    # for tkey, color in zip(TMP_plotvals['n-ancestors'], (None, 'darkgreen')):
    #     ax.scatter(TMP_plotvals['n-ancestors'][tkey], TMP_plotvals[lb_metric][tkey], label=tkey, alpha=0.4, color=color)
    # mpl_finish(ax, plotdir, plotname, title='%s on true tree' % lb_metric.upper(), xlabel='N ancestors since affinity increase', ylabel=lb_label)
    plot_2d_scatter(plotname, plotdir, n_ancestor_vals, lb_metric, lb_label, '%s (true tree)' % lb_metric.upper(), xvar='n-ancestors', xlabel='N ancestors since affinity increase', log='y' if lb_metric == 'lbr' else '')
    fnames.append('%s/%s.svg' % (plotdir, plotname))

    # then plot potential lb cut thresholds with percentiles
    if debug:
        print '    ptile   %s     mean N ancestors' % lb_metric
    lb_ptile_vals = {'lb_ptiles' : [], 'mean_n_ancestors' : []}
    lbvals = n_ancestor_vals[lb_metric]
    n_anc_vals = n_ancestor_vals['n-ancestors']
    for percentile in numpy.arange(*ptile_range_tuple):
        lb_ptile_val = numpy.percentile(lbvals, percentile)  # lb value corresponding to <percentile>
        corresponding_n_anc_vals = [n_anc for lb, n_anc in zip(lbvals, n_anc_vals) if lb >= lb_ptile_val]  # n-anc vals corresponding to lb greater than <lb_ptile_val> (i.e. the n-anc vals that you'd get if you took all the lb values greater than that)
        if len(corresponding_n_anc_vals) == 0:
            if debug:
                print '   %5.0f    no vals' % percentile
            continue
        lb_ptile_vals['lb_ptiles'].append(percentile)
        lb_ptile_vals['mean_n_ancestors'].append(numpy.mean(corresponding_n_anc_vals))
        if debug:
            print '   %5.0f   %5.2f   %8.4f' % (percentile, lb_ptile_val, lb_ptile_vals['mean_n_ancestors'][-1])

    # then make the "perfect" line
    if debug:
        print '    ptile  n_taken   %s     mean N ancestors' % 'affy'
    perfect_ptile_vals = {'ptiles' : [], 'mean_n_ancestors' : []}
    sorted_n_anc_vals = sorted(n_ancestor_vals['n-ancestors'])
    for percentile in numpy.arange(*ptile_range_tuple):
        n_to_take = int((1. - percentile / 100) * len(sorted_n_anc_vals))
        corresponding_n_anc_vals = sorted_n_anc_vals[:n_to_take]
        if len(corresponding_n_anc_vals) == 0:
            if debug:
                print '   %5.0f   %5d    no vals' % (percentile, n_to_take)
            continue
        perfect_ptile_vals['ptiles'].append(percentile)
        perfect_ptile_vals['mean_n_ancestors'].append(numpy.mean(corresponding_n_anc_vals))
        if debug:
            print '   %5.0f   %5d   %8.4f' % (percentile, n_to_take, perfect_ptile_vals['mean_n_ancestors'][-1])

    fig, ax = mpl_init()
    ax.plot(lb_ptile_vals['lb_ptiles'], lb_ptile_vals['mean_n_ancestors'], linewidth=3, alpha=0.7)
    mean_n_anc = numpy.mean(n_anc_vals)
    ax.plot(ax.get_xlim(), (mean_n_anc, mean_n_anc), linewidth=3, alpha=0.7, color='darkred', linestyle='--', label='no correlation')
    ax.plot(perfect_ptile_vals['ptiles'], perfect_ptile_vals['mean_n_ancestors'], linewidth=3, alpha=0.7, color='darkgreen', linestyle='--', label='perfect correlation')
    plotname = '%s-vs-n-ancestors-%s-tree-ptiles' % (lb_metric, plot_str)
    ymax = max([mean_n_anc] + lb_ptile_vals['mean_n_ancestors'] + perfect_ptile_vals['mean_n_ancestors'])
    mpl_finish(ax, plotdir, plotname, xbounds=(ptile_range_tuple[0], ptile_range_tuple[1]), ybounds=(0, 1.1 * ymax), leg_loc=(0.035, 0.05), title='potential %s thresholds (%s tree)' % (lb_metric.upper(), plot_str), xlabel='%s threshold (percentile)' % lb_metric.upper(), ylabel='mean N ancestors since affinity increase')
    fnames.append('%s/%s.svg' % (plotdir, plotname))

    return [fnames]

# ----------------------------------------------------------------------------------------
def plot_true_vs_inferred_lb(plotdir, true_lines, inf_lines, lb_metric, lb_label, debug=False):
    plotvals = {val_type : {uid : l['tree-info']['lb'][lb_metric][uid] for l in lines for uid in l['unique_ids']}
                for val_type, lines in (('true', true_lines), ('inf', inf_lines))}
    common_uids = set(plotvals['true']) & set(plotvals['inf'])  # there should/may be a bunch of internal nodes in the simulation lines but not in the inferred lines, but otherwise they should have the same uids
    plotvals = {val_type : [plotvals[val_type][uid] for uid in common_uids] for val_type in plotvals}
    plotname = '%s-true-vs-inferred' % lb_metric
    plot_2d_scatter(plotname, plotdir, plotvals, 'inf', '%s on inferred tree' % lb_metric.upper(), 'true vs inferred %s' % lb_metric.upper(), xvar='true', xlabel='%s on true tree' % lb_metric.upper())
    return ['%s/%s.svg' % (plotdir, plotname)]

# ----------------------------------------------------------------------------------------
def plot_per_mutation_lonr(plotdir, lines_to_use, reco_info):
    fig, ax = mpl_init()

    plotvals = {'lonr' : [], 'affinity_change' : []}
    for line in lines_to_use:
        true_affinities = {uid : reco_info[uid]['affinities'][0] for uid in line['unique_ids']}
        nodefos = line['tree-info']['lonr']['nodes']
        for lfo in line['tree-info']['lonr']['values']:
            if lfo['parent'] not in true_affinities:
                print '    %s parent \'%s\' not in true affinities, skipping lonr values' % (utils.color('red', 'warning'), lfo['parent'])
                continue
            if lfo['child'] not in true_affinities:
                print '    %s child \'%s\' not in true affinities, skipping lonr values' % (utils.color('red', 'warning'), lfo['child'])
                continue

            plotvals['lonr'].append(lfo['lonr'])
            plotvals['affinity_change'].append(true_affinities[lfo['child']] - true_affinities[lfo['parent']])

    ax.scatter(plotvals['affinity_change'], plotvals['lonr'], alpha=0.7) #, info['ccf_under'][meth], label='clonal fraction', color='#cc0000', linewidth=4)
    plotname = 'lonr-per-mut-vs-affinity'
    mpl_finish(ax, plotdir, plotname, xlabel='change in affinity', ylabel='LONR') #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

# ----------------------------------------------------------------------------------------
def plot_aggregate_lonr(plotdir, lines_to_use, reco_info, debug=False):
    fig, ax = mpl_init()
    plotvals = {'S' : [], 'NS' : []}
    for line in lines_to_use:
        for lfo in line['tree-info']['lonr']['values']:
            if lfo['synonymous']:
                plotvals['S'].append(lfo['lonr'])
            else:
                plotvals['NS'].append(lfo['lonr'])
    # ax.plot(plotvals['S'], label='S', linewidth=3, alpha=0.7)
    # ax.plot(plotvals['NS'], label='NS', linewidth=3, alpha=0.7)
    xmin, xmax = [mfcn([x for mtlist in plotvals.values() for x in mtlist]) for mfcn in (min, max)]
    hists = {mt : Hist(30, xmin, xmax, value_list=plotvals[mt], title=mt, xtitle='LONR', ytitle='mutations') for mt in plotvals}
    plotname = 'lonr-ns-vs-s'

    lonr_score = hists['NS'].get_mean() - hists['S'].get_mean()
    draw_no_root(hists['NS'], more_hists=[hists['S']], plotname=plotname, plotdir=plotdir, alphas=[0.7, 0.7], plottitle='NS - S: %.2f' % lonr_score, errors=True, remove_empty_bins=True)

    # for mt, hist in hists.items():
    #     hist.mpl_plot(ax, label=mt, remove_empty_bins=True)
    # mpl_finish(ax, plotdir, plotname, xlabel='LONR', ylabel='mutations') #, xbounds=(minfrac*xmin, maxfrac*xmax), ybounds=(-0.05, 1.05), log='x', xticks=xticks, xticklabels=[('%d' % x) for x in xticks], leg_loc=(0.8, 0.55 + 0.05*(4 - len(plotvals))), leg_title=leg_title, title=title)

# ----------------------------------------------------------------------------------------
def plot_laplacian_spectra(plotdir, plotname, eigenvalues, title):
    hist = Hist(30, min(eigenvalues), max(eigenvalues), value_list=eigenvalues)
    fig, ax = mpl_init()
    hist.mpl_plot(ax)
    mpl_finish
    mpl_finish(ax, plotdir, plotname, xlabel='eigenvalues', ylabel='count', title=title)
