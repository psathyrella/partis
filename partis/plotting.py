from __future__ import unicode_literals

from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import copy
import matplotlib as mpl
from io import open
mpl.use('Agg')
mpl.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d
import os
import glob
import sys
import csv
import numpy
import operator
import itertools
import collections
import six

from . import utils
from . import plotconfig
from .hist import Hist
from . import treeutils
from . import hutils
from .clusterpath import ClusterPath

#                   green    dark red  light blue  light red  sky blue  pink/purple   grey
default_colors = ['#006600', '#990012', '#2b65ec', '#cc0000', '#3399ff', '#a821c7', '#808080']
default_linewidths = ['5', '3', '2', '2', '2']
default_markersizes = ['20', '15', '8', '5', '5', '5']
pltcolors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # pyplot/matplotlib default colors
#                     blue      orange     green
frozen_pltcolors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']  # default colors from version 2.2.4 (so we don't get different colors on different machines/installs)
legend_groups = {  # replace <key> by <value> when making legend (i.e. group several together). Colors shoud be the same if they share a legend entry, otherwise it'll be a random one of the colors
    None : 'None',
    'ASC_1-gs' : 'GS ASC',
    'ASC_2-gs' : 'GS ASC',
    'ASC_3-gs' : 'GS ASC',
    'ASC_4-gs' : 'GS ASC',
    'ASC_5-gs' : 'GS ASC',
    'ASC_T-gs' : 'GS ASC',
    'ASC/T-gs' : 'GS ASC',
    'B_1-gs' : 'GS B',
    'B_2-gs' : 'GS B',
    'B_3-gs' : 'GS B',
    'B_4-gs' : 'GS B',
    'Bcell_2-gs' : 'GS B',
    'B/T-gs' : 'GS B',
    'None-gs' : 'GS n/a',
    '0-gs' : 'GS n/a',
    '8wph-gs' : 'containment',
    '16wph-gs' : 'containment',
    'entry-gs' : 'containment',
    'ctmt-gs' : 'containment',
    'lesion-gs' : 'lesion',
    'None-pbmc' : 'PBMC',
    'asc-pbmc' : 'PBMC',
    'mbc-pbmc' : 'PBMC',
    'naive-pbmc' : 'PBMC',
    'total-pbmc' : 'PBMC',
}

legend_labels = {
    'pbmc' : 'PBMC',
    'gs' : 'GS',
    'timepoint-tissues' : 'timepoint',
    'cell-types' : 'cell type',
    'cell-type-tissues' : 'cell type',
    'None' : 'inf. ancest.',
}

asc_col = '#d62728'
na_col = '#006600'
pbmc_col = '#2b65ec'
ctmt_col = na_col

hard_meta_colors = {
    # # leslie goo:
    # 'IGHM' : '#9467bd',  # purple
    # 'IGHD' : 'black',
    # 'IGHG1' : '#1f77b4', 'IGHG2' : '#6eb7e8', 'IGHG3' : '#6088a2', 'IGHG4' : '#1c47bb',  # shades of blue
    # 'IGHA1' : '#d62728', 'IGHA2' : '#ea7979',  # shades of red
    # 'memory' : '#1f77b4',  # blue
    # 'naive' : '#ff7f0e',  # orange
    # 'pb' : '#2ca02c',  # green
    # 'prepb' : '#90e690',
    # 'w2' : '#1f77b4',  # blue
    # 'w23' : '#2ca02c',  # green
    # 'w25' : '#ff7f0e',  # orange
    # 'd385' : '#2b65ec',
    # 'd765' : '#ff7f0e',  # orange
    # 'dm351' : '#006600',  # green
    # 'dm379' : '#9467bd',
    # anton
    # TODO update subjects
    # 'al' : '#1f77b4',
    # 'aw' : '#006600',
    # 'cf' : '#ff7f0e',
    # 'kh' : '#be4f48', #'#cc0000',
    # 'lm' : '#9467bd',
    # 'gs' : '#1f77b4',  # blue
    # # 'pbmc' : '#006600',  # green
    # # 'lesion' : '#be4f48',  # red
    'pbmc' : '#1f77b4',
    'PBMC' : '#2b65ec',
    'gs' : '#ff7f0e',
    'entry' : '#006600',  # green
    'None-gs' : na_col,
    'GS n/a' : na_col,
    '0-gs' : na_col,
    'GS ASC' : asc_col,
    'ASC_1-gs' : asc_col,
    'ASC_2-gs' : asc_col,
    'ASC_3-gs' : asc_col,
    'ASC_4-gs' : asc_col,
    'ASC_5-gs' : asc_col,
    'ASC_T-gs' : asc_col,
    'ASC/T-gs' : asc_col,
    'GS B' : '#ff7f0e',
    'B_1-gs' : '#ff7f0e',
    'B_2-gs' : '#ff7f0e',
    'Bcell_2-gs' : '#ff7f0e',
    'B_3-gs' : '#ff7f0e',
    'B_4-gs' : '#ff7f0e',
    'B/T-gs' : '#ff7f0e',
    'PBMC-gs' : pbmc_col,
    'asc-pbmc' : pbmc_col,
    'mbc-pbmc' : pbmc_col,
    'naive-pbmc' : pbmc_col,
    'total-pbmc' : pbmc_col,
    'None-pbmc' : pbmc_col,
    '8wph-gs' : ctmt_col,
    '16wph-gs' : ctmt_col,
    'entry-gs' : ctmt_col,
    'ctmt-gs' : ctmt_col,
    '8wph-gs' : ctmt_col,
    'lesion-gs' : '#ff7f0e',
}

# check for inconsistencies in the legend groups:
tlist = [(legend_groups.get(l, l), l, c) for l, c in hard_meta_colors.items()]
def kfn(p): return p[0]
for lgrp, y in utils.group_seqs_by_value(tlist, kfn, return_values=True):
    # print('  %20s     %s' % (lgrp, y))
    lgroups, labels, colors = zip(*y)
    if len(set(colors)) > 1:
        raise Exception('multiple colors specified for legend group %s:\n    %s\n    %s' % (lgrp, labels, colors))

# ----------------------------------------------------------------------------------------
def get_cluster_size_xticks(xmin=None, xmax=None, hlist=None):  # pass in either xmin and xmax, or hlist NOTE pretty similar to get_auto_y_ticks() (for search: log_bins log bins)
    if xmin is None or xmax is None:
        assert xmin is None and xmax is None  # would have to implement it if you want to be able to set just one
        assert hlist is not None
        minlist, maxlist = zip(*[h.get_filled_bin_xbounds() for h in hlist])
        xmin, xmax = [mfcn(mlist) for mfcn, mlist in zip((min, max), (minlist, maxlist))]
    def tstr(xt): return ('%.0f'%xt) if xt < 500 else '%.0e'%xt
    default_xticks = [1, 2, 3, 10, 30, 75, 200, 1000, 5000, 10000]
    xticks = [xt for xt in default_xticks if xt >= xmin and xt <= xmax]
    if len(xticks) < 3:
        xticks = [int(xmin) + 1, int((xmin + xmax)/2.), int(xmax)]
    if xmax > 2*xticks[-1]:  # just big enough that they don't overlap
        xticks.append(xmax)
    # this was another way of getting x ticks, if you still have the partition, and it's kind of nice:
    # csizes = sorted([len(c) for c in partition])
    # xticks = [x for x in numpy.logspace(math.log(csizes[0], 10), math.log(csizes[-1], 10), num=5)]
    return xticks, [tstr(xt) for xt in xticks]

# ----------------------------------------------------------------------------------------
def make_csize_hist(partition, n_bins=10, xbins=None, xtitle=None, weight_by_n_seqs=False, debug=False):  # would be nice to combine with plot_cluster_size_hists() but I'm not sure it really makes sense
    if len(partition) == 0:
        return Hist(n_bins=n_bins, xmin=0, xmax=100)  # ugh, have to set some values
    cslist = [len(c) for c in partition]
    if xbins is None:
        xbins, n_bins = hutils.auto_volume_bins(cslist, n_bins, int_bins=True, debug=debug)
    else:
        if debug:
            print('   using user specified bins %s:' % xbins)
            hutils.binprint(xbins, cslist)
        if any(x==int(x) for x in xbins):
            raise Exception('bin boundaries should be non-integers (e.g. 3.5) to avoid integer values seeming to fall on the boundary')
        if max(cslist) >= xbins[-1]:  # last entry in xbins is the low edge of overflow bin, and we *really* want to know if anything's overflowing (well maybe in future i'll want to allow this, but not now)
            raise Exception('cluster size %d would fall in overflow bin (see previous lines)' % max(cslist))
        if min(cslist) < xbins[0]:  # same for underflow
            raise Exception('cluster size %d would fall in underflow bin (see previous lines)' % min(cslist))
        n_bins = len(xbins) - 1
    thist = Hist(n_bins=n_bins, xmin=xbins[0], xmax=xbins[-1], xbins=xbins, value_list=cslist, weight_list=cslist if weight_by_n_seqs else None, sumw2=weight_by_n_seqs)
    for ib in range(1, thist.n_bins + 1):  # set bin labels
        lo, hi = thist.low_edges[ib], thist.low_edges[ib+1]
        ivals = [i for i in range(int(math.ceil(lo)), int(math.floor(hi)) + 1)]  # integers that fall in this bin
        if len(ivals) > 0:  # if there aren't any values in this bin, the bins really suck, but I just want it to not crash here
            thist.bin_labels[ib] = str(ivals[0]) if len(ivals)==1 else '%d-%d'%(ivals[0], ivals[-1])  # set bin label to min/max of integers in this bin
    if xtitle is not None:
        thist.xtitle = xtitle
    return thist

# ----------------------------------------------------------------------------------------
plot_ratios = {
    'v' : (30, 3),
    'd' : (8, 4),
    'j' : (8, 3)
}

# ----------------------------------------------------------------------------------------
# specify either all_emph_vals or clusters and antn_dict
def meta_emph_init(meta_info_key_to_color, clusters=None, antn_dict=None, all_emph_vals=None, formats=None, lgroups=None):
    # tme_colors = alt_colors + [c for c in frozen_pltcolors if c not in alt_colors]
    def get_tcols(): return [c for c in frozen_pltcolors if c not in ['#d62728', '#7f7f7f']]  # can't use red or grey
    tme_colors = get_tcols()
    if all_emph_vals is None:
        all_emph_vals = set(utils.meta_emph_str(meta_info_key_to_color, v, formats=formats, legend_groups=lgroups) for c in clusters for v in antn_dict.get(':'.join(c), {}).get(meta_info_key_to_color, [None for _ in c]))  # set of all possible values that this meta info key takes on in any cluster
    else:  # NOTE all_emph_vals needs to be a set if you pass it in
        assert clusters is None and antn_dict is None
    def cfcn(i, v): return 'grey' if v in [None, 'None'] else hard_meta_colors.get(v, tme_colors[i%len(tme_colors)])
    emph_colors = []
    for iv, val in enumerate(sorted(all_emph_vals - set([None, 'None']), key=str)):
        tcol = cfcn(iv, val)
        emph_colors.append((val, cfcn(iv, val)))
        if tcol in tme_colors:
            tme_colors.remove(tcol)
        if len(tme_colors) == 0:
            print('  %s ran out of colors, re-adding all of them' % utils.wrnstr())
            tme_colors = get_tcols()
    emph_colors += [('None', 'grey'), (None, 'grey')]  # want to make sure None is last, so it's at the bottom of the legend
    return all_emph_vals, emph_colors

# ----------------------------------------------------------------------------------------
def make_meta_info_legend(plotdir, plotname, meta_info_key_to_color, emph_colors, all_emph_vals, meta_emph_formats=None, alpha=None):
    title = meta_info_key_to_color
    if meta_emph_formats is not None and meta_emph_formats.get(meta_info_key_to_color) not in ['len', None]:
        title = meta_emph_formats[meta_info_key_to_color].replace('@', ' ')
    elif title in legend_labels:
        title = legend_labels[title]
    emph_colors = [(v, c) for v, c in emph_colors if v in all_emph_vals]  # remove 'None' if there weren't any in the actual annotations
    if any(c==title for c, _ in emph_colors):  # if it's actually a color (i.e. probably a bool) no point in adding title)
        title = None
    lfn = plotname + '-legend'
    plot_legend_only({legend_groups.get(l, legend_labels.get(l, l)) : {'color' : c, 'alpha' : alpha} for l, c in emph_colors}, plotdir, lfn, title=title)
    return lfn

# # ----------------------------------------------------------------------------------------
# def _hls2hex(rgb_tuple):
#     h, l, s, alpha = rgb_tuple
#     hexstr = '#%02x%02x%02x' %tuple(map(lambda x: int(x*255), colorsys.hls_to_rgb(h, l, s)))
#     print rgb_tuple, hexstr
#     return hexstr

# ----------------------------------------------------------------------------------------
def getgrey(gtype='medium'):
    if gtype == 'medium':
        return '#929292'
    elif gtype == 'light-medium':
        return '#cdcdcd'
    elif gtype == 'light':
        return '#d3d3d3'
    elif gtype == 'white':
        return '#ffffff'
    else:
        assert False

# ----------------------------------------------------------------------------------------
def rgb_to_hex(rgb_tuple):
    assert len(rgb_tuple) == 3
    return '#%02x%02x%02x' %tuple([int(x*255) for x in rgb_tuple[:3]])

# ----------------------------------------------------------------------------------------
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

# ----------------------------------------------------------------------------------------
def get_color_norm(vals, remove_top_end=False, hard_min=None):
    if len(vals) == 0:
        print('  %s zero values passed to get_color_norm' % utils.color('yellow', 'warning'))
        vals = [0.]
    sorted_vals = sorted(vals)
    vmin = sorted_vals[0] - 0.2 * (sorted_vals[-1] - sorted_vals[0]) if hard_min is None else hard_min  # don't want anybody to be white, so set <vmin> to a bit less than the actual min value (i.e. so white corresponds to a value that's a bit less than any of our values)
    vmax = sorted_vals[-1]
    if remove_top_end:  # remove the top end of the color spectrum (at least when I'm adding this, it's because the Reds and Blues top color is way too close to black [and I can't set opacity on lines in ete3, which would also fix it])
        vmax = vmax + 0.3 * (vmax - vmin)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    return norm

# ----------------------------------------------------------------------------------------
def get_normalized_scalar_map(vals, cmap, remove_top_end=False, hard_min=None):
    # if cmap is None:
    #     cmap = plt.cm.Blues  # 'Blues'
    norm = get_color_norm(vals, remove_top_end=remove_top_end, hard_min=hard_min)
    scalarMap = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    return scalarMap

# ----------------------------------------------------------------------------------------
def get_smap_color(smap, info, key=None, val=None):  # specify *either* <key> or <val> (don't need <info> if you're passing <val>)
    if val is None:
        assert key is not None
        if key not in info or info[key] is None:
            return getgrey()
        val = info[key]
    rgb_code = smap.to_rgba(val)[:3]
    return rgb_to_hex(rgb_code)

# ----------------------------------------------------------------------------------------
def get_leg_entries(n_entries=5, vals=None, min_val=None, max_val=None, colorfcn=None):
    if min_val is None:
        min_val = min(vals)
    if max_val is None:
        max_val = max(vals)
    if min_val == max_val:
        max_val = min_val + (1 if min_val == 0 else 0.1 * min_val)
    max_diff = max(utils.eps, (max_val - min_val) / float(n_entries - 1))
    leg_vals = list(numpy.arange(min_val, max_val + utils.eps, max_diff))  # first value is exactly <min_val>, last value is exactly <max_val> (eps is to keep it from missing the last one)
    if colorfcn is None:  # just return the values, let the calling fcn work out the colors
        return leg_vals
    else:
        leg_entries = [(v, {'color' : colorfcn(v)}) for v in leg_vals]
        return collections.OrderedDict(leg_entries)

# ----------------------------------------------------------------------------------------
def expand_bounds(bounds, only_down=False):
    assert len(bounds) == 2
    return hutils.get_expanded_bounds(bounds, abs(bounds[1] - bounds[0]), only_down=only_down)

# ----------------------------------------------------------------------------------------
# returns modified copy of input list
def add_jitter(xvals, delta=None, frac=0.02):
    if delta is None:
        delta =  max(xvals) - min(xvals)
    jvals = numpy.random.uniform(-frac * delta, frac * delta, size=len(xvals))
    # jvals = [j for j in numpy.arange(-frac * delta, frac * delta, 2 * delta / float(len(xvals) - 1))]  # doesn't work yet, but maybe it'd make more sense to make it actual uniform, rather than uniform random?
    return [x + j for x, j in zip(xvals, jvals)]

# ----------------------------------------------------------------------------------------
def make_bool_hist(n_true, n_false, hist_label):
    """ fill a two-bin histogram with the fraction false in the first bin and the fraction true in the second """
    if 'partis.fraction_uncertainty' not in sys.modules:
        from . import fraction_uncertainty

    hist = Hist(2, -0.5, 1.5, ytitle='freq')

    def set_bin(numer, denom, ibin, label):
        frac = float(numer) / denom
        bounds = sys.modules['partis.fraction_uncertainty'].err(numer, denom)
        err = max(abs(frac - bounds[0]), abs(frac - bounds[1]))
        hist.set_ibin(ibin, frac, error=err, label=label)

    set_bin(n_true, n_true + n_false, 1, 'right')
    set_bin(n_false, n_true + n_false, 2, 'wrong')

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
# returns new hists with under/overflow bins shifted to adjacent bins (within plot range)
def shift_hist_overflows(hists, xmin, xmax):
    new_hists = []
    for ihist in range(len(hists)):
        htmp = copy.deepcopy(hists[ihist])
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

        if first_shown_bin == -1 or last_shown_bin == -1:
            raise RuntimeError(f"Can\'t shift overflows for hist '{htmp.title}' since it has no bins within shift range {xmin}, {xmax} (hist range: {htmp.xmin}, {htmp.xmax})")
        htmp.set_ibin(first_shown_bin,
                      underflows + htmp.bin_contents[first_shown_bin],
                      error=math.sqrt(under_err2 + htmp.errors[first_shown_bin]**2))
        htmp.set_ibin(last_shown_bin,
                      overflows + htmp.bin_contents[last_shown_bin],
                      error=math.sqrt(over_err2 + htmp.errors[last_shown_bin]**2))
        new_hists.append(htmp)
    return new_hists

# ----------------------------------------------------------------------------------------
# NOTE now you should set <hist> to None if you have more than one hist
def draw_no_root(hist, log='', plotdir=None, plotname='', more_hists=None, scale_errors=None, normalize=False, bounds=None, ybounds=None,
                 figsize=None, shift_overflows=False, colors=None, errors=False, write_csv=False, xline=None, yline=None, xyline=None, linestyles=None,
                 linewidths=None, plottitle=None, csv_fname=None, stats='', print_stats=False, translegend=(0., 0.), rebin=None,
                 xtitle=None, ytitle=None, markersizes=None, no_labels=False, only_csv=False, alphas=None, remove_empty_bins=False,
                 square_bins=False, xticks=None, xticklabels=None, yticks=None, yticklabels=None, leg_title=None, no_legend=False, hfile_labels=None,
                 text_dict=None, adjust=None, make_legend_only_plot=False, rotation=None):
    assert os.path.exists(plotdir)

    hists = [hist,] if hist is not None else []  # use <hist> if it's set (i.e. backwards compatibility for old calls), otherwise <hist> should be None if <more_hists> is set
    if more_hists is not None:
        hists = hists + more_hists
    if hist is None:  # gets used below for some label stuff, and yes this sucks and is ugly
        hist = hists[0]
    if sum(h.integral(True) for h in hists) == 0:
        print('  %s total integral of %d hists is zero, so not plotting anything' % (utils.wrnstr(), len(hists)))
        return None

    multiply_by_bin_width = False
    if normalize and len(set((h.n_bins, h.xmin, h.xmax) for h in hists)) > 1:
        print('    %s (%s) normalizing hists with different bins, which will *not* work/look right if there\'s empty bins (turning on square_bins should better show this)' % (utils.wrnstr(), plotname))
        multiply_by_bin_width = True

    xmin, xmax, ymin, ymax = None, None, None, None
    for htmp in hists:
        if htmp.title == 'null':  # empty hists
            continue
        if scale_errors is not None:
            factor = float(scale_errors[0]) if len(scale_errors) == 1 else float(scale_errors[hists.index(htmp)])
            for ibin in range(htmp.n_bins + 2):
                htmp.errors[ibin] *= factor
        if normalize:  # NOTE removed <normalization_bounds> option, hopefully I'm not using it any more
            htmp.normalize(multiply_by_bin_width=multiply_by_bin_width)
        if htmp.integral(True)>0 and (ymin is None or htmp.get_minimum(xbounds=bounds, exclude_empty='y' in log) < ymin):  # adding this afterwards, so might screw up something below
            ymin = htmp.get_minimum(xbounds=bounds, exclude_empty='y' in log)
        if ymax is None or htmp.get_maximum(xbounds=bounds) > ymax:
            ymax = htmp.get_maximum(xbounds=bounds)
        if htmp.integral(True) > 0:
            fbmin, fbmax = htmp.get_filled_bin_xbounds()
            if xmin is None or fbmin < xmin:  # overridden by <bounds> below
                xmin = fbmin
            if xmax is None or fbmax > xmax:
                xmax = fbmax

    if bounds is not None:
        xmin, xmax = bounds
    if ybounds is not None:  # ugly, but adding it long after the rest of the fcn
        ymin, ymax = [utils.non_none([yb, ym]) for ym, yb in zip((ymin, ymax), ybounds)]  # allows to set only one of them in the args

    if shift_overflows:
        if '_vs_per_gene_support' in plotname or '_fraction_correct_vs_mute_freq' in plotname or plotname in [r + '_gene' for r in utils.regions]:
            print(        '%s overriding overflow shifting for %s' % (utils.color('yellow', 'warning'), plotname))
        else:
            hists = shift_hist_overflows(hists, xmin, xmax)
        # assert '_vs_per_gene_support' not in plotname and '_fraction_correct_vs_mute_freq' not in plotname and plotname.find('_gene') != 1  # really, really, really don't want to shift overflows for these

    if write_csv:
        for ih, htmp in enumerate(hists):
            fn = utils.non_none([csv_fname, plotdir + '/' + plotname + '.csv'])
            if len(hists) > 1:
                fn = utils.insert_before_suffix('-'+(hfile_labels[ih] if hfile_labels is not None else htmp.title), fn)
            htmp.write(fn)

    if only_csv:
        return None

    # this is the slow part of plotting (well, writing the svg is also slow)
    fig, ax = mpl_init(figsize=figsize)
    mpl.rcParams.update({'legend.fontsize' : 15})

    tmpcolors = copy.deepcopy(colors)  # don't want to modify the arguments
    if tmpcolors is None:  # fiddle here http://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
        tmpcolors = ['royalblue', 'darkred', 'green', 'darkorange']
    n_different_colors = len(tmpcolors)
    while len(tmpcolors) < len(hists):
        tmpcolors += tmpcolors

    tmplinestyles = [] if linestyles is None else copy.deepcopy(linestyles)
    itmp = 0
    availstyles = ['-', '--', '-.', ':']
    while len(tmplinestyles) < len(hists):
        tmplinestyles += [availstyles[itmp % len(availstyles)] for _ in range(n_different_colors)]
        itmp += 1

    def floatstr(val):
        ndig = 1 if val > 1 else 3
        return ('%.'+str(ndig)+'f') % val
    for ih in range(len(hists)):
        htmp = hists[ih]
        statstr = None
        if stats == 'mean':
            statstr = ' (mean %s)' % floatstr(htmp.get_mean())
        elif stats == 'absmean':
            statstr = ' (abs av %s)' % floatstr(htmp.get_mean(absval=True))
        elif stats == '0-bin':
            statstr = ' (%s %s)' % ('0-bin' if htmp.bin_labels[1]=='' else htmp.bin_labels[1], floatstr(htmp.bin_contents[1]))
        elif stats is not None and stats != '':  # damnit, I ended up with both of the damn things as possible defaults
            raise Exception('unexpected stats str \'%s\'' % stats)
        markersize = None
        if markersizes is not None:
            imark = ih if len(markersizes) > 1 else 0
            # if imark > len(markersizes) - 1:
            #     print '  %s wrapping imark %d to match len of markersizes %d (%s)' % (utils.wrnstr(), imark, len(markersizes), markersizes)
            markersize = markersizes[imark % len(markersizes)]
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
        # i'm not sure why the linewidths get to here as strings, I guess that used to work, but now it kicks this really opaque error TypeError: Cannot cast array data from dtype('<U1') to dtype('float64') according to the rule 'safe'
        htmp.mpl_plot(ax, color=tmpcolors[ih], linewidth=linewidth, linestyle=tmplinestyles[ih], ignore_overflows=True, errors=errors, alpha=alpha, markersize=markersize, remove_empty_bins=remove_empty_bins, square_bins=square_bins)

    # NOTE it would be nice to combine xline, yline, and xyline (I don't want to go find everwhere that calls this right now)
    if xline is not None:
        ax.plot([xline, xline], [-0.1*ymax, 0.5*ymax], color='black', linestyle='--', linewidth=3)
    if yline is not None:
        print('%s fix y line' % utils.color('red', 'error'))
    if xyline is not None:
        assert len(xyline) == 2
        assert len(xyline[0]) == 2 and len(xyline[1]) == 2
        ax.plot([xyline[0][0], xyline[1][0]], [xyline[0][1], xyline[1][1]], color='black', linestyle='--', linewidth=3)
    # if yline is not None:
    #     # if yline < hframe.GetYaxis().GetXmin() or xline > hframe.GetYaxis().GetXmax():  # make sure we got valid a x position for the line
    #     #     print 'WARNING plotting y line at %f out of bounds (%f, %f)' % (float(ymin), hframe.GetYaxis().GetXmin(), hframe.GetYaxis().GetXmax())
    #     yl = TLine(hframe.GetXaxis().GetXmin(), yline, hframe.GetXaxis().GetXmax(), yline)
    #     yl.Draw()

    if xticklabels is not None and 'y' in log:  # if xticklabels is set we need to also set the y ones so the fonts match up
        if yticks is not None:
            print('  %s resetting yticks' % utils.color('yellow', 'warning'))
        yticks, yticklabels = get_auto_y_ticks(ymin, ymax, log=log)
        ymin = min(yticks + [ymin])
        ymax = max(yticks + [ymax])
    if xticks is None:
        if not no_labels and hist.bin_labels.count('') != len(hist.bin_labels):
            xticks = hist.get_bin_centers()
            xticklabels = hist.bin_labels

    if plottitle is not None:
        tmptitle = plottitle
    elif plotname in plotconfig.plot_titles:
        tmptitle = plotconfig.plot_titles[plotname]
    else:
        tmptitle = hist.title  # hm, maybe shouldn't be hist.title? I think that's usually supposed to be the legend
    if statstr is not None:
        tmptitle += statstr
        if print_stats:
            print('    %s %s' % (plotname, statstr))

    if xtitle is not None:
        tmpxtitle = xtitle
    elif plotname in plotconfig.xtitles:
        tmpxtitle = plotconfig.xtitles[plotname]
    else:
        tmpxtitle = hist.xtitle  # hm, maybe shouldn't be hist.title? I think that's usually supposed to be the legend

    if text_dict is not None:
        fig.text(text_dict['x'], text_dict['y'], text_dict['text'], fontdict={'weight' : 'bold'}) #, fontsize=)
    ymin = 0.8 * ymin if 'y' in log else ymin  # why tf was this here? -0.03*ymax
    default_adjust = {'left' : 0.2}
    if adjust is not None:
        default_adjust.update(adjust)
    fn = mpl_finish(ax, plotdir, plotname,
                    title=tmptitle,
                    xlabel=tmpxtitle,
                    ylabel=hist.ytitle if ytitle is None else ytitle,
                    xbounds=[xmin, xmax],
                    ybounds=[ymin, 1.15*ymax],
                    leg_loc=(0.72 + translegend[0], 0.7 + translegend[1]),
                    log=log, xticks=xticks, xticklabels=xticklabels, yticks=yticks, yticklabels=yticklabels,
                    no_legend=(no_legend or len(hists) <= 1), adjust=default_adjust, leg_title=leg_title, rotation=rotation)

    if make_legend_only_plot:
        plot_legend_only(None, plotdir, plotname+'-legend', ax=ax, title=tmpxtitle)

    return fn

# ----------------------------------------------------------------------------------------
# I'm not sure that the name of this function makes sense, for instance
#   - I don't think the violin and swarm options are really stacking anything
# need to describe what form <plotvals> should have
# i think you're supposed to either 1) set template_hist or 2) set no_hist and either swarm_plots or violin_plots
def stack_meta_hists(plotname, plotdir, mekey, plotvals, template_hist=None, colors=None, formats=None, only_csv=False, xtitle=None, normalize=False, swarm_plots=False, violin_plots=False, no_hist=False, xticks=None, is_bin_contents=False,
                     log='', figsize=None, ytitle=None, rotation=None, plottitle='', remove_empty_bins=True):
    all_emph_vals, emph_colors = meta_emph_init(mekey, formats=formats, all_emph_vals=set(plotvals))
    all_emph_vals = sorted(all_emph_vals, key=str)

    if not no_hist:
        if template_hist is None:
            raise Exception('need to specify template_hist if no_hist isn\'t set')
        hdict = collections.OrderedDict()
        for meval in all_emph_vals:
            hdict[meval] = Hist(value_list=None if is_bin_contents else plotvals[meval], template_hist=template_hist, title=str(meval), xtitle=xtitle)
            if is_bin_contents:
                hdict[meval].bin_contents = plotvals[meval]
        mcolors = {v : c for v, c in emph_colors}
        hist_list, hist_colors = zip(*[(h, mcolors[m]) for m, h in hdict.items()])
        hfile_labels = [h.title for h in hist_list]
        for htmp in hist_list:
            htmp.title = '%s (%s)' % (htmp.title, utils.round_to_n_digits(htmp.get_mean(), 2))
        if ytitle is None:
            ytitle = 'fraction of total' if normalize else 'counts'
        fn = draw_no_root(None, more_hists=list(hist_list), plotname=plotname, colors=list(hist_colors), plotdir=plotdir, write_csv=True, only_csv=only_csv, hfile_labels=hfile_labels,
                          shift_overflows=True, leg_title='%s (mean)'%mekey.rstrip('s'), alphas=[0.7 for _ in hist_list], linewidths=[5, 3, 3], markersizes=[15, 10, 8], errors=True, #square_bins=True, #errors=True,
                          remove_empty_bins=remove_empty_bins, ytitle=ytitle, log=log, figsize=figsize, rotation=rotation, plottitle=plottitle) # , normalize=True

    if violin_plots or swarm_plots:
        mvals, vlnvals = zip(*[(v, plotvals[v]) for v in all_emph_vals if v in plotvals])  # get order right
        if 'seaborn' not in sys.modules:
            import seaborn
        sns = sys.modules['seaborn']

    if violin_plots:
        ax = sns.violinplot(data=vlnvals)
        xticklabels = [str(v) for v in mvals]
        fn = mpl_finish(ax, plotdir, 'violin-'+plotname, xticklabels=xticklabels, xlabel=mekey.rstrip('s'), ylabel=xtitle)

    if swarm_plots:
        ax = sns.swarmplot(data=vlnvals)
        if colors is not None:  # should really apply these to the violin and hist plots as well
            assert len(colors) == len(ax.collections)
            for swrm, mval in zip(ax.collections, sorted(colors)):
                swrm.set_facecolors(colors[mval])
        xticklabels = [str(v) for v in mvals]
        ax.set_xticklabels(xticklabels) #, rotation='vertical', size=8 if xticklabelsize is None else xticklabelsize)
        fn = mpl_finish(ax, plotdir, 'swarm-'+plotname, xlabel=mekey.rstrip('s'), ylabel=xtitle, yticks=xticks)  # , xticklabels=xticklabels

    return fn

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
def make_mean_hist(hists, ignore_empty_bins=False):
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
        if ignore_empty_bins:
            vlist = [v for v in vlist if v > 0]
        if len(vlist) == 0:
            continue
        meanhist.set_ibin(ib, numpy.mean(vlist), error=0 if len(vlist)==1 else (numpy.std(vlist, ddof=1) / math.sqrt(len(vlist))))
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
timeticks = [0.1, 1, 10, 60, 600, 3600, 36000, 86400, 604800]  # seconds
timeticklabels = ['0.1 sec', '1 sec', '10 sec', '1 min', '10 min', '1 hour', '10 hours', '1 day', '1 week']

# would be nicer to call this 'titles' or something now that i'm adding a separate one for axes, but i don't want to change it
legends = {'vollmers-0.9' : 'VJ CDR3 0.9',
           'vjcdr3-0.9' : 'VJ CDR3 0.9',
           'vjcdr3-0.8' : 'VJ CDR3 0.8',
           # 'partition partis' : 'full partis',
           'partition' : 'full partis',
           # 'naive-hamming-partition partis' : 'point partis',
           'naive-hamming-partition' : 'point partis',
           # 'vsearch-partition partis' : 'vsearch partis',
           'vsearch-partition' : 'vsearch partis',
           'subset-partition' : 'subset partis',
           'vsearch-subset-partition' : 'vsearch subset partis',
           'star-partition' : 'star partis',
           'phylo-naive-iqtree' : 'IQ-TREE (CDR3 N)',
           'phylo-naive-iqtree-fuzz-0' : 'IQ-TREE (fuzz 0)',
           'phylo-naive-iqtree-fuzz-1' : 'IQ-TREE (fuzz 1)',
           'phylo-naive-iqtree-fuzz-2' : 'IQ-TREE (fuzz 2)',
           'single-chain-partis' : 'single chain partis',
           'annotate' : '1-seq. partis',
           'seed-partition' : 'full partis (seed)',
           'seed-naive-hamming-partition' : 'point partis (seed)',
           'changeo' : 'IMGT + Change-O',
           'scoper' : 'SCOPer',
           'mobille' : 'MobiLLe',
           'igblast' : 'IgBLAST',
           'linearham' : 'linearham',
           'iqtree-tree-perf' : 'IQ-TREE',  #  1.6.12
           'iqtree-1.6.beta3-tree-perf' : 'IQ-TREE 1.6.beta3',
           'iqtree-2.3.1-tree-perf' : 'IQ-TREE 2.3.1',
           'raxml-tree-perf' : 'RAxML',
           'igphyml-tree-perf' : 'IgPhyML',
           'gctree-tree-perf' : 'GCtree',
           'gctree-base-tree-perf' : 'GCtree (base)',
           'gctree-mut-mult-tree-perf' : 'GCtree (context)',
           'gctree-no-dag-tree-perf' : 'GCtree (no DAG)',
           # '0.1-true-singletons' : '10% random singletons',
           # '0.1-true-reassign' : '10% random reassign',
           'misassign-0.60-singletons' : 'synth. 60%\nsingleton',
           'synth-singletons-0.20' : 'synth. 20%\nsingleton',
           'misassign-0.10-reassign' : 'synth. 10%\nreassign',
           'misassign-distance-0.03' : 'synth.\nneighbor 0.03',
           'synth-distance-0.03' : 'synth.\nneighbor 0.03',
           'mixcr' : 'MiXCR',
           'adj_mi' : 'adj MI',
           'ccf_under' : 'precision',
           'ccf_over' : 'sensitivity',
           'ccf_product' : 'F1 score',
           'f1' : 'F1 score',
           'coar' : 'COAR',
           'mrca' : 'MRCA',
           'rf' : 'RF (weighted)',
           'cln-frac' : 'collision frac.',
           'n-clusters' : 'N clusters',
           'n-leaves' : 'family size',
           'constant-number-of-leaves' : '',
           'n-sim-events' : 'N families',
           'scratch-mute-freq' : 'SHM fraction (nuc)',
           'obs-times' : 'N generations',
           'mfreq' : 'SHM fraction (nuc)',
           'mfreq-pct' : 'SHM % (nuc)',
           'time-reqd' : 'time required',
           'pcfrac-correct' : 'correctly paired',
           'pcfrac-mispaired' : 'mispaired',
           'pcfrac-unpaired' : 'unpaired',
           'pcfrac-correct-family' : 'paired with correct family',
           'pcfrac-near-family' : 'paired with similar family',
           'pcfrac-correct-ns' : 'correctly paired (non-singleton)',
           'pcfrac-mispaired-ns' : 'mispaired (non-singleton)',
           'pcfrac-unpaired-ns' : 'unpaired (non-singleton)',
           'pcfrac-correct-family-ns' : 'paired with correct family (non-singleton)',
           'naive-hdist' : 'mean N incorrect bases', #'ham. dist. to true naive',  # NOTE duplicates/similar to entries in plotconfig.py
           'n-seqs' : 'N seqs',
           'biggest-logprob-cluster-to-calculate' : 'max calc\'d cluster size',
           'dataset-in' : 'input dataset',
           'bulk-data-fraction' : 'bulk data frac',
           'max-abundances' : 'max abundance\nwgtd. simu - data',
           'distr-abundances' : 'abundance distr.\nwgtd. simu - data',
           'distr-hdists' : 'SHM distr\nsimu - data',
           'xscale-test-mae' : 'MAE\nxscale test',
           'xscale-train-vs-test' : 'xscale\nabs(train - test) width',
           'xshift' : 'xshift\nmean abs(true - inf)',
           'xshift-test-mae' : 'MAE\nxshift test',
           'xshift-train-vs-test' : 'xshift\nabs(train - test) width',
           'n-trials' : 'N trees',
           'n-trees-per-expt' : 'N trees/expt.',
           'vrc01-muts' : 'N VRC01 muts',
           'gct-data' : 'gctree',
           'gct-data-d15' : 'gctree d15',
           'gct-data-d20' : 'gctree d20',
           'gct-data-w10' : 'gctree w10',
           'bst-data-d20' : 'beast d20',
           'bst-data' : 'data (beast)',
           'iqt-data' : 'data (iqtree)',
           'iqt-data-d15' : 'iqtree d15',
           'iqt-data-d20' : 'iqtree d20',
           'iqt-data-w10' : 'iqtree w10',
           'simu-iqtree' : 'simu (iqtree)',
           }

axis_labels = {
    'pcfrac-correct' : 'frac. correctly paired',
    'pcfrac-mispaired' : 'frac. mispaired',
    'pcfrac-unpaired' : 'frac. unpaired',
    'pcfrac-correct-family' : 'frac. correct family',
    'pcfrac-near-family' : 'frac. similar family',
    'pcfrac-correct-ns' : 'frac. correctly paired',
    'pcfrac-mispaired-ns' : 'frac. mispaired',
    'pcfrac-unpaired-ns' : 'frac. unpaired',
    'pcfrac-correct-family-ns' : 'frac. correct family',
}

val_cfgs = {
    'legends' : {
        'constant-number-of-leaves' : {'0' : {'n-leaves' : {'default' : 'geom.', 'hist' : 'data'}}, '1' : 'const.'},
        'n-leaves' : {'hist' : 'distr.'},
    },
    'colors' : {
        'default' : 'black',
        'n-leaves' : {'hist' : 'darkred', '1' : 'black', '2' : '#006600', '3' : '#1f77b4', '10' : '#ff7f0e'},
    },
    'linestyles' : {
        'default' : '-',
        'constant-number-of-leaves' : {'0' : 'solid', '1' : 'dashed'},
    }
}

dmp = plt.get_cmap("Dark2").colors
colors = {'true' : '#006600',
          'partition' : '#cc0000',  # 8c001a',
          'vsearch-partition' : '#990012',  #c04000',
          'subset-partition' : '#ff7f0e',  #c04000',
          'vsearch-subset-partition' : '#ff7f0e',  #c04000',
          'star-partition' :  '#ff7f0e',
          'phylo-naive-iqtree' :  'grey',
          'phylo-naive-iqtree-fuzz-0' :  'grey',
          'phylo-naive-iqtree-fuzz-1' :  'grey',
          'phylo-naive-iqtree-fuzz-2' :  'grey',
          'single-chain-partis' : '#006600',
          'single-chain-scoper' : '#006600',
          'annotate' : '#1f77b4',
          'naive-hamming-partition' : '#990012',
          'seed-partition' : '#990012',
          'seed-naive-hamming-partition' : '#990012',
          'vollmers-0.5' : '#3333ff',
          'vollmers-0.9' : '#3399ff',
          'vjcdr3-0.9' : '#3399ff',
          'vjcdr3-0.8' : '#3399ff',
          'changeo' :  '#2b65ec',
          'scoper' : '#2b65ec',
          'mobille' : '#a821c7',
          'igblast' : '#a821c7',
          'linearham' : '#006600',
          'iqtree-tree-perf' : dmp[0],
          'iqtree-1.6.beta3-tree-perf' : dmp[0],
          'iqtree-2.3.1-tree-perf' : dmp[0],
          'raxml-tree-perf' : dmp[1],
          'igphyml-tree-perf' : dmp[2],
          'gctree-tree-perf' : dmp[3],
          'gctree-base-tree-perf' : dmp[3],
          'gctree-mut-mult-tree-perf' : dmp[3],
          'gctree-no-dag-tree-perf' : dmp[3],
          'enclone' : 'green',
          'mixcr' : '#2b65ec',
          'misassign-0.60-singletons' : '#808080',
          'synth-singletons-0.20' : '#808080',
          'misassign-0.10-reassign' : '#808080',
          'misassign-distance-0.03' : '#808080',
          'synth-distance-0.03' : '#808080',
}

linewidths = {'true' : 15,
              'vsearch-partition' : 3,
              'vsearch-subset-partition' : 3,
              'star-partition' :  2,
              'single-chain-partis' : 2,
              'single-chain-scoper' : 2,
              'annotate' : 2,
              'naive-hamming-partition' : 3,
              'seed-partition' : 2,
              'seed-naive-hamming-partition' : 2,
              'partition' : 6,
              'vollmers-0.5' : 4,
              'vollmers-0.9' : 6,
              'changeo' : 3,
              'scoper' : 3,
              'mobille' : 3,
              'igblast' : 3,
              'mixcr' : 6,
              'misassign-0.60-singletons' : 4,
              'synth-singletons-0.20' : 4,
              'misassign-0.10-reassign' : 3,
              'misassign-distance-0.03' : 2,
              'synth-distance-0.03' : 2,
}

linestyles = {'naive-hamming-partition' : 'dashed',
              'vsearch-partition' : 'dotted',
              'vsearch-subset-partition' : 'dotted',
              'star-partition' :  'dashed',
              'single-chain-partis' : 'dashed',
              'single-chain-scoper' : 'dotted',
              'annotate' : 'dashed',
              'changeo' : 'dashed',
              'scoper' : 'dashed',
              'mobille' : 'dashed',
              'phylo-naive-iqtree-fuzz-0' : 'dotted',
              'phylo-naive-iqtree-fuzz-1' : 'dashdot',
              'phylo-naive-iqtree-fuzz-2' : 'dashed',
              'igblast' : 'dotted',
              'mixcr' : 'dotted',
              'misassign-distance-0.03' : 'dashed',
              'synth-distance-0.03' : 'dashed',
              'gctree-base-tree-perf' : 'dashed',
              'gctree-no-dag-tree-perf' : 'dashed',
              'gctree-mut-mult-tree-perf' : 'dotted',
}

alphas = {'true' : 0.6,
          'vollmers-0.9' : 0.6,
          'misassign-0.60-singletons' : 0.5,
          'synth-singletons-0.20' : 0.5,
          'misassign-distance-0.03' : 0.8,
          'synth-distance-0.03' : 0.8,
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
# pass in <fnames> instead of <hists> if you want the bins to match
def plot_cluster_size_hists(plotdir, plotname, hists, title='', xmin=None, xmax=None, log='xy', normalize=False, hcolors=None, ytitle=None, fnames=None, translegend=None, alphas=None,
                            stacked_bars=False, no_legend=False):
    if fnames is not None:
        # NOTE this hasn't really been tested (with new hutils fcn)
        assert hists is None
        cslists = collections.OrderedDict()
        for label, fn in fnames.items():
            cslists[label] = [len(c) for c in ClusterPath(fname=fn).best()]
        hists = hutils.make_hists_from_lists_of_values(cslists, 'cluster-size', is_log_x=True)

    hist_list, tmpcolors = [], []
    for ih, (name, hist) in enumerate(hists.items()):
        if 'misassign' in name:
            continue
        if 'vollmers' in name:
            if '0.7' in name or '0.8' in name or '0.95' in name or '0.5' in name:
                continue
        if hist.integral(True) == 0:
            continue

        if hcolors is None:
            tmpcolors.append(colors.get(name, 'grey' if len(hists)==1 else default_colors[ih%len(default_colors)]))
        else:
            tmpcolors.append(hcolors.get(name, 'grey'))
        hxmin, hxmax = hist.get_filled_bin_xbounds()
        if xmin is None or hxmin < xmin:
            xmin = hxmin
        if xmax is None or hxmax > xmax:
            xmax = hxmax
        hist.title = legends.get(name, name)
        hist_list.append(hist)
        xticks = hist.get_bin_centers()
        xticklabels = hist.bin_labels

    rotation = None
    if any(len(t) > 5 for t in xticklabels):
        rotation = 'vertical'

    if len(hist_list) == 0:
        return

    if alphas is None:
        alphas = [0.7 for _ in hist_list]
    if 'x' in log:
        if xmin < 1:  # the above gives us the bin low edge, which with log x scale is way too far left of the lowest point
            xmin = 0.9
        xmax *= 1.025
    # xticks, xticklabels = get_cluster_size_xticks(xmin, xmax)  # NOTE could also pass list of hists in here to get xmin, xmax
    if ytitle is None:
        ytitle = '%s of clusters' % ('fraction' if normalize else 'number')
    if translegend is None:
        translegend = (0, 0) if len(hists)==1 else (-0.7, -0.65)
    ybounds = None
    if 'y' not in log:
        ybounds = [-0.05, None]

    if stacked_bars:
        fig, ax = mpl_init()
        base_vals = None
        for ih, thist in enumerate(hist_list):
            ax.bar(thist.bin_labels, thist.bin_contents, label=thist.title, bottom=base_vals, alpha=alphas[ih], color=tmpcolors[ih], linewidth=0)
            base_vals = [(base_vals[i] if base_vals is not None else 0) + c for i, c in enumerate(thist.bin_contents)]
        return mpl_finish(ax, plotdir, plotname, xbounds=(0.5, thist.n_bins + 0.5), xlabel='cluster size', ylabel=ytitle, no_legend=no_legend, log=log, rotation=rotation) #, leg_loc=(0.1, 0.2), xbounds=(xticks[0], xticks[-1]), ybounds=(ymin, 1.01), title=title, xlabel=xlabel, ylabel='metric value')
    else:
        return draw_no_root(None, more_hists=hist_list, plotdir=plotdir, plotname=plotname, log=log, normalize=normalize, remove_empty_bins=True, colors=tmpcolors, xticks=xticks, xticklabels=xticklabels,
                            bounds=(xmin, xmax), ybounds=ybounds, plottitle=title, xtitle='cluster size', ytitle=ytitle, alphas=alphas, translegend=translegend, linewidths=[5, 2], markersizes=[20, 8],
                            no_legend=no_legend, errors=True, rotation=rotation, figsize=(6, 4))

# ----------------------------------------------------------------------------------------
def plot_tree_mut_stats(args, plotdir, antn_list, is_simu, only_leaves=False, only_csv=False, fnames=None):
    # ----------------------------------------------------------------------------------------
    def add_to_distr_dict(ucounts, udistr):
        for scount in ucounts.values():
            if scount not in udistr:
                udistr[scount] = 0
            udistr[scount] += 1
    # ----------------------------------------------------------------------------------------
    def add_antn(iln, line):
        unique_seqs, unique_muts = {}, {}
        for uid, mseq in zip(line['unique_ids'], line['seqs']):
            if only_leaves:
                node = treefos[iln]['tree'].find_node_with_taxon_label(uid)
                if not node.is_leaf():
                    continue
            if mseq not in unique_seqs:
                unique_seqs[mseq] = 0
            unique_seqs[mseq] += 1
            umuts = utils.get_mut_codes(line['naive_seq'], mseq) #, amino_acid=False debug=True
            for mcd in umuts:
                if mcd['str'] not in unique_muts:
                    unique_muts[mcd['str']] = 0
                unique_muts[mcd['str']] += 1
            n_muts = utils.per_seq_val(line, 'n_mutations', uid)
            if n_muts not in n_mut_dict:
                n_mut_dict[n_muts] = 0
            n_mut_dict[n_muts] += 1
        add_to_distr_dict(unique_seqs, useq_distr)
        add_to_distr_dict(unique_muts, umut_distr)
    # ----------------------------------------------------------------------------------------
    def finalize(udistr, plotname, title):
        hist = hutils.make_hist_from_dict_of_counts(udistr, 'int', '')
        fn = hist.fullplot(plotdir, plotname, pargs={'remove_empty_bins' : True}, fargs={'xlabel' : 'N observations', 'ylabel' : 'counts', 'title' : title, 'log' : 'xy'}, only_csv=only_csv) #, texts=[(0.7, 0.8, 'N gtypes %d'%len(unique_seqs))])
        if fnames is not None:
            fnames[-1].append(fn)
    # ----------------------------------------------------------------------------------------
    print('    plotting tree mutation stats %s' % ('using only leaves' if only_leaves else 'with all seqs'))
    from . import lbplotting
    utils.prep_dir(plotdir, wildlings=['*.csv', '*.svg'])
    if fnames is not None:
        fnames.append([])
    useq_distr, umut_distr, n_mut_dict = {}, {}, {}
    treefos = treeutils.get_treefos(args, antn_list) if only_leaves else None  # missing some args to this fcn here since i don't want to propagate them, but i don't think it matters
    for iln, line in enumerate(antn_list):
        add_antn(iln, line)
    for ulabel, udistr in zip(['seq', 'mut'], (useq_distr, umut_distr)):
        finalize(udistr, 'u%s_distr'%ulabel, 'unique %s distr'%ulabel)
    finalize(n_mut_dict, 'n_muts', 'distance to root')

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
# NOTE can use get_leg_entries()
def plot_legend_only(leg_entries, plotdir, plotname, ax=None, title=None, n_digits=None):
    lfn = plotdir+'/'+plotname+'.svg'
    if leg_entries is None:
        assert ax is not None
        handles, labels = ax.get_legend_handles_labels()
        leg_entries = collections.OrderedDict((l, None) for l in labels)
    if len(leg_entries) == 0:
        return lfn
    max_label_len = max(len(str(l)) for l in leg_entries)
    figlegend = plt.figure(figsize=(2 + max_label_len / 12., 2 + len(leg_entries) / 4.))
    if ax is None:
        fig, ax = mpl_init()
        for tlab, lfo in sorted((str(l), fo) for l, fo in leg_entries.items()):  # have to convert labels to str in case e.g. one of them's None
            if n_digits is not None and tlab is not None:
                tlab = utils.round_to_n_digits(float(tlab), 2)
            ax.plot([None], [None], label=str(tlab), color=lfo['color'], linewidth=lfo.get('linewidth', 5), linestyle=lfo.get('linestyle', '-'), alpha=lfo.get('alpha', 0.6))  # str() is to convert None to 'None', otherwise it doesn't show up
        handles, labels = ax.get_legend_handles_labels()
    figlegend.legend(handles, labels, loc='center', title=title)
    figlegend.savefig(lfn)
    return lfn

# ----------------------------------------------------------------------------------------
def mpl_init(figsize=None, fsize=20, label_fsize=15, fweight='bold'):
    if 'seaborn' not in sys.modules:
        import seaborn  # really #$*$$*!ing slow to import, but only importing part of it doesn't seem to help
    sys.modules['seaborn'].set_style('ticks')
    mpl.rcParams.update({
        'font.size': fsize, 'axes.titlesize': fsize, 'axes.labelsize': fsize,
        'xtick.labelsize': label_fsize, 'ytick.labelsize': label_fsize,  # NOTE this gets (maybe always?) overriden by xticklabelsize/yticklabelsize in mpl_finis()
        'legend.fontsize': fsize,
        'font.family': 'Lato',
        'font.weight': fweight,
        'axes.labelweight': fweight,
        'axes.titleweight': fweight,
        'figure.autolayout': True})
    fig, ax = plt.subplots(figsize=figsize)
    fig.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.92)

    return fig, ax

# ----------------------------------------------------------------------------------------
# initially copied from https://matplotlib.org/3.5.1/gallery/lines_bars_and_markers/scatter_piecharts.html
def plot_pie_chart_marker(ax, xpos, ypos, radius, fracfos, alpha=None, debug=False):
    np = numpy
    total = sum(f['fraction'] for f in fracfos)
    if not utils.is_normed(total):
        print('  %s fractions add to %f (should add to 1): %s' % (utils.wrnstr(), total, '  '.join('%12s %-.3f'%(f['label'], f['fraction']) for f in sorted(fracfos, key=lambda x: x['fraction'], reverse=True))))

    if debug:
        print('   frac    label    N   min   max')
        lwd = max(len(f['label']) for f in fracfos)
    total = 0
    for ifo, ffo in enumerate(fracfos):
        if ffo['fraction'] == 0:
            continue
        lnsp = np.linspace(total, total + ffo['fraction'] if ifo < len(fracfos) - 1 else 1)  # evenly spaced values from <total> to <total + frac> for this slice's <frac>
        xvals = np.cos(2 * np.pi * lnsp)
        yvals = np.sin(2 * np.pi * lnsp)
        xyvals = np.row_stack([[0, 0], np.column_stack([xvals, yvals])])
        s1 = np.abs(xyvals).max()  # max x or y val (i guess s= arg in ax.scatter() is based on max x or y size?)
        ax.scatter([xpos], [ypos], marker=xyvals, s=(270*radius*s1)**2, facecolor=ffo['color'], alpha=alpha, linewidth=0)  # s= is in "points squared", but radius is in axis/fig coords ([0, 1], or maybe [-1, 1]?), and I can't figure out how to convert and I'm tired of googling so using 275 which seems about right, hopefully it keeps working
        total += ffo['fraction']
        if debug:
            print('   %.3f  %s  %3d  %5.3f %5.3f %.3f %.3f' % (ffo['fraction'], utils.wfmt(ffo['label'], lwd), len(lnsp), min(lnsp), max(lnsp), s1, max(max(x, y) for x, y in zip(xvals, yvals) ))) #, [math.sqrt(x*x + y*y) for x, y in zip(xvals, yvals)])
    if debug:
        print('')

# ----------------------------------------------------------------------------------------
def mpl_finish(ax, plotdir, plotname, title='', xlabel='', ylabel='', xbounds=None, ybounds=None, leg_loc=(0.04, 0.6), leg_prop=None, log='',
               xticks=None, xticklabels=None, xticklabelsize=None, yticklabelsize=None, yticks=None, yticklabels=None, no_legend=False, adjust=None,
               suffix='svg', leg_title=None, legend_fontsize=None, fig=None, right_y_axis=False, rotation=None):
    if 'seaborn' not in sys.modules:
        import seaborn  # really #$*$$*!ing slow to import, but only importing part of it doesn't seem to help
    if not no_legend:
        handles, labels = ax.get_legend_handles_labels()
        if len(handles) > 0:
            if len(handles) > 5:
                leg_loc = leg_loc[0], leg_loc[1] - 0.15
            legend = ax.legend(handles, labels, loc=leg_loc, prop=leg_prop, title=leg_title, fontsize=legend_fontsize)
    default_adjust = {'bottom' : 0.20, 'left' : 0.18, 'right' : 0.95, 'top' : 0.92}
    if adjust is not None:
        default_adjust.update(adjust)  # this makes things still work if the calling fcn has {} as the default rather than None
    plt.gcf().subplots_adjust(**default_adjust)  # ok, it's not default any more, but I don't want to change the name
    sys.modules['seaborn'].despine(right=not right_y_axis)  #trim=True, bottom=True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if 'x' in log:
        ax.set_xscale('symlog')  # 'log' used to work, but now it screws up the x axis labels
    if 'y' in log:
        ax.set_yscale('log')
        if yticks is None:
            yticks, yticklabels = get_auto_y_ticks(ax.get_ylim()[0], ax.get_ylim()[1], log=log)
    if xticks is not None:  # if these are after the xbound setting they override it
        plt.xticks(xticks)
    if yticks is not None:
        plt.yticks(yticks)
    if xticklabels is not None:
        # xticklabels = [str(i) for i, _ in enumerate(xticklabels)]
        ax.set_xticklabels(xticklabels, size=xticklabelsize, rotation=rotation)
        if fig is None:
            fig = ax.get_figure()
        renderer = fig.canvas.get_renderer()
        bboxes = []
        for xlab in ax.get_xticklabels():
            if xlab.get_text():  # skip empty labels
                # print(xlab.get_text(), xlab.get_window_extent(renderer))
                bbox = xlab.get_window_extent(renderer)
                bboxes.append(bbox)
        # print([bboxes[i].overlaps(bboxes[i + 1]) for i in range(len(bboxes) - 1)])
        if any(bboxes[i].overlaps(bboxes[i + 1]) for i in range(len(bboxes) - 1)):
            ax.set_xticklabels(xticklabels, rotation='vertical', size=14 if xticklabelsize is None else xticklabelsize)
    if xbounds is not None and xbounds[0] != xbounds[1]:
        plt.xlim(xbounds[0], xbounds[1])
    if ybounds is not None and ybounds[0] != ybounds[1]:
        plt.ylim(ybounds[0], ybounds[1])
    if yticklabels is not None:  # NOTE there's some complicated auto log y tick stuff in draw_no_root()
        ax.set_yticklabels(yticklabels, size=yticklabelsize)
    plt.title(title, fontweight='bold', fontsize=20 if len(title) < 25 else 15)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    fullname = plotdir + '/' + plotname + '.' + suffix
    plt.savefig(fullname)
    plt.close('all')
    # subprocess.check_call(['chmod', '664', fullname])
    return fullname  # this return is being added long after this fcn was written, so it'd be nice to go through all the places where it's called and take advantage of the return value

# ----------------------------------------------------------------------------------------
def get_auto_y_ticks(ymin, ymax, log=''):  # NOTE pretty similar to get_cluster_size_xticks() (for search: log_bins log bins)
    # ----------------------------------------------------------------------------------------
    def tstr(y):
        y = utils.round_to_n_digits(y, 3)
        if 'y' in log:
            if y > 1000 or y < 0.01:
                return '%.0e' % y
            else:
                return str(y) #'%.1f' % y
        else:
            return str(y)
    # ----------------------------------------------------------------------------------------
    if ymin == 0:
        ymin = 1e-10  # not sure what to use here
    tstart, tstop = math.floor(math.log(ymin, 10)), math.ceil(math.log(ymax, 10))
    n_ticks = tstop - tstart + 1
    if n_ticks == 2:  # i don't think it can be 0 or 1, but not sure
        n_ticks += 1
    yticks = [y for y in numpy.logspace(tstart, tstop, num=n_ticks)]
    # ymax = yticks[-1]
    if ymax > 1.4*yticks[-1]:
        yticks.append(ymax)
    return yticks, [tstr(t) for t in yticks]

# ----------------------------------------------------------------------------------------
def plot_csim_matrix_from_files(plotdir, plotname, meth1, ofn1, meth2, ofn2, n_biggest_clusters, title='', debug=False):
    # fpath = 'partitions/sizes/cluster-sizes'
    partitions = {}
    for mstr, ofn in zip((meth1, meth2), (ofn1, ofn2)):
        _, _, cpath = utils.read_output(ofn, skip_annotations=True, is_partition_file=True)
        partitions[mstr] = cpath.best()
    plot_cluster_similarity_matrix(plotdir, plotname, meth1, partitions[meth1], meth2, partitions[meth2], n_biggest_clusters, title=title, debug=debug)

# this script snippet may be useful:
# import plotting
# bd = '/fh/fast/matsen_e/dralph/partis/paired-loci/lcdr3/v1/seed-0/allowed-cdr3-lengths-igh,12-13'
# inf_fn = '%s/partis/partition-igh.yaml' % bd
# tru_fn = '%s/simu/igh.yaml' % bd
# def getptn(fn):
#     _, _, cpath = utils.read_output(fn)
#     return cpath.best()
# tru_ptn, inf_ptn = [getptn(f) for f in [tru_fn, inf_fn]]
# plotdir = '%s/csim-plots' % bd
# print plotting.plot_cluster_similarity_matrix(plotdir, 'csim-matrix', 'true', tru_ptn, 'inferred', inf_ptn, 30, debug=True) #, debug=True)
# print plotting.plot_cluster_similarity_matrix(plotdir, 'csim-matrix', 'inferred', inf_ptn, 'true', tru_ptn, 30, debug=True) #, debug=True)
# sys.exit()

# ----------------------------------------------------------------------------------------
# iscn_denominator: 'max' if you want to look at a method that's oversplitting since it will show how much of the larger cluster is split among various clusters in the other partition
#                   'min'                                        overmerging, since... eh, maybe not? not sure
def plot_cluster_similarity_matrix(plotdir, plotname, meth1, partition1, meth2, partition2, n_biggest_clusters, iscn_denominator='max', title='', debug=False):
    # ----------------------------------------------------------------------------------------
    def get_ticks(clens, n_big):
        modulo = 2
        # axis_max = min(len(partition1), min(n_biggest_clusters, len(partition2)))
        axis_max = min(n_big, len(clens))
        if axis_max > 20:
            modulo = 3
        if n_big > 40:
            modulo = int(n_big / 15.)
        return modulo, axis_max, [n - 0.5 for n in range(1, axis_max + 1, modulo)]
    # ----------------------------------------------------------------------------------------
    # partition1 = [['4'], ['7', '8'], ['6', '5'], ['99', '3', '1']]
    # # partition2 = [['1', '2', '3'], ['4'], ['5', '6'], ['7', '8']]
    # partition2 = [['3'], ['5'], ['6'], ['7'], ['8'], ['99', '3', '4']]

    a_cluster_lengths, b_cluster_lengths, smatrix = utils.partition_similarity_matrix(partition1, partition2, n_biggest_clusters, iscn_denominator=iscn_denominator, a_label=meth1, b_label=meth2, debug=debug)

    fig, ax = mpl_init()
    data = numpy.array(smatrix)
    cmap = plt.cm.get_cmap('viridis') #Blues  #cm.get_cmap('jet')
    cmap.set_under('w')
    heatmap = ax.pcolormesh(data, cmap=cmap, vmin=0., vmax=1.)
    cbar = plt.colorbar(heatmap, label='overlap / %s size' % iscn_denominator.replace('min', 'smaller').replace('max', 'larger'), pad=0.09)
    
    n_big_a, n_big_b = n_biggest_clusters if hasattr(n_biggest_clusters, '__iter__') else (n_biggest_clusters, n_biggest_clusters)
    xmod, x_axis_max, xticks = get_ticks(b_cluster_lengths, n_big_b)
    ymod, y_axis_max, yticks = get_ticks(a_cluster_lengths, n_big_a)
    xticklabels = [str(b_cluster_lengths[it]) for it in range(0, len(b_cluster_lengths), xmod)]
    yticklabels = [str(a_cluster_lengths[it]) for it in range(0, len(a_cluster_lengths), ymod)]
    return mpl_finish(ax, plotdir, plotname, title=title, xlabel='%s cluster size'%legends.get(meth2, meth2), ylabel='%s cluster size'%legends.get(meth1, meth1),
                     xticks=xticks, yticks=yticks, xticklabels=xticklabels, yticklabels=yticklabels, xticklabelsize=15, yticklabelsize=15,
                     xbounds=(0, x_axis_max), ybounds=(0, y_axis_max), rotation='vertical')

# ----------------------------------------------------------------------------------------
# NOTE set unset/empty values to float('nan') to keep them transparent
def plot_smatrix(plotdir, plotname, xydicts=None, xylists=None, kfcn=None, n_max_bins=None, smatrix=None, xybins=None, float_vals=False, lfcn=lambda x: str(x), y_lfcn=None, xlabel='', ylabel='', title='', blabel='counts', tdbg=False):
    # ----------------------------------------------------------------------------------------
    def truncate_bins(xbins, ybins, n_max_bins):
        print('    truncating x bins %d --> %d and ybins %d --> %d' % (len(xbins), n_max_bins, len(ybins), n_max_bins))
        return xbins[:n_max_bins], ybins[:n_max_bins]
    # ----------------------------------------------------------------------------------------
    def get_smatrix_from_xy_dicts(xvals, yvals, kfcn=None, n_max_bins=None):  # NOTE lots of duplication with next fcn (but I think it's cleaner to have them separate)
        xbins, ybins = [[-1] + sorted(sorted(set(td.values())), key=kfcn) for td in (xvals, yvals)]  # sort first with default order (e.g. alphabetical) then with the supplied kfcn (e.g. len) so that e.g. length ties have a defined order
        if n_max_bins is not None:  # it would be nice to skip bins based on the number of counts, but then we have to resize <smatrix> and the bins afterwards, and this bin sorting at least atm will usually give us the highest count bins
            xbins, ybins = truncate_bins(xbins, ybins, n_max_bins)
        smatrix = [[0 for _ in xbins] for _ in ybins]
        if tdbg > 1:
            uid_matrix = [[[] for _ in xbins] for _ in ybins]
        n_skipped = 0
        for uid in sorted(set(yvals) | set(xvals)):
            yv, xv = yvals.get(uid, -1), xvals.get(uid, -1)
            if yv not in ybins or xv not in xbins:
                n_skipped += 1
                continue
            smatrix[ybins.index(yv)][xbins.index(xv)] += 1
            if tdbg > 1:
                uid_matrix[ybins.index(yv)][xbins.index(xv)].append(uid)
        if tdbg > 1:
            lb = str(max(len(str(b)) for b in ybins + xbins))  # max length (when converted to str) of any bin label
            print('  uids in smatrix')
            for iy, yb in enumerate(ybins):
                for ix, xb in enumerate(xbins):
                    print(('  %'+lb+'s   %'+lb+'s    %4d   %s') % (yb, xb, len(uid_matrix[iy][ix]), ':'.join(uid_matrix[iy][ix])))

        return smatrix, xbins, ybins, n_skipped
    # ----------------------------------------------------------------------------------------
    def get_smatrix_from_xy_lists(xvals, yvals, kfcn=None, n_max_bins=None):  # NOTE lots of duplication with previous fcn (but I think it's cleaner to have them separate)
        assert len(xvals) == len(yvals)
        xbins, ybins = [sorted(set(tl), key=kfcn) for tl in (xvals, yvals)]
        if n_max_bins is not None:  # it would be nice to skip bins based on the number of counts, but then we have to resize <smatrix> and the bins afterwards, and this bin sorting at least atm will usually give us the highest count bins
            xbins, ybins = truncate_bins(xbins, ybins, n_max_bins)
        smatrix = [[0 for _ in xbins] for _ in ybins]
        for xv, yv in zip(xvals, yvals):
            smatrix[ybins.index(yv)][xbins.index(xv)] += 1

        return smatrix, xbins, ybins, 0
    # ----------------------------------------------------------------------------------------
    assert xydicts is not None or xylists is not None or smatrix is not None  # specify exactly one of them
    if y_lfcn is None:
        y_lfcn = lfcn
    if smatrix is None:
        sfcn = get_smatrix_from_xy_dicts if xydicts is not None else get_smatrix_from_xy_lists
        xvals, yvals = utils.non_none([xydicts, xylists])
        smatrix, xbins, ybins, n_skipped = sfcn(xvals, yvals, kfcn=kfcn, n_max_bins=n_max_bins)
    else:
        n_skipped = 0
        if xybins is None:  # you probably want to pass in the bins, but if not they just get set to the 0-based row/column index
            xbins, ybins = list(range(len(smatrix[0]))), list(range(len(smatrix)))
        else:
            xbins, ybins = xybins
    if tdbg:
        def vstr(v):
            if v is None or numpy.isnan(v): return ''
            return ('%.2f' if float_vals else '%d') % v
        lb = str(max(len(str(b)) for b in ybins + xbins))  # max length (when converted to str) of any bin label
        print('  detailed smatrix')
        print('        %s' % ''.join((('  %'+lb+'s')%ib for ib in [''] + ybins)))
        for iff, fb in enumerate(xbins):
            print('       ', end=' ')
            for ii, ib in enumerate(ybins):
                print(('%s %'+lb+'s') % ((('%'+lb+'s')%fb) if ii==0 else '', vstr(smatrix[ii][iff])), end=' ')
            print('')

    fig, ax = mpl_init()
    cmap = plt.cm.get_cmap('viridis') #Blues  #cm.get_cmap('jet')
    cmap.set_under('w')
    smtx_min = min([v for r in smatrix for v in r])
    vmin = min(0., smtx_min) if float_vals else 0.5
    # smatrix = [[utils.non_none([v, vmin]) for v in vl] for vl in smatrix]  # we *want* the Nones, since that's what makes them blank (rather than all freaking purple)
    if float_vals and any(v is not None and v < vmin for vl in smatrix for v in vl):  # would be easy to fiddle with this but i don't want to right now, and I'm only plotting positive values atm
        raise Exception('value(s) %s less than vmin %.2f in plot_smatrix()' % ([v for vl in smatrix for v in vl if v < vmin], vmin))
    heatmap = ax.pcolormesh(numpy.ma.array(smatrix, mask=numpy.isnan(smatrix)), cmap='viridis', vmin=vmin) #, vmax=1.) #, norm=mpl.colors.LogNorm()
    cbar = plt.colorbar(heatmap, label=('%s (skipped %d)'%(blabel, n_skipped)) if n_skipped > 0 else blabel, pad=0.12)
    def ltsize(n): return 15 if n < 15 else 8
    return mpl_finish(ax, plotdir, plotname, title=title, xlabel=xlabel, ylabel=ylabel,
                      xticks=[n - 0.5 for n in range(1, len(xbins) + 1)], yticks=[n - 0.5 for n in range(1, len(ybins) + 1)],
                      xticklabels=[lfcn(b) for b in xbins], yticklabels=[y_lfcn(b) for b in ybins],
                      xticklabelsize=ltsize(len(xbins)), yticklabelsize=ltsize(len(ybins)))

# ----------------------------------------------------------------------------------------
# NOTE pt 'header' to line to add extra headers below the main title
def make_html(plotdir, n_columns=3, extension='svg', fnames=None, title='', bgcolor='000000', new_table_each_row=False, htmlfname=None, extra_links=None, subdname=None):
    if fnames is not None:  # make sure it's formatted properly
        for rowfnames in fnames:
            if not isinstance(rowfnames, list):
                raise Exception('each entry in fnames should be a list of strings, but got a %s: %s' % (type(rowfnames), rowfnames))
            for fn in rowfnames:
                if not isinstance(fn, six.string_types):
                    raise Exception('each entry in each row should be a string (file name), but got a %s: %s' % (type(fn), fn))
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
             '<body bgcolor="%s">' % bgcolor,
             '<h3 style="text-align:left; color:DD6600;">' + title + '</h3>',
             extra_link_str,
             '<table>',
             '<tr>']

    def add_newline(lines, header=None):
        if new_table_each_row:
            endlines, startlines = ['</tr>', '</table>'], ['<table>', '<tr>']
        else:
            endlines, startlines = ['</tr>'], ['<tr>']
        lines += endlines
        if header is not None:
            lines += ['<h3 style="text-align:left; color:DD6600;">' + header + '</h3>']
        lines += startlines
    def add_fname(lines, fullfname):  # NOTE <fullname> may, or may not, be a base name (i.e. it might have a subdir tacked on the left side)
        fname = fullfname.replace(plotdir, '').lstrip('/')
        if htmlfname is None or subdname is not None:
            fname = '%s/%s' % (utils.non_none([subdname, dirname]) , fname)  # if specifying htmlfname, the calling function needs to either add the subdirs/<dirname> to each fname, or set subdname)
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

        # and group insertion lengths and contents together
        found_bound_fnames = [fn for fn in fnamelist if any(b + '_insertion' == utils.getprefix(os.path.basename(fn)) for b in utils.all_boundaries)]
        if len(found_bound_fnames) == len(utils.all_boundaries):  # if we don't have all of them, we can't make a complete row, so just let them get grouped below
            fnames.append(found_bound_fnames)
            for fn in found_bound_fnames:
                fnamelist.remove(fn)
        found_content_fnames = [fn for fn in fnamelist if '_content' in fn]
        if len(found_content_fnames) > 0:
            fnames.append(found_content_fnames)
            for fn in found_content_fnames:
                fnamelist.remove(fn)

        # then do the rest in groups of <n_columns>
        while len(fnamelist) > 0:
            fnames.append(fnamelist[:n_columns])
            fnamelist = fnamelist[n_columns:]

    # write the meat of the html
    for rowlist in fnames:
        if 'header' in rowlist:
            if len(rowlist) != 2:
                raise Exception('malformed header row list in fnames (should be len 2 but got %d): %s' % (len(rowlist), rowlist))
            add_newline(lines, header=rowlist[1])
            continue
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
    if 'partis.fraction_uncertainty' not in sys.modules:
        from . import fraction_uncertainty

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

    tmphilos = [sys.modules['partis.fraction_uncertainty'].err(r, r + w) for r, w in zip(right, wrong)]
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
    if 'partis.fraction_uncertainty' not in sys.modules:
        from . import fraction_uncertainty
    fraction_uncertainty = sys.modules['partis.fraction_uncertainty']

    def get_single_vals(pv):
        yvals = [float(c) / t for c, t in zip(pv['ycounts'], pv['ytotals'])]  # total shouldn't be able to be zero
        tmphilos = [fraction_uncertainty.err(c, t) for c, t in zip(pv['ycounts'], pv['ytotals'])]
        yerrs = [err[1] - err[0] for err in tmphilos]
        print('  %s                    %s' % (xlabel, ylabel))
        for iv in range(len(pv['xvals'])):
            print('   %8.0f     %5.0f / %-5.0f  =  %5.2f   +/-  %.3f' % (pv['xvals'][iv], pv['ycounts'][iv], pv['ytotals'][iv], yvals[iv], yerrs[iv]))
        return pv['xvals'], yvals, yerrs

    fig, ax = mpl_init()
    mpl.rcParams.update({'legend.fontsize' : 15})

    xmin, xmax, xticks = None, None, None
    for ii in range(len(labels)):
        print(labels[ii])
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
def plot_laplacian_spectra(plotdir, plotname, eigenvalues, title):
    hist = Hist(30, min(eigenvalues), max(eigenvalues), value_list=eigenvalues)
    fig, ax = mpl_init()
    hist.mpl_plot(ax)
    mpl_finish
    mpl_finish(ax, plotdir, plotname, xlabel='eigenvalues', ylabel='count', title=title)

# ----------------------------------------------------------------------------------------
# if <high_x_val> is set, clusters with median x above <high_x_val> get skipped by default and returned, the idea being that you call this fcn again at the end with <plot_high_x> set just on the thereby-returned high-x clusters
def make_single_joyplot(sorted_clusters, annotations, repertoire_size, plotdir, plotname, x1key='n_mutations', x1label='N mutations', x2key=None, x2label=None, high_x_val=None, plot_high_x=False,
                        cluster_indices=None, title=None, queries_to_include=None, meta_info_to_emphasize=None, meta_info_key_to_color=None, meta_emph_formats=None, all_emph_vals=None, emph_colors=None, global_max_vals=None,
                        make_legend=False, remove_none_vals=False, sortlabel='?', add_clone_id=False, add_index=False, dont_label_queries_to_include=False, debug=False):
    from . import lbplotting
    smetrics = treeutils.affy_metrics + treeutils.daffy_metrics  # treeutils.lb_metrics.keys() + treeutils.dtr_metrics
    # NOTE <xvals> must be sorted
    # ----------------------------------------------------------------------------------------
    def offcolor(offset):
        offcolors = {'up' : '#386cc2', 'down' : '#bf1328'}  # greyish blue, mild red
        return offcolors[offset]
    # ----------------------------------------------------------------------------------------
    def bexpand(bpair, fuzz=0.02):
        if bpair[0] == bpair[1]:
            return expand_bounds(bpair)  # ok this sucks, i should really use this for all cases here, but i don't want to run through all possibilities and all calls right now, so just using it for the edge case where things are breaking right now
        diff = max(fuzz, bpair[-1] - bpair[0])  # if max and min are the same, use <fuzz> for <diff>
        return [bpair[0] - fuzz * diff, bpair[1] + fuzz * diff]
    # ----------------------------------------------------------------------------------------
    def get_xval_list(cluster, xkey):  # NOTE this *has* to return values in the same order they're in line['unique_ids']
        line = annotations[':'.join(cluster)]
        if xkey in line:
            return line[xkey]
        else:
            return treeutils.smvals(line, xkey)
    # ----------------------------------------------------------------------------------------
    def get_xval_dict(uids, xkey):
        line = annotations[':'.join(cluster)]
        if xkey in smetrics:
            return {u : line['tree-info']['lb'][xkey][u] for u in uids}
        else:
            return {u : utils.per_seq_val(line, xkey, u) for u in uids}
    # ----------------------------------------------------------------------------------------
    def getbounds(xkey):
        all_xvals = [x for c in sorted_clusters for x in get_xval_list(c, xkey) if x is not None]  # NOTE can't ignore/skip None vals in the list/dict getter fcn above, since order has to match line['unique_ids']
        if len(all_xvals) == 0:
            print('    %s no (non-None) xvals for %s in single joyplot fcn' % (utils.wrnstr(), xkey))
            return None
        bounds = [f(all_xvals) for f in [min, max]]
        if bounds[0] == bounds[1]:  # if min and max are the same (i.e. all vals are the same), just use the value +/- 10%
            bounds = expand_bounds(bounds)  # NOTE old version doesn't work e.g. for negative vals: [(1 + s * 0.1) * bounds[0] for s in [-1, 1]]
        if global_max_vals is not None and xkey in global_max_vals:
            bounds[1] = global_max_vals[xkey]
        return bounds
    # ----------------------------------------------------------------------------------------
    def uselog(xkey):  # the low end (zero bin) of these distributions always dominates, but we're actually interested in the upper tail, so logify it
        return xkey in smetrics or xkey == 'affinities'
    # ----------------------------------------------------------------------------------------
    def add_hist(xkey, sorted_xvals, yval, iclust, cluster, median_x1, fixed_x1max, base_alpha, repfracstr, offset=None):
        if None in sorted_xvals:
            if remove_none_vals:
                sorted_xvals = [v for v in sorted_xvals if v is not None]
            else:
                raise Exception('None type value[s] for %s: %s' % (xkey, sorted_xvals))
        qti_x_vals = {}
        tqtis = []  # queries to emphasize in this cluster, as pairs of (uid, label)
        if queries_to_include is not None and not dont_label_queries_to_include:
            tqtis += [(u, u) for u in set(cluster) & set(queries_to_include)]
        if meta_info_to_emphasize is not None or meta_info_key_to_color is not None:
            antn = annotations[':'.join(cluster)]
            if meta_info_to_emphasize is not None and meta_emph_key in antn:
                def eqfcn(v): return utils.meta_info_equal(meta_emph_key, meta_emph_val, v, formats=meta_emph_formats)
                estr = utils.meta_emph_str(meta_emph_key, meta_emph_val, formats=meta_emph_formats)
                emphids = [u for u, v in zip(cluster, antn[meta_emph_key]) if eqfcn(v)]
                tqtis += [(u, estr) for u in emphids]
        if len(tqtis) > 0:
            qti_x_vals = get_xval_dict([u for u, _ in tqtis], xkey)  # add a red line for each of 'em (i.e. color that hist bin red)
            if any(v > fixed_x1max for v in qti_x_vals.values()):
                fixed_x1max = 1.05 + max(qti_x_vals.values())
            if plot_high_x:
                xfac = 1.1
            elif float(median_x1) / fixed_x1max < 0.5:  # if seqs are mostly on the left, put text on right
                xfac = min(0.75, max(sorted_xvals) / float(fixed_x1max) + 0.1)
            else:  # vice versa
                xfac = 0.1
            qtistrs = [l for u, l in sorted(tqtis, key=lambda x: qti_x_vals[x[0]])]  # sort by x value, then label with value from tqtis
            if any(qtistrs.count(s)>1 for s in set(qtistrs)):
                for qstr in [s for s in set(qtistrs) if qtistrs.count(s)>1]:
                    qtistrs = [s for s in qtistrs if s!=qstr] + ['%s(x%d)' % (qstr, qtistrs.count(qstr))]
            ax.text(xfac * fixed_x1max, yval, ', '.join(qtistrs), color='red', fontsize=8)

        if debug:
            fstr = '6.1f' if xkey == 'n_mutations' else '6.4f'
            print(('     %5s  %-10s  %4.1f  %'+fstr+'  %'+fstr) % ('%d' % csize if iclust == 0 else '', repfracstr if iclust == 0 else '', yval, numpy.median(sorted_xvals), numpy.mean(sorted_xvals)), end=' ')
            if len(tqtis) > 0:
                print('   ' + ' '.join(qtistrs), end=' ')
            print('')

        if xkey == 'n_mutations':
            nbins = sorted_xvals[-1] - sorted_xvals[0] + 1
            hist = Hist(nbins, sorted_xvals[0] - 0.5, sorted_xvals[-1] + 0.5)
        else:
            nbins = 30 if xkey in smetrics else 15
            hist = Hist(nbins, *bexpand(xbounds[xkey], fuzz=0.01))
        hist.list_fill(sorted_xvals)
        if uselog(xkey):
            hist.logify(0.3)  # the factor is kind of arbitrary, but we need to set the scale for the smallest bin contents (in this case 1)
        assert hist.overflow_contents() == 0.  # includes underflows
        max_contents = max(hist.bin_contents)
        for ibin in range(1, hist.n_bins + 1):
            barheight = 0. if max_contents==0 else utils.intexterpolate(0., min_bar_height, max_contents, max_bar_height, hist.bin_contents[ibin])
            if meta_info_key_to_color is not None:
                bin_ids = [u for u, x in zip(antn['unique_ids'], get_xval_list(cluster, xkey)) if hist.find_bin(x)==ibin]  # uids in this bin
                def psfcn(u): return utils.meta_emph_str(meta_info_key_to_color, utils.per_seq_val(antn, meta_info_key_to_color, u, use_default=True), formats=meta_emph_formats)
                me_vals = [psfcn(u) for u in bin_ids]  # meta info values for the uids in this bin
                me_color_fracs = [(c, me_vals.count(v) / float(len(me_vals))) for v, c in emph_colors if v in me_vals]
            bin_color = base_color
            if offset is None:  # default: bar extends equally above + below center
                y_lower, y_upper = yval - barheight/2., yval + barheight/2.
            else:  # this bar is only up or down (and presumably a different bar is being drawn the other direction)
                y_lower, y_upper = yval, yval
                bin_color = offcolor(offset)
                if offset == 'up':
                    y_upper += barheight / 2.
                elif offset == 'down':
                    y_lower -= barheight / 2.
                else:
                    assert False
            alpha = base_alpha
            # alpha = utils.intexterpolate(0, min_alpha, max_contents, max_alpha, hist.bin_contents[ibin])
            if hist.bin_contents[ibin] == 0.:
                bin_color = 'grey'
                alpha = 0.4
            xlo, xhi = hist.low_edges[ibin], hist.low_edges[ibin+1]
            if xkey == x2key:  # if it's the second one, we need to rescale the x vals to correspond to the existing x1key x axis
                xlo, xhi = [utils.intexterpolate(xbounds[x2key][0], xbounds[x1key][0], xbounds[x2key][1], xbounds[x1key][1], x) for x in [xlo, xhi]]
            if meta_info_key_to_color is None or hist.bin_contents[ibin] == 0.:  # normal/default: one bin color
                ax.fill_between([xlo, xhi], [y_lower, y_lower], [y_upper, y_upper], color=bin_color, alpha=alpha)
            else:  # color different fractions of the bar according to input meta info
                t_y_lo = y_lower
                for tcol, tfrac in me_color_fracs:
                    t_y_hi = t_y_lo + tfrac * (y_upper - y_lower)
                    ax.fill_between([xlo, xhi], [t_y_lo, t_y_lo], [t_y_hi, t_y_hi], color=tcol, alpha=alpha)
                    t_y_lo = t_y_hi
            if any(hist.find_bin(x)==ibin for x in qti_x_vals.values()):
                xmid, delta_x = xlo + 0.5 * (xhi - xlo), (xhi - xlo) / 12.
                ymin, ymax = yval - max_bar_height/4., yval + max_bar_height/4.
                ax.fill_between([xmid - delta_x, xmid + delta_x], [ymin, ymin], [ymax, ymax], color='red', alpha=alpha)
        return fixed_x1max  # ick ick ick

    # ----------------------------------------------------------------------------------------
    alt_colors = ['#006600', '#3399ff', '#ffa500']
    # goldenrod '#daa520'
    # red '#cc0000',
    # dark red '#990012'
    # purple '#a821c7'
    # grey '#808080'

    dpi = 80
    xpixels = 550 #450
    min_ypixels = 500 #400
    total_delta_y = len(sorted_clusters)
    y_bar_pixels = 12 if x2key is None else 25
    if meta_info_key_to_color is not None:
        y_bar_pixels = 20
    min_bar_height, max_bar_height = 0.3 / min_ypixels * total_delta_y, float(y_bar_pixels) / min_ypixels * total_delta_y
    ypixels = max(min_ypixels, y_bar_pixels * total_delta_y)
    fig, ax = mpl_init(figsize=(xpixels // dpi, ypixels // dpi))
    # min_alpha, max_alpha = 0.1, 1.
    base_alpha = 0.55

    ymin, ymax = 9999, 0
    iclust_global = 0  # index within this plot
    yticks, yticklabels = [], []

    high_x_clusters = []
    xbounds = {x1key : getbounds(x1key)}  # these are the smallest/largest x values in any of <sorted_clusters>, whereas <high_x_val> is a fixed calling-fcn-specified value that may be more or less (kind of wasteful to get all the x vals here and then also in the main loop)
    if x2key is not None:
        xbounds[x2key] = getbounds(x2key)
    if any(xbounds[xk] is None for xk in xbounds):
        print('    %s None type xbounds in single joyplot: %s' % (utils.wrnstr(), xbounds))
        return 'no values' if high_x_val is None else high_x_clusters  # 'no values' isn't really a file name, it just shows up as a dead link in the html
    fixed_xmax = high_x_val if high_x_val is not None else xbounds[x1key][1]  # xmax to use for the plotting (ok now there's three max x values, this is getting confusing)
    if meta_info_to_emphasize is not None:
        meta_emph_key, meta_emph_val = list(meta_info_to_emphasize.items())[0]
        if all(meta_emph_key not in l for l in annotations.values()):
            print('  %s emphasis key \'%s\' not found in any of %d annotations' % (utils.color('yellow', 'warning'), meta_emph_key, len(annotations)))

    if debug:
        print('  %s   %d x %d   %s' % (plotname, xpixels, ypixels, utils.color('red', 'high %s'%x1key) if plot_high_x else ''))
        print('      size   frac      yval    median   mean')

    if x2key is None:
        cgroup_iter = itertools.groupby(sorted_clusters, key=lambda c: len(c))  # this doesn't re-sort anything, it just creates groups by size (like |sort|uniq)
    else:
        cgroup_iter = [(len(c), [c]) for c in sorted_clusters]  # creates a structure similar to the previous clause, but with just trivial groups (one for each cluster), since in this case the clusters aren't sorted by size (and are instead sorted with continuous-valued variables) so we don't need/want the groupby stuff to a get decent y axis
    for csize, cluster_group in cgroup_iter:
        if x2key is None:
            cluster_group = sorted(list(cluster_group), key=lambda c: numpy.median(get_xval_list(c, x1key)))  # sort ties in the default sorting by median <x1key> (has no effect if x2key is set since the groups are always length 1)
        repfracstr = utils.get_repfracstr(csize, repertoire_size)
        for iclust, cluster in enumerate(cluster_group):
            x1vals = sorted([x for x in get_xval_list(cluster, x1key) if x is not None])
            median_x1 = numpy.median([v for v in x1vals if v is not None])  # maybe should use mean instead of median?

            if high_x_val is not None and median_x1 > high_x_val and not plot_high_x:  # if <high_x_val> is not set, we don't skip any clusters
                high_x_clusters.append(cluster)
                continue

            yval = len(sorted_clusters) - iclust_global
            if yval < ymin:
                ymin = yval
            if yval > ymax:
                ymax = yval
            yticks.append(yval)
            yticklabels.append(repfracstr if x2key is None else '%d'%csize)

            base_color = alt_colors[iclust_global % len(alt_colors)] if meta_info_key_to_color is None else 'black'

            fixed_xmax = add_hist(x1key, x1vals, yval, iclust, cluster, median_x1, fixed_xmax, base_alpha, repfracstr, offset=None if x2key is None else 'up')
            if x2key is not None:
                tmpval = lbplotting.mean_of_top_quintile([v for v in x1vals if v is not None])  # NOTE presumably this needs to match sortlabel
                ax.plot([tmpval, tmpval], [yval, yval + 1./4], linewidth=2.5, alpha=0.55, color='green')
                x2vals = sorted(get_xval_list(cluster, x2key))
                fixed_xmax = add_hist(x2key, x2vals, yval, iclust, cluster, median_x1, fixed_xmax, base_alpha, repfracstr, offset='down')

            if cluster_indices is not None:  # add the (global) cluster index (i.e. 1 - rank) and cluster size as text on the right side of the plot
                xtext = x1vals[-1] if plot_high_x else fixed_xmax  # would be nice to use fig coords and fig.text(), but we want these to be equal to the y values of the slugs
                xtext -= 3
                xwidth = ax.get_xlim()[1] - ax.get_xlim()[0] if plot_high_x else fixed_xmax
                if iclust_global == 0:
                    if add_clone_id:
                        ax.text(0.02 * xwidth + xtext - 7, yval + 0.55, 'clone id', color=base_color, fontsize=6, alpha=base_alpha, fontdict={'weight' : 'bold'})
                    if add_index:
                        ax.text(0.05 * xwidth + xtext - 3, yval + 0.55, 'index', color=base_color, fontsize=6, alpha=base_alpha, fontdict={'weight' : 'bold'})
                    ax.text(0.12 * xwidth + xtext - 2, yval + 0.55, 'size', color=base_color, fontsize=6, alpha=base_alpha, fontdict={'weight' : 'bold'})
                if add_clone_id:
                    ax.text(0.02 * xwidth + xtext - 7, yval, utils.get_clone_id(cluster), color=base_color, fontsize=6, alpha=base_alpha, fontdict={'weight' : 'bold'})
                if add_index:
                    ax.text(0.05 * xwidth + xtext, yval, str(cluster_indices[':'.join(cluster)]), color=base_color, fontsize=6, alpha=base_alpha, fontdict={'weight' : 'bold'})
                ax.text(0.10 * xwidth + xtext, yval, str(csize), color=base_color, fontsize=6, alpha=base_alpha, fontdict={'weight' : 'bold'})

            iclust_global += 1

    fsize = 12
    if x2key is not None:
        fig.text(0.8, 0.25, x1label if not uselog(x1key) else '%s (log)'%x1label, color=offcolor('up'), alpha=base_alpha, fontdict={'weight' : 'bold'}, fontsize=fsize)
        fig.text(0.8, 0.215, x2label if not uselog(x2key) else '%s (log)'%x2label, color=offcolor('down'), alpha=base_alpha, fontdict={'weight' : 'bold'}, fontsize=fsize)

    plot_x_bounds = [high_x_val, xbounds[x1key][1]] if plot_high_x else bexpand((xbounds[x1key][0], fixed_xmax))
    n_x_ticks, xlabel, xticks, xticklabels = 4, x1label, None, None
    if x2key is not None:
        xlabel = x2label
        xticks = [x for x in numpy.arange(xbounds[x1key][0], xbounds[x1key][1], (xbounds[x1key][1] - xbounds[x1key][0]) / float(n_x_ticks-1))] + [xbounds[x1key][1]]
        xticklabels = ['%.1f' % utils.intexterpolate(xbounds[x1key][0], xbounds[x2key][0], xbounds[x1key][1], xbounds[x2key][1], x) for x in xticks]  # translate x1 tick positions to x2 tick labels
        fig.text(0.13, 0.07, '%.3f'%xbounds[x1key][0], color=offcolor('up'), alpha=base_alpha, fontdict={'weight' : 'bold'}, fontsize=fsize)
        fig.text(0.89, 0.07, '%.3f'%xbounds[x1key][1], color=offcolor('up'), alpha=base_alpha, fontdict={'weight' : 'bold'}, fontsize=fsize)
        fig.text(0.52, 0.03, x1label, color=offcolor('up'), alpha=base_alpha, fontdict={'weight' : 'bold'}, fontsize=fsize)
        fig.text(0.05, 0.9, 'sorted by\n%s'%sortlabel, color=offcolor('up'), alpha=base_alpha, fontdict={'weight' : 'bold'}, fontsize=fsize)
    n_y_ticks = 5
    if x2key is None and len(yticks) > n_y_ticks:
        yticks = [yticks[i] for i in range(0, len(yticks), int(len(yticks) / float(n_y_ticks - 1)))]
        yticklabels = [yticklabels[i] for i in range(0, len(yticklabels), int(len(yticklabels) / float(n_y_ticks - 1)))]
    fn = mpl_finish(ax, plotdir, plotname, xlabel=xlabel, ylabel=('family size (frac. of %d)' % repertoire_size) if x2key is None else 'clonal family size', title=title,
                    xbounds=plot_x_bounds, ybounds=bexpand((ymin, ymax), fuzz=0.03 if x2key is None else 0.07), xticks=xticks, xticklabels=xticklabels, yticks=yticks, yticklabels=yticklabels, yticklabelsize=11, adjust={'left' : 0.2, 'right' : 0.85})

    if meta_info_key_to_color is not None and make_legend:
        make_meta_info_legend(plotdir, plotname, meta_info_key_to_color, emph_colors, all_emph_vals, meta_emph_formats=meta_emph_formats, alpha=base_alpha)

    if high_x_val is None:
        return fn
    else:
        return high_x_clusters

# ----------------------------------------------------------------------------------------
def bubble_plot(plotname, plotdir, bubfos, title='', xtra_text=None, alpha=0.4):
    # ----------------------------------------------------------------------------------------
    def check_missing_extra():  # this shouldn't happen any more now that i'm using strings for 'id' key, but eh can't hurt
        iids, fids = set(d['id'] for d in bubfos), set(d['id'] for d in bubble_positions)
        if len(fids - iids) > 0:
            print('    %s extra bubble ids read from file: %s' % (utils.wrnstr(), ' '.join(fids - iids)))
        if len(iids - fids) > 0:
            print('    %s missing bubble ids from file: %s' % (utils.wrnstr(), ' '.join(iids - fids)))
        assert len(bubble_positions) == len(bubfos)
    # ----------------------------------------------------------------------------------------
    rfn = '%s/csize-radii.csv' % plotdir
    bpfn = '%s/bubble-positions.csv' % plotdir
    workfnames = [rfn, bpfn]
    with open(rfn, utils.csv_wmode()) as rfile:
        writer = csv.DictWriter(rfile, ['id', 'radius'])
        writer.writeheader()
        for bfo in bubfos:
            writer.writerow({k : bfo[k] for k in ['id', 'radius']})
    cmd = '%s/bin/circle-plots.py %s %s' % (utils.get_partis_dir(), rfn, bpfn)
    utils.simplerun(cmd, extra_str='        ')
    with open(bpfn) as bpfile:
        def cfn(k, v): return v if k=='id' else float(v)
        reader = csv.DictReader(bpfile)
        bubble_positions = []  # they come back out of order, so need to sort
        for line in reader:
            bubble_positions.append({k : cfn(k, v) for k, v in line.items()})
        check_missing_extra()
        bdict = {b['id'] : b for b in bubfos}
        for posfo in bubble_positions:
            bdict[posfo['id']].update(posfo)
    fig, ax = mpl_init()
    if len(bubfos) == 0:
        print('      %s no bubble positions, returning' % utils.wrnstr())
        return None
    lim = max(max(abs(bfo['x']) + bfo['r'], abs(bfo['y']) + bfo['r']) for bfo in bubfos)
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    ax.axis('off')
    plt.gca().set_aspect('equal')

    for bfo in bubfos:
        if 'texts' in bfo:
            for tfo in bfo['texts']:
                ax.text(bfo['x']+tfo.get('dx', 0), bfo['y']+tfo.get('dy', 0), tfo['tstr'], fontsize=tfo['fsize'], alpha=1, color=tfo['tcol'])
        if bfo['fracs'] is None or len(bfo['fracs']) == 0:
            ax.add_patch(plt.Circle((bfo['x'], bfo['y']), bfo['r'], alpha=alpha, linewidth=2, fill=True))  # plain circle
        else:
            plot_pie_chart_marker(ax, bfo['x'], bfo['y'], bfo['r'], bfo['fracs'], alpha=alpha)
    if xtra_text is not None:
        fig.text(xtra_text['x'], xtra_text['y'], xtra_text['text'], fontsize=xtra_text.get('fontsize', 8), color=xtra_text.get('color', 'black'))
    mpl_finish(ax, plotdir, plotname, title=title)

    for wfn in workfnames:
        os.remove(wfn)
    return plotname

# ----------------------------------------------------------------------------------------
def plot_vrc01_class_muts(plotdir, plotname, annotation, mekey=None, formats=None, only_csv=False, debug=False):
    # ----------------------------------------------------------------------------------------
    def get_vmut_counts(infos, tdbg=False):
        obs_muts = vrc01.glvrc01_mutation_set(infos, debug=tdbg)
        vclass_muts = vrc01.vrc01_class_mutation_set(debug=tdbg)
        obs_vclass_muts = [m for m in obs_muts if m in vclass_muts]
        vmuts, vtot = len(obs_vclass_muts), len(obs_muts)
        return vmuts, vtot
    # ----------------------------------------------------------------------------------------
    if only_csv:
        print('  only_csv not implemented for vrc01 mut plots, returning')
        return [[]]
    import partis.vrc01 as vrc01
    fnames = [[]]

    vrc01.check_naive_seq(annotation)

    utils.add_seqs_aa(annotation)
    mutvals, totvals = {}, {}
    for meval in set(annotation[mekey]):
        infos = [{'name' : u, 'seq' : s} for u, s, m in zip(annotation['unique_ids'], annotation['seqs_aa'], annotation[mekey]) if m==meval]
        vmuts, vtot = get_vmut_counts(infos)
        mutvals[meval], totvals[meval] = vmuts, vtot

    all_emph_vals, emph_colors = meta_emph_init(mekey, formats=formats, all_emph_vals=set(annotation[mekey]))
    all_emph_vals = sorted(all_emph_vals, key=str)

    n_bins = len(all_emph_vals)
    mcolors = {v : c for v, c in emph_colors}
    hlist, colors = [], []
    for ibin, meval in enumerate(all_emph_vals, 1):
        hist = Hist(n_bins, -0.5, n_bins - 0.5)
        hist.bin_labels = [''] + all_emph_vals + ['']
        hist.set_ibin(ibin, 0 if totvals[meval]==0 else mutvals[meval] / totvals[meval], 0, label=meval)  # it's already labeled, but this make me feel warmer + fuzzier
        hlist.append(hist)
        colors.append(mcolors[meval])

    draw_no_root(None, plotdir=plotdir, plotname=plotname, more_hists=hlist, square_bins=True, linewidths=[3 for _ in hlist], alphas=[0.7 for _ in hlist], colors=colors, plottitle='vrc01 mut frac', ybounds=[0, 1], ytitle='fraction of mutations')

    fnames[-1].append(plotname)
    return fnames
