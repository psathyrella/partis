from __future__ import unicode_literals
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')

import math
from scipy.interpolate import interp1d
import os
import glob
import sys
import stat
import copy
import csv
import numpy
from scipy import stats
from array import array
from subprocess import check_call
import re
from collections import OrderedDict

import utils
import fraction_uncertainty
import plotconfig
from hist import Hist

from opener import opener

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
    with opener('w')(fname) as histfile:
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
    hist = Hist(2, -0.5, 1.5)

    def set_bin(numer, denom, ibin, label):
        frac = float(numer) / denom
        bounds = fraction_uncertainty.err(numer, denom, use_beta=True)
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
def draw_no_root(hist, log='', plotdir=None, plotname='foop', more_hists=None, scale_errors=None, normalize=False, bounds=None,
                 figsize=None, shift_overflows=False, colors=None, errors=False, write_csv=False, xline=None, yline=None, linestyles=None,
                 linewidths=None, plottitle=None, csv_fname=None, stats='', translegend=(0., 0.), rebin=None,
                 xtitle=None, ytitle=None, markersizes=None, no_labels=False, only_csv=False):
    assert os.path.exists(plotdir)

    fig, ax = mpl_init(figsize=figsize)
    mpl.rcParams.update({'legend.fontsize' : 15})

    hists = [hist,]
    if more_hists is not None:
        hists = hists + more_hists

    xmin, xmax, ymax = None, None, None
    ih = 0
    for htmp in hists:
        if scale_errors is not None:
            factor = float(scale_errors[0]) if len(scale_errors) == 1 else float(scale_errors[ih])
            for ibin in range(htmp.n_bins+2):
                htmp.errors[ibin] *= factor
        if normalize:  # NOTE removed <normalization_bounds> option, hopefully I'm not using it any more
            htmp.normalize()
        if ymax is None or htmp.get_maximum(xbounds=bounds) > ymax:
            ymax = htmp.get_maximum(xbounds=bounds)
        if xmin is None or htmp.xmin < xmin:  # overridden by <bounds> below
            xmin = htmp.xmin
        if xmax is None or htmp.xmax > xmax:
            xmax = htmp.xmax

        ih += 1

    if bounds is not None:
        xmin, xmax = bounds

    assert not shift_overflows  # TODO implement this inside of Hist
    # if shift_overflows:
    #     for htmp in hists:
    #         if htmp == None:
    #             continue
    #         underflows, overflows = 0.0, 0.0
    #         first_shown_bin, last_shown_bin = -1, -1
    #         for ib in range(0, htmp.GetXaxis().GetNbins()+2):
    #             if htmp.GetXaxis().GetBinCenter(ib) <= xmin:
    #                 underflows += htmp.GetBinContent(ib)
    #                 htmp.SetBinContent(ib, 0.0)
    #             elif first_shown_bin == -1:
    #                 first_shown_bin = ib
    #             else:
    #                 break
    #         for ib in reversed(range(0, htmp.GetXaxis().GetNbins()+2)):
    #             if htmp.GetXaxis().GetBinCenter(ib) >= xmax:
    #                 overflows += htmp.GetBinContent(ib)
    #                 htmp.SetBinContent(ib, 0.0)
    #             elif last_shown_bin == -1:
    #                 last_shown_bin = ib
    #             else:
    #                 break

    #         if 'd_hamming' in plotname:
    #             print htmp.GetTitle()
    #             print '  underflow', underflows, htmp.GetBinContent(first_shown_bin)
    #             print '  overflow', overflows, htmp.GetBinContent(last_shown_bin)
    #             print '  first', htmp.GetXaxis().GetBinCenter(first_shown_bin)
    #             print '  last', htmp.GetXaxis().GetBinCenter(last_shown_bin)
    #         htmp.SetBinContent(first_shown_bin, underflows + htmp.GetBinContent(first_shown_bin))
    #         htmp.SetBinContent(last_shown_bin, overflows + htmp.GetBinContent(last_shown_bin))

    if colors is None:  # fiddle here http://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
        assert len(hists) < 5
        colors = ('royalblue', 'darkred', 'green', 'darkorange')
    else:
        assert len(hists) <= len(colors)
    if linestyles is None:
        linestyles = ['-' for _ in range(len(hists))]
    else:
        assert len(hists) <= len(linestyles)

    for ih in range(len(hists)):
        htmp = hists[ih]
        # if 'rms' in stats:
        #     htmp.SetTitle(htmp.GetTitle() + (' (%.2f)' % htmp.GetRMS()))
        if 'mean' in stats:
            htmp.title += ' (%.2f)' % htmp.get_mean()
        if '0-bin' in stats:
            htmp.title += ' (%.2f)' % htmp.bin_contents[1]
        markersize = None
        if markersizes is not None:
            imark = ih if len(markersizes) > 1 else 0
            markersize = markersizes[imark]
        linewidth = None
        if linewidths is None:
            if ih < 6 and len(hists) > 1:
                linewidth = 6-ih
        else:
            ilw = ih if len(linewidths) > 1 else 0
            linewidth = linewidths[ilw]
        if no_labels:
            htmp.bin_labels = ['' for _ in htmp.bin_labels]
        if rebin is not None:
            htmp.rebin(rebin)
        htmp.mpl_plot(ax, color=colors[ih], linewidth=linewidth, linestyle=linestyles[ih], ignore_overflows=True, errors=errors)

    if xline is not None:
        ax.plot([xline, xline], [-0.1*ymax, 0.5*ymax], color='black', linestyle='--', linewidth=3)
    if yline is not None:
        print 'TODO fix y line'
    # if yline is not None:
    #     # if yline < hframe.GetYaxis().GetXmin() or xline > hframe.GetYaxis().GetXmax():  # make sure we got valid a x position for the line
    #     #     print 'WARNING plotting y line at %f out of bounds (%f, %f)' % (float(ymin), hframe.GetYaxis().GetXmin(), hframe.GetYaxis().GetXmax())
    #     yl = TLine(hframe.GetXaxis().GetXmin(), yline, hframe.GetXaxis().GetXmax(), yline)
    #     yl.Draw()

    xticks, xticklabels = None, None
    if hist.bin_labels.count('') != len(hist.bin_labels):
        xticks = hist.get_bin_centers()
        xticklabels = hist.bin_labels

    if not os.path.exists(plotdir + '/plots'):
        raise Exception('ERROR dir \'' + plotdir + '/plots\' d.n.e.')

    if not only_csv:
        mpl_finish(ax, plotdir, plotname,
                   title=plotname if plottitle is None else plottitle,
                   xlabel=hist.xtitle if xtitle is None else xtitle,
                   ylabel=hist.ytitle if ytitle is None else ytitle,
                   xbounds=[xmin, xmax],
                   ybounds=[-0.03*ymax, 1.15*ymax],
                   leg_loc=(0.72 + translegend[0], 0.7 + translegend[1]),
                   log=log, xticks=xticks, xticklabels=xticklabels)

    if write_csv:
        assert more_hists is None
        if csv_fname is None:
            hist.write(plotdir + '/plots/' + plotname + '.csv')
        else:
            hist.write(csv_fname)

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
def get_mean_info(hists):
    raise Exception('needs to be converted off of root')
    # assert len(hists) > 0
    # means, sems, normalized_means = [], [], []
    # sum_total, total_entries = 0.0, 0.0
    # bin_values = {}  # map from bin centers to list (over hists) of entries (corresponds to <ibin> in unihist)
    # bin_labels = {}
    # unihist = get_unified_bin_hist(hists)  # empty hist with logical OR of bins in <hists> (with some assumptions... shit is complicated yo)
    # for hist in hists:
    #     # if hist.Integral() == 0.0:
    #     #     continue
    #     means.append(hist.GetMean())
    #     sems.append(hist.GetMeanError())  # NOTE is this actually right? depends if root uses .Integral() or .GetEntries() in the GetMeanError() call
    #     sum_total += hist.GetMean() * hist.Integral()
    #     total_entries += hist.Integral()

    #     for ib in range(1, hist.GetNbinsX()+1):  # NOTE ignoring under/overflows
    #         ibin = unihist.find_bin(hist.GetBinCenter(ib))
    #         if ibin not in bin_values:
    #             bin_values[ibin] = []
    #             bin_labels[ibin] = hist.GetXaxis().GetBinLabel(ib)
    #         bin_values[ibin].append(hist.GetBinContent(ib))
    #         assert bin_labels[ibin] == hist.GetXaxis().GetBinLabel(ib)

    # # find the mean over hists
    # mean_of_means = 0.0 if total_entries == 0 else sum_total / total_entries
    # # then "normalize" each human's mean by this mean over humans, and that human's variance
    # normalized_means = []
    # for im in range(len(means)):
    #     if sems[im] > 0.0:
    #         normalized_means.append((means[im] - mean_of_means) / sems[im])
    #     else:
    #         normalized_means.append(0)

    # binlist = sorted(bin_values)  # list of all bin indices that were filled in any hist
    # for ibin in binlist:  # NOTE this is *not* a weighted mean, i.e. you better have made sure the subsets all have the same sample size
    #     unihist.set_ibin(ibin, numpy.mean(bin_values[ibin]), error=numpy.std(bin_values[ibin]), label=bin_labels[ibin])

    # return { 'means':means, 'sems':sems, 'normalized_means':normalized_means, 'mean_bin_hist':unihist }

# ----------------------------------------------------------------------------------------
def add_gene_calls_vs_mute_freq_plots(args, hists, rebin=1., debug=False):
    print 'TODO what\'s up with rebin rescaling below?'
    for idir in range(len(args.names)):
        name = args.names[idir]
        for region in utils.regions:
            hright = hists[idir][region + '_gene_right_vs_mute_freq']
            hwrong = hists[idir][region + '_gene_wrong_vs_mute_freq']
            hdenom = copy.deepcopy(hright)
            hdenom.add(hwrong)
            hfrac = copy.deepcopy(hright)
            hfrac.divide_by(hdenom)
            # hfrac.Scale(1. / rebin)  
            if debug:
                print name, region
            for ib in range(hfrac.n_bins + 2):
                lo, hi, cached = fraction_uncertainty.err(hright.bin_contents[ib], hdenom.bin_contents[ib], for_paper=True)
                hfrac.errors[ib] = (hi - lo) / 2.
                if debug:
                    print '%5d %5d   %.2f   %.3f' % (hright.bin_contents[ib], hdenom.bin_contents[ib], hfrac.bin_contents[ib], (hi - lo) / 2.)
            hists[idir][region + '_gene_fraction_vs_mute_freq'] = hfrac

# ----------------------------------------------------------------------------------------
def get_hists_from_dir(dirname, histname, string_to_ignore=None):
    hists = {}
    for fname in glob.glob(dirname + '/*.csv'):
        varname = os.path.basename(fname).replace('.csv', '')
        if string_to_ignore is not None:
            varname = varname.replace(string_to_ignore, '')
        try:
            hists[varname] = Hist(fname=fname, title=histname)
        except KeyError:  # probably not a histogram csv
            pass
    if len(hists) == 0:
        raise Exception('ERROR no csvs in %s' % dirname)
    return hists

# ----------------------------------------------------------------------------------------
def compare_directories(args, xtitle='', use_hard_bounds=''):
    """ 
    Read all the histograms stored as .csv files in <args.plotdirs>, and overlay them on a new plot.
    If there's a <varname> that's missing from any dir, we skip that plot entirely and print a warning message.
    """
    utils.prep_dir(args.outdir + '/plots', multilings=['*.png', '*.svg', '*.csv'])
    if args.leaves_per_tree is not None:
        assert len(args.leaves_per_tree) == len(args.plotdirs)

    # read hists from <args.plotdirs>
    hists = []
    for idir in range(len(args.plotdirs)):
        string_to_ignore = None if args.strings_to_ignore is None else args.strings_to_ignore[idir]
        hist_list = get_hists_from_dir(args.plotdirs[idir] + '/plots', args.names[idir], string_to_ignore=string_to_ignore)
        hists.append(hist_list)

    # then loop over all the <varname>s we found
    all_names, all_means, all_sems, all_normalized_means = [], [], [], []
    # ----------------------------------------------------------------------------------------
    # vs_rebin = 2
    vs_rebin = 1
    if 'v_gene_right_vs_mute_freq' in hists[0].keys():
        add_gene_calls_vs_mute_freq_plots(args, hists, rebin=vs_rebin)
    # ----------------------------------------------------------------------------------------
    for varname, hist in hists[0].items():
        # add the hists
        all_hists = [hist,]
        missing_hist = False
        for idir in range(1, len(args.plotdirs)):
            try:  # add the hist
                all_hists.append(hists[idir][varname])
            except KeyError:  # oops, didn't find it in this dir, so skip this variable entirely
                print args.names[idir], varname
                all_hists.append(Hist(1, 0, 1))

        if '_gene' in varname and '_vs_' not in varname:  # for the gene usage frequencies we need to make sure all the plots have the genes in the same order
            all_hists = add_bin_labels_not_in_all_hists(all_hists)

        if args.calculate_mean_info:
            meaninfo = get_mean_info(all_hists)
            all_names.append(varname)
            all_means.append(meaninfo['means'])
            all_sems.append(meaninfo['sems'])
            all_normalized_means.append(meaninfo['normalized_means'])
            meaninfo['mean_bin_hist'].write(args.outdir + '/plots/' + varname + '-mean-bins.csv')

        # bullshit complicated config stuff
        bounds, no_labels, figsize = None, False, None
        translegend = (0.0, -0.2)
        extrastats, log = '', ''
        xtitle, ytitle, xline, normalization_bounds = hist.xtitle, hist.ytitle, None, None
        simplevarname = varname.replace('-mean-bins', '')
        plottitle = plotconfig.plot_titles[simplevarname] if simplevarname in plotconfig.plot_titles else simplevarname

        if args.normalize:
            ytitle = 'frequency'

        if 'mute-freqs/v' in args.plotdirs[0] or 'mute-freqs/d' in args.plotdirs[0] or 'mute-freqs/j' in args.plotdirs[0]:
            assert not args.normalize
            ytitle = 'mutation freq'

        if '_gene' in varname and '_vs_' not in varname:
            xtitle = 'allele'
            if hist.n_bins == 2:
                extrastats = ' 0-bin'  # print the fraction of entries in the zero bin into the legend (i.e. the fraction correct)
        else:
            xtitle = 'bases'

        line_width_override = None
        rebin = args.rebin
        errors = not args.no_errors
        if args.plot_performance:
            if 'hamming_to_true_naive' in varname:
                xtitle = 'hamming distance'
                if '_normed' in varname:
                    xtitle = 'fractional ' + xtitle
            elif '_vs_mute_freq' in varname:
                xtitle = 'mutation freq'
                ytitle = 'fraction correct'
                if varname[0] == 'v' or varname[0] == 'j':
                    translegend = (-0.4, -0.4)
                # errors = True
                rebin = vs_rebin
            else:
                xtitle = 'inferred - true'
            bounds = plotconfig.true_vs_inferred_hard_bounds.setdefault(varname, None)
        else:
            bounds = plotconfig.default_hard_bounds.setdefault(varname.replace('-mean-bins', ''), None)
            if bounds is None and 'insertion' in varname:
                bounds = plotconfig.default_hard_bounds.setdefault('all_insertions', None)
            if '_gene' in varname and '_vs_' not in varname:
                no_labels = True
                if 'j_' not in varname:
                    figsize = (3, 1.5)  #1000, 500
                line_width_override = 1
            elif 'mute-freqs/v' in args.plotdirs[0] or 'mute-freqs/j' in args.plotdirs[0]:
                figsize = (3, 1.5)  #1000, 500
                bounds = plotconfig.default_hard_bounds.setdefault(utils.unsanitize_name(varname.replace('-mean-bins', '')), None)

        if 'IGH' in varname:
            if 'mute-freqs' in args.plotdirs[0]:
                gene = utils.unsanitize_name(simplevarname)
                plottitle = gene  # + ' -- mutation frequency'
                xtitle = 'position'
                if utils.get_region(gene) == 'j':
                    translegend = (0.1, 0.)  #(-0.35, -0.02)
                else:
                    translegend = (0.15, -0.02)
                xline = None
                if utils.get_region(gene) == 'v' and args.cyst_positions is not None:
                    xline = args.cyst_positions[gene]
                    # normalization_bounds = (int(cyst_positions[gene]) - 70, None)
                elif utils.get_region(gene) == 'j' and args.tryp_positions is not None:
                    xline = args.tryp_positions[gene]
                    # normalization_bounds = (None, int(tryp_positions[gene]) + 5)
            else:
                ilastdash = simplevarname.rfind('-')
                gene = utils.unsanitize_name(simplevarname[:ilastdash])
                base_varname = simplevarname[ilastdash + 1 :]
                base_plottitle = plotconfig.plot_titles[base_varname] if base_varname in plotconfig.plot_titles else ''
                plottitle = gene + ' -- ' + base_plottitle

        # draw that little #$*(!
        linewidths = [line_width_override, ] if line_width_override is not None else args.linewidths
        assert args.leaves_per_tree is None
        # scale_errors = math.sqrt(args.leaves_per_tree[idir]) if args.leaves_per_tree is not None else args.scale_errors
        draw_no_root(all_hists[0], plotname=varname, plotdir=args.outdir, more_hists=all_hists[1:], write_csv=False, stats=args.stats + ' ' + extrastats, bounds=bounds,
                     shift_overflows=False, errors=errors, scale_errors=args.scale_errors, rebin=rebin, plottitle=plottitle, colors=args.colors, linestyles=args.linestyles,
                     xtitle=xtitle, ytitle=ytitle, xline=xline, normalize=(args.normalize and '_vs_mute_freq' not in varname),
                     linewidths=linewidths, markersizes=args.markersizes, figsize=figsize, no_labels=no_labels, log=log, translegend=translegend)

    if args.calculate_mean_info:
        # write mean info
        with opener('w')(args.outdir + '/plots/means.csv') as meanfile:
            writer = csv.DictWriter(meanfile, ('name', 'means', 'sems', 'normalized-means'))
            writer.writeheader()
            for ivar in range(len(all_means)):
                writer.writerow({
                    'name':all_names[ivar],
                    'means':':'.join([str(m) for m in all_means[ivar]]),
                    'sems':':'.join([str(s) for s in all_sems[ivar]]),
                    'normalized-means':':'.join([str(nm) for nm in all_normalized_means[ivar]])
                })

    if not args.only_csv_plots:
        check_call(['./bin/permissify-www', args.outdir])  # NOTE this should really permissify starting a few directories higher up
        check_call(['./bin/makeHtml', args.outdir, '3', 'null', 'svg'])

# ----------------------------------------------------------------------------------------
def make_mean_plots(plotdir, subdirs, outdir):
    meanlist, variancelist = [], []
    normalized_means = []
    for sd in subdirs:
        with opener('r')(plotdir + '/' + sd + '/plots/means.csv') as meanfile:
            reader = csv.DictReader(meanfile)
            for line in reader:
                means = [ float(m) for m in line['means'].split(':') ]
                meanlist.append(numpy.mean(means))
                variancelist.append(numpy.var(means))
                nmvals = [ float(nm) for nm in line['normalized-means'].split(':') ]
                normalized_means += nmvals

    # ----------------------------------------------------------------------------------------
    # first make hexbin plot
    plt.subplot(111)
    plt.hexbin(meanlist, variancelist, gridsize=20, cmap=matplotlib.cm.gist_yarg, bins=None)
    # plt.axis([0, 5, 0, 2])
    plt.xlabel('mean')
    plt.ylabel('variance')

    cb = plt.colorbar()
    cb.set_label('mean value')
    utils.prep_dir(outdir + '/plots', multilings=['*.png', '*.svg', '*.csv'])
    plt.savefig(outdir + '/plots/hexmeans.png')
    plt.close()
    plt.clf()

    # ----------------------------------------------------------------------------------------
    # then make normalized mean plot
    n, bins, patches = plt.hist(normalized_means, 50)
    plt.xlabel(r'$(x_i - \mu) / \sigma_i$')
    plt.title(r'$\sigma=' + str(math.sqrt(numpy.var(normalized_means))) + '$')
    # plt.axis([-10, 10, 0, 220])

    plt.savefig(outdir + '/plots/means.png')
    plt.close()


    check_call(['./permissify-www', outdir])  # NOTE this should really permissify starting a few directories higher up

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
    # hist.all_data = sizes
    return hist

# ----------------------------------------------------------------------------------------
def make_mean_hist(hists, debug=False):
    """ return the hist with bin contents the mean over <hists> of each bin """
    binvals = {}
    all_data = None
    for hist in hists:
        if debug:
            print '    sub',
        for ib in range(0, hist.n_bins + 2):
            low_edge = hist.low_edges[ib]
            if low_edge not in binvals:
                binvals[low_edge] = 0.
            binvals[low_edge] += hist.bin_contents[ib]
            if debug:
                print '   ', low_edge, hist.bin_contents[ib],
        if all_data is not None and hist.all_data is None:
            raise Exception('tried to average hists with and without all_data set')
        if hist.all_data is not None:
            if all_data is None:
                all_data = []
            all_data += hist.all_data
        if debug:
            print ''
    binlist = sorted(binvals.keys())
    meanhist = Hist(len(binlist) - 2, binlist[1], binlist[-1], binlist[1 : -1])
    meanhist.all_data = all_data
    if debug:
        print '   mean',
    for ib in range(len(binlist)):
        meanhist.set_ibin(ib, binvals[binlist[ib]])
        if debug:
            print '   ', meanhist.low_edges[ib], meanhist.bin_contents[ib],
    if debug:
        print ''

    meanhist.normalize()
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
legends = {'vollmers-0.9' : 'VJ CDR3 0.9',
           'partition partis' : 'full partis',
           'partition' : 'full partis',
           'naive-hamming-partition partis' : 'point partis',
           'naive-hamming-partition' : 'point partis',
           'vsearch-partition partis' : 'vsearch partis',
           'vsearch-partition' : 'vsearch partis',
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

# linewidths['v-true'] = 10
# linewidths['cdr3-true'] = 10
# colors['v-true'] = '#006600'
# colors['cdr3-true'] = '#006600'
# colors['v-indels'] = '#cc0000'
# colors['cdr3-indels'] = '#cc0000'

# ----------------------------------------------------------------------------------------
def plot_cluster_size_hists(outfname, hists, title, xmax=None, log='x'):

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

        # was using this:
        hist.normalize()
        plots[name] = ax.plot(hist.get_bin_centers(), hist.bin_contents_no_zeros(1e-8), linewidth=linewidths.get(name, 4), label=legends.get(name, name), color=colors.get(name, 'grey'), linestyle=linestyles.get(name, 'solid'), alpha=alphas.get(name, 1.))

        # # or maybe try a hist?
        # if n_queries is None:
        #     n_queries = sum(hist.all_data)
        # if tmpmax is None or max(hist.all_data) > tmpmax:
        #     tmpmax = max(hist.all_data)
        # plots[name] = ax.hist(hist.all_data, bins=[0.5, 1.5, 2.5, 3.5, 5.5, 10, 25, 50, 80, 110, 200, 350, 500, 1000], histtype='step', normed=True, linewidth=linewidths.get(name, 4), label=legends.get(name, name), color=colors.get(name, 'grey'), linestyle=linestyle, alpha=alpha)


    legend = ax.legend()
    sns.despine()  #trim=True, bottom=True)

    # xmax = tmpmax
    if xmax is None:
        xmax = plt.gca().get_xlim()[1]
    else:
        ax.set_xlim(0.9, xmax)

    if 'stanford' in title:
        ymin = 5e-4
    else:
        ymin = 5e-5
    plt.ylim(ymin, 1)
    plt.title(title)
    plt.xlabel('cluster size')
    plt.ylabel('fraction of clusters')
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
def plot_adj_mi_and_co(plotvals, mut_mult, plotdir, valname, xvar, title=''):
    fig, ax = mpl_init()
    mpl.rcParams.update({
        'legend.fontsize': 15,})
    plots = {}
    for meth, xyvals in plotvals.items():
        xvals = xyvals.keys()
        yvals = [ve[0] for ve in xyvals.values()]
        yerrs = [ve[1] for ve in xyvals.values()]
        kwargs = {'linewidth' : linewidths.get(meth, 4),
                  'label' : legends.get(meth, meth),
                  'color' : colors.get(meth, 'grey'),
                  'linestyle' : linestyles.get(meth, 'solid'),
                  'alpha' : alphas.get(meth, 1.)
                  }
        if xvar == 'n_leaves':
            kwargs['fmt'] = '-o'
            plots[meth] = ax.errorbar(xvals, yvals, yerr=yerrs, **kwargs)
        else:  # darn it, the order in the legend gets messed up if I do some as .plot and some as .errorbar
            plots[meth] = ax.plot(xvals, yvals, **kwargs)
    
    lx, ly = 1.6, 0.7
    legend = ax.legend(bbox_to_anchor=(lx, ly))
    # legend.get_frame().set_facecolor('white')
    ymin = -0.01
    ax.set_ylim(ymin, 1.01)
    sns.despine()  #trim=True, bottom=True)
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
        ax.set_xlim(xvals[0], xvals[-1])
    elif xvar == 'nseqs':
        # xticks = [xticks[i] for i in range(0, len(xticks), 2)]
        # if 750 in xticks:
        #     xticks.remove(750)
        # xticks += xvals[-1:]
        # xticks = [100, 5000, 10000, 15000]
        xticks = [100, 1000, 10000, 200000]
        ax.set_xscale('log')
        ax.set_xlim(0.9 * xvals[0], 1.05 * xvals[-1])

    xticklabels = xticks if xvar == 'n_leaves' else ['%.0e' % xt for xt in xticks]
    plt.xticks(xticks, xticklabels)

    yticks = [yt for yt in [0., .2, .4, .6, .8, 1.] if yt >= ymin]
    yticklabels = [str(yt) for yt in yticks]
    plt.yticks(yticks, yticklabels)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)
    plotname =  valname + '-%d-mutation.svg' % mut_mult
    plt.savefig(plotdir + '/' + plotname)
    plt.close()

# ----------------------------------------------------------------------------------------
def mpl_init(figsize=None, fontsize=20):
    sns.set_style('ticks')
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
def mpl_finish(ax, plotdir, plotname, title='', xlabel='', ylabel='', xbounds=None, ybounds=None, leg_loc=(0.04, 0.6), log='', xticks=None, xticklabels=None):
    # xticks[0] = 0.000001
    legend = ax.legend(loc=leg_loc)
    plt.gcf().subplots_adjust(bottom=0.14, left=0.18, right=0.95, top=0.92)
    sns.despine()  #trim=True, bottom=True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if 'x' in log:
        ax.set_xscale('log')
    if 'y' in log:
        ax.set_yscale('log')
    if xbounds is not None:
        plt.xlim(xbounds[0], xbounds[1])
    if ybounds is not None:
        plt.ylim(ybounds[0], ybounds[1])
    if xticks is not None:
        plt.xticks(xticks)
    if xticklabels is not None:
        ax.set_xticklabels(xticklabels)
    plt.title(title)
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    plt.savefig(plotdir + '/' + plotname + '.svg')
    plt.close()

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
    
    ticks = [n - 0.5 for n in range(1, n_biggest_clusters + 1, 2)]
    xticklabels = [str(int(n + 0.5)) for n in ticks]
    yticklabels = xticklabels
    if n_biggest_clusters > 20:
        modulo = 3
        ticks = [ticks[it] for it in range(0, len(ticks), modulo)]
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
def make_html(plotdir, n_columns=3, extension='svg'):
    if plotdir[-1] == '/':  # remove trailings slash, if present
        plotdir = plotdir[:-1]
    if not os.path.exists(plotdir):
        raise Exception('plotdir %s d.n.e.' % plotdir)
    dirname = os.path.basename(plotdir)
    lines = ['<!DOCTYPE html', 
             '    PUBLIC "-//W3C//DTD HTML 3.2//EN">', 
             '<html>', 
             '<head><title>foop</title></head>', 
             '<body bgcolor="000000">', 
             '<h3 style="text-align:left; color:DD6600;">foop</h3>', 
             '', 
             '<table border="0" cellspacing="5" width="100%">', 
             '<tr>']

    fnames = sorted(glob.glob(plotdir + '/*.' + extension))
    for ifn in range(len(fnames)):
        if ifn > 0 and ifn % n_columns == 0:
            lines += ['</tr>', '<tr>']
        fname = os.path.basename(fnames[ifn])
        line = '<td width="25%"><a target="_blank" href="' + dirname + '/' + fname + '"><img src="' + dirname + '/' + fname + '" alt="' + dirname + '/' + fname + '" width="100%"></a></td>"'
        lines.append(line)

    lines += ['</tr>',
              '</table>',
              '</body>',
              '</html>']

    htmlfname = os.path.dirname(plotdir) + '/' + dirname + '.html'  # more verbose than necessary
    with open(htmlfname, 'w') as htmlfile:
        htmlfile.write('\n'.join(lines))
    check_call(['chmod', '664', htmlfname])
