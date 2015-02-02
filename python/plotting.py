import math
import os
import glob
import sys
import csv
import numpy
from scipy import stats
from array import array
from subprocess import check_call

import utils
import fraction_uncertainty
import plotconfig
from hist import Hist

def check_root():
    try:
        from ROOT import kBlue
        return True
    except ImportError:
        return False

sys.argv.append('-b')  # root just loves its stupid little splashes
has_root = check_root()
if has_root:
    from ROOT import gStyle, TH1D, TCanvas, kRed, gROOT, TLine, TH2Poly, TLegend, kBlue, kGreen, kCyan, kOrange
    gROOT.Macro("plotting/MitStyleRemix.cc+")
else:
    print ' ROOT not found, proceeding without plotting'

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
    assert False  # I *think* (well, hope) I'm not using this any more
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
def make_hist_from_bin_entry_file(fname, hist_label='', log='', normalize=False, rescale_entries=None):
    """ 
    Return root histogram with each bin low edge and bin content read from <fname> 
    E.g. from the results of hist.Hist.write()
    """
    low_edges, contents, bin_labels, bin_errors, sum_weights_squared = [], [], [], [], []
    xtitle = ''
    with opener('r')(fname) as infile:
        reader = csv.DictReader(infile)
        for line in reader:
            low_edges.append(float(line['bin_low_edge']))
            refactor = 1.0 if rescale_entries is None else rescale_entries
            contents.append(float(line['contents']) / refactor)
            if 'sum-weights-squared' in line:
                sum_weights_squared.append(float(line['sum-weights-squared']) / (refactor*refactor))
            if 'error' in line:
                assert 'sum-weights-squared' not in line
                bin_errors.append(float(line['error']) * math.sqrt(refactor))
            if 'binlabel' in line:
                bin_labels.append(line['binlabel'])
            else:
                bin_labels.append('')
            if 'xtitle' in line:
                xtitle = line['xtitle']

    n_bins = len(low_edges) - 2  # file should have a line for the under- and overflow bins
    xbins = array('f', [0.0 for i in range(n_bins+1)])  # NOTE has to be n bins *plus* 1
    low_edges = sorted(low_edges)
    for ib in range(n_bins+1):
        xbins[ib] = low_edges[ib+1]  # low_edges[1] is the lower edge of the first bin, i.e. the first bin after the underflow bin, and this will set the last entry in xbins to lower[n_bins+1], i.e. the lower edge of the overflow bin. Which, I bloody well think, is correct
    hist = TH1D(hist_label, '', n_bins, xbins)  # this will barf if the csv file wasn't sorted by bin low edge
    hist.GetXaxis().SetTitle(xtitle)
    for ib in range(n_bins+2):
        hist.SetBinContent(ib, contents[ib])
        if len(sum_weights_squared) > 0:
            hist.SetBinError(ib, math.sqrt(sum_weights_squared[ib]))
        elif len(bin_errors) > 0:
            hist.SetBinError(ib, bin_errors[ib])
        else:
            hist.SetBinError(ib, math.sqrt(contents[ib]))
        if bin_labels[ib] != '':
            hist.GetXaxis().SetBinLabel(ib, bin_labels[ib])

    if normalize and hist.Integral() > 0.0:
        hist.Scale(1./hist.Integral())
        hist.GetYaxis().SetTitle('freq')

    return hist
    
# ----------------------------------------------------------------------------------------
def make_hist_from_observation_file(fname, column, hist_label='', n_bins=30, log=''):
    """ return root histogram filled with each value from <column> in csv file <fname> """
    if not has_root:
        return None
    values = []
    with opener('r')(fname) as infile:
        reader = csv.DictReader(infile)
        for line in reader:
            values.append(float(line[column]))

    values = sorted(values)
    xbins = array('f', [0 for i in range(n_bins+1)])  # NOTE has to be n_bins *plus* 1
    set_bins(values, n_bins, 'x' in log, xbins, var_type='float')
    hist = TH1D(hist_label, '', n_bins, xbins)
    for value in values:
        hist.Fill(value)

    return hist

# ----------------------------------------------------------------------------------------
def make_bool_hist(n_true, n_false, hist_label):
    """ fill a two-bin histogram with the fraction false in the first bin and the fraction true in the second """
    hist = TH1D(hist_label, '', 2, -0.5, 1.5)
    hist.Sumw2()

    true_frac = float(n_true) / (n_true + n_false)
    hist.SetBinContent(1, true_frac)
    true_bounds = fraction_uncertainty.err(n_true, n_true + n_false, use_beta=True)
    hist.SetBinError(1, max(abs(true_frac - true_bounds[0]), abs(true_bounds[1] - true_bounds[1])))
    false_frac = float(n_false) / (n_true + n_false)
    hist.SetBinContent(2, false_frac)
    false_bounds = fraction_uncertainty.err(n_false, n_true + n_false, use_beta=True)
    hist.SetBinError(2, max(abs(false_frac - false_bounds[0]), abs(false_bounds[1] - false_bounds[1])))

    hist.GetXaxis().SetNdivisions(0)
    hist.GetXaxis().SetBinLabel(1, 'right')
    hist.GetXaxis().SetBinLabel(2, 'wrong')
    hist.GetXaxis().SetLabelSize(0.1)

    return hist

# ----------------------------------------------------------------------------------------
def make_hist_from_list(values, hist_label, n_bins=30):
    """ Fill a histogram with float values in a list """
    if len(values) == 0:
        print 'WARNING no values for %s in make_hist' % hist_label
        return TH1D(hist_label, '', 1, 0, 1)
    xbins = array('f', [0 for i in range(n_bins+1)])  # NOTE has to be n_bins *plus* 1
    set_bins(values, n_bins, is_log_x=False, xbins=xbins, var_type='float')
    hist = TH1D(hist_label, '', n_bins, xbins)
    for val in values:
        hist.Fill(val)
    return hist

# ----------------------------------------------------------------------------------------
# <values> is of form {<bin 1>:<counts 1>, <bin 2>:<counts 2>, ...}
def make_hist_from_dict_of_counts(values, var_type, hist_label, log='', xmin_force=0.0, xmax_force=0.0, normalize=False, sort=False):
    """ Fill a histogram with values from a dictionary (each key will correspond to one bin) """
    assert var_type == 'int' or var_type == 'string'  # floats should be handled by Hist class in hist.py

    if not has_root:
        return
    if len(values) == 0:
        print 'WARNING no values for %s in make_hist' % hist_label
        return TH1D(hist_label, '', 1, 0, 1)

    bin_labels = sorted(values)
    if not sort and var_type == 'string':  # for strings, sort so most common value is to left side
        bin_labels = sorted(values, key=values.get, reverse=True)

    if var_type == 'string':
        n_bins = len(values)
    else:
        n_bins = bin_labels[-1] - bin_labels[0] + 1

    hist = None
    xbins = array('f', [0 for i in range(n_bins+1)])  # NOTE has to be n_bins *plus* 1
    if xmin_force == xmax_force:  # if boundaries aren't set explicitly, work out what they should be
        if var_type == 'string':
            set_bins(bin_labels, n_bins, 'x' in log, xbins, var_type)
            hist = TH1D(hist_label, '', n_bins, xbins)
        else:
            hist = TH1D(hist_label, '', n_bins, bin_labels[0] - 0.5, bin_labels[-1] + 0.5)  # for integers, just go from the first to the lat bin label (they're sorted)
    else:
      hist = TH1D(hist_label, '', n_bins, xmin_force, xmax_force)
    hist.Sumw2()

    for ival in range(len(values)):
        if var_type == 'string':
            label = bin_labels[ival]
            ibin = ival + 1
            hist.GetXaxis().SetBinLabel(ibin, label)
        else:
            ibin = hist.FindBin(bin_labels[ival])
        hist.SetBinContent(ibin, values[bin_labels[ival]])
        hist.SetBinError(ibin, math.sqrt(values[bin_labels[ival]]))
  
    # make sure there's no overflows
    if hist.GetBinContent(0) != 0.0 or hist.GetBinContent(hist.GetNbinsX()+1) != 0.0:
        print 'overflows in ',hist.GetName(),'exiting'
        for ibin in range(0, hist.GetNbinsX()+2):
            print '%d %f %f %f' % (ibin, hist.GetXaxis().GetBinLowEdge(ibin), hist.GetXaxis().GetBinUpEdge(ibin), hist.GetBinContent(ibin))
        sys.exit()

    hist.GetYaxis().SetTitle('counts')
    if normalize:
        hist.Scale(1./hist.Integral())
        hist.GetYaxis().SetTitle('freq')
    
    return hist

# ----------------------------------------------------------------------------------------
def make_hist_from_my_hist_class(myhist, name):
    roothist = TH1D(name, '', myhist.n_bins, myhist.xmin, myhist.xmax)
    for ibin in range(myhist.n_bins + 2):
        roothist.SetBinContent(ibin, myhist.bin_contents[ibin])
        roothist.SetBinError(ibin, math.sqrt(myhist.sum_weights_squared[ibin]))
        # if 'all-mean-freq' in name:
        #     print '  ', myhist.bin_contents[ibin], math.sqrt(myhist.sum_weights_squared[ibin])
    return roothist

# ----------------------------------------------------------------------------------------
def add_bin_labels_not_in_all_hists(hists):
    """ find the OR of all bin labels present in <hists>, and remake each hist in <hists> to have zero bins for any that weren't there already """
    # first convert each hist to a map from bin label to entries
    all_labels = []
    histmaps = []
    for hist in hists:
        histmaps.append({})
        for ibin in range(1, hist.GetNbinsX() + 1):  # ignore under/over flows, they're kinda useless for bin-labelled hists
            label = hist.GetXaxis().GetBinLabel(ibin)
            histmaps[-1][label] = (hist.GetBinContent(ibin), hist.GetBinError(ibin))  # 2-tuple with (content, error)
            if label not in all_labels:
                all_labels.append(label)

    all_labels = sorted(all_labels)

    # then go through and make new histograms for everybody
    finalhists = []
    for ih in range(len(histmaps)):
        original_hist = hists[ih]
        hmap = histmaps[ih]
        finalhists.append(TH1D('uni-bin-label-' + original_hist.GetName(), original_hist.GetTitle(), len(all_labels), 0.5, len(all_labels) + 0.5))
        for ilabel in range(len(all_labels)):
            label = all_labels[ilabel]
            ibin = ilabel + 1  # root conventions
            finalhists[-1].GetXaxis().SetBinLabel(ibin, label)
            if label in hmap:
                finalhists[-1].SetBinContent(ibin, hmap[label][0])
                finalhists[-1].SetBinError(ibin, hmap[label][1])
            else:
                finalhists[-1].SetBinContent(ibin, 0.0)
                finalhists[-1].SetBinError(ibin, 0.0)

    return finalhists

# ----------------------------------------------------------------------------------------
def draw(hist, var_type, log='', plotdir=None, plotname='foop', more_hists=None, write_csv=False, stats=None, bounds=None,
         errors=False, shift_overflows=False, csv_fname=None, scale_errors=None, rebin=None, plottitle='',
         colors=None, linestyles=None, cwidth=700, cheight=600, imagetype='svg', xtitle=None, ytitle=None,
         xline=None, yline=None, draw_str=None):
    assert os.path.exists(plotdir)
    if not has_root:
        return

    if hist.GetNbinsX() > 2:         #and ('_gene' in plotname or
        if 'mute-freqs/v' in plotdir:
            cwidth, cheight = 1200, 600
        elif 'mute-freqs/j' in plotdir:
            cwidth, cheight = 1000, 600
    cvn = TCanvas('cvn-'+plotname, '', cwidth, cheight)

    # if 'd_gene' in plotname:
    #     for ib in range(hist.GetNbinsX()+1):
    #         print hist.GetXaxis().GetBinLowEdge(ib), hist.GetBinContent(ib)
    hists = [hist,]
    if more_hists != None:
        hists = hists + more_hists

    xmin, xmax, ymax = None, None, None
    for htmp in hists:
        if rebin is not None:
            htmp.Rebin(rebin)
        if scale_errors is not None:
            for ibin in range(htmp.GetNbinsX()+2):
                htmp.SetBinError(ibin, htmp.GetBinError(ibin) * scale_errors)

        if xmin is None or htmp.GetBinLowEdge(1) < xmin:
            xmin = htmp.GetBinLowEdge(1)
        if xmax is None or htmp.GetXaxis().GetBinUpEdge(htmp.GetNbinsX()) > xmax:
            xmax = htmp.GetXaxis().GetBinUpEdge(htmp.GetNbinsX())
        if ymax is None or htmp.GetMaximum() > ymax:
            ymax = htmp.GetMaximum()
    if bounds != None:
        xmin, xmax = bounds
    hframe = TH1D('hframe', '', hist.GetNbinsX(), xmin, xmax)
    if var_type == 'string' or var_type == 'bool':
        for ib in range(1, hframe.GetNbinsX()+1):
            hframe.GetXaxis().SetBinLabel(ib, hist.GetXaxis().GetBinLabel(ib))

    hframe.SetMaximum(1.2*ymax)
    if var_type == 'bool':
        hframe.GetXaxis().SetLabelSize(0.1)
    if '_gene' in plotname:
        hframe.GetXaxis().SetLabelSize(0.03)
    if plottitle == '':
        plottitle = plotname
    if xtitle is None:
        xtitle = hist.GetXaxis().GetTitle()
        if xtitle == '' and '_gene' not in plotname:
            xtitle = 'bases'
    if ytitle is None:
        ytitle = hframe.GetYaxis().GetTitle()
        if ytitle == '':
            ytitle = 'freq'
    hframe.SetTitle(plottitle + ';' + xtitle + ';' + ytitle)
    hframe.Draw('txt')

    if shift_overflows:
        for htmp in hists:
            if htmp == None:
                continue
            underflows, overflows = 0.0, 0.0
            first_shown_bin, last_shown_bin = -1, -1
            for ib in range(0, htmp.GetXaxis().GetNbins()+2):
                if htmp.GetXaxis().GetBinCenter(ib) <= xmin:
                    underflows += htmp.GetBinContent(ib)
                    htmp.SetBinContent(ib, 0.0)
                elif first_shown_bin == -1:
                    first_shown_bin = ib
                else:
                    break
            for ib in reversed(range(0, htmp.GetXaxis().GetNbins()+2)):
                if htmp.GetXaxis().GetBinCenter(ib) >= xmax:
                    overflows += htmp.GetBinContent(ib)
                    htmp.SetBinContent(ib, 0.0)
                elif last_shown_bin == -1:
                    last_shown_bin = ib
                else:
                    break

            if 'd_hamming' in plotname:
                print htmp.GetTitle()
                print '  underflow', underflows, htmp.GetBinContent(first_shown_bin)
                print '  overflow', overflows, htmp.GetBinContent(last_shown_bin)
                print '  first', htmp.GetXaxis().GetBinCenter(first_shown_bin)
                print '  last', htmp.GetXaxis().GetBinCenter(last_shown_bin)
            htmp.SetBinContent(first_shown_bin, underflows + htmp.GetBinContent(first_shown_bin))
            htmp.SetBinContent(last_shown_bin, overflows + htmp.GetBinContent(last_shown_bin))

    if colors is None:
        assert len(hists) < 5
        colors = (kRed, kBlue-4, kGreen+2, kOrange+1)  # 632, 596, 418, 801
    else:
        assert len(hists) <= len(colors)
    if linestyles is None:
        # assert len(hists) < 5
        linestyles = [1 for _ in range(len(hists))]
    else:
        assert len(hists) <= len(linestyles)

    if draw_str is None:
        draw_str = 'hist same'
    else:
        draw_str += ' same'
    if errors:  # not working!
        draw_str = 'e ' + draw_str
    for ih in range(len(hists)):
        htmp = hists[ih]
        htmp.SetLineColor(colors[ih])
        # if ih == 0:
        htmp.SetMarkerSize(0)
        htmp.SetMarkerColor(colors[ih])
        htmp.SetLineStyle(linestyles[ih])
        if ih < 6 and len(hists) < 5:
            htmp.SetLineWidth(6-ih)
        htmp.Draw(draw_str)

    if len(hists) < 5:
        leg = TLegend(0.57, 0.72, 0.99, 0.9)
    else:
        leg = TLegend(0.57, 0.6, 0.99, 0.9)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    for htmp in hists:
        if stats is not None:
            if 'rms' in stats:
                htmp.SetTitle(htmp.GetTitle() + (' (%.2f)' % htmp.GetRMS()))
            if 'mean' in stats:
                htmp.SetTitle(htmp.GetTitle() + (' (%.2f)' % htmp.GetMean()))
            if '0-bin' in stats:
                htmp.SetTitle(htmp.GetTitle() + (' (%.2f)' % htmp.GetBinContent(1)))
        leg.AddEntry(htmp, htmp.GetTitle() , 'l')
    leg.Draw()

    if xline is not None:
        if xline <= hframe.GetXaxis().GetXmin() or xline >= hframe.GetXaxis().GetXmax():  # make sure we got valid a x position for the line
            print 'WARNING plotting x line at %f out of bounds (%f, %f)' % (float(xmin), hframe.GetXaxis().GetXmin(), hframe.GetXaxis().GetXmax())
        # xl = TLine(xline, hframe.GetYaxis().GetXmin(), xline, 0.5*ymax)
        xl = TLine(xline, -0.1*ymax, xline, 0.5*ymax)
        xl.SetLineStyle(2)
        xl.Draw()
    if yline is not None:
        if yline <= hframe.GetYaxis().GetXmin() or xline >= hframe.GetYaxis().GetXmax():  # make sure we got valid a x position for the line
            print 'WARNING plotting y line at %f out of bounds (%f, %f)' % (float(ymin), hframe.GetYaxis().GetXmin(), hframe.GetYaxis().GetXmax())
        yl = TLine(hframe.GetXaxis().GetXmin(), yline, hframe.GetXaxis().GetXmax(), yline)
        yl.Draw()

    cvn.SetLogx('x' in log)
    cvn.SetLogy('y' in log)
    if not os.path.exists(plotdir + '/plots'):
        print 'ERROR dir \'' + plotdir + '/plots\' d.n.e.'
        assert False

    if write_csv:
        assert more_hists == None
        if csv_fname == None:
            write_hist_to_file(plotdir + '/plots/' + plotname + '.csv', hist)
        else:
            write_hist_to_file(csv_fname, hist)
    cvn.SaveAs(plotdir + '/plots/' + plotname + '.' + imagetype)
    # if '_gene' in plotname:
    #     cvn.SaveAs(plotdir + '/plots/' + plotname + '.png')

# ----------------------------------------------------------------------------------------
def get_hists_from_dir(dirname, histname, rescale_entries=None, normalize=False):
    hists = {}
    for fname in glob.glob(dirname + '/*.csv'):
        varname = os.path.basename(fname).replace('.csv', '')
        try:
            hists[varname] = make_hist_from_bin_entry_file(fname, histname + '-csv-' + varname, normalize=normalize, rescale_entries=rescale_entries)
            hists[varname].SetTitle(histname)
        except KeyError:  # probably not a histogram csv
            pass
    if len(hists) == 0:
        print 'ERROR no csvs in',dirname
        sys.exit()
    return hists

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
    assert len(hists) > 0
    means, sems, normalized_means = [], [], []
    sum_total, total_entries = 0.0, 0.0
    bin_values = {}  # map from bin centers to list (over hists) of entries (corresponds to <ibin> in unihist)
    bin_labels = {}
    unihist = get_unified_bin_hist(hists)
    for hist in hists:
        means.append(hist.GetMean())
        sems.append(hist.GetMeanError())  # NOTE is this actually right? depends if root uses .Integral() or .GetEntries() in the GetMeanError() call
        sum_total += hist.GetMean() * hist.Integral()
        total_entries += hist.Integral()

        for ib in range(1, hist.GetNbinsX()+1):  # NOTE ignoring under/overflows
            ibin = unihist.find_bin(hist.GetBinCenter(ib))
            if ibin not in bin_values:
                bin_values[ibin] = []
                bin_labels[ibin] = hist.GetXaxis().GetBinLabel(ib)
            bin_values[ibin].append(hist.GetBinContent(ib))
            assert bin_labels[ibin] == hist.GetXaxis().GetBinLabel(ib)

    # find the mean over hists
    mean_of_means = sum_total / total_entries
    # then "normalize" each human's mean by this mean over humans, and that human's variance
    normalized_means = []
    for im in range(len(means)):
        if sems[im] > 0.0:
            normalized_means.append((means[im] - mean_of_means) / sems[im])
        else:
            normalized_means.append(0)

    binlist = sorted(bin_values)  # list of all bin indices that were filled in any hist
    # for k,v in bin_values.items():
    #     print k, v
    # for ib in binlist:
    #     print ib
    for ibin in binlist:  # NOTE this is *not* a weighted mean, i.e. you better have made sure the subsets all have the same sample size
        unihist.set_ibin(ibin, numpy.mean(bin_values[ibin]), error=numpy.std(bin_values[ibin]), label=bin_labels[ibin])

    return { 'means':means, 'sems':sems, 'normalized_means':normalized_means, 'mean_bin_hist':unihist }

# ----------------------------------------------------------------------------------------
def compare_directories(outdir, dirs, names, xtitle='', use_hard_bounds='', stats='', errors=True, scale_errors=None, rebin=None, colors=None, linestyles=None, plot_performance=False,
                        cyst_positions=None, tryp_positions=None, leaves_per_tree=None, calculate_mean_info=True):
    """ 
    Read all the histograms stored as .csv files in <dirs>, and overlay them on a new plot.
    If there's a <varname> that's missing from any dir, we skip that plot entirely and print a warning message.
    """
    utils.prep_dir(outdir + '/plots', multilings=['*.png', '*.svg', '*.csv'])
    if leaves_per_tree is not None:
        assert len(leaves_per_tree) == len(dirs)

    # read hists from <dirs>
    hists = []
    for idir in range(len(dirs)):
        rescale_entries = leaves_per_tree[idir] if leaves_per_tree is not None else None
        hists.append(get_hists_from_dir(dirs[idir] + '/plots', names[idir], rescale_entries=rescale_entries, normalize=True))

    # then loop over all the <varname>s we found
    histmisses = []
    all_names, all_means, all_sems, all_normalized_means = [], [], [], []
    for varname, hist in hists[0].iteritems():
        # add the hists
        all_hists = [hist,]
        missing_hist = False
        for idir in range(1, len(dirs)):
            try:  # add the hist
                all_hists.append(hists[idir][varname])
            except KeyError:  # oops, didn't find it in this dir, so skip this variable entirely
                histmisses.append(varname)
                missing_hist = True
                break
        if missing_hist:
            continue        

        if '_gene' in varname:  # for the gene usage frequencies we need to make sure all the plots have the genes in the same order
            all_hists = add_bin_labels_not_in_all_hists(all_hists)

        if calculate_mean_info:
            meaninfo = get_mean_info(all_hists)
            all_names.append(varname)
            all_means.append(meaninfo['means'])
            all_sems.append(meaninfo['sems'])
            all_normalized_means.append(meaninfo['normalized_means'])
            meaninfo['mean_bin_hist'].write(outdir + '/plots/' + varname + '-mean-bins.csv')

        # bullshit complicated config stuff
        var_type = 'int' if hist.GetXaxis().GetBinLabel(1) == '' else 'bool'
        plottitle = plotconfig.plot_titles[varname] if varname in plotconfig.plot_titles else ''
        bounds = None
        if plot_performance:
            bounds = plotconfig.true_vs_inferred_hard_bounds.setdefault(varname, None)
        else:
            bounds = None
            # bounds = default_hard_bounds.setdefault(varname, None)
            # if '_insertion' in varname and 'content' not in varname:  # subsetting by gene now, so the above line doesn't always work
            #     tmpname = varname[ varname.find('_insertion') - 2 : ]
            #     bounds = default_hard_bounds[tmpname]
        extrastats = ''
        if '_gene' in varname:
            if hist.GetNbinsX() == 2:
                extrastats = ' 0-bin'  # print the fraction of entries in the zero bin into the legend (i.e. the fraction correct)
        xtitle, xline, draw_str = None, None, None
        if 'IGH' in varname:
            if 'mute-freqs' in dirs[0]:
                gene = utils.unsanitize_name(varname)
                plottitle = gene + ' -- mutation frequency'
                xtitle = 'position'
                xline = None
                if utils.get_region(gene) == 'v' and cyst_positions is not None:
                    xline = cyst_positions[gene]['cysteine-position']
                elif utils.get_region(gene) == 'j' and tryp_positions is not None:
                    xline = int(tryp_positions[gene])
                draw_str = 'e'
            else:
                ilastdash = varname.rfind('-')
                gene = utils.unsanitize_name(varname[:ilastdash])
                base_varname = varname[ilastdash + 1 :]
                base_plottitle = plotconfig.plot_titles[base_varname] if base_varname in plotconfig.plot_titles else ''
                plottitle = gene + ' -- ' + base_plottitle

        # draw that little #$*(!
        draw(all_hists[0], var_type, plotname=varname, plotdir=outdir, more_hists=all_hists[1:], write_csv=False, stats=stats + ' ' + extrastats, bounds=bounds,
             shift_overflows=False, errors=errors, scale_errors=scale_errors, rebin=rebin, plottitle=plottitle, colors=colors, linestyles=linestyles,
             xtitle=xtitle, xline=xline, draw_str=draw_str)

    if len(histmisses) > 0:
        print 'WARNING: missing hists for %s' % ' '.join(histmisses)

    if calculate_mean_info:
        # write mean info
        with opener('w')(outdir + '/plots/means.csv') as meanfile:
            writer = csv.DictWriter(meanfile, ('name', 'means', 'sems', 'normalized-means'))
            writer.writeheader()
            for ivar in range(len(all_means)):
                writer.writerow({
                    'name':all_names[ivar],
                    'means':':'.join([str(m) for m in all_means[ivar]]),
                    'sems':':'.join([str(s) for s in all_sems[ivar]]),
                    'normalized-means':':'.join([str(nm) for nm in all_normalized_means[ivar]])
                })

    check_call(['./permissify-www', outdir])  # NOTE this should really permissify starting a few directories higher up
    check_call(['./makeHtml', outdir, '3', 'null', 'svg'])

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

    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot

    # ----------------------------------------------------------------------------------------
    # first make hexbin plot
    pyplot.subplot(111)
    pyplot.hexbin(meanlist, variancelist, gridsize=20, cmap=matplotlib.cm.gist_yarg, bins=None)
    # pyplot.axis([0, 5, 0, 2])
    pyplot.xlabel('mean')
    pyplot.ylabel('variance')

    cb = pyplot.colorbar()
    cb.set_label('mean value')
    utils.prep_dir(outdir + '/plots', multilings=['*.png', '*.svg', '*.csv'])
    pyplot.savefig(outdir + '/plots/hexmeans.png')
    pyplot.clf()

    # ----------------------------------------------------------------------------------------
    # then make normalized mean plot
    n, bins, patches = pyplot.hist(normalized_means, 50)
    pyplot.xlabel(r'$(x_i - \mu) / \sigma_i$')
    pyplot.title(r'$\sigma=' + str(math.sqrt(numpy.var(normalized_means))) + '$')
    # pyplot.axis([-10, 10, 0, 220])

    pyplot.savefig(outdir + '/plots/means.png')


    check_call(['./permissify-www', outdir])  # NOTE this should really permissify starting a few directories higher up

# # ----------------------------------------------------------------------------------------
# def make_html(plotdir, filetype):
#     outfname = plotdir + '/plots.html'
#     imagefiles = glob.glob(plotdir + )
