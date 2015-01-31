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

plot_titles = {
    'v_gene' : 'V gene',
    'd_gene' : 'D gene',
    'j_gene' : 'J gene',
    'hamming_to_true_naive' : 'HTTN',
    'v_hamming_to_true_naive' : 'V HTTN',
    'd_hamming_to_true_naive' : 'D HTTN',
    'j_hamming_to_true_naive' : 'J HTTN',
    'v_hamming_to_true_naive_normed' : 'V HTTN',
    'd_hamming_to_true_naive_normed' : 'D HTTN',
    'j_hamming_to_true_naive_normed' : 'J HTTN',
    'd_3p_del' : 'D 3p del',
    'd_5p_del' : 'D 5p del',
    'dj_insertion' : 'DJ insert length',
    'dj_insertion_content' : 'DJ insert base content',
    'j_5p_del' : 'J 5p del',
    'mute_freqs' : 'mute freq',
    'v_3p_del' : 'V 3p del',
    'vd_insertion' : 'VD insert length',
    'vd_insertion_content' : 'VD insert base content',
    'all-mean-freq' : 'sequence mutation freq',
    'v-mean-freq' : 'V mutation freq',
    'd-mean-freq' : 'D mutation freq',
    'j-mean-freq' : 'J mutation freq',
}

true_vs_inferred_hard_bounds = {
    'hamming_to_true_naive' : (-0.5, 19.5),
    'v_hamming_to_true_naive' : (-0.5, 8.5),
    'd_hamming_to_true_naive' : (-0.5, 10.5),
    'j_hamming_to_true_naive' : (-0.5, 12.5),
    'v_hamming_to_true_naive_normed' : (-0.5, 8.5),
    'd_hamming_to_true_naive_normed' : (-0.5, 50.5),
    'j_hamming_to_true_naive_normed' : (-0.5, 20),
    'd_3p_del' : (-8.5, 8.5),
    'd_5p_del' : (-8.5, 8.5),
    'dj_insertion' : (-10.5, 15.5),
    'j_5p_del' : (-10.5, 15.5),
    'mute_freqs' : (-.05, .05),  # NOTE make sure you know where the decimal place is here!
    'v_3p_del' : (-3.5, 3.5),
    'vd_insertion' : (-8.5, 8.5)}

default_hard_bounds = {
    'hamming_to_true_naive' : (-0.5, 19.5),
    'v_hamming_to_true_naive' : (-0.5, 8.5),
    'd_hamming_to_true_naive' : (-0.5, 10.5),
    'j_hamming_to_true_naive' : (-0.5, 12.5),
    'v_hamming_to_true_naive_normed' : (-0.5, 8.5),
    'd_hamming_to_true_naive_normed' : (-0.5, 50.5),
    'j_hamming_to_true_naive_normed' : (-0.5, 20),
    'd_3p_del' : (-1, 15),
    'd_5p_del' : (-1, 18),
    'dj_insertion' : (-1, 13),
    'jf_insertion' : (-1, 13),
    'fv_insertion' : (-1, 13),
    'j_5p_del' : (-1, 15),
    'all-mean-freq' : (0.0, 0.25),  # NOTE make sure you know where the decimal place is here!
    'v-mean-freq' : (0.0, 0.25),  # NOTE make sure you know where the decimal place is here!
    'd-mean-freq' : (0.0, 0.4),  # NOTE make sure you know where the decimal place is here!
    'j-mean-freq' : (0.0, 0.3),  # NOTE make sure you know where the decimal place is here!
    'v_3p_del' : (-1, 6),
    'vd_insertion' : (-1, 15)}

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
            if 'binerror' in line:
                bin_errors.append(float(line['binerror']) * math.sqrt(refactor))
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
        finalhists.append(TH1D(original_hist.GetName(), original_hist.GetTitle(), len(all_labels), 0.5, len(all_labels) + 0.5))
        for ilabel in range(len(all_labels)):
            label = all_labels[ilabel]
            ibin = ilabel + 1  # root conventions
            if label in hmap:
                finalhists[-1].SetBinContent(ibin, hmap[label][0])
                finalhists[-1].SetBinError(ibin, hmap[label][1])
            else:
                finalhists[-1].SetBinContent(ibin, 0.0)
                finalhists[-1].SetBinError(ibin, 0.0)

    return finalhists

# ----------------------------------------------------------------------------------------
def draw(hist, var_type, log='', plotdir=None, plotname='foop', more_hists=None, write_csv=False, stats='', bounds=None,
         errors=False, shift_overflows=False, csv_fname=None, scale_errors=None, rebin=None, plottitle='',
         colors=None, linestyles=None, unify_bin_labels=False, cwidth=700, cheight=600, imagetype='svg', xtitle=None, ytitle=None,
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

    hists = [hist,]
    if more_hists != None:
        hists = hists + more_hists

    if unify_bin_labels:
        hists = add_bin_labels_not_in_all_hists(hists)

    xmin, xmax, ymax = None, None, None
    for htmp in hists:
        if rebin is not None:
            htmp.Rebin(rebin)
        if scale_errors is not None:
            for ibin in range(htmp.GetNbinsX()+2):
                htmp.SetBinError(ibin, htmp.GetBinError(ibin) * scale_errors)

        if xmin == None or htmp.GetBinLowEdge(1) < xmin:
            xmin = htmp.GetBinLowEdge(1)
        if xmax == None or htmp.GetXaxis().GetBinUpEdge(htmp.GetNbinsX()) > xmax:
            xmax = htmp.GetXaxis().GetBinUpEdge(htmp.GetNbinsX())
        if ymax == None or htmp.GetMaximum() > ymax:
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
        hists[varname] = make_hist_from_bin_entry_file(fname, histname + '-csv-' + varname, normalize=normalize, rescale_entries=rescale_entries)
        hists[varname].SetTitle(histname)
    if len(hists) == 0:
        print 'ERROR no csvs in',dirname
        sys.exit()
    return hists

# ----------------------------------------------------------------------------------------
def compare_directories(outdir, dirs, names, xtitle='', use_hard_bounds='', stats='', errors=True, scale_errors=None, rebin=None, colors=None, linestyles=None, plot_performance=False,
                        cyst_positions=None, tryp_positions=None, leaves_per_tree=None):
    """ read all the histograms stored as .csv files in dir1 and dir2, and for those with counterparts overlay them on a new plot """
    utils.prep_dir(outdir + '/plots', multilings=['*.png', '*.svg', '*.csv'])
    if leaves_per_tree is not None:
        assert len(leaves_per_tree) == len(dirs)
    hists = []
    for idir in range(len(dirs)):
        rescale_entries = leaves_per_tree[idir] if leaves_per_tree is not None else None
        hists.append(get_hists_from_dir(dirs[idir] + '/plots', names[idir], rescale_entries=rescale_entries, normalize=True))


    histmisses = []
    names, means, sems, normalized_means = [], [], [], []
    for varname, hist in hists[0].iteritems():
        if 'hamming_to_true_naive' in varname:
            hist.GetXaxis().SetTitle('hamming')
        if 'hamming_to_true_naive_normed' in varname:
            hist.GetXaxis().SetTitle('% hamming')

        meanvals, semvals = [hist.GetMean(),], [hist.GetMeanError(),]  # means for each human for <varname> (and standard error on each mean)
        sum_total, total_entries = hist.GetMean() * hist.GetEntries(), hist.GetEntries()  # sums of weights for taking the weighted average

        more_hists = []
        missing_hist = False
        for idir in range(1, len(dirs)):
            try:
                more_hists.append(hists[idir][varname])

                meanvals.append(hists[idir][varname].GetMean())
                semvals.append(hists[idir][varname].GetMeanError())
                sum_total += hists[idir][varname].GetMean() * hists[idir][varname].Integral()
                total_entries += hists[idir][varname].Integral()
            except:
                # print 'WERRING skipping plot with missing hist (%s in %s but not %s)' % (varname, dirs[0], dirs[idir])
                histmisses.append(varname)
                missing_hist = True
                break
        if missing_hist:
            continue        

        # first find the mean over humans
        mean_of_means = sum_total / total_entries
        # then "normalize" each human's mean by this mean over humans, and that human's variance
        normalized_meanvals = []
        for im in range(len(meanvals)):
            if semvals[im] > 0:
                normalized_meanvals.append((meanvals[im] - mean_of_means) / semvals[im])
            else:
                normalized_meanvals.append(0)
        names.append(varname)
        means.append(meanvals)
        sems.append(semvals)
        normalized_means.append(normalized_meanvals)

        # print '%s %.2f' % (varname, mean_of_means)
        # for i in range(len(dirs)):
        #     print '   %.2f %.2f %.2f %.2f' % (hists[i][varname].GetMean(), hists[i][varname].GetMeanError(), hists[i][varname].Integral(), meanvals[i]),
        # print ''

        var_type = 'int'
        if hist.GetXaxis().GetBinLabel(1) != '':
            var_type = 'bool'
        bounds = None
        if plot_performance:
            # use_hard_bounds == 'true_vs_inferred':
            if varname in true_vs_inferred_hard_bounds:
                bounds = true_vs_inferred_hard_bounds[varname]
        else:
            if varname in default_hard_bounds:
                bounds = default_hard_bounds[varname]
            if '_insertion' in varname and 'content' not in varname:  # subsetting by gene now, so the above line doesn't always work
                tmpname = varname[ varname.find('_insertion') - 2 : ]
                bounds = default_hard_bounds[tmpname]
        unify_bin_labels = False
        extrastats = ''
        if '_gene' in varname:
            if hist.GetNbinsX() == 2:
                extrastats = ' 0-bin'
            else:
                unify_bin_labels = True
        plottitle = plot_titles[varname] if varname in plot_titles else ''
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
                base_plottitle = plot_titles[base_varname] if base_varname in plot_titles else ''
                plottitle = gene + ' -- ' + base_plottitle
        draw(hist, var_type, plotname=varname, plotdir=outdir, more_hists=more_hists, write_csv=False, stats=stats + ' ' + extrastats, bounds=bounds,
             shift_overflows=False, errors=errors, scale_errors=scale_errors, rebin=rebin, plottitle=plottitle, colors=colors, linestyles=linestyles, unify_bin_labels=unify_bin_labels,
             xtitle=xtitle, xline=xline, draw_str=draw_str)

    if len(histmisses) > 0:
        print 'WARNING: missing hists for %s' % ' '.join(histmisses)

    with opener('w')(outdir + '/plots/means.csv') as meanfile:
        writer = csv.DictWriter(meanfile, ('name', 'means', 'sems', 'normalized-means'))
        writer.writeheader()
        for ivar in range(len(means)):
            writer.writerow({
                'name':names[ivar],
                'means':':'.join([str(m) for m in means[ivar]]),
                'sems':':'.join([str(s) for s in sems[ivar]]),
                'normalized-means':':'.join([str(nm) for nm in normalized_means[ivar]])
            })
    # meanlist, variancelist = [], []
    # for ibin in range(0, len(means)):
    #     # print '%5.2f %5.2f' % (numpy.mean(means[ibin]), numpy.var(means[ibin]))
    #     meanlist.append(numpy.mean(means[ibin]))
    #     variancelist.append(numpy.var(means[ibin]))

    # # ----------------------------------------------------------------------------------------
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
    # gridsize=30
    # pyplot.subplot(111)
    # pyplot.hexbin(meanlist, variancelist, gridsize=gridsize, cmap=matplotlib.cm.jet, bins=None)
    # cb = pyplot.colorbar()
    # cb.set_label('mean value')

    flatmeans = []
    for mv in means:
        flatmeans += mv
    n, bins, patches = pyplot.hist(flatmeans)
    pyplot.xlabel(r'$(x - \mu) / \sigma$')
    # pyplot.ylabel('')
    # pyplot.title('Histogram of IQ')
    pyplot.text(-5, 15, r'$\sigma=' + str(math.sqrt(numpy.var(flatmeans))) + '$')
    # pyplot.axis([-5, 5, 0, 25])

    pyplot.savefig(outdir + '/plots/means.png')
    # ----------------------------------------------------------------------------------------

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
