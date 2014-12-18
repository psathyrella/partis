import math
import os
import glob
import sys
import csv
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
    from ROOT import TH1F, TCanvas, kRed, gROOT, TLine, TLegend, kBlue, kGreen, kCyan, kOrange
    gROOT.Macro("plotting/MitStyleRemix.cc+")
else:
    print ' ROOT not found, proceeding without plotting'

from opener import opener

hard_bounds = {
    'hamming_to_true_naive' : (-0.5, 25.5),
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
    'mute_freqs' : (-5.5, 5.5),
    'v_3p_del' : (-3.5, 3.5),
    'vd_insertion' : (-8.5, 8.5)
}

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
    with opener('w')(fname) as histfile:
        writer = csv.DictWriter(histfile, ('bin_low_edge', 'contents', 'xtitle', 'binlabel'))  # this is a really crummy way of writing style information, but root files *suck*, so this is what I do for now
        writer.writeheader()
        for ibin in range(hist.GetNbinsX() + 2):
            writer.writerow({'bin_low_edge' : hist.GetXaxis().GetBinLowEdge(ibin),
                             'contents' : hist.GetBinContent(ibin),
                             'xtitle' : hist.GetXaxis().GetTitle(),
                             'binlabel' : hist.GetXaxis().GetBinLabel(ibin)})

# ----------------------------------------------------------------------------------------
def make_hist_from_bin_entry_file(fname, hist_label='', log=''):
    """ return root histogram with each bin low edge and bin content read from <fname> """
    low_edges, contents, bin_labels = [], [], []
    xtitle = ''
    with opener('r')(fname) as infile:
        reader = csv.DictReader(infile)
        for line in reader:
            low_edges.append(float(line['bin_low_edge']))
            contents.append(float(line['contents']))
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
    hist = TH1F(hist_label, '', n_bins, xbins)  # this will barf if the csv file wasn't sorted by bin low edge
    hist.GetXaxis().SetTitle(xtitle)
    for ib in range(n_bins+2):
        hist.SetBinContent(ib, contents[ib])
        if bin_labels[ib] != '':
            hist.GetXaxis().SetBinLabel(ib, bin_labels[ib])

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
    hist = TH1F(hist_label, '', n_bins, xbins)
    for value in values:
        hist.Fill(value)

    return hist

# ----------------------------------------------------------------------------------------
def make_bool_hist(n_true, n_false, hist_label):
    """ fill a two-bin histogram with the fraction false in the first bin and the fraction true in the second """
    hist = TH1F(hist_label, '', 2, -0.5, 1.5)
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
        return TH1F(hist_label, '', 1, 0, 1)
    xbins = array('f', [0 for i in range(n_bins+1)])  # NOTE has to be n_bins *plus* 1
    set_bins(values, n_bins, is_log_x=False, xbins=xbins, var_type='float')
    hist = TH1F(hist_label, '', n_bins, xbins)
    for val in values:
        hist.Fill(val)
    return hist

# ----------------------------------------------------------------------------------------
def make_hist(values, var_type, hist_label, log='', xmin_force=0.0, xmax_force=0.0, normalize=False, sort=False):
    """ Fill a histogram with values from a dictionary (each key will correspond to one bin) """
    if not has_root:
        return
    if len(values) == 0:
        print 'WARNING no values for %s in make_hist' % hist_label
        return TH1F(hist_label, '', 1, 0, 1)

    bin_labels = sorted(values)
    if not sort and var_type == 'string':  # for strings, sort so most common value is to left side
        bin_labels = sorted(values, key=values.get, reverse=True)

    if var_type == 'string':
        n_bins = len(values)
    elif var_type == 'int':
        n_bins = bin_labels[-1] - bin_labels[0] + 1
    else:
        assert False

    hist = None
    xbins = array('f', [0 for i in range(n_bins+1)])  # NOTE has to be n_bins *plus* 1
    if xmin_force == xmax_force:  # if boundaries aren't set explicitly, work them out dynamically
        if var_type == 'int':
            hist = TH1F(hist_label, '', n_bins, bin_labels[0] - 0.5, bin_labels[-1] + 0.5)
        else:
            set_bins(bin_labels, n_bins, 'x' in log, xbins, var_type)
            hist = TH1F(hist_label, '', n_bins, xbins)
    else:
      hist = TH1F(hist_label, '', n_bins, xmin_force, xmax_force)
    hist.Sumw2()
        

    for ival in range(len(values)):
        if var_type == 'string':
            label = bin_labels[ival]
            hist.GetXaxis().SetBinLabel(ival+1, label)
            hist.SetBinContent(ival+1, values[bin_labels[ival]])
            hist.SetBinError(ival+1, math.sqrt(values[bin_labels[ival]]))
        else:
            hist.Fill(bin_labels[ival], values[bin_labels[ival]])
  
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
def draw(hist, var_type, log='', plotdir=os.getenv('www'), plotname='foop', more_hists=None, write_csv=False, stats='', bounds=None, errors=False, shift_overflows=False):
    if not has_root:
        return
    cvn = TCanvas('cvn', '', 700, 600)

    hists = [hist,]
    if more_hists != None:
        hists = hists + more_hists

    xmin, xmax, ymax = None, None, None
    for htmp in hists:
        if xmin == None or htmp.GetBinLowEdge(1) < xmin:
            xmin = htmp.GetBinLowEdge(1)
        if xmax == None or htmp.GetXaxis().GetBinUpEdge(htmp.GetNbinsX()) > xmax:
            xmax = htmp.GetXaxis().GetBinUpEdge(htmp.GetNbinsX())
        if ymax == None or htmp.GetMaximum() > xmax:
            ymax = htmp.GetMaximum()
    if bounds != None:
        xmin, xmax = bounds
    hframe = TH1F('hframe', '', hist.GetNbinsX(), xmin, xmax)
    if var_type == 'string' or var_type == 'bool':
        for ib in range(1, hframe.GetNbinsX()+1):
            hframe.GetXaxis().SetBinLabel(ib, hist.GetXaxis().GetBinLabel(ib))

    hframe.SetMaximum(1.35*ymax)
    if var_type == 'bool':
        hframe.GetXaxis().SetLabelSize(0.1)
    hframe.SetTitle(plotname + ';' + hist.GetXaxis().GetTitle() + ';')
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

    colors = (kRed, kBlue-4, kGreen+2, kOrange+1)
    draw_str = 'hist same'
    if errors:  # not working!
        draw_str = 'e ' + draw_str
    for ih in range(len(hists)):
        htmp = hists[ih]
        htmp.SetLineColor(colors[ih])
        if ih == 0:
            htmp.SetMarkerSize(0)
        assert ih < 6
        htmp.SetLineWidth(6-ih)
        htmp.Draw(draw_str)

    leg = TLegend(0.57, 0.72, 0.99, 0.9)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    for htmp in hists:
        if 'rms' in stats:
            htmp.SetTitle(htmp.GetTitle() + (' (%.2f)' % htmp.GetRMS()))
        leg.AddEntry(htmp, htmp.GetTitle() , 'l')
    leg.Draw()

    cvn.SetLogx('x' in log)
    cvn.SetLogy('y' in log)
    if not os.path.exists(plotdir + '/plots'):
        print 'ERROR dir \'' + plotdir + '/plots\' d.n.e.'
        assert False

    if write_csv:
        assert more_hists == None
        write_hist_to_file(plotdir + '/plots/' + plotname + '.csv', hist)
    cvn.SaveAs(plotdir + '/plots/' + plotname + '.svg')

# ----------------------------------------------------------------------------------------
def get_hists_from_dir(dirname, histname):
    hists = {}
    for fname in glob.glob(dirname + '/*.csv'):
        varname = os.path.basename(fname).replace('.csv', '')
        hists[varname] = make_hist_from_bin_entry_file(fname, histname + '-csv-' + varname)
        hists[varname].SetTitle(histname)
    if len(hists) == 0:
        print 'ERROR no hists in',dirname
        sys.exit()
    return hists

# ----------------------------------------------------------------------------------------
def compare_directories(outdir, dirs, names, xtitle='', stats=''):
    """ read all the histograms stored as .csv files in dir1 and dir2, and for those with counterparts overlay them on a new plot """
    utils.prep_dir(outdir + '/plots', '*.svg')
    hists = []
    for idir in range(len(dirs)):
        hists.append(get_hists_from_dir(dirs[idir] + '/plots', names[idir]))
    for varname, hist in hists[0].iteritems():
        if 'hamming_to_true_naive' in varname:
            hist.GetXaxis().SetTitle('hamming')
        if 'hamming_to_true_naive_normed' in varname:
            hist.GetXaxis().SetTitle('% hamming')
        # deal with log axes
        # if varname.find('hamming_to_true_naive') > 0:
        #     log = 'y'
        # else:
        log = ''

        more_hists = []
        for idir in range(1, len(dirs)):
            more_hists.append(hists[idir][varname])
            
        var_type = 'int'
        if hist.GetXaxis().GetBinLabel(1) != '':
            var_type = 'bool'
        bounds = None
        if varname in hard_bounds:
            bounds = hard_bounds[varname]
        draw(hist, var_type, plotname=varname, plotdir=outdir, more_hists=more_hists, write_csv=False, stats=stats, bounds=bounds, log=log, shift_overflows=False)
    check_call(['./permissify-www', outdir])  # NOTE this should really permissify starting a few directories higher up
    check_call(['makeHtml', outdir, '3', 'null', 'svg'])
        
