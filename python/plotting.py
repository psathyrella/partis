import math
import os
import glob
import sys
import csv
from array import array
from subprocess import check_call

import utils

def check_root():
    try:
        from ROOT import kBlue
        return True
    except ImportError:
        return False

has_root = check_root()
if has_root:
    sys.argv.append('-b')
    from ROOT import TH1F, TCanvas, kRed, gROOT, TLine, TLegend, kBlue
    gROOT.Macro("plotting/MitStyleRemix.cc+")
else:
    print ' ROOT not found, proceeding without plotting'

from opener import opener

# ----------------------------------------------------------------------------------------
def set_bins(values, n_bins, is_log_x, xbins, var_type='float'):
    """ NOTE <values> should be sorted """
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
        writer = csv.DictWriter(histfile, ('bin_low_edge', 'contents'))
        writer.writeheader()
        for ibin in range(hist.GetNbinsX() + 2):
            writer.writerow({'bin_low_edge':hist.GetXaxis().GetBinLowEdge(ibin), 'contents':hist.GetBinContent(ibin)})

# ----------------------------------------------------------------------------------------
def make_hist_from_bin_entry_file(fname, hist_label='', log=''):
    """ return root histogram with each bin low edge and bin content read from <fname> """
    low_edges = []
    contents = []
    with opener('r')(fname) as infile:
        reader = csv.DictReader(infile)
        for line in reader:
            low_edges.append(float(line['bin_low_edge']))
            contents.append(float(line['contents']))

    n_bins = len(low_edges) - 2  # file should have a line for the under- and overflow bins
    xbins = array('f', [0.0 for i in range(n_bins+1)])  # NOTE has to be n bins *plus* 1
    low_edges = sorted(low_edges)
    for ib in range(n_bins+1):
        xbins[ib] = low_edges[ib+1]  # low_edges[1] is the lower edge of the first bin, i.e. the first bin after the underflow bin, and this will set the last entry in xbins to lower[n_bins+1], i.e. the lower edge of the overflow bin. Which, I bloody well think, is correct
    hist = TH1F(hist_label, '', n_bins, xbins)  # this will barf if the csv file wasn't sorted by bin low edge
    for ib in range(n_bins+2):
        hist.SetBinContent(ib, contents[ib])

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
def make_hist(values, var_type, hist_label, log='', xmin_force=0.0, xmax_force=0.0, normalize=False):
    """ fill a histogram with values from a dictionary """
    if not has_root:
        return

    bin_labels = sorted(values)
    if var_type == 'string':  # for strings, sort so most common value is to left side
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
        set_bins(bin_labels, n_bins, 'x' in log, xbins, var_type)
        hist = TH1F(hist_label, '', n_bins, xbins)
    else:
      hist = TH1F(hist_label, '', n_bins, xmin_force, xmax_force)
    hist.Sumw2()

    for ival in range(len(values)):
        if var_type == 'string':
            label = bin_labels[ival]
            if 'IGH' in label:
                label = label.replace('IGH','').replace('*',' ').lower()
            hist.GetXaxis().SetBinLabel(ival+1, label)
            hist.SetBinContent(ival+1, values[bin_labels[ival]])
        else:
            hist.Fill(bin_labels[ival], values[bin_labels[ival]])
        
    # print '  %d (mean %f) values for %s' % (len(values), hist.GetMean(), 'what?')
  
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
def draw(hist, var_type, log='', plotdir=os.getenv('www'), plotname='foop', hist2=None, write_csv=True):
    if not has_root:
        return
    cvn = TCanvas('cvn', '', 700, 600)
    xmin = hist.GetBinLowEdge(1)
    xmax = hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX())
    if hist2 != None:
        xmin = min(xmin, hist2.GetBinLowEdge(1))
        xmax = max(xmax, hist2.GetXaxis().GetBinUpEdge(hist2.GetNbinsX()))
    hframe = TH1F('hframe', '', hist.GetNbinsX(), xmin, xmax)
    if var_type == 'string':
        for ib in range(1, hframe.GetNbinsX()+1):
            hframe.GetXaxis().SetBinLabel(ib, hist.GetXaxis().GetBinLabel(ib))

    ymax = hist.GetMaximum()
    if hist2 != None:
        ymax = max(ymax, hist2.GetMaximum())
    hframe.SetMaximum(1.35*ymax)
    hframe.SetTitle(plotname + ';' + hist.GetXaxis().GetTitle() + ';')
    hframe.Draw("txt")

    hist.SetLineColor(kBlue);
    # hist.SetMarkerSize(0);
    hist.SetLineWidth(4);
    hist.Draw("hist same");
    if hist2 != None:
        hist2.SetLineColor(kRed+2)  #419);
        # hist2.SetMarkerSize(0);
        hist2.SetLineWidth(4);
        hist2.Draw("hist same");

    leg = TLegend(0.6, 0.75, 0.9, 0.9)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.AddEntry(hist, hist.GetTitle() , 'l')
    if hist2 != None:
        leg.AddEntry(hist2, hist2.GetTitle() , 'l')
    leg.Draw()

    cvn.SetLogx('x' in log)
    cvn.SetLogy('y' in log)
    if not os.path.exists(plotdir + '/plots'):
        print 'ERROR dir \'' + plotdir + '/plots\' d.n.e.'
        assert False

    if write_csv:
        write_hist_to_file(plotdir + '/plots/' + plotname + '.csv', hist)
    cvn.SaveAs(plotdir + '/plots/' + plotname + '.svg')

# ----------------------------------------------------------------------------------------
def compare_directories(dir1, name1, dir2, name2, xtitle=''):
    """ read all the histograms stored as .csv files in dir1 and dir2, and for those with counterparts overlay them on a new plot """
    plotdir = os.getenv('www') + '/partis/compare-performance'
    utils.prep_dir(plotdir + '/plots', '*.svg')
    dir1_hists = {}
    for fname in glob.glob(dir1 + '/*.csv'):
        varname = os.path.basename(fname).replace('.csv', '')
        dir1_hists[varname] = make_hist_from_bin_entry_file(fname, name1 + '-csv-' + varname)
        dir1_hists[varname].SetTitle(name1)
    if len(dir1_hists) == 0:
        print 'ERROR no hists in',dir1
        sys.exit()
    dir2_hists = {}
    for fname in glob.glob(dir2 + '/*.csv'):
        varname = os.path.basename(fname).replace('.csv', '')
        dir2_hists[varname] = make_hist_from_bin_entry_file(fname, name2 + '-csv-' + varname)
        dir2_hists[varname].SetTitle(name2)
    if len(dir2_hists) == 0:
        print 'ERROR no hists in',dir2
        sys.exit()
    for varname, hist in dir1_hists.iteritems():
        hist.GetXaxis().SetTitle(xtitle)
        if varname not in dir2_hists:
            print 'WARNING %s not found in %s' % (varname, dir2)
            continue
        hist2 = dir2_hists[varname]
        draw(hist, 'int', plotname=varname, plotdir=plotdir, hist2=hist2, write_csv=False)
    check_call(['./permissify-www', plotdir])  # NOTE this should really permissify starting a few directories higher up
    check_call(['makeHtml', plotdir, '3', 'null', 'svg'])
        
