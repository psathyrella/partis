import math
import os
import sys
from array import array
sys.argv.append('-b')
from ROOT import TH1F, TCanvas, kRed, gROOT, TLine, TLegend, kBlue
gROOT.Macro("plotting/MitStyleRemix.cc+")

from opener import opener

# ----------------------------------------------------------------------------------------
def set_bins(values, n_bins, is_log_x, var_type, xbins):
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
def make_hist(values, var_type, hist_label, log='', xmin_force=0.0, xmax_force=0.0, normalize=False):
    """ fill a histogram with values from a dictionary """

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
        set_bins(bin_labels, n_bins, 'x' in log, var_type, xbins)
        hist = TH1F('h'+hist_label, '', n_bins, xbins)
    else:
      hist = TH1F('h'+hist_label, '', n_bins, xmin_force, xmax_force)
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
    
    return hist;

# ----------------------------------------------------------------------------------------
def draw(hist, var_type, log='', plotdir=os.getenv('www'), plotname='foop'):
    cvn = TCanvas('cvn', '', 700, 600)
    xmin = hist.GetBinLowEdge(1)
    xmax = hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX())
    hframe = TH1F('hframe', '', hist.GetNbinsX(), xmin, xmax)
    if var_type == 'string':
        for ib in range(1, hframe.GetNbinsX()+1):
            hframe.GetXaxis().SetBinLabel(ib, hist.GetXaxis().GetBinLabel(ib))

    hframe.SetMaximum(1.35*hist.GetMaximum())
    hframe.SetTitle(plotname + ';' + hist.GetXaxis().GetTitle() + ';')
    hframe.Draw("txt")
    leg = TLegend(0.65, 0.75, 0.94, 0.9)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.AddEntry(hist, hist.GetName() , 'lp')
    # leg.Draw()

    hist.SetLineColor(kBlue);
    # hist.SetMarkerSize(0);
    hist.SetLineWidth(4);
    hist.Draw("hist same");
    cvn.SetLogx('x' in log)
    cvn.SetLogy('y' in log)
    cvn.SaveAs(plotdir + '/plots/' + plotname + '.svg')
