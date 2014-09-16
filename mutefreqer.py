#!/usr/bin/env python

import sys
import os
import csv
import math
sys.argv.append('-b')
from ROOT import TH1F, TCanvas, kRed, gROOT, TLine
gROOT.Macro("plotting/root/MitStyleRemix.cc+")
from scipy.stats import beta
from opener import opener
import utils

# TODO read versions from reference file to check

def fraction_uncertainty(obs, total):
    """ Return uncertainty on the ratio n / total """
    assert obs <= total
    if total == 0.0:
        return 0.0
    lo = beta.ppf(1./6, 1 + obs, 1 + total - obs)
    hi = beta.ppf(1. - 1./6, 1 + obs, 1 + total - obs)
    if float(obs) / total < lo:  # if k/n very small (probably zero), take a one-sided c.i. with 2/3 the mass
        lo = 0.
        hi = beta.ppf(2./3, 1 + obs, 1 + total - obs)
    if float(obs) / total > hi:  # same deal if k/n very large (probably one)
        lo = beta.ppf(1./3, 1 + obs, 1 + total - obs)
        hi = 1.
    return (lo,hi)

class MuteFreqReader(object):
    def __init__(self, inputdir, human, naivety, imax=-1):
        self.human = human
        self.naivety = naivety
        self.freqs = {}
        infname = inputdir + '/' + self.human + '/' + self.naivety + '/mute-counts.csv.bz2'
        print ' opening ',infname
        with opener('r')(infname) as infile:
            reader = csv.DictReader(infile)
            il = 0
            for line in reader:
                il += 1
                assert line['subject'] == self.human
                gene_name = line['reference']
                if gene_name not in self.freqs:
                    self.freqs[gene_name] = {}
                assert utils.maturity_to_naivety(line['subset']) == self.naivety
                position = int(line['position'])
                assert position not in self.freqs[gene_name]
                self.freqs[gene_name][position] = {}
                self.freqs[gene_name][position]['ref'] = line['ref_base']
                self.freqs[gene_name][position]['n_reads'] = int(line['n_reads'])
                # assert line['N'] == ''
                for nuke in utils.nukes:
                    self.freqs[gene_name][position][nuke] = float(line[nuke]) / int(line['n_reads'])
                if imax > 0 and il > imax:
                    break
    
    def write_freqs(self, baseplotdir, baseoutdir, total_frequency=True, only_gene_name='', calculate_uncertainty=True):
        cvn = TCanvas("cvn", "", 1700, 600)
        for gene_name in self.freqs:
            if only_gene_name != '' and gene_name != only_gene_name:
                continue
            print '  %-20s' % (gene_name)
            mute_freqs = self.freqs[gene_name]
            sorted_positions = sorted(mute_freqs)

            # calculate mute freq and its uncertainty
            for position in sorted_positions:
                n_conserved, n_mutated = 0, 0
                total = mute_freqs[position]['n_reads']
                for nuke in utils.nukes:
                    obs = int(round(mute_freqs[position][nuke] * total))
                    if nuke == mute_freqs[position]['ref']:
                        n_conserved += obs
                    else:
                        n_mutated += obs
                    # uncert = fraction_uncertainty(obs, total)  # uncertainty for each nuke
                assert n_mutated + n_conserved == total
                mute_freqs[position]['mute_freq'] = float(n_mutated) / total
                mutated_fraction_err = (0.0, 0.0)
                if calculate_uncertainty:  # it's kinda slow
                    mutated_fraction_err = fraction_uncertainty(n_mutated, total)
                mute_freqs[position]['mute_freq_lo_err'] = mutated_fraction_err[0]
                mute_freqs[position]['mute_freq_hi_err'] = mutated_fraction_err[1]

            # write to csv
            outdir = baseoutdir + '/' + self.human + '/' + self.naivety + '/mute-freqs'
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            outfname = outdir +  '/' + utils.sanitize_name(gene_name) + '.csv'
            # TODO there's kind of starting to be a lot of differenct scripts producing inputs for recombinator. I should unify them
            with opener('w')(outfname) as outfile:  # write out mutation freqs for use by recombinator
                outfile.write('position,mute_freq,lo_err,hi_err\n')
                for position in sorted_positions:
                    outfile.write('%d,%f,%f,%f\n' % (position, mute_freqs[position]['mute_freq'], mute_freqs[position]['mute_freq_lo_err'],mute_freqs[position]['mute_freq_hi_err']))
                
            # and make a plot
            hist = TH1F('hist_' + utils.sanitize_name(gene_name), '',
                        sorted_positions[-1] - sorted_positions[0] + 1,
                        sorted_positions[0] - 0.5, sorted_positions[-1] + 0.5)
            lo_err_hist = TH1F(hist)
            hi_err_hist = TH1F(hist)
            for position in sorted_positions:
                hist.SetBinContent(hist.FindBin(position), mute_freqs[position]['mute_freq'])
                lo_err_hist.SetBinContent(hist.FindBin(position), mute_freqs[position]['mute_freq_lo_err'])
                hi_err_hist.SetBinContent(hist.FindBin(position), mute_freqs[position]['mute_freq_hi_err'])
            hframe = TH1F(hist)
            hframe.SetTitle(gene_name + ';;')
            hframe.Reset()
            hframe.SetMinimum(lo_err_hist.GetMinimum() - 0.03)
            hframe.SetMaximum(1.1*hi_err_hist.GetMaximum())
            hframe.Draw('')
            line = TLine(hist.GetXaxis().GetXmin(), 0., hist.GetXaxis().GetXmax(), 0.)
            line.SetLineColor(0)
            line.Draw()  # can't figure out how to convince hframe not to draw a horizontal line at y=0, so... cover it up
            hist.SetLineColor(419)
            hist.Draw('same')
            lo_err_hist.SetLineColor(kRed+2)
            hi_err_hist.SetLineColor(kRed+2)
            lo_err_hist.SetMarkerColor(kRed+2)
            hi_err_hist.SetMarkerColor(kRed+2)
            lo_err_hist.SetMarkerStyle(22)
            hi_err_hist.SetMarkerStyle(23)
            lo_err_hist.Draw('p same')
            hi_err_hist.Draw('p same')
            plotdir = baseplotdir + '/' + self.human + '/' + self.naivety + '/plots'
            if not os.path.exists(plotdir):
                os.makedirs(plotdir)
            outfname = plotdir + '/' + utils.sanitize_name(gene_name) + '.png'
            cvn.SaveAs(outfname)

if __name__ == "__main__":
    recombinator_dir = '/home/dralph/Dropbox/work/recombinator'
    human = 'C'
    naivety = 'M'
    mfr = MuteFreqReader(recombinator_dir+ '/data/human-beings', human, naivety)
    mfr.write_freqs(os.getenv('PWD').replace('/home/dralph/Dropbox', os.getenv('www')),
                    recombinator_dir + '/data/human-beings')
