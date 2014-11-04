#!/usr/bin/env python

import sys
import os
from subprocess import check_call
import csv

import plotting
has_root = plotting.has_root
if has_root:
    from ROOT import TCanvas, TH1F, TLine, kRed

import utils
from opener import opener

# ----------------------------------------------------------------------------------------
class MuteFreqer(object):
    def __init__(self, base_outdir, base_plotdir, germline_seqs):
        self.outdir = base_outdir + '/mute-freqs'
        self.base_plotdir = base_plotdir
        if self.base_plotdir != '':
            self.base_plotdir += '/mute-freqs'
            for region in utils.regions:
                utils.prep_dir(self.base_plotdir + '/' + region + '/plots', '*.svg')
        utils.prep_dir(self.outdir, '*.csv')
        self.germline_seqs = germline_seqs
        self.counts = {}
    
    # ----------------------------------------------------------------------------------------
    def increment(self, info):
        for region in utils.regions:
            if info[region + '_gene'] not in self.counts:
                self.counts[info[region + '_gene']] = {}
            mute_counts = self.counts[info[region + '_gene']]  # temporary variable to avoid long dict access
            germline_seq = info[region + '_gl_seq']
            query_seq = info[region + '_qr_seq']
            # utils.color_mutants(germline_seq, query_seq, print_result=True, extra_str='  ')
            assert len(germline_seq) == len(query_seq)
            for inuke in range(len(germline_seq)):
                i_germline = inuke + int(info[region + '_5p_del'])  # account for left-side deletions in the indexing
                if i_germline not in mute_counts:  # if we have not yet observed this position in a query sequence, initialize it
                    mute_counts[i_germline] = {'A':0, 'C':0, 'G':0, 'T':0, 'total':0, 'gl_nuke':germline_seq[inuke]}
                mute_counts[i_germline]['total'] += 1
                mute_counts[i_germline][query_seq[inuke]] += 1

    # ----------------------------------------------------------------------------------------
    def write(self, calculate_uncertainty=True):
        cvn = None
        if has_root:
            cvn = TCanvas("cvn", "", 1700, 600)
        for gene in self.counts:
            mute_counts = self.counts[gene]
            sorted_positions = sorted(mute_counts)

            # calculate mute freq and its uncertainty
            for position in sorted_positions:
                # sum over A,C,G,T (TODO don't sum over them)
                n_conserved, n_mutated = 0, 0
                for nuke in utils.nukes:
                    if nuke == mute_counts[position]['gl_nuke']:
                        n_conserved += mute_counts[position][nuke]
                    else:
                        n_mutated += mute_counts[position][nuke]
                    # uncert = utils.fraction_uncertainty(obs, total)  # uncertainty for each nuke
                mute_counts[position]['freq'] = float(n_mutated) / mute_counts[position]['total']
                mutated_fraction_err = (0.0, 0.0)
                if calculate_uncertainty:  # it's kinda slow
                    mutated_fraction_err = utils.fraction_uncertainty(n_mutated, mute_counts[position]['total'])
                mute_counts[position]['freq_lo_err'] = mutated_fraction_err[0]
                mute_counts[position]['freq_hi_err'] = mutated_fraction_err[1]


            # TODO there's kind of starting to be a lot of differenct scripts producing inputs for recombinator. I should unify them

            # write to csv
            outfname = self.outdir + '/' + utils.sanitize_name(gene) + '.csv'
            with opener('w')(outfname) as outfile:
                writer = csv.DictWriter(outfile, ('position', 'mute_freq', 'lo_err', 'hi_err'))
                writer.writeheader()
                for position in sorted_positions:
                    row = {'position':position, 'mute_freq':mute_counts[position]['freq'], 'lo_err':mute_counts[position]['freq_lo_err'], 'hi_err':mute_counts[position]['freq_hi_err']}
                    writer.writerow(row)
                
            if has_root: # make a plot
                hist = TH1F('hist_' + utils.sanitize_name(gene), '',
                            sorted_positions[-1] - sorted_positions[0] + 1,
                            sorted_positions[0] - 0.5, sorted_positions[-1] + 0.5)
                lo_err_hist = TH1F(hist)
                hi_err_hist = TH1F(hist)
                for position in sorted_positions:
                    hist.SetBinContent(hist.FindBin(position), mute_counts[position]['freq'])
                    lo_err_hist.SetBinContent(hist.FindBin(position), mute_counts[position]['freq_lo_err'])
                    hi_err_hist.SetBinContent(hist.FindBin(position), mute_counts[position]['freq_hi_err'])
                hframe = TH1F(hist)
                hframe.SetTitle(gene + ';;')
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
                if self.base_plotdir != '':
                    plotfname = self.base_plotdir + '/' + utils.get_region(gene) + '/plots/' + utils.sanitize_name(gene) + '.svg'
                    cvn.SaveAs(plotfname)

        if has_root:
            if self.base_plotdir != '':
                check_call(['./permissify-www', self.base_plotdir])  # NOTE this should really permissify starting a few directories higher up
                for region in utils.regions:
                    check_call(['makeHtml', self.base_plotdir + '/' + region, '2', 'null', 'svg'])

    # ----------------------------------------------------------------------------------------
    def clean(self):
        """ remove all the parameter files """
        for gene in self.counts:
            outfname = self.outdir + '/' + utils.sanitize_name(gene) + '.csv'
            os.remove(outfname)
        os.rmdir(self.outdir)
