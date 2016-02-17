#!/usr/bin/env python

import sys
import os
from subprocess import check_call
import csv

import plotting

import utils
import fraction_uncertainty
import paramutils
from hist import Hist
from opener import opener

# ----------------------------------------------------------------------------------------
class MuteFreqer(object):
    def __init__(self, germline_seqs):  #, base_outdir='', base_plotdir='', write_parameters=True, plot_parameters=True):
        self.germline_seqs = germline_seqs
        self.counts, self.freqs, self.plotting_info = {}, {}, {}
        n_bins, xmin, xmax = 100, 0.0, 0.5
        self.mean_rates = {'all':Hist(n_bins, xmin, xmax)}
        for region in utils.regions:
            self.mean_rates[region] = Hist(n_bins, xmin, xmax)
        self.finalized = False
        
    # ----------------------------------------------------------------------------------------
    def increment(self, info):
        # first do overall mute freqs
        freq = utils.get_mutation_rate(self.germline_seqs, info)
        self.mean_rates['all'].fill(freq)

        # then per-region stuff
        for region in utils.regions:
            # per-region mean freqs
            freq = utils.get_mutation_rate(self.germline_seqs, info, restrict_to_region=region)
            self.mean_rates[region].fill(freq)

            # per-gene per-position
            if info[region + '_gene'] not in self.counts:
                self.counts[info[region + '_gene']] = {}
            mute_counts = self.counts[info[region + '_gene']]  # temporary variable to avoid long dict access
            germline_seq = info[region + '_gl_seq']
            query_seq = info[region + '_qr_seq']
            assert len(germline_seq) == len(query_seq)
            for inuke in range(len(germline_seq)):
                i_germline = inuke + int(info[region + '_5p_del'])  # account for left-side deletions in the indexing
                if germline_seq[inuke] in utils.ambiguous_bases or query_seq[inuke] in utils.ambiguous_bases:
                    continue
                if i_germline not in mute_counts:  # if we have not yet observed this position in a query sequence, initialize it
                    mute_counts[i_germline] = {n : 0 for n in utils.nukes + ['total', ]}
                    mute_counts[i_germline]['gl_nuke'] = germline_seq[inuke]
                mute_counts[i_germline]['total'] += 1
                mute_counts[i_germline][query_seq[inuke]] += 1

                # if i_germline not in tmpmute_counts:  # if we have not yet observed this position in a query sequence, initialize it
                #     tmpmute_counts[i_germline] = {}
                #     tmpmute_counts[i_germline]['muted'] = Hist(20, 0.0, 0.3, True)
                #     tmpmute_counts[i_germline]['total'] = Hist(20, 0.0, 0.3, True)
                # freq = utils.get_mutation_rate(self.germline_seqs, info)
                # if query_seq[inuke] != germline_seq[inuke]:
                #     tmpmute_counts[i_germline]['muted'].fill(freq)
                # tmpmute_counts[i_germline]['total'].fill(freq)

    # ----------------------------------------------------------------------------------------
    def finalize(self, calculate_uncertainty=True):
        """ convert from counts to mut freqs """
        assert not self.finalized

        self.n_cached, self.n_not_cached = 0, 0
        for gene in self.counts:
            self.freqs[gene], self.plotting_info[gene] = {}, []
            # NOTE <counts> hold the overall (not per-base) frequencies, while <freqs> holds the per-base frequencies
            counts, freqs, plotting_info = self.counts[gene], self.freqs[gene], self.plotting_info[gene]
            sorted_positions = sorted(counts)
            for position in sorted_positions:
                freqs[position] = {}
                plotting_info.append({})
                plotting_info[-1]['name'] = utils.sanitize_name(gene) + '_' + str(position)
                plotting_info[-1]['nuke_freqs'] = {}
                n_conserved, n_mutated = 0, 0
                for nuke in utils.nukes:
                    nuke_freq = float(counts[position][nuke]) / counts[position]['total']
                    freqs[position][nuke] = nuke_freq
                    plotting_info[-1]['nuke_freqs'][nuke] = nuke_freq
                    if calculate_uncertainty:  # it's kinda slow
                        errs = fraction_uncertainty.err(counts[position][nuke], counts[position]['total'])
                        if errs[2]:
                            self.n_cached += 1
                        else:
                            self.n_not_cached += 1
                        # print nuke_freq, errs[0], errs[1], '(', counts[position][nuke], ',', counts[position]['total'], ')'
                        assert errs[0] <= nuke_freq  # these checks are probably unnecessary. EDIT and totally saved my ass about ten minutes after writing the previous statement
                        assert nuke_freq <= errs[1]
                        freqs[position][nuke + '_lo_err'] = errs[0]
                        freqs[position][nuke + '_hi_err'] = errs[1]

                    if nuke == counts[position]['gl_nuke']:
                        n_conserved += counts[position][nuke]
                    else:
                        n_mutated += counts[position][nuke]  # sum over A,C,G,T
                    # uncert = fraction_uncertainty.err(obs, total)  # uncertainty for each nuke
                counts[position]['freq'] = float(n_mutated) / counts[position]['total']
                mutated_fraction_err = (0.0, 0.0)
                if calculate_uncertainty:  # it's kinda slow
                    mutated_fraction_err = fraction_uncertainty.err(n_mutated, counts[position]['total'])
                    if mutated_fraction_err[2]:
                        self.n_cached += 1
                    else:
                        self.n_not_cached += 1
                counts[position]['freq_lo_err'] = mutated_fraction_err[0]
                counts[position]['freq_hi_err'] = mutated_fraction_err[1]

        self.mean_rates['all'].normalize()  # we expect overflows in mute freq hists, so no need to warn us
        for region in utils.regions:
            self.mean_rates[region].normalize()

        self.finalized = True

    # ----------------------------------------------------------------------------------------
    def write(self, base_outdir, mean_freq_outfname):
        if not self.finalized:
            self.finalize()

        outdir = base_outdir + '/mute-freqs'
        utils.prep_dir(outdir, '*.csv')

        for gene in self.counts:
            counts, freqs, plotting_info = self.counts[gene], self.freqs[gene], self.plotting_info[gene]
            sorted_positions = sorted(counts)
            outfname = outdir + '/' + utils.sanitize_name(gene) + '.csv'
            with opener('w')(outfname) as outfile:
                nuke_header = []
                for nuke in utils.nukes:
                    nuke_header.append(nuke)
                    nuke_header.append(nuke + '_obs')
                    nuke_header.append(nuke + '_lo_err')
                    nuke_header.append(nuke + '_hi_err')
                writer = csv.DictWriter(outfile, ('position', 'mute_freq', 'lo_err', 'hi_err') + tuple(nuke_header))
                writer.writeheader()
                for position in sorted_positions:
                    row = {'position':position,
                           'mute_freq':counts[position]['freq'],
                           'lo_err':counts[position]['freq_lo_err'],
                           'hi_err':counts[position]['freq_hi_err']}
                    for nuke in utils.nukes:
                        row[nuke] = freqs[position][nuke]
                        row[nuke + '_obs'] = counts[position][nuke]
                        row[nuke + '_lo_err'] = freqs[position][nuke + '_lo_err']
                        row[nuke + '_hi_err'] = freqs[position][nuke + '_hi_err']
                    writer.writerow(row)

        assert 'REGION' in mean_freq_outfname
        self.mean_rates['all'].write(mean_freq_outfname.replace('REGION', 'all'))  # hackey hackey hackey replacement... *sigh*
        for region in utils.regions:
            self.mean_rates[region].write(mean_freq_outfname.replace('REGION', region))

    # ----------------------------------------------------------------------------------------
    def plot(self, base_plotdir, cyst_positions=None, tryp_positions=None, only_csv=False):
        if not self.finalized:
            self.finalize()

        plotdir = base_plotdir + '/mute-freqs'
        overall_plotdir = plotdir + '/overall'
        utils.prep_dir(overall_plotdir, multilings=('*.csv', '*.svg'))
        for region in utils.regions:
            utils.prep_dir(plotdir + '/' + region, multilings=('*.csv', '*.svg'))
            # utils.prep_dir(plotdir + '/' + region + '-per-base/plots', multilings=('*.csv', '*.png'))

        for gene in self.counts:
            counts, plotting_info = self.counts[gene], self.plotting_info[gene]
            sorted_positions = sorted(counts)
            genehist = Hist(sorted_positions[-1] - sorted_positions[0] + 1, sorted_positions[0] - 0.5, sorted_positions[-1] + 0.5, xtitle='fixme', ytitle='fixme')  #, title=utils.sanitize_name(gene))
            for position in sorted_positions:
                hi_diff = abs(counts[position]['freq'] - counts[position]['freq_hi_err'])
                lo_diff = abs(counts[position]['freq'] - counts[position]['freq_lo_err'])
                err = 0.5*(hi_diff + lo_diff)
                genehist.set_ibin(genehist.find_bin(position), counts[position]['freq'], error=err)
            xline = None
            figsize = [3, 3]
            if utils.get_region(gene) == 'v' and cyst_positions is not None:
                xline = cyst_positions[gene]
                figsize[0] *= 3.5
            elif utils.get_region(gene) == 'j' and tryp_positions is not None:
                xline = tryp_positions[gene]
                figsize[0] *= 2
            plotting.draw_no_root(genehist, plotdir=plotdir + '/' + utils.get_region(gene), plotname=utils.sanitize_name(gene), errors=True, write_csv=True, xline=xline, figsize=figsize, only_csv=only_csv)
            # paramutils.make_mutefreq_plot(plotdir + '/' + utils.get_region(gene) + '-per-base', utils.sanitize_name(gene), plotting_info)  # needs translation to mpl

        # make mean mute freq hists
        plotting.draw_no_root(self.mean_rates['all'], plotname='all-mean-freq', plotdir=overall_plotdir, stats='mean', bounds=(0.0, 0.4), write_csv=True, only_csv=only_csv)
        for region in utils.regions:
            plotting.draw_no_root(self.mean_rates[region], plotname=region+'-mean-freq', plotdir=overall_plotdir, stats='mean', bounds=(0.0, 0.4), write_csv=True, only_csv=only_csv)

        if not only_csv:  # write html file and fix permissiions
            plotting.make_html(overall_plotdir)
            for region in utils.regions:
                plotting.make_html(plotdir + '/' + region, n_columns=1)

    # ----------------------------------------------------------------------------------------
    def clean(self):
        """ remove all the parameter files """
        for gene in self.counts:
            outfname = self.outdir + '/' + utils.sanitize_name(gene) + '.csv'
            os.remove(outfname)
        os.rmdir(self.outdir)
