#!/usr/bin/env python

import sys
import os
import math
from subprocess import check_call
import numpy
import csv

import plotting

import utils
import fraction_uncertainty
import paramutils
from hist import Hist
from opener import opener

# ----------------------------------------------------------------------------------------
class MuteFreqer(object):
    def __init__(self, germline_seqs, calculate_uncertainty=True):
        self.germline_seqs = germline_seqs
        self.calculate_uncertainty = calculate_uncertainty
        self.counts, self.freqs = {}, {}
        n_bins, xmin, xmax = 200, 0., 1.
        self.mean_rates = {'all':Hist(n_bins, xmin, xmax)}
        for region in utils.regions:
            self.mean_rates[region] = Hist(n_bins, xmin, xmax)
        self.finalized = False

        # tigger stuff
        self.tigger = False
        self.n_max_mutes = 20
        self.n_obs_min = 10
        self.min_y_intercept = 1./8

    # ----------------------------------------------------------------------------------------
    def increment(self, info):
        self.mean_rates['all'].fill(utils.get_mutation_rate(self.germline_seqs, info))  # mean freq over whole sequence (excluding insertions)

        for region in utils.regions:
            regional_freq, len_excluding_ambig = utils.get_mutation_rate(self.germline_seqs, info, restrict_to_region=region, return_len_excluding_ambig=True)
            n_mutes = regional_freq * len_excluding_ambig  # total number of mutations in the region (for tigger stuff)
            if abs(n_mutes - int(n_mutes)) > 1e6:
                raise Exception('n mutated %f not an integer' % n_mutes)
            n_mutes = int(n_mutes)
            self.mean_rates[region].fill(regional_freq)  # per-region mean freq

            # per-gene per-position freqs
            gene = info[region + '_gene']
            if gene not in self.counts:
                self.counts[gene] = {}
            gcounts = self.counts[gene]  # temporary variable to avoid long dict access
            germline_seq = info[region + '_gl_seq']
            query_seq = info[region + '_qr_seq']
            assert len(germline_seq) == len(query_seq)
            for ipos in range(len(germline_seq)):
                igl = ipos + int(info[region + '_5p_del'])  # account for left-side deletions in the indexing
                if germline_seq[ipos] in utils.ambiguous_bases or query_seq[ipos] in utils.ambiguous_bases:
                    continue
                if igl not in gcounts:  # if we have not yet observed this position in a query sequence, initialize it
                    gcounts[igl] = {n : 0 for n in utils.nukes + ['total', ]}
                    gcounts[igl]['gl_nuke'] = germline_seq[ipos]
                    gcounts[igl]['tigger'] = {}
                gcounts[igl]['total'] += 1
                gcounts[igl][query_seq[ipos]] += 1  # note that if <query_seq[ipos]> isn't among <utils.nukes>, this will toss a key error

                if self.tigger:
                    if igl not in gcounts:
                        gcounts[igl]['tigger'] = {}
                    if utils.get_region(gene) == 'v':
                        if n_mutes not in gcounts[igl]['tigger']:
                            gcounts[igl]['tigger'][n_mutes] = {'muted' : 0, 'total' : 0}
                        gcounts[igl]['tigger'][n_mutes]['total'] += 1
                        if query_seq[ipos] != germline_seq[ipos]:  # if this position is mutated
                            gcounts[igl]['tigger'][n_mutes]['muted'] += 1  # mark that we saw this germline position mutated once in a sequence with <n_mutes> regional mutation frequency

    # ----------------------------------------------------------------------------------------
    def get_uncertainty(self, obs, total):
        if self.calculate_uncertainty:  # it's kinda slow
            errs = fraction_uncertainty.err(obs, total)
            if errs[2]:
                self.n_cached += 1
            else:
                self.n_not_cached += 1
        else:
            errs = 0., 1.

        return errs[0], errs[1]

    # ----------------------------------------------------------------------------------------
    def tigger_calcs(self, position, gcounts, mean_x_icpt):
        iterinfo = gcounts['tigger'].items()

        obs = [d['muted'] for nm, d in iterinfo if nm < self.n_max_mutes]
        if sum(obs) < self.n_obs_min:  # ignore positions with only a few observed mutations
            return

        lohis = [fraction_uncertainty.err(d['muted'], d['total'], use_beta=True) for nm, d in iterinfo if nm < self.n_max_mutes]
        errs = [(hi - lo) / 2 for lo, hi, _ in lohis]
        weights = [1./(e*e) for e in errs]

        freqs = [float(d['muted']) / d['total'] for nm, d in iterinfo if nm < self.n_max_mutes]
        total = [d['total'] for nm, d in iterinfo if nm < self.n_max_mutes]

        # for i in range(len(freqs)):
        #     print '  %3d / %3d = %6.2f    %6.2f    %6.2f' % (obs[i], total[i], freqs[i], errs[i], weights[i])
    
        n_mutelist = [nm for nm in gcounts['tigger'].keys() if nm < self.n_max_mutes]

        params, cov = numpy.polyfit(n_mutelist, freqs, 1, w=weights, cov=True)
        slope, slope_err = params[0], math.sqrt(cov[0][0])
        y_icpt, y_icpt_err = params[1], math.sqrt(cov[1][1])
    
        interesting = False
        if y_icpt + y_icpt_err < 1./8:
            x_icpt, x_icpt_err = -y_icpt / slope, abs(y_icpt / slope) * math.sqrt((y_icpt_err/y_icpt)**2 + (slope_err/slope)**2)
            mean_x_icpt['sum'] += x_icpt / x_icpt_err
            mean_x_icpt['total'] += 1. / x_icpt_err
        else:
            x_icpt, x_icpt_err = 0, 0
            interesting = True
        print_str = '   %3d   %9.3f +/- %-9.3f   %9.3f +/- %-9.3f   %7.4f +/- %7.4f      %3d / %3d' % (position, x_icpt, x_icpt_err, y_icpt, y_icpt_err, slope, slope_err, sum(obs), sum(total))
        if interesting:
            print_str = utils.color('red', print_str)
        print print_str

        # for testing it's easier to make plots here:
        # plotinfo = {'n_muted' : n_mutelist, 'freqs' : freqs, 'errs' : errs, 'slope' : slope, 'intercept' : y_icpt}
        # plotting.make_tigger_plot('IGHVX', position, plotinfo)
        return plotinfo

    # ----------------------------------------------------------------------------------------
    def finalize_tigger(self):
        utils.prep_dir(os.getenv('www') + '/partis/tmp', wildling='*.svg')
        for gene in self.counts:
            if utils.get_region(gene) != 'v':
                continue
            print '\n%s' % gene
            print ' position         x-icpt                  y-icpt                   slope              mut / total'
            mean_x_icpt = {'sum' : 0., 'total' : 0.}
            for position in sorted(self.counts[gene].keys()):
                self.freqs[gene][position]['tigger'] = self.tigger_calcs(position, self.counts[gene][position], mean_x_icpt)
            print mean_x_icpt
            if mean_x_icpt['total'] > 0.:
                print mean_x_icpt['sum'] / mean_x_icpt['total']
        assert False
        for gene in self.freqs:
            if utils.get_region(gene) != 'v':
                continue
            info = {p : self.freqs[gene][p]['tigger-fits'] for p in self.freqs[gene]}
            x_intercepts = [-v['intercept'] / v['slope'] for k, v in info.items() if v['intercept'] is not None and v['intercept'] < 0.3]
            print sorted(x_intercepts)
            print sum(x_intercepts) / float(len(x_intercepts))
            print numpy.median(x_intercepts)

    # ----------------------------------------------------------------------------------------
    def finalize(self):
        """ convert from counts to mut freqs """
        assert not self.finalized

        self.n_cached, self.n_not_cached = 0, 0
        for gene in self.counts:
            gcounts = self.counts[gene]
            freqs = {position : {} for position in gcounts}
            for position in sorted(gcounts.keys()):
                n_conserved, n_mutated = 0, 0
                for nuke in utils.nukes:
                    ncount, total = gcounts[position][nuke], gcounts[position]['total']
                    nuke_freq = float(ncount) / total
                    freqs[position][nuke] = nuke_freq
                    freqs[position][nuke + '_lo_err'], freqs[position][nuke + '_hi_err'] = self.get_uncertainty(ncount, total)
                    if nuke == gcounts[position]['gl_nuke']:
                        n_conserved += ncount
                    else:
                        n_mutated += ncount  # sum over A,C,G,T
                freqs[position]['freq'] = float(n_mutated) / total
                freqs[position]['freq_lo_err'], freqs[position]['freq_hi_err'] = self.get_uncertainty(n_mutated, total)

            self.freqs[gene] = freqs

        for hist in self.mean_rates.values():
            hist.normalize()

        if self.tigger:
            self.finalize_tigger()

        self.finalized = True

    # ----------------------------------------------------------------------------------------
    def write(self, base_outdir, mean_freq_outfname):
        if not self.finalized:
            self.finalize()

        outdir = base_outdir + '/mute-freqs'
        utils.prep_dir(outdir, '*.csv')

        for gene in self.counts:
            gcounts, freqs = self.counts[gene], self.freqs[gene]
            outfname = outdir + '/' + utils.sanitize_name(gene) + '.csv'
            with opener('w')(outfname) as outfile:
                nuke_header = [n + xtra for n in utils.nukes for xtra in ('', '_obs', '_lo_err', '_hi_err')]
                writer = csv.DictWriter(outfile, ('position', 'mute_freq', 'lo_err', 'hi_err') + tuple(nuke_header))
                writer.writeheader()
                for position in sorted(gcounts.keys()):
                    row = {'position':position,
                           'mute_freq':freqs[position]['freq'],
                           'lo_err':freqs[position]['freq_lo_err'],
                           'hi_err':freqs[position]['freq_hi_err']}
                    for nuke in utils.nukes:
                        row[nuke] = freqs[position][nuke]
                        row[nuke + '_obs'] = gcounts[position][nuke]
                        row[nuke + '_lo_err'] = freqs[position][nuke + '_lo_err']
                        row[nuke + '_hi_err'] = freqs[position][nuke + '_hi_err']
                    writer.writerow(row)

        assert 'REGION' in mean_freq_outfname
        self.mean_rates['all'].write(mean_freq_outfname.replace('REGION', 'all'))  # hackey hackey hackey replacement... *sigh*
        for region in utils.regions:
            self.mean_rates[region].write(mean_freq_outfname.replace('REGION', region))

    # ----------------------------------------------------------------------------------------
    def tigger_plot(self, only_csv=False):
        if only_csv:  # not implemented
            return
        for gene in self.freqs:
            for position in self.freqs[gene]:
                plotting.make_tigger_plot(gene, position, self.freqs[gene][position]['tigger'])

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
        if self.tigger:
            utils.prep_dir(plotdir + '/tigger', multilings=('*.csv', '*.svg'))

        for gene in self.freqs:
            freqs = self.freqs[gene]
            sorted_positions = sorted(freqs.keys())
            genehist = Hist(sorted_positions[-1] - sorted_positions[0] + 1, sorted_positions[0] - 0.5, sorted_positions[-1] + 0.5, xtitle='fixme', ytitle='fixme')  #, title=utils.sanitize_name(gene))
            for position in sorted_positions:
                hi_diff = abs(freqs[position]['freq'] - freqs[position]['freq_hi_err'])
                lo_diff = abs(freqs[position]['freq'] - freqs[position]['freq_lo_err'])
                err = 0.5*(hi_diff + lo_diff)
                genehist.set_ibin(genehist.find_bin(position), freqs[position]['freq'], error=err)
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

        if self.tigger:
            self.tigger_plot(only_csv)

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
