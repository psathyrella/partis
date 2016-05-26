#!/usr/bin/env python

import sys
import os
import math
from subprocess import check_call
import numpy
import csv
from collections import OrderedDict
import operator
import scipy
import time

import plotting

import utils
import fraction_uncertainty
import paramutils
from hist import Hist
from opener import opener

# ----------------------------------------------------------------------------------------
class MuteFreqer(object):
    def __init__(self, germline_seqs, find_new_alleles=False, calculate_uncertainty=True):
        self.germline_seqs = germline_seqs
        self.find_new_alleles = find_new_alleles
        self.calculate_uncertainty = calculate_uncertainty
        self.counts, self.freqs = {}, {}
        n_bins, xmin, xmax = 200, 0., 1.
        self.mean_rates = {'all':Hist(n_bins, xmin, xmax)}
        for region in utils.regions:
            self.mean_rates[region] = Hist(n_bins, xmin, xmax)
        self.finalized = False

        self.gene_obs_counts = {}

        # allele finding stuff
        self.small_number = 1e-5
        self.n_max_mutations_per_segment = 20  # don't look a v segments that have more than this many mutations
        self.n_max_snps = self.n_max_mutations_per_segment - 8  # try excluding up to this many bins (on the left) when doing the fit (leaves at least 8 points for fit)
        self.n_obs_min = 10  # don't fit positions that have fewer observed mutations than this
        self.min_non_candidate_positions_to_fit = 5  # always fit at least a few non-candidate positions
        self.min_y_intercept = 0.3  # roughly speaking, use this as the boundary between snp and non-snp positions
        self.default_slope_bounds = (-0.2, 0.2)  # fitting function needs some reasonable bounds from which to start
        self.big_y_icpt_bounds = (self.min_y_intercept, 1.5)  # snp-candidate positions should fit well when forced to use these bounds, but non-snp positions should fit like &*@!*
        self.min_score = 2  # (mean ratio over snp candidates) - (first non-candidate ratio) must be greater than this
        self.min_candidate_ratio = 2  # every candidate ratio must be greater than this

    # ----------------------------------------------------------------------------------------
    def increment(self, info):
        self.mean_rates['all'].fill(utils.get_mutation_rate(self.germline_seqs, info))  # mean freq over whole sequence (excluding insertions)

        for region in utils.regions:
            regional_freq, len_excluding_ambig = utils.get_mutation_rate(self.germline_seqs, info, restrict_to_region=region, return_len_excluding_ambig=True)
            n_mutes = regional_freq * len_excluding_ambig  # total number of mutations in the region (for allele finding stuff)
            if abs(n_mutes - int(n_mutes)) > 1e6:
                raise Exception('n mutated %f not an integer' % n_mutes)
            n_mutes = int(n_mutes)
            self.mean_rates[region].fill(regional_freq)  # per-region mean freq

            # per-gene per-position freqs
            gene = info[region + '_gene']
            if gene not in self.counts:
                self.counts[gene] = {}
                self.gene_obs_counts[gene] = 0
            self.gene_obs_counts[gene] += 1
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
                    gcounts[igl]['allele-finding'] = {}
                gcounts[igl]['total'] += 1
                gcounts[igl][query_seq[ipos]] += 1  # note that if <query_seq[ipos]> isn't among <utils.nukes>, this will toss a key error

                if self.find_new_alleles:
                    if igl not in gcounts:
                        gcounts[igl]['allele-finding'] = {}
                    if utils.get_region(gene) == 'v':
                        if n_mutes not in gcounts[igl]['allele-finding']:
                            gcounts[igl]['allele-finding'][n_mutes] = {'muted' : 0, 'total' : 0}
                        gcounts[igl]['allele-finding'][n_mutes]['total'] += 1
                        if query_seq[ipos] != germline_seq[ipos]:  # if this position is mutated
                            gcounts[igl]['allele-finding'][n_mutes]['muted'] += 1  # mark that we saw this germline position mutated once in a sequence with <n_mutes> regional mutation frequency

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
    def get_residual_sum(self, xvals, yvals, errs, slope, intercept):
        def expected(x):
            return slope * x + intercept
        residual_sum = sum([(y - expected(x))**2 / err**2 for x, y, err in zip(xvals, yvals, errs)])
        return residual_sum

    # ----------------------------------------------------------------------------------------
    def get_curvefit(self, pos, n_mutelist, freqs, errs, y_icpt_bounds=None):
        def func(x, slope, y_icpt):
            return slope*x + y_icpt

        bounds = (-float('inf'), float('inf'))
        if y_icpt_bounds is not None:
            bounds = [[s, y] for s, y in zip(self.default_slope_bounds, y_icpt_bounds)]
        params, cov = scipy.optimize.curve_fit(func, n_mutelist, freqs, sigma=errs, bounds=bounds)
        slope, slope_err = params[0], math.sqrt(cov[0][0])
        y_icpt, y_icpt_err = params[1], math.sqrt(cov[1][1])
        residual_sum = self.get_residual_sum(n_mutelist, freqs, errs, slope, y_icpt)
        ndof = len(n_mutelist) - 1
        fitfo = {
            'slope'  : slope,
            'y_icpt' : y_icpt,
            'slope_err'  : slope_err,
            'y_icpt_err' : y_icpt_err,
            'residuals_over_ndof' : float(residual_sum) / ndof,
            'print_str' : '    %9.3f +/- %-9.3f   %7.4f +/- %7.4f    %7.4f' % (y_icpt, y_icpt_err, slope, slope_err, float(residual_sum) / ndof),
            'n_mutelist' : n_mutelist,
            'freqs' : freqs,
            'errs' : errs
        }

        return fitfo

    # ----------------------------------------------------------------------------------------
    def get_allele_finding_xyvals(self, position, gpcounts):
        iterinfo = gpcounts['allele-finding'].items()

        obs = [d['muted'] for nm, d in iterinfo if nm < self.n_max_mutations_per_segment]

        lohis = [fraction_uncertainty.err(d['muted'], d['total'], use_beta=True) for nm, d in iterinfo if nm < self.n_max_mutations_per_segment]
        errs = [(hi - lo) / 2 for lo, hi, _ in lohis]
        weights = [1./(e*e) for e in errs]

        freqs = [float(d['muted']) / d['total'] if d['total'] > 0 else 0. for nm, d in iterinfo if nm < self.n_max_mutations_per_segment]
        total = [d['total'] for nm, d in iterinfo if nm < self.n_max_mutations_per_segment]
   
        n_mutelist = [nm for nm in gpcounts['allele-finding'].keys() if nm < self.n_max_mutations_per_segment]

        return {'obs' : obs, 'total' : total, 'n_mutelist' : n_mutelist, 'freqs' : freqs, 'errs' : errs, 'weights' : weights}

    # ----------------------------------------------------------------------------------------
    def finalize_allele_finding(self, debug=False):
        start = time.time()
        for gene in self.counts:
            if utils.get_region(gene) != 'v':
                continue
            if debug:
                print '\n%s (observed %d times)' % (utils.color_gene(gene), self.gene_obs_counts[gene])

            positions = sorted(self.counts[gene].keys())
            xyvals = {pos : self.get_allele_finding_xyvals(pos, self.counts[gene][pos]) for pos in positions}
            positions_to_fit = [pos for pos in positions if sum(xyvals[pos]['obs']) > self.n_obs_min]  # ignore positions with only a few observed mutations
            if debug:
                print '          skipping %d / %d positions (with fewer than %d observed mutations)' % (len(positions) - len(positions_to_fit), len(positions), self.n_obs_min)

            for pos in positions_to_fit:
                self.freqs[gene][pos]['allele-finding'] = xyvals[pos]

            scores, min_ratios, candidates = {}, {}, {}
            for istart in range(1, self.n_max_snps):
                if debug:
                    if istart == 1:
                        print 'snps       position   ratio   (m=0 / m>%5.2f)       muted / obs ' % self.big_y_icpt_bounds[0]
                    print '%-3d' % istart

                subxyvals = {pos : {k : v[istart:] for k, v in xyvals[pos].items()} for pos in positions_to_fit}

                residuals = {}
                for pos in positions_to_fit:
                    # as long as we already have a few non-candidate positions, skip positions that have no frequencies greater than the min y intercept (note that they could in principle still have a large y intercept, but we don't really care)
                    if len(residuals) > istart + self.min_non_candidate_positions_to_fit and len([f for f in subxyvals[pos]['freqs'] if f > self.min_y_intercept]) == 0:
                        continue

                    zero_icpt_fit = self.get_curvefit(pos, subxyvals[pos]['n_mutelist'], subxyvals[pos]['freqs'], subxyvals[pos]['errs'], y_icpt_bounds=(0. - self.small_number, 0. + self.small_number))
                    big_icpt_fit = self.get_curvefit(pos, subxyvals[pos]['n_mutelist'], subxyvals[pos]['freqs'], subxyvals[pos]['errs'], y_icpt_bounds=self.big_y_icpt_bounds)

                    residuals[pos] = {'zero_icpt' : zero_icpt_fit['residuals_over_ndof'], 'big_icpt' : big_icpt_fit['residuals_over_ndof']}

                residual_ratios = {pos : r['zero_icpt'] / r['big_icpt'] for pos, r in residuals.items()}
                sorted_ratios = sorted(residual_ratios.items(), key=operator.itemgetter(1), reverse=True)  # sort the positions in decreasing order of residual ratio
                candidate_snps = [pos for pos, _ in sorted_ratios[:istart]]  # the first <istart> positions are the candidate snps
                first_non_snp, _ = sorted_ratios[istart]

                mean_of_candidates = numpy.mean([residual_ratios[cs] for cs in candidate_snps])
                first_non_snp_ratio = residual_ratios[first_non_snp]
                scores[istart] = mean_of_candidates - first_non_snp_ratio
                min_ratios[istart] = min([residual_ratios[cs] for cs in candidate_snps])
                candidates[istart] = {cp : residual_ratios[cp] for cp in candidate_snps}

                if debug:
                    for pos in candidate_snps + [first_non_snp, ]:
                        xtrastrs = ('[', ']') if pos == first_non_snp else (' ', ' ')
                        print '             %s %3d    %5.2f   (%5.2f / %-5.2f)       %3d / %3d %s' % (xtrastrs[0], pos, residual_ratios[pos],
                                                                                                      residuals[pos]['zero_icpt'], residuals[pos]['big_icpt'], sum(subxyvals[pos]['obs']), sum(subxyvals[pos]['total']), xtrastrs[1])
                    print ' %7.2f   (%5.2f - %-5.2f)       (min %7.2f)' % (scores[istart], mean_of_candidates, first_non_snp_ratio, min_ratios[istart])

            n_candidate_snps = None
            for istart, ratio in sorted(scores.items(), key=operator.itemgetter(1), reverse=True):
                if n_candidate_snps is None and ratio > self.min_score and min_ratios[istart] > self.min_candidate_ratio:  # take the biggest score that also has a big enough min candidate ratio
                    n_candidate_snps = istart

            if debug:
                print '    snps      score     min'
                for istart, score in sorted(scores.items(), key=operator.itemgetter(1), reverse=True):
                    print_str = '    %2d     %7.2f   %7.2f' % (istart, score, min_ratios[istart])
                    if istart == n_candidate_snps:
                        print_str = utils.color('red', print_str)
                    print print_str

            if n_candidate_snps is not None:
                print '\n    %s found a new allele candidate separated from %s by %d SNP%s at position%s %s' % (utils.color('red', 'note'), utils.color_gene(gene), n_candidate_snps,
                                                                                                                's' if n_candidate_snps > 1 else '', 's' if n_candidate_snps > 1 else '',
                                                                                                                ', '.join(['%d' % p for p in candidates[n_candidate_snps]]))
            else:
                print '\n      no new alleles found'

        if debug:
            print '      allele finding time: %.1f' % (time.time()-start)

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

        if self.find_new_alleles:
            self.finalize_allele_finding()

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
    def allele_finding_plot(self, gene, plotdir, only_csv=False):
        if only_csv:  # not implemented
            return
        if utils.get_region(gene) != 'v':
            return
        utils.prep_dir(plotdir, multilings=('*.csv', '*.svg'))
        for position in self.freqs[gene]:
            if 'allele-finding' in self.freqs[gene][position] and self.freqs[gene][position]['allele-finding'] is not None:
                plotting.make_allele_finding_plot(plotdir, gene, position, self.freqs[gene][position]['allele-finding'])

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
            # per-base plots:
            # paramutils.make_mutefreq_plot(plotdir + '/' + utils.get_region(gene) + '-per-base', utils.sanitize_name(gene), plotting_info)  # needs translation to mpl

            if self.find_new_alleles:
                self.allele_finding_plot(gene, plotdir + '/allele-finding/' + utils.sanitize_name(gene), only_csv)

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
