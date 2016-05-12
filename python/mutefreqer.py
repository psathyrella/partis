#!/usr/bin/env python

import sys
import os
import math
from subprocess import check_call
import numpy
import csv
from collections import OrderedDict

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

        self.gene_obs_counts = {}

        # tigger stuff
        self.tigger = True
        self.n_max_mutes = 20
        self.n_max_bins_to_exclude = self.n_max_mutes - 5  # try excluding up to this many bins (on the left) when doing the fit
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
    def get_residual_sum(self, xvals, yvals, errs, slope, intercept):
        def expected(x):
            return slope * x + intercept
        residual_sum = sum([(y - expected(x))**2 / err**2 for x, y, err in zip(xvals, yvals, errs)])
        return residual_sum

    # ----------------------------------------------------------------------------------------
    def get_chi_squared(self, xvals, yvals, slope, intercept):
        def expected(x):
            return slope * x + intercept
        chi_squared = sum([(y - expected(x))**2 / expected(x) for x, y in zip(xvals, yvals)])
        return chi_squared

    # # ----------------------------------------------------------------------------------------
    # def get_r_squared(self, xvals, yvals, slope, intercept):
    #     ybar = float(sum(yvals)) / len(yvals)
    #     numer = sum([(slope * x + intercept - ybar)**2 for x in xvals])
    #     denom = sum([(y - ybar)**2 for y in yvals])
    #     # print '    %8.3f / %8.3f = %8.3f' % (numer, denom, numer / denom)
    #     return numer / denom

    # ----------------------------------------------------------------------------------------
    def get_fit(self, n_mutelist, freqs, weights, errs):
        if len(n_mutelist) < 3:
            return {'print_str' : '    %9s +/- %-9s   %7s +/- %7s    %7s' % ('', '', '', '', '')}
        
        # if freqs.count(0.) == len(freqs):  # no need to fit if it's all zeros
        #     slope = 0.

        params, cov = numpy.polyfit(n_mutelist, freqs, 1, w=weights, cov=True)
        slope, slope_err = params[0], math.sqrt(cov[0][0])
        y_icpt, y_icpt_err = params[1], math.sqrt(cov[1][1])

        # r_squared = self.get_r_squared(n_mutelist, freqs, slope, y_icpt)
        # chi_squared = self.get_chi_squared(n_mutelist, freqs, slope, y_icpt)
        residual_sum = self.get_residual_sum(n_mutelist, freqs, errs, slope, y_icpt)
        ndof = len(n_mutelist) - 1
    
        x_icpt, x_icpt_err = 0, 0
        # interesting = y_icpt - 3*y_icpt_err > self.min_y_intercept  # if y_icpt minus 3 std devs don't get you below the min, it's probably a snp
        # if interesting:
        #     x_icpt, x_icpt_err = 0, 0
        # else:
        #     x_icpt = -y_icpt / slope
        #     x_icpt_err = abs(y_icpt / slope) * math.sqrt((y_icpt_err/y_icpt)**2 + (slope_err/slope)**2)
        #     # err_squared = x_icpt_err * x_icpt_err
        #     # mean_x_icpt['sum'] += x_icpt / err_squared
        #     # mean_x_icpt['total'] += 1. / err_squared

        fitfo = {
            'slope'  : slope,
            'x_icpt' : x_icpt,
            'y_icpt' : y_icpt,
            'slope_err'  : slope_err,
            'x_icpt_err' : x_icpt_err,
            'y_icpt_err' : y_icpt_err,
            # 'r_squared' : r_squared,
            # 'chi_squared' : chi_squared,
            'residual_sum' : residual_sum,
            'ndof' : ndof,
            'print_str' : '    %9.3f +/- %-9.3f   %7.4f +/- %7.4f    %7.4f' % (y_icpt, y_icpt_err, slope, slope_err, float(residual_sum) / ndof),
        }

        return fitfo

    # ----------------------------------------------------------------------------------------
    def get_tigger_xyvals(self, position, gpcounts):
        iterinfo = gpcounts['tigger'].items()

        obs = [d['muted'] for nm, d in iterinfo if nm < self.n_max_mutes]

        lohis = [fraction_uncertainty.err(d['muted'], d['total'], use_beta=True) for nm, d in iterinfo if nm < self.n_max_mutes]
        errs = [(hi - lo) / 2 for lo, hi, _ in lohis]
        weights = [1./(e*e) for e in errs]

        freqs = [float(d['muted']) / d['total'] if d['total'] > 0 else 0. for nm, d in iterinfo if nm < self.n_max_mutes]
        total = [d['total'] for nm, d in iterinfo if nm < self.n_max_mutes]
   
        n_mutelist = [nm for nm in gpcounts['tigger'].keys() if nm < self.n_max_mutes]

        return {'obs' : obs, 'total' : total, 'n_mutelist' : n_mutelist, 'freqs' : freqs, 'errs' : errs, 'weights' : weights}

    # # ----------------------------------------------------------------------------------------
    # def tigger_calcs(self, position, n_mutelist, freqs, weights):
    #     fitfo = self.get_fit(n_mutelist, freqs, weights)
    #     # print_str = '       %3d   %9.3f    %9.3f +/- %-9.3f   %9.3f +/- %-9.3f   %7.4f +/- %7.4f      %3d / %3d' % (istart, freqs[istart], fitfo['x_icpt'], fitfo['x_icpt_err'], fitfo['y_icpt'], fitfo['y_icpt_err'], fitfo['slope'], fitfo['slope_err'], sum(obs), sum(total))
    #     print_str = '       %3d' % istart
    #     print_str += fitfo['print_str']
    #     print_str += '    %3d / %3d' % (sum(obs), sum(total))
    #     if self.is_interesting(fitfo):
    #         print_str = utils.color('red', print_str)
    #     print print_str

    #     plotinfo = {'n_muted' : n_mutelist[1:], 'freqs' : freqs[1:], 'errs' : errs[1:], 'slope' : fitfo['slope'], 'intercept' : fitfo['y_icpt']}
    #     return plotinfo

    # ----------------------------------------------------------------------------------------
    def big_y_intercept(self, y_icpt, y_icpt_err):
        return y_icpt - 2*y_icpt_err > self.min_y_intercept  # if y_icpt minus two std devs don't get you below the min, it's probably a snp

    # ----------------------------------------------------------------------------------------
    def find_istart_of_interest(self, residuals, y_icpts, y_icpt_errs, debug=False):
        min_residual_ratio = 5.
        fraction_defining_most = 4. / 5

        def interesting_residuals(istart):
            # require that the fits gets dramatically better at this point
            if residuals[istart-1] / residuals[istart] < min_residual_ratio:
                if debug:
                    print '      residual ratio too small %.2f / %.2f = %.2f < %.2f' % (residuals[istart-1], residuals[istart], residuals[istart-1] / residuals[istart], min_residual_ratio)
                return False

            # require that fits are uniformly crappy to left of here
            crappy_fit_val = 5.
            if float(sum(residuals[:istart])) / istart < crappy_fit_val:
                if debug:
                    print '      fits to left aren\'t bad enough: %.2f / %d = %.2f < %.2f' % (float(sum(residuals[:istart])), istart, float(sum(residuals[:istart])) / istart, crappy_fit_val)
                return False

            # and uniformly good above here
            other_crappy_fit_val = 2.5
            if float(sum(residuals[istart + 1 :])) / (len(residuals) - istart - 1) > other_crappy_fit_val:
                if debug:
                    print istart
                    print residuals
                    print '      fits to right aren\'t good enough: %.2f / %d = %.2f > %.2f' % (float(sum(residuals[istart + 1 :])), (len(residuals) - istart - 1), float(sum(residuals[istart + 1 :])) / (len(residuals) - istart - 1), other_crappy_fit_val)
                return False

            return True

        def interesting_y_icpts(istart):
            if not self.big_y_intercept(y_icpts[istart], y_icpt_errs[istart]):
                if debug:
                    print '      y intercept not big enough %.2f - 2*%.2f = %.2f' % (y_icpts[istart], y_icpt_errs[istart], y_icpts[istart] - 2*y_icpt_errs[istart])
                return False

            # require that most of the istarts above <istart> also have large y intercepts
            big_y_list = [self.big_y_intercept(y_icpts[i], y_icpt_errs[i]) for i in range(istart + 1, len(y_icpts))]
            if float(big_y_list.count(True)) / len(big_y_list) < fraction_defining_most:
                if debug:
                    print '      too many small y intercepts to right of here'
                return False

            return True
            
        istart_of_interest = None
        for istart in range(1, len(residuals)):
            if debug:
                print '  %2d' % istart
            if residuals[istart] is None:
                if debug:
                    print '      none-type residual'
                continue
            if interesting_residuals(istart) and interesting_y_icpts(istart):
                if istart_of_interest is not None:
                    raise Exception('  multiple istarts of interest')
                istart_of_interest = istart

        return istart_of_interest

    # ----------------------------------------------------------------------------------------
    def find_positions_of_interest(self, fitfo, positions, debug=False):
        for position in positions:
            residuals, y_icpts, y_icpt_errs = [], [], []
            for istart in fitfo:
                # residuals = [f['residual_sum'] for f in fitfo[istart][position] for istart in fitfo]
                if position in fitfo[istart]:
                    residuals.append(float(fitfo[istart][position]['residual_sum']) / fitfo[istart][position]['ndof'])
                    y_icpts.append(fitfo[istart][position]['y_icpt'])
                    y_icpt_errs.append(fitfo[istart][position]['y_icpt_err'])
                else:
                    residuals.append(None)
                    y_icpts.append(None)
                    y_icpt_errs.append(None)

            if debug:
                print position
            istart_of_interest = self.find_istart_of_interest(residuals, y_icpts, y_icpt_errs, debug=debug)
            if istart_of_interest is not None:
                print '  interesting: %d for position %d' % (istart_of_interest, position)

    # ----------------------------------------------------------------------------------------
    def finalize_tigger(self):

        def is_printably_interesting(tmpfitfo):
            # return tmpfitfo['y_icpt'] - 3*tmpfitfo['y_icpt_err'] > self.min_y_intercept  # if y_icpt minus 3 std devs don't get you below the min, it's probably a snp
            return self.big_y_intercept(tmpfitfo['y_icpt'], tmpfitfo['y_icpt_err'])

        for gene in self.counts:
            if utils.get_region(gene) != 'v':
                continue
            print '\n%s (observed %d times)' % (gene, self.gene_obs_counts[gene])

            xyvals = {}
            positions = sorted(self.counts[gene].keys())
            for position in positions:
                xyvals[position] = self.get_tigger_xyvals(position, self.counts[gene][position])

            fitfo = OrderedDict()
            for istart in range(self.n_max_bins_to_exclude):
                print 'istart %3d' % istart
                print '      position            y-icpt                   slope           residuals       mut / total'
                n_calcd_positions, n_skipped_positions = 0, 0  # positions that didn't have enough observed mutations
                fitfo[istart] = OrderedDict()
                for position in positions:
                    subxyvals = {k : v[istart:] for k, v in xyvals[position].items()}
                    if sum(subxyvals['obs']) < self.n_obs_min:  # ignore positions with only a few observed mutations
                        n_skipped_positions += 1
                        continue

                    n_calcd_positions += 1
                    fitfo[istart][position] = self.get_fit(subxyvals['n_mutelist'], subxyvals['freqs'], subxyvals['weights'], subxyvals['errs'])

                    print_str = '        %3d   %s       %3d / %3d' % (position, fitfo[istart][position]['print_str'], sum(subxyvals['obs']), sum(subxyvals['total']))
                    if is_printably_interesting(fitfo[istart][position]):
                        print_str = utils.color('red', print_str)
                    print print_str

                    # self.freqs[gene][position]['tigger'] =

                print '          %d / %d positions with more than %d observed mutations (i.e. skipped %d positions)' % (n_calcd_positions, n_calcd_positions + n_skipped_positions, self.n_obs_min, n_skipped_positions)

            self.find_positions_of_interest(fitfo, positions)

        sys.exit(1)

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
    def tigger_plot(self, gene, plotdir, only_csv=False):
        if only_csv:  # not implemented
            return
        if utils.get_region(gene) != 'v':
            return
        utils.prep_dir(plotdir, multilings=('*.csv', '*.svg'))
        for position in self.freqs[gene]:
            if self.freqs[gene][position]['tigger'] is not None:
                plotting.make_tigger_plot(plotdir, gene, position, self.freqs[gene][position]['tigger'])

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

            if self.tigger:
                self.tigger_plot(gene, plotdir + '/tigger/' + utils.sanitize_name(gene), only_csv)

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
