import itertools
import time
import math
import sys
import os
import operator
from subprocess import check_call
import scipy
import glob
import numpy

from hist import Hist
import fraction_uncertainty
from mutefreqer import MuteFreqer
import plotting
import utils
import glutils

# ----------------------------------------------------------------------------------------
def fstr(fval):
    if fval is None:
        return 'none'
    elif fval < 10:
        return '%.2f' % fval
    elif fval < 1e4:
        return '%.0f.' % fval
    else:
        return '%.0e' % fval

# ----------------------------------------------------------------------------------------
class AlleleFinder(object):
    def __init__(self, glfo, args, itry, cpath=None):
        self.glfo = glfo
        self.args = args
        self.itry = itry
        self.cpath = cpath

        self.fraction_of_seqs_to_exclude = 0.01  # exclude the fraction of sequences with largest v_{5,3}p deletions whose counts add up to this fraction of total sequences NOTE you don't want to make this too big, because although you'll be removing all the seqs with large 4p deletions, this number also gets used when you're deciding whether your new allele is in the default glfo
        self.n_bases_to_exclude = {'5p' : {}, '3p' : {}}  # i.e. on all the seqs we keep, we exclude this many bases; and any sequences that have larger deletions than this are not kept
        self.genes_to_exclude = set()  # genes that, with the above restrictions, are entirely excluded

        self.max_fit_length = 10

        self.n_snps_to_switch_to_two_piece_method = 5

        self.n_muted_min = 30  # don't fit positions that have fewer total mutations than this (i.e. summed over bins)
        self.n_total_min = 150  # ...or fewer total observations than this
        self.n_muted_min_per_bin = 8  # <istart>th bin has to have at least this many mutated sequences (i.e. 2-3 sigma from zero)
        self.min_fraction_per_bin = 0.005  # require that every bin (i.e. n_muted) from 0 through <self.args.max_n_snps> has at least 1% of the total (unless overridden)

        self.min_min_candidate_ratio = 2.25  # every candidate ratio must be greater than this
        self.min_mean_candidate_ratio = 2.75  # mean of candidate ratios must be greater than this
        self.min_bad_fit_residual = 2.
        self.max_good_fit_residual = 2.5
        self.max_ok_fit_residual = 10.

        self.min_min_candidate_ratio_to_plot = 1.5  # don't plot positions that're below this (for all <istart>)

        self.default_slope_bounds = (-0.1, 1.)  # fitting function needs some reasonable bounds from which to start (I could at some point make slope part of the criteria for candidacy, but it wouldn't add much sensitivity)
        self.unbounded_y_icpt_bounds = (-1., 1.5)

        self.counts, self.fitfos = {}, {}
        self.new_allele_info = []
        self.positions_to_plot = {}
        self.n_seqs_too_highly_mutated = {}  # sequences (per-gene) that had more than <self.args.n_max_mutations_per_segment> mutations
        self.gene_obs_counts = {}
        self.overall_mute_counts = Hist(self.args.n_max_mutations_per_segment - 1, 0.5, self.args.n_max_mutations_per_segment - 0.5)  # i.e. 0th (underflow) bin corresponds to zero mutations
        self.per_gene_mute_counts = {}  # crappy name -- this is the denominator for each position in <self.counts>. Which is usually the same for most positions, but it's cleaner to keep it separate than choose one of 'em.
        self.n_big_del_skipped = {s : {} for s in self.n_bases_to_exclude}

        self.n_fits = 0

        self.default_initial_glfo = None
        if self.args.default_initial_germline_dir is not None:  # if this is set, we want to take any new allele names from this directory's glfo if they're in there
            self.default_initial_glfo = glutils.read_glfo(self.args.default_initial_germline_dir, glfo['chain'])

        self.finalized = False

        self.reflengths = {}
        self.alleles_with_evidence = set()

    # ----------------------------------------------------------------------------------------
    def init_gene(self, gene):
        self.counts[gene] = {}
        for igl in range(len(self.glfo['seqs'][utils.get_region(gene)][gene])):
            self.counts[gene][igl] = {}
            for istart in range(self.args.n_max_mutations_per_segment + 1):  # istart and n_mutes are equivalent
                self.counts[gene][igl][istart] = {n : 0 for n in ['muted', 'total'] + utils.nukes}
        self.gene_obs_counts[gene] = 0
        self.per_gene_mute_counts[gene] = Hist(self.args.n_max_mutations_per_segment - 1, 0.5, self.args.n_max_mutations_per_segment - 0.5)  # i.e. 0th (underflow) bin corresponds to zero mutations
        for side in self.n_big_del_skipped:
            self.n_big_del_skipped[side][gene] = 0
        self.n_seqs_too_highly_mutated[gene] = 0

    # ----------------------------------------------------------------------------------------
    def set_excluded_bases(self, swfo, debug=False):
        region = 'v'
        for side in self.n_bases_to_exclude:
            dcounts = {}
            for query in swfo['queries']:
                gene = swfo[query][region + '_gene']
                if gene not in dcounts:
                    dcounts[gene] = {}
                dlen = swfo[query][region + '_' + side + '_del']
                if dlen not in dcounts[gene]:
                    dcounts[gene][dlen] = 0
                dcounts[gene][dlen] += 1

            for gene in dcounts:
                observed_deletions = sorted(dcounts[gene].keys())
                total_obs = sum(dcounts[gene].values())
                running_sum = 0
                if debug > 1:
                    print gene
                    print '  observed %s deletions: %s (counts %s)' % (side, ' '.join([str(d) for d in observed_deletions]), ' '.join([str(c) for c in dcounts[gene].values()]))
                    print '     len   fraction'
                for dlen in observed_deletions:
                    self.n_bases_to_exclude[side][gene] = dlen  # setting this before the "if" means that if we fall through (e.g. if there aren't enough sequences to get above the threshold) we'll still have a reasonable default
                    running_sum += dcounts[gene][dlen]
                    if debug > 1:
                        print '    %4d    %5.3f' % (dlen, float(running_sum) / total_obs)
                    if float(running_sum) / total_obs > 1. - self.fraction_of_seqs_to_exclude:  # if we've already added deletion lengths accounting for most of the sequences, ignore the rest of 'em
                        break
                if debug > 1:
                    print '     choose', self.n_bases_to_exclude[side][gene]

        # print choices and check consistency
        if debug:
            print '    exclusions:  5p   3p'
        for gene in dcounts:
            if debug:
                print '                %3d  %3d  %s' % (self.n_bases_to_exclude['5p'][gene], self.n_bases_to_exclude['3p'][gene], utils.color_gene(gene, width=15))
            if self.n_bases_to_exclude['5p'][gene] + self.n_bases_to_exclude['3p'][gene] >= len(self.glfo['seqs'][utils.get_region(gene)][gene]):
                self.genes_to_exclude.add(gene)
                if debug:
                    print '%s excluding from analysis' % utils.color('red', 'too long:')
            # print ''

    # ----------------------------------------------------------------------------------------
    def get_seqs(self, info, region, gene):
        germline_seq = info[region + '_gl_seq']
        assert len(info['seqs']) == 1
        query_seq = info[region + '_qr_seqs'][0]
        assert len(germline_seq) == len(query_seq)

        left_exclusion = self.n_bases_to_exclude['5p'][gene] - info[region + '_5p_del']
        right_exclusion = self.n_bases_to_exclude['3p'][gene] - info[region + '_3p_del']
        assert left_exclusion >= 0  # internal consistency check -- we should've already removed all the sequences with bigger deletions
        assert right_exclusion >= 0
        if left_exclusion + right_exclusion >= len(germline_seq):
            raise Exception('excluded all bases for %s: %d + %d >= %d' % (info['unique_ids'], left_exclusion, right_exclusion, len(germline_seq)))
        germline_seq = germline_seq[left_exclusion : len(germline_seq) - right_exclusion]
        query_seq = query_seq[left_exclusion : len(query_seq) - right_exclusion]
        # NOTE <germline_seq> and <query_seq> no longer correspond to <info>, but that should be ok
        if gene not in self.reflengths:
            self.reflengths[gene] = len(query_seq)
        assert self.reflengths[gene] == len(query_seq)  # just an internal consistency check now -- they should all be identical

        n_mutes = utils.hamming_distance(germline_seq, query_seq)
        return n_mutes, germline_seq, query_seq

    # ----------------------------------------------------------------------------------------
    def increment(self, info):
        for region in ['v', ]:
            gene = info[region + '_gene']
            if gene in self.genes_to_exclude:
                continue
            if gene not in self.counts:
                self.init_gene(gene)
            gcts = self.counts[gene]  # shorthand name

            self.gene_obs_counts[gene] += 1

            skip_this = False
            for side in self.n_bases_to_exclude:  # NOTE this is important, because if there is a snp to the right of the cysteine, all the sequences in which it is removed by the v_3p deletion will be shifted one bin leftward, screwing everything up (same goes now that we're doing the same thing on the left side)
                if info[region + '_' + side + '_del'] > self.n_bases_to_exclude[side][gene]:
                    self.n_big_del_skipped[side][gene] += 1
                    skip_this = True
                    break  # don't really have to break, but it makes it so the counters add up more nicely
            if skip_this:
                continue

            n_mutes, germline_seq, query_seq = self.get_seqs(info, region, gene)  # NOTE no longer necessarily correspond to <info>
            self.overall_mute_counts.fill(n_mutes)  # NOTE this is almost the same as the hists in mutefreqer.py, except those divide by sequence length, and more importantly, they don't make sure all the sequences are the same length (i.e. the base exclusions stuff)
            self.per_gene_mute_counts[gene].fill(n_mutes)  # NOTE this is almost the same as the hists in mutefreqer.py, except those divide by sequence length, and more importantly, they don't make sure all the sequences are the same length (i.e. the base exclusions stuff)

            # i.e. do *not* use <info> after this point

            if n_mutes > self.args.n_max_mutations_per_segment:
                self.n_seqs_too_highly_mutated[gene] += 1
                continue

            assert len(germline_seq) == len(self.glfo['seqs'][region][gene]) - self.n_bases_to_exclude['5p'][gene] - self.n_bases_to_exclude['3p'][gene]
            for ipos in range(len(germline_seq)):
                igl = ipos + self.n_bases_to_exclude['5p'][gene]  # position in original (i.e. complete) germline gene

                if germline_seq[ipos] in utils.ambiguous_bases or query_seq[ipos] in utils.ambiguous_bases:  # skip if either germline or query sequence is ambiguous at this position
                    continue

                gcts[igl][n_mutes]['total'] += 1
                if query_seq[ipos] != germline_seq[ipos]:  # if this position is mutated
                    gcts[igl][n_mutes]['muted'] += 1  # mark that we saw this germline position mutated once in a sequence with <n_mutes> regional mutation frequency
                gcts[igl][n_mutes][query_seq[ipos]] += 1  # if there's a new allele, we need this to work out what the snp'd base is

    # ----------------------------------------------------------------------------------------
    def get_residual_sum(self, xvals, yvals, errs, slope, intercept, debug=False):
        def expected(x):
            return slope * x + intercept
        residuals = [(y - expected(x))**2 / err**2 for x, y, err in zip(xvals, yvals, errs)]
        if debug:
            print 'resid: ' + ' '.join(['%5.3f' % r for r in residuals])
        return sum(residuals)

    # ----------------------------------------------------------------------------------------
    def dbgstr(self, fitfo, extra_str='', pvals=None):
        return_strs = []
        if pvals is not None:
            return_strs.append(' '.join(['%5d' % n for n in pvals['n_mutelist']]))
            return_strs.append(' '.join(['%5.3f' % f for f in pvals['freqs']]))
            return_strs.append(' '.join(['%5.3f' % e for e in pvals['errs']]))
        return_strs.append('        %s  m: %5.3f +/- %5.3f   b: %5.3f +/- %5.3f     %5.3f / %-3d = %5.3f' % (extra_str, fitfo['slope'], fitfo['slope_err'], fitfo['y_icpt'], fitfo['y_icpt_err'],
                                                                                                             fitfo['residual_sum'], fitfo['ndof'], fitfo['residuals_over_ndof']))
        return '\n'.join(return_strs)

    # ----------------------------------------------------------------------------------------
    def default_fitfo(self, ndof=0, y_icpt_bounds=None, xvals=None, yvals=None, errs=None):
        fitfo = {
            'slope'  : 0.,
            'y_icpt' : 0.,
            'slope_err'  : float('inf'),
            'y_icpt_err' : float('inf'),
            'residual_sum' : 0.,
            'ndof' : ndof,
            'residuals_over_ndof' : 0.,  # NOTE not necessarily just the one over the other
            'y_icpt_bounds' : y_icpt_bounds if y_icpt_bounds is not None else self.unbounded_y_icpt_bounds,
            'xvals' : xvals,
            'yvals' : yvals,
            'errs' : errs
        }
        return fitfo

    # ----------------------------------------------------------------------------------------
    def get_tmp_fitvals(self, pvals):
        n_mutelist, freqs, errs = [], [], []
        for im in range(len(pvals['n_mutelist'])):
            if pvals['total'][im] > 0:
                n_mutelist.append(pvals['n_mutelist'][im])
                freqs.append(pvals['freqs'][im])
                errs.append(pvals['errs'][im])
        return n_mutelist, freqs, errs

    # ----------------------------------------------------------------------------------------
    def get_curvefit(self, pvals, y_icpt_bounds, debug=False):
        n_mutelist, freqs, errs = self.get_tmp_fitvals(pvals)  # this is probably kind of slow
        if y_icpt_bounds[0] == y_icpt_bounds[1]:  # fixed y-icpt
            fitfo = self.default_fitfo(len(n_mutelist) - 1, y_icpt_bounds, xvals=n_mutelist, yvals=freqs, errs=errs)
            fitfo['y_icpt'] = y_icpt_bounds[0]
            if fitfo['ndof'] > 0:
                def linefunc(x, slope):
                    return slope*x + fitfo['y_icpt']
                bounds = self.default_slope_bounds
                params, cov = scipy.optimize.curve_fit(linefunc, n_mutelist, freqs, sigma=errs, bounds=bounds)
                self.n_fits += 1
                fitfo['slope'], fitfo['slope_err'] = params[0], math.sqrt(cov[0][0])
            elif fitfo['ndof'] == 0:
                fitfo = self.approx_fit_vals(pvals, fixed_y_icpt=fitfo['y_icpt'], debug=debug)
        else:  # floating y-icpt
            fitfo = self.default_fitfo(len(n_mutelist) - 2, y_icpt_bounds, xvals=n_mutelist, yvals=freqs, errs=errs)
            if fitfo['ndof'] > 0:
                def linefunc(x, slope, y_icpt):
                    return slope*x + y_icpt
                bounds = [[s, y] for s, y in zip(self.default_slope_bounds, y_icpt_bounds)]
                params, cov = scipy.optimize.curve_fit(linefunc, n_mutelist, freqs, sigma=errs, bounds=bounds)
                self.n_fits += 1
                fitfo['slope'], fitfo['slope_err'] = params[0], math.sqrt(cov[0][0])
                fitfo['y_icpt'], fitfo['y_icpt_err'] = params[1], math.sqrt(cov[1][1])
            elif fitfo['ndof'] == 0:
                fitfo = self.approx_fit_vals(pvals, fixed_y_icpt=None, debug=debug)

        if fitfo['ndof'] > 0:
            fitfo['residual_sum'] = self.get_residual_sum(n_mutelist, freqs, errs, fitfo['slope'], fitfo['y_icpt'])
            fitfo['residuals_over_ndof'] = float(fitfo['residual_sum']) / fitfo['ndof']
        else:
            fitfo['residual_sum'] = 1.
            fitfo['residuals_over_ndof'] = 0.

        if debug:
            print self.dbgstr(fitfo, extra_str='fit', pvals={'n_mutelist' : n_mutelist, 'freqs' : freqs, 'errs': errs})  # not necessarily the same as <pvals>

        return fitfo

    # ----------------------------------------------------------------------------------------
    def get_reweights(self, gene, position, debug=False):
        """
        Reweight mutation frequencies to correct for bin-to-bin variations in allele prevalence (due to non-flat mutation frequency distributions).
        I.e. reweight to maintain linearity in the presence of super-non-continuous mutation frequency distributions.
        I don't think it's right yet, although it kind of works.
        In any case, shit works fine without reweighting, so I don't feel like dealing with the complication a.t.m. (especially don't want to deal with correcting the uncertainties for reweighting).
        """
        gcts = self.counts[gene][position]  # shorthand name
        overall_nuke_totals = {n : 0 for n in utils.nukes}
        per_bin_nuke_totals = {n_muted : {n : 0 for n in utils.nukes} for n_muted in gcts}
        for n_muted in gcts:
            for nuke in utils.nukes:
                overall_nuke_totals[nuke] += gcts[n_muted][nuke]
                per_bin_nuke_totals[n_muted][nuke] += gcts[n_muted][nuke]
        freqs = [float(d['muted']) / d['total'] if d['total'] > 0 else 0. for d in gcts.values()]
        if debug:
            print ' ', ' '.join(['%5d' % n for n in gcts])
            for nuke in utils.nukes:
                print nuke, ' '.join(['%5d' % gcts[n][nuke] for n in gcts]), overall_nuke_totals[nuke]
            print ' ', ' '.join(['%5.3f' % f for f in freqs])

        reweighted_freqs = []
        for n_muted in gcts:
            freq, total = 0., 0.
            for nuke in utils.nukes:
                reweight = 0.
                if per_bin_nuke_totals[n_muted][nuke] > 0:
                    reweight = float(overall_nuke_totals[nuke]) / per_bin_nuke_totals[n_muted][nuke]
                total += reweight * gcts[n_muted][nuke]
                if nuke != self.glfo['seqs'][utils.get_region(gene)][gene][position]:
                    freq += reweight * gcts[n_muted][nuke]
            reweighted_freqs.append(freq / total if total > 0. else 0.)
        if debug:
            print ' ', ' '.join(['%5.3f' % f for f in reweighted_freqs])
        return reweighted_freqs

    # ----------------------------------------------------------------------------------------
    def get_allele_finding_xyvals(self, gene, position):
        gcts = self.counts[gene][position]  # shorthand name

        obs = [d['muted'] for d in gcts.values()]

        lohis = [fraction_uncertainty.err(d['muted'], d['total'], use_beta=True) if d['total'] > 0. else (0., 1., None) for d in gcts.values()]  # set uncertainty bounds to (0., 1.) for zero-denominator bins
        errs = [(hi - lo) / 2 for lo, hi, _ in lohis]
        weights = [1./(e*e) for e in errs]

        freqs = [float(d['muted']) / d['total'] if d['total'] > 0 else 0. for d in gcts.values()]
        total = [d['total'] for d in gcts.values()]

        n_mutelist = gcts.keys()

        return {'obs' : obs, 'total' : total, 'n_mutelist' : n_mutelist, 'freqs' : freqs, 'errs' : errs, 'weights' : weights}

    # ----------------------------------------------------------------------------------------
    def is_a_candidate(self, gene, fitfo, istart, debug=False):
        if fitfo['min_snp_ratios'][istart] < self.min_min_candidate_ratio:  # worst snp candidate has to be pretty good on its own
            if debug:
                print '    min snp ratio %s too small (less than %s)' % (fstr(fitfo['min_snp_ratios'][istart]), fstr(self.min_min_candidate_ratio)),
            return False
        if fitfo['mean_snp_ratios'][istart] < self.min_mean_candidate_ratio:  # mean of snp candidate ratios has to be even better
            if debug:
                print '    mean snp ratio %s too small (less than %s)' % (fstr(fitfo['mean_snp_ratios'][istart]), fstr(self.min_mean_candidate_ratio)),
            return False

        # return false if any of the candidate positions don't have enough mutated counts in the <istart>th bin NOTE this is particularly important because it's the handle that tells us it's *this* <istart> that's correct, rather than <istart> + 1
        for candidate_pos in fitfo['candidates'][istart]:
            n_istart_muted = self.counts[gene][candidate_pos][istart]['muted']
            if n_istart_muted < self.n_muted_min_per_bin:
                if debug:
                    print '    not enough mutated counts at candidate position %d with %d %s (%d < %d)' % (candidate_pos, istart, utils.plural_str('mutations', n_istart_muted), n_istart_muted, self.n_muted_min_per_bin),
                return False

        # return false if any of the candidate positions don't have enough mutated counts in the <istart>th bin NOTE this is particularly important because it's the handle that tells us it's *this* <istart> that's correct, rather than <istart> + 1
        if istart >= self.n_snps_to_switch_to_two_piece_method:
            for pos_1, pos_2 in itertools.combinations(fitfo['candidates'][istart], 2):
                # make sure all the snp positions have similar fits
                fitfo_1, fitfo_2 = self.fitfos[gene]['fitfos'][istart][pos_1], self.fitfos[gene]['fitfos'][istart][pos_2]
                if not self.consistent_slope_and_y_icpt(2.5, fitfo_1['postfo'], fitfo_2['postfo']):  # NOTE this has to be very permissive, since with multiple new alleles at the same position the y-icpt (at least) is expected to be quite different
                    if debug:
                        print '    positions %d and %d have inconsistent post-istart fits' % (pos_1, pos_2)
                    return False
                if not self.consistent_slope_and_y_icpt(2.5, fitfo_1['prefo'], fitfo_2['prefo']):
                    if debug:
                        print '    positions %d and %d have inconsistent pre-istart fits' % (pos_1, pos_2)
                    return False

        if debug:
            print '    candidate',
        return True

    # ----------------------------------------------------------------------------------------
    def approx_fit_vals(self, pvals, fixed_y_icpt=None, debug=False):
        # NOTE uncertainties are kinda complicated if you do the weighted mean, so screw it, it works fine with the plain mean
        fitfo = self.default_fitfo()

        def getslope(i1, i2, shift=False):
            if not shift:
                tmp_y = yv
            else:
                tmp_y = [yv[i] + (-1) ** (i%2) * ev[i] for i in range(len(yv))]  # shift alternately up or down by the uncertainty
            return (tmp_y[i2] - tmp_y[i1]) / (xv[i2] - xv[i1])

        xv, yv, ev = pvals['n_mutelist'], pvals['freqs'], pvals['errs']  # tmp shorthand
        assert len(xv) > 0
        if len(xv) == 1:
            if fixed_y_icpt is None:
                return fitfo  # just leave the default values
                # raise Exception('can\'t handle single point with floating y-icpt')
            if xv[0] > 0.:  # <xv[0]> can't be able to be negative, but if it's zero, then the fit values aren't well-defined
                fitfo['slope'] = (yv[0] - fixed_y_icpt) / (xv[0] - 0.)
        else:
            slopes = [getslope(i-1, i) for i in range(1, len(xv))]  # only uses adjacent points, and double-counts interior points, but we don't care (we don't use steps of two, because then we'd the last one if it's odd-length)
            if len(xv) == 2:  # add two points, shifting each direction by each point's uncertainty
                slopes.append(getslope(0, 1, shift=True))
            fitfo['slope'] = numpy.average(slopes)
            fitfo['slope_err'] = numpy.std(slopes, ddof=1) / math.sqrt(len(xv))

        if fixed_y_icpt is None:
            y_icpts = [yv[i] - getslope(i-1, i) * xv[i] for i in range(1, len(xv))]
            fitfo['y_icpt'] = numpy.average(y_icpts)
            if len(xv) == 2:
                y_icpts.append(yv[1] - getslope(0, 1, shift=True) * xv[1])
            fitfo['y_icpt_err'] = numpy.std(y_icpts, ddof=1) / math.sqrt(len(xv))
        else:
            fitfo['y_icpt'] = fixed_y_icpt

        if debug:
            print self.dbgstr(fitfo, extra_str='apr', pvals=pvals)

        return fitfo

    # ----------------------------------------------------------------------------------------
    def consistent(self, factor, v1, v1err, v2, v2err, debug=False):
        # i.e. if both slope and intercept are within <factor> std deviations of each other, don't bother fitting, because the fit isn't going to say they're wildly inconsistent
        lo, hi = sorted([v1, v2])
        joint_err = max(v1err, v2err)
        if debug:
            print '      %6.3f +/- %6.3f   %6.3f +/- %6.3f   -->   %6.3f + %3.1f * %6.3f = %6.3f >? %6.3f   %s' % (v1, v1err, v2, v2err, lo, factor, joint_err, lo + factor * joint_err, hi, lo + factor * joint_err > hi)
        return lo + factor * joint_err > hi

    # ----------------------------------------------------------------------------------------
    def consistent_slope_and_y_icpt(self, factor, vals1, vals2, debug=False):
        consistent_slopes = self.consistent(factor, vals1['slope'], vals1['slope_err'], vals2['slope'], vals2['slope_err'], debug=debug)
        consistent_y_icpts = self.consistent(factor, vals1['y_icpt'], vals1['y_icpt_err'], vals2['y_icpt'], vals2['y_icpt_err'], debug=debug)
        return consistent_slopes and consistent_y_icpts

    # ----------------------------------------------------------------------------------------
    def empty_pre_bins(self, gene, istart, positions_to_try_to_fit, debug=False):
        """ return true if fewer than <istart> positions have enough entries in the bins before <istart> """
        prexyvals = {pos : {k : v[:istart] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}  # arrays up to, but not including, <istart>
        good_positions = [pos for pos in positions_to_try_to_fit if sum(prexyvals[pos]['total']) > self.n_total_min]  # almost the same as the line where we get <positions_to_try_to_fit>, except now it's only the bins before <istart>
        return len(good_positions) < istart

    # ----------------------------------------------------------------------------------------
    def get_big_y(self, pvals):
        big_y_icpt = numpy.average(pvals['freqs'], weights=pvals['weights'])  # corresponds, roughly, to the expression level of the least common allele to which we have sensitivity NOTE <istart> is at index 0
        big_y_icpt_err = numpy.std(pvals['freqs'], ddof=1)  # NOTE this "err" is from the variance over bins, and ignores the sample statistics of each bin. This is a little weird, but good: it captures cases where the points aren't well-fit by a line, either because of multiple alleles with very different prevalences, or because the sequences aren't very independent NOTE the former case is very important)
        return big_y_icpt, big_y_icpt_err

    # ----------------------------------------------------------------------------------------
    def very_different_bin_totals(self, pvals, istart, debug=False):
        # i.e. if there's a homozygous new allele at <istart> + 1
        factor = 2.  # i.e. check everything that's more than <factor> sigma away UPDATE this will have to change -- totals per bin vary too much
        joint_total_err = max(math.sqrt(pvals['total'][istart - 1]), math.sqrt(pvals['total'][istart]))
        last_total = pvals['total'][istart - 1]
        istart_total = pvals['total'][istart]
        if debug:
            print '    different bin totals: %.0f - %.0f = %.0f ?> %.1f * %5.3f = %5.3f'  % (istart_total, last_total, istart_total - last_total, factor, joint_total_err, factor * joint_total_err)
        return istart_total - last_total > factor * joint_total_err  # it the total (denominator) is very different between the two bins

    # ----------------------------------------------------------------------------------------
    def big_discontinuity(self, pvals, istart, debug=False):
        factor = 2.  # i.e. check everything that's more than <factor> sigma away
        joint_freq_err = max(pvals['errs'][istart - 1], pvals['errs'][istart])
        last_freq = pvals['freqs'][istart - 1]
        istart_freq = pvals['freqs'][istart]
        if debug:
            print '    discontinuity: %5.3f - %5.3f = %5.3f ?> %.1f * %5.3f = %5.3f'  % (istart_freq, last_freq, istart_freq - last_freq, factor, joint_freq_err, factor * joint_freq_err)
        return istart_freq - last_freq > factor * joint_freq_err

    # ----------------------------------------------------------------------------------------
    def fit_two_piece_istart(self, gene, istart, positions_to_try_to_fit, print_dbg_header=False, debug=False):
        if debug and print_dbg_header:
            print '             position   ratio       (one piece / two pieces)  ',
            print '%0s %s' % ('', ''.join(['%11d' % nm for nm in range(self.args.n_max_mutations_per_segment + 1)]))  # NOTE *has* to correspond to line at bottom of fcn below

        # NOTE I'm including the zero bin here -- do I really want to do that? UPDATE yes, I think so -- we know there will be zero mutations in that bin, but the number of sequences in it still contains information (uh, I think)
        min_ibin = max(0, istart - self.max_fit_length)
        max_ibin = min(self.args.n_max_mutations_per_segment, istart + self.max_fit_length)
        prexyvals = {pos : {k : v[min_ibin : istart] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}  # arrays up to, but not including, <istart>
        postxyvals = {pos : {k : v[istart : max_ibin] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}  # arrays from <istart> onwards
        bothxyvals = {pos : {k : v[min_ibin : max_ibin] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}

        candidate_ratios, residfo = {}, {}  # NOTE <residfo> is really just for dbg printing... but we have to sort before we print, so we need to keep the info around
        for pos in positions_to_try_to_fit:
            dbg = False #(pos in [10, 11] and istart==2)
            if dbg:
                print 'pos %d' % pos
            prevals = prexyvals[pos]
            postvals = postxyvals[pos]
            bothvals = bothxyvals[pos]
            big_y_icpt, big_y_icpt_err = self.get_big_y(postvals)
            big_y_icpt_bounds = (big_y_icpt - 1.5*big_y_icpt_err, big_y_icpt + 1.5*big_y_icpt_err)  # we want the bounds to be lenient enough to accomodate non-zero slopes (in the future, we could do something cleverer like extrapolating with the slope of the line to x=0)

            if sum(postvals['obs']) < self.n_muted_min or sum(postvals['total']) < self.n_total_min:
                continue

            # if the discontinuity is less than <factor> sigma, and the bin totals are closer than <factor> sigma, skip it
            if not self.big_discontinuity(bothvals, istart) and not self.very_different_bin_totals(bothvals, istart):
                continue

            # if both slope and intercept are quite close to each other, the fits aren't going to say they're wildly inconsistent
            if istart < self.n_snps_to_switch_to_two_piece_method:
                # if the bounds include zero, there won't be much difference between the two fits
                if big_y_icpt_bounds[0] <= 0.:
                    continue

                # if there's only two points in <prevals>, we can't use the bad fit there to tell us this isn't a candidate, so we check and skip if the <istart - 1>th freq isn't really low
                if istart == 2 and bothvals['freqs'][istart - 1] > big_y_icpt - 1.5 * big_y_icpt_err:
                    continue

                # if a rough estimate of the y-icpt is less than zero, the zero-icpt fit is probably going to be pretty good
                approx_fitfo = self.approx_fit_vals(postvals)
                if approx_fitfo['y_icpt'] < 0.:
                    continue
            else:
                pre_approx = self.approx_fit_vals(prevals)
                post_approx = self.approx_fit_vals(postvals)
                if pre_approx['slope'] > post_approx['slope']:  #  or self.consistent_slope_and_y_icpt(pre_approx, post_approx):  # UPDATE really don't want to require inconsistent slope and y-icpt
                    continue

            onefit = self.get_curvefit(bothvals, y_icpt_bounds=(0., 0.), debug=dbg)  # self.unbounded_y_icpt_bounds)

            # don't bother with the two-piece fit if the one-piece fit is pretty good
            if onefit['residuals_over_ndof'] < self.min_bad_fit_residual:
                continue

            prefit = self.get_curvefit(prevals, y_icpt_bounds=(0., 0.), debug=dbg)  # self.unbounded_y_icpt_bounds)
            postfit = self.get_curvefit(postvals, y_icpt_bounds=big_y_icpt_bounds, debug=dbg)  # self.unbounded_y_icpt_bounds)
            twofit_residuals = prefit['residuals_over_ndof'] * prefit['ndof'] + postfit['residuals_over_ndof'] * postfit['ndof']
            twofit_ndof = prefit['ndof'] + postfit['ndof']
            twofit_residuals_over_ndof = twofit_residuals / twofit_ndof

            # at least for large <istart> pre-slope should be smaller than post-slope
            if istart >= self.n_snps_to_switch_to_two_piece_method and prefit['slope'] > postfit['slope']:
                continue

            # pre-<istart> should actually be a good line, at least for small <istart>
            if istart < self.n_snps_to_switch_to_two_piece_method and prefit['residuals_over_ndof'] > self.max_good_fit_residual:
                continue

            # <istart> and above should be a kinda-sort-ok line
            if postfit['residuals_over_ndof'] > self.max_ok_fit_residual:
                continue

            candidate_ratios[pos] = onefit['residuals_over_ndof'] / twofit_residuals_over_ndof if twofit_residuals_over_ndof > 0. else float('inf')
            residfo[pos] = {'onefo' : onefit, 'prefo' : prefit, 'postfo' : postfit, 'twofo' : {'residuals_over_ndof' : twofit_residuals_over_ndof}}
            if dbg or candidate_ratios[pos] > self.min_min_candidate_ratio_to_plot:
                self.positions_to_plot[gene].add(pos)  # if we already decided to plot it for another <istart>, it'll already be in there

        candidates = [pos for pos, _ in sorted(candidate_ratios.items(), key=operator.itemgetter(1), reverse=True)]  # sort the candidate positions in decreasing order of residual ratio
        candidates = candidates[:istart]  # remove any extra ones
        candidate_ratios = {pos : ratio for pos, ratio in candidate_ratios.items() if pos in candidates}

        if debug:
            if len(candidates) > 0:
                print '    %d %s' % (istart, utils.plural_str('snp', istart))
            for pos in candidates:
                pos_str = '%3s' % str(pos)
                if candidate_ratios[pos] > self.min_min_candidate_ratio and residfo[pos]['onefo']['residuals_over_ndof'] > self.min_bad_fit_residual:
                    pos_str = utils.color('yellow', pos_str)
                print_str = ['                 %s    %5s            %5s / %-5s               ' % (pos_str, fstr(candidate_ratios[pos]),
                                                                                                  fstr(residfo[pos]['onefo']['residuals_over_ndof']), fstr(residfo[pos]['twofo']['residuals_over_ndof']))]
                for n_mutes in range(self.args.n_max_mutations_per_segment + 1):
                    if n_mutes in bothxyvals[pos]['n_mutelist']:
                        inm = bothxyvals[pos]['n_mutelist'].index(n_mutes)
                        print_str.append('%4d / %-4d' % (bothxyvals[pos]['obs'][inm], bothxyvals[pos]['total'][inm]))
                    else:
                        print_str.append('           ')
                print ''.join(print_str)

        if len(candidates) < istart:
            return

        self.fitfos[gene]['min_snp_ratios'][istart] = min(candidate_ratios.values())
        self.fitfos[gene]['mean_snp_ratios'][istart] = numpy.mean(candidate_ratios.values())
        self.fitfos[gene]['candidates'][istart] = candidate_ratios
        self.fitfos[gene]['fitfos'][istart] = residfo

    # ----------------------------------------------------------------------------------------
    def fit_istart(self, gene, istart, positions_to_try_to_fit, fitfo, print_dbg_header=False, debug=False):
        if debug and print_dbg_header:
            print '             position   ratio    (m=0 / m=big)      big bounds',
            print '%0s %s' % ('', ''.join(['%11d' % nm for nm in range(self.args.n_max_mutations_per_segment + 1)]))  # NOTE *has* to correspond to line at bottom of fcn below

        subxyvals = {pos : {k : v[istart : istart + self.max_fit_length] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}  # arrays from <istart> onwards

        candidate_ratios, residfo = {}, {}  # NOTE <residfo> is really just for dbg printing... but we have to sort before we print, so we need to keep the info around
        for pos in positions_to_try_to_fit:
            pvals = subxyvals[pos]

            if sum(pvals['obs']) < self.n_muted_min or sum(pvals['total']) < self.n_total_min:
                continue

            big_y_icpt, big_y_icpt_err = self.get_big_y(pvals)
            big_y_icpt_bounds = (big_y_icpt - 1.5*big_y_icpt_err, big_y_icpt + 1.5*big_y_icpt_err)  # we want the bounds to be lenient enough to accomodate non-zero slopes (in the future, we could do something cleverer like extrapolating with the slope of the line to x=0)

            # if the bounds include zero, there won't be much difference between the two fits
            if big_y_icpt_bounds[0] <= 0.:
                continue

            # if a rough estimate of the x-icpt is less than zero, the zero-icpt fit is probably going to be pretty good
            approx_fitfo = self.approx_fit_vals(pvals)
            if approx_fitfo['y_icpt'] < 0.:
                continue

            zero_icpt_fit = self.get_curvefit(pvals, y_icpt_bounds=(0., 0.))

            # don't bother with the big-icpt fit if the zero-icpt fit is pretty good
            if zero_icpt_fit['residuals_over_ndof'] < self.min_bad_fit_residual:
                continue

            big_icpt_fit = self.get_curvefit(pvals, y_icpt_bounds=big_y_icpt_bounds)

            # # big-icpt fit should actually be at least ok
            # if big_icpt_fit['residuals_over_ndof'] > self.max_good_fit_residual:
            #     continue

            candidate_ratios[pos] = zero_icpt_fit['residuals_over_ndof'] / big_icpt_fit['residuals_over_ndof'] if big_icpt_fit['residuals_over_ndof'] > 0. else float('inf')
            residfo[pos] = {'zerofo' : zero_icpt_fit, 'bigfo' : big_icpt_fit}
            if candidate_ratios[pos] > self.min_min_candidate_ratio_to_plot:
                self.positions_to_plot[gene].add(pos)  # if we already decided to plot it for another <istart>, it'll already be in there

        candidates = [pos for pos, _ in sorted(candidate_ratios.items(), key=operator.itemgetter(1), reverse=True)]  # sort the candidate positions in decreasing order of residual ratio
        candidates = candidates[:istart]  # remove any extra ones
        candidate_ratios = {pos : ratio for pos, ratio in candidate_ratios.items() if pos in candidates}

        if debug:
            if len(candidates) > 0:
                print '    %d %s' % (istart, utils.plural_str('snp', istart))
            for pos in candidates:
                pos_str = '%3s' % str(pos)
                if candidate_ratios[pos] > self.min_min_candidate_ratio and residfo[pos]['zerofo']['residuals_over_ndof'] > self.min_bad_fit_residual:
                    pos_str = utils.color('yellow', pos_str)
                print_str = ['                 %s    %5s    %5s / %-5s   [%5.3f - %5.3f] ' % (pos_str, fstr(candidate_ratios[pos]),
                                                                                              fstr(residfo[pos]['zerofo']['residuals_over_ndof']), fstr(residfo[pos]['bigfo']['residuals_over_ndof']),
                                                                                              residfo[pos]['bigfo']['y_icpt_bounds'][0], residfo[pos]['bigfo']['y_icpt_bounds'][1])]
                print_str += '    '
                for n_mutes in range(self.args.n_max_mutations_per_segment + 1):
                    if n_mutes in subxyvals[pos]['n_mutelist']:
                        inm = subxyvals[pos]['n_mutelist'].index(n_mutes)
                        print_str.append('%4d / %-4d' % (subxyvals[pos]['obs'][inm], subxyvals[pos]['total'][inm]))
                    else:
                        print_str.append('           ')
                print ''.join(print_str)

        if len(candidates) < istart:
            return

        fitfo['min_snp_ratios'][istart] = min(candidate_ratios.values())
        fitfo['mean_snp_ratios'][istart] = numpy.mean(candidate_ratios.values())
        fitfo['candidates'][istart] = candidate_ratios

    # ----------------------------------------------------------------------------------------
    def see_if_new_allele_is_in_default_initial_glfo(self, new_name, new_seq, template_gene, debug=False):
        region = utils.get_region(template_gene)
        if new_name in self.default_initial_glfo['seqs'][region]:  # if we removed an existing allele and then re-added it, it'll already be in the default glfo, so there's nothing for us to do in this fcn
            return new_name, new_seq
        chain = self.glfo['chain']
        assert region == 'v'  # conserved codon stuff below will have to be changed for j
        newpos = self.glfo[utils.conserved_codons[chain][region] + '-positions'][template_gene]  # codon position for template gene should be ok
        for oldname_gene, oldname_seq in self.default_initial_glfo['seqs'][region].items():  # NOTE <oldname_{gene,seq}> is the old *name* corresponding to the new (snp'd) allele, whereas <old_seq> is the allele from which we inferred the new (snp'd) allele
            # first see if they match up through the cysteine
            oldpos = self.default_initial_glfo[utils.conserved_codons[chain][region] + '-positions'][oldname_gene]
            if oldname_seq[self.n_bases_to_exclude['5p'][template_gene] : oldpos + 3] != new_seq[self.n_bases_to_exclude['5p'][template_gene] : newpos + 3]:  # uh, I think we want to use the ones for the template gene
                continue

            # then require that any bases in common to the right of the cysteine in the new allele match the ones in the old one (where "in common" means either of them can be longer, since this just changes the insertion length)
            bases_to_right_of_cysteine = min(len(oldname_seq) - (oldpos + 3), len(new_seq) - self.n_bases_to_exclude['3p'][template_gene] - (newpos + 3))

            if bases_to_right_of_cysteine > 0 and oldname_seq[oldpos + 3 : oldpos + 3 + bases_to_right_of_cysteine] != new_seq[newpos + 3 : newpos + 3 + bases_to_right_of_cysteine]:
                continue

            print '        using old name %s for new allele %s (blue bases are not considered):' % (utils.color_gene(oldname_gene), utils.color_gene(new_name))
            def print_sequence_chunks(seq, cpos, name):
                print '            %s%s%s%s%s   %s' % (utils.color('blue', seq[:self.n_bases_to_exclude['5p'][template_gene]]),
                                                       seq[self.n_bases_to_exclude['5p'][template_gene] : cpos],
                                                       utils.color('reverse_video', seq[cpos : cpos + 3]),
                                                       seq[cpos + 3 : cpos + 3 + bases_to_right_of_cysteine],
                                                       utils.color('blue', seq[cpos + 3 + bases_to_right_of_cysteine:]),
                                                       utils.color_gene(name))
            print_sequence_chunks(oldname_seq, oldpos, oldname_gene)
            print_sequence_chunks(new_seq, newpos, new_name)

            new_name = oldname_gene
            new_seq = oldname_seq  # *very* important
            break

        return new_name, new_seq

    # ----------------------------------------------------------------------------------------
    def add_new_allele(self, template_gene, fitfo, n_candidate_snps, debug=False):
        # figure out what the new nukes are
        old_seq = self.glfo['seqs'][utils.get_region(template_gene)][template_gene]
        new_seq = old_seq
        mutfo = {}
        for pos in sorted(fitfo['candidates'][n_candidate_snps]):
            obs_counts = {nuke : self.counts[template_gene][pos][n_candidate_snps][nuke] for nuke in utils.nukes}  # NOTE it's super important to only use the counts from sequences with <n_candidate_snps> total mutations
            sorted_obs_counts = sorted(obs_counts.items(), key=operator.itemgetter(1), reverse=True)
            original_nuke = self.glfo['seqs'][utils.get_region(template_gene)][template_gene][pos]
            new_nuke = None
            for nuke, _ in sorted_obs_counts:  # take the most common one that isn't the existing gl nuke
                if nuke != original_nuke:
                    new_nuke = nuke
                    break
            assert old_seq[pos] == original_nuke
            mutfo[pos] = {'original' : original_nuke, 'new' : new_nuke}
            new_seq = new_seq[:pos] + new_nuke + new_seq[pos+1:]

        new_name, mutfo = glutils.get_new_allele_name_and_change_mutfo(template_gene, mutfo)

        if self.default_initial_glfo is not None:  # if this is set, we want to take the names from this directory's glfo (i.e. see if there's already a name for <new_name>'s sequence)
            new_name, new_seq = self.see_if_new_allele_is_in_default_initial_glfo(new_name, new_seq, template_gene, debug=debug)

        if new_name in self.glfo['seqs'][utils.get_region(template_gene)]:
            print '    new gene %s already in glfo (probably due to 3p end length issues), so skipping it' % utils.color_gene(new_name)
            return

        print '    found a new allele candidate separated from %s by %d %s at %s:' % (utils.color_gene(template_gene), n_candidate_snps,
                                                                                      utils.plural_str('snp', n_candidate_snps), utils.plural_str('position', n_candidate_snps)),
        for pos in sorted(mutfo):
            print '   %d (%s --> %s)' % (pos, mutfo[pos]['original'], mutfo[pos]['new']),
        print ''

        old_len_str, new_len_str = '', ''
        old_seq_for_cf, new_seq_for_cf = old_seq, new_seq
        if len(new_seq) > len(old_seq):  # i.e if <old_seq> (the template gene) is shorter than the sequence corresponding to the original name for the new allele that we found from it
            new_len_str = utils.color('blue', new_seq[len(old_seq):])
            new_seq_for_cf = new_seq[:len(old_seq)]
            print '         %d extra (blue) bases in new seq were not considered' % (len(new_seq) - len(old_seq))
        elif len(old_seq) > len(new_seq):
            old_len_str = utils.color('blue', old_seq[len(new_seq):])
            old_seq_for_cf = old_seq[:len(new_seq)]
            print '         %d extra (blue) bases in old seq were not considered' % (len(old_seq) - len(new_seq))
        print '          %s%s   %s' % (old_seq_for_cf, old_len_str, utils.color_gene(template_gene))
        print '          %s%s   %s' % (utils.color_mutants(old_seq_for_cf, new_seq_for_cf), new_len_str, utils.color_gene(new_name))

        # and add it to the set of new alleles for this gene
        self.new_allele_info.append({
            'template-gene' : template_gene,
            'gene' : new_name,
            'seq' : new_seq,
            'snp-positions' : mutfo.keys(),
            'aligned-seq' : None
        })

    # ----------------------------------------------------------------------------------------
    def print_skip_debug(self, gene, positions, positions_to_try_to_fit):

        def get_skip_str(skip_positions):
            skip_str = ''
            if len(skip_positions) > 0 and len(skip_positions) < 20:
                skip_str = ' (' + ' '.join([str(p) for p in skip_positions]) + ')'
            return skip_str

        glseq = self.glfo['seqs'][utils.get_region(gene)][gene]
        # new_allele_names = [i['gene'] for i in self.new_allele_info]  # should make this less hackey
        # glseq = self.new_allele_info[gene]['seq'] if gene in new_allele_names else self.glfo['seqs'][utils.get_region(gene)][gene]
        too_close_to_ends = range(self.n_bases_to_exclude['5p'][gene]) + range(len(glseq) - self.n_bases_to_exclude['3p'][gene], len(glseq))
        not_enough_counts = set(positions) - set(positions_to_try_to_fit) - set(too_close_to_ends)  # well, not enough counts, *and* not too close to the ends

        print '          skipping',
        print '%d / %d positions:' % (len(positions) - len(positions_to_try_to_fit), len(positions)),
        print '%d were too close to the ends%s' % (len(too_close_to_ends), get_skip_str(too_close_to_ends)),
        print 'and %d had fewer than %d mutations and fewer than %d observations%s' % (len(not_enough_counts), self.n_muted_min, self.n_total_min, get_skip_str(not_enough_counts))

    # ----------------------------------------------------------------------------------------
    def finalize(self, debug=False):
        assert not self.finalized

        region = 'v'
        print '%s: looking for new alleles' % utils.color('blue', 'try ' + str(self.itry))
        if not self.args.always_find_new_alleles:  # NOTE this is (on purpose) summed over all genes -- genes with homozygous unknown alleles would always fail this criterion
            binline, contents_line = self.overall_mute_counts.horizontal_print(bin_centers=True, bin_decimals=0, contents_decimals=0)
            print '             n muted in v' + binline + '(and up)'
            print '                   counts' + contents_line
            total = int(self.overall_mute_counts.integral(include_overflows=True))  # underflow bin is zero mutations, and we want overflow, too NOTE why the fuck isn't this quite equal to sum(self.gene_obs_counts.values())?
            for n_mutes in range(self.args.n_max_snps):
                if self.overall_mute_counts.bin_contents[n_mutes] / float(total) < self.min_fraction_per_bin:
                    print '        not looking for new alleles: not enough counts (%d / %d = %.3f < %.3f)' % (self.overall_mute_counts.bin_contents[n_mutes], total, self.overall_mute_counts.bin_contents[n_mutes] / float(total), self.min_fraction_per_bin),
                    print 'with %d %s mutations (override this with --always-find-new-alleles, or by reducing --n-max-snps)' % (n_mutes, region)
                    self.finalized = True
                    return

        start = time.time()
        self.xyvals = {}
        self.positions_to_plot = {gene : set() for gene in self.counts}
        for gene in sorted(self.counts):
            if debug:
                sys.stdout.flush()
                print ' %s: %d %s (%d unmutated)' % (utils.color_gene(gene), self.gene_obs_counts[gene], utils.plural_str('observation', self.gene_obs_counts[gene]), self.per_gene_mute_counts[gene].bin_contents[0])
                print '          skipping',
                print '%d seqs that are too highly mutated,' % self.n_seqs_too_highly_mutated[gene],
                print '%d that had 5p deletions larger than %d,' % (self.n_big_del_skipped['5p'][gene], self.n_bases_to_exclude['5p'][gene]),
                print 'and %d that had 3p deletions larger than %d)' % (self.n_big_del_skipped['3p'][gene], self.n_bases_to_exclude['3p'][gene])

            if self.gene_obs_counts[gene] < self.n_total_min:
                continue

            positions = sorted(self.counts[gene])
            self.xyvals[gene] = {pos : self.get_allele_finding_xyvals(gene, pos) for pos in positions}
            positions_to_try_to_fit = [pos for pos in positions if sum(self.xyvals[gene][pos]['obs']) > self.n_muted_min or sum(self.xyvals[gene][pos]['total']) > self.n_total_min]  # ignore positions with neither enough mutations nor total observations

            if debug and len(positions) > len(positions_to_try_to_fit):
                self.print_skip_debug(gene, positions, positions_to_try_to_fit)

            if len(positions_to_try_to_fit) < self.args.n_max_snps:
                if debug:
                    print '          not enough positions with enough observations to fit %s' % utils.color_gene(gene)
                continue

            if self.per_gene_mute_counts[gene].bin_contents[0] > self.n_total_min:  # UPDATE this will have to change to accomodate repertoires with very few unmutated sequences
                self.alleles_with_evidence.add(gene)

            # loop over each snp hypothesis
            self.fitfos[gene] = {n : {} for n in ('min_snp_ratios', 'mean_snp_ratios', 'candidates', 'fitfos')}
            not_enough_candidates = []  # just for dbg printing
            for istart in range(1, self.args.n_max_snps + 1):
                self.fit_two_piece_istart(gene, istart, positions_to_try_to_fit, print_dbg_header=(istart==self.n_snps_to_switch_to_two_piece_method), debug=debug)

                if istart not in self.fitfos[gene]['candidates']:  # just for dbg printing
                    not_enough_candidates.append(istart)

            if debug and len(not_enough_candidates) > 0:
                print '      not enough candidates for istarts: %s' % ' '.join([str(i) for i in not_enough_candidates])

            if debug and len(self.fitfos[gene]['candidates']) > 0:
                print '  evaluating each snp hypothesis'
                print '    snps       min ratio'
            istart_candidates = []
            for istart in self.fitfos[gene]['candidates']:  # note that not all <istart>s get added to self.fitfos[gene]
                if debug:
                    print '    %2d     %9s' % (istart, fstr(self.fitfos[gene]['min_snp_ratios'][istart])),
                if self.is_a_candidate(gene, self.fitfos[gene], istart, debug=debug):
                    istart_candidates.append(istart)
                if debug:
                    print ''

            # first take the biggest one, then if there's any others that have entirely non-overlapping positions, we don't need to re-run
            already_used_positions = set()
            n_new_alleles_for_this_gene = 0  # kinda messy way to implement this
            for istart in sorted(istart_candidates, reverse=True):
                if n_new_alleles_for_this_gene >= self.args.n_max_new_alleles_per_gene_per_iteration:
                    print '    skipping any additional new alleles for this gene (already have %d)' % n_new_alleles_for_this_gene
                    break
                these_positions = set(self.fitfos[gene]['candidates'][istart])
                if len(these_positions & already_used_positions) > 0:
                    continue
                already_used_positions |= these_positions
                n_new_alleles_for_this_gene += 1
                self.add_new_allele(gene, self.fitfos[gene], istart, debug=debug)

        if debug:
            if len(self.new_allele_info) > 0:
                print '  found %d new %s: %s' % (len(self.new_allele_info), utils.plural_str('allele', len(self.new_allele_info)), ' '.join([utils.color_gene(nfo['gene']) for nfo in self.new_allele_info]))
            else:
                print '    no new alleles'
            print '  allele finding: %d fits in %.1f sec' % (self.n_fits, time.time()-start)

        self.finalized = True

    # ----------------------------------------------------------------------------------------
    def plot(self, base_plotdir, only_csv=False):
        if not self.finalized:
            self.finalize(debug=debug)

        plotdir = base_plotdir + '/allele-finding'
        if self.itry is not None:
            plotdir = plotdir + '/try-' + str(self.itry)

        print '    plotting allele finding',
        sys.stdout.flush()

        for old_gene_dir in glob.glob(plotdir + '/*'):  # has to be a bit more hackey than elsewhere, since we have no way of knowing what genes might have had their own directories written last time we wrote to this dir
            if not os.path.isdir(old_gene_dir):
                raise Exception('not a directory: %s' % old_gene_dir)
            utils.prep_dir(old_gene_dir, wildlings=('*.csv', '*.svg'))
            os.rmdir(old_gene_dir)
        utils.prep_dir(plotdir, wildlings=('*.csv', '*.svg'))

        if only_csv:  # not implemented
            print '    <only_csv> not yet implemented in allelefinder'
            return

        start = time.time()
        template_genes = [newfo['template-gene'] for newfo in self.new_allele_info]
        for gene in self.positions_to_plot:  # we can make plots for the positions we didn't fit, but there's a *lot* of them and they're slow
            snp_positions = [] if gene not in template_genes else self.new_allele_info[template_genes.index(gene)]['snp-positions']
            for position in self.positions_to_plot[gene]:
                fitfos = None
                # snp_positions = [10, 11, 12, 13, 14]
                if position in snp_positions and len(snp_positions) in self.fitfos[gene]['fitfos'] and position in self.fitfos[gene]['fitfos'][len(snp_positions)]:  # not sure why I need the second and third one
                    fitfos = self.fitfos[gene]['fitfos'][len(snp_positions)][position]
                plotting.make_allele_finding_plot(plotdir + '/' + utils.sanitize_name(gene), gene, position, self.xyvals[gene][position], xmax=self.args.n_max_mutations_per_segment, fitfos=fitfos)

        check_call(['./bin/permissify-www', plotdir])
        print '(%.1f sec)' % (time.time()-start)
