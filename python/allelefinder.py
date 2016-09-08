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
    def __init__(self, glfo, args, itry):
        self.glfo = glfo
        self.args = args
        self.itry = itry

        self.fraction_of_seqs_to_exclude = 0.01  # exclude the fraction of sequences with largest v_3p deletions whose counts add up to this fraction of total sequences NOTE you don't want to make this too big, because although you'll be removing all the seqs with large 4p deletions, this number also gets used when you're deciding whether your new allele is in the default glfo
        self.n_bases_to_exclude = {'5p' : {}, '3p' : {}}  # i.e. on all the seqs we keep, we exclude this many bases; and any sequences that have larger deletions than this are not kept

        self.n_max_snps = 20  # max number of snps, i.e. try excluding up to this many bins on the left
        self.n_max_mutations_per_segment = 30  # don't look at sequences whose v segments have more than this many mutations
        self.max_fit_length = 99999  # UPDATE nevermind, I no longer think there's a reason not to fit the whole thing OLD: don't fit more than this many bins for each <istart> (the first few positions in the fit are the most important, and if we fit too far to the right these important positions get diluted) UPDATE I'm no longer so sure that I shouldn't fit the whole shebang 

        self.n_snps_to_switch_to_two_piece_method = 5

        self.n_muted_min = 30  # don't fit positions that have fewer total mutations than this (i.e. summed over bins)
        self.n_total_min = 150  # ...or fewer total observations than this
        self.n_muted_min_per_bin = 8 # <istart>th bin has to have at least this many mutated sequences (i.e. 2-3 sigma from zero)

        self.min_min_candidate_ratio = 2.25  # every candidate ratio must be greater than this
        self.min_mean_candidate_ratio = 2.75  # mean of candidate ratios must be greater than this
        self.min_zero_icpt_residual = 2.  # TODO rename this

        self.min_min_candidate_ratio_to_plot = 1.5  # don't plot positions that're below this (for all <istart>)

        self.default_slope_bounds = (-0.1, 0.2)  # fitting function needs some reasonable bounds from which to start (I could at some point make slope part of the criteria for candidacy, but it wouldn't add much sensitivity)

        self.counts = {}
        self.new_allele_info = []
        self.positions_to_plot = {}
        self.n_seqs_too_highly_mutated = {}  # sequences (per-gene) that had more than <self.n_max_mutations_per_segment> mutations
        self.gene_obs_counts = {}
        self.n_big_del_skipped = {s : {} for s in self.n_bases_to_exclude}

        self.n_fits = 0

        self.default_initial_glfo = None
        if self.args.default_initial_germline_dir is not None:  # if this is set, we want to take any new allele names from this directory's glfo if they're in there
            self.default_initial_glfo = glutils.read_glfo(self.args.default_initial_germline_dir, glfo['chain'])

        self.finalized = False

        self.reflengths = {}

    # ----------------------------------------------------------------------------------------
    def init_gene(self, gene):
        self.counts[gene] = {}
        for igl in range(len(self.glfo['seqs'][utils.get_region(gene)][gene])):
            self.counts[gene][igl] = {}
            for istart in range(self.n_max_mutations_per_segment + 1):  # istart and n_mutes are equivalent
                self.counts[gene][igl][istart] = {n : 0 for n in ['muted', 'total'] + utils.nukes}
        self.gene_obs_counts[gene] = 0
        for side in self.n_big_del_skipped:
            self.n_big_del_skipped[side][gene] = 0
        self.n_seqs_too_highly_mutated[gene] = 0

    # ----------------------------------------------------------------------------------------
    def set_excluded_bases(self, swfo, debug=True):
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
                if debug:
                    print gene
                    print '  observed %s deletions %s' % (side, ' '.join([str(d) for d in observed_deletions]))
                total_obs = sum(dcounts[gene].values())
                running_sum = 0
                for dlen in observed_deletions:
                    if debug:
                        print '    %4d %3d / %3d = %5.3f' % (dlen, float(running_sum), total_obs, float(running_sum) / total_obs)
                    self.n_bases_to_exclude[side][gene] = dlen  # setting this before the "if" means that if we fall through (e.g. if there aren't enough sequences to get above the threshold) we'll still have a reasonable default
                    if float(running_sum) / total_obs > 1. - self.fraction_of_seqs_to_exclude:  # if we've already added deletion lengths accounting for most of the sequences, ignore the rest of 'em
                        if debug:
                            print '                 choose', dlen
                        break
                    running_sum += dcounts[gene][dlen]

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

            # i.e. do *not* use <info> after this point

            if n_mutes > self.n_max_mutations_per_segment:
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
    def get_residual_sum(self, xvals, yvals, errs, slope, intercept):
        def expected(x):
            return slope * x + intercept
        residual_sum = sum([(y - expected(x))**2 / err**2 for x, y, err in zip(xvals, yvals, errs)])
        return residual_sum

    # ----------------------------------------------------------------------------------------
    def dbgstr(self, slope, slope_err, y_icpt, y_icpt_err, extra_str=''):
        return '        %s  m: %5.3f +/- %5.3f   b: %5.3f +/- %5.3f' % (extra_str, slope, slope_err, y_icpt, y_icpt_err)

    # ----------------------------------------------------------------------------------------
    def get_curvefit(self, n_mutelist, freqs, errs, y_icpt_bounds, debug=False):
        # bounds = (-float('inf'), float('inf'))
        if y_icpt_bounds == (0., 0.):
            def linefunc(x, slope):
                return slope*x
            bounds = self.default_slope_bounds
            params, cov = scipy.optimize.curve_fit(linefunc, n_mutelist, freqs, sigma=errs, bounds=bounds)
            slope, slope_err = params[0], math.sqrt(cov[0][0])
            y_icpt, y_icpt_err = 0., 0.
            ndof = len(n_mutelist) - 1
        else:
            def linefunc(x, slope, y_icpt):
                return slope*x + y_icpt
            bounds = [[s, y] for s, y in zip(self.default_slope_bounds, y_icpt_bounds)]
            params, cov = scipy.optimize.curve_fit(linefunc, n_mutelist, freqs, sigma=errs, bounds=bounds)
            slope, slope_err = params[0], math.sqrt(cov[0][0])
            y_icpt, y_icpt_err = params[1], math.sqrt(cov[1][1])
            ndof = len(n_mutelist) - 2

        residual_sum = self.get_residual_sum(n_mutelist, freqs, errs, slope, y_icpt)

        if debug:
            print self.dbgstr(slope, slope_err, y_icpt, y_icpt_err, extra_str='fit')

        fitfo = {
            'slope'  : slope,
            'y_icpt' : y_icpt,
            'slope_err'  : slope_err,
            'y_icpt_err' : y_icpt_err,
            'residuals_over_ndof' : float(residual_sum) / ndof,
            'ndof' : ndof,
            'y_icpt_bounds' : y_icpt_bounds
        }
        self.n_fits += 1
        return fitfo

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

        for candidate_pos in fitfo['candidates'][istart]:  # return false if any of the candidate positions don't have enough mutated counts in the <istart>th bin NOTE this is particularly important because it's the handle that tells us it's *this* <istart> that's correct, rather than <istart> + 1
            n_istart_muted = self.counts[gene][candidate_pos][istart]['muted']
            if n_istart_muted < self.n_muted_min_per_bin:
                if debug:
                    print '    not enough mutated counts at candidate position %d with %d %s (%d < %d)' % (candidate_pos, istart, utils.plural_str('mutations', n_istart_muted), n_istart_muted, self.n_muted_min_per_bin),
                return False

        if debug:
            print '    candidate',
        return True

    # ----------------------------------------------------------------------------------------
    def approx_fit_vals(self, pvals, debug=False):
        def getslope(i1, i2):
            return (y[i2] - y[i1]) / (x[i2] - x[i1])

        # NOTE uncertainties are kinda complicated if you do the weighted mean, so screw it, it works fine with the plain mean
        x, y = pvals['n_mutelist'], pvals['freqs']  # tmp shorthand
        slopes = [getslope(i-1, i) for i in range(1, len(x))]  # only uses adjacent points, and double-counts interior points, but we don't care (we don't use steps of two, because then we'd the last one if it's odd-length)
        slope = numpy.average(slopes)
        slope_err = numpy.std(slopes, ddof=1) / math.sqrt(len(x))

        y_icpts = [y[i] - getslope(i-1, i) * x[i] for i in range(1, len(x))]
        y_icpt = numpy.average(y_icpts)
        y_icpt_err = numpy.std(y_icpts, ddof=1) / math.sqrt(len(x))

        if debug:
            # self.get_curvefit(pvals['n_mutelist'], pvals['freqs'], pvals['errs'], y_icpt_bounds=(-1., 1.), debug=True)
            print self.dbgstr(slope, slope_err, y_icpt, y_icpt_err, extra_str='apr')

        return {'m' : slope, 'm_err' : slope_err, 'b' : y_icpt, 'b_err' : y_icpt_err}

    # ----------------------------------------------------------------------------------------
    def consistent(self, v1, v1err, v2, v2err, debug=False):
        factor = 1.  # i.e. if both slope and intercept are within <factor> std deviations of each other, don't bother fitting, because the fit isn't going to say they're wildly inconsistent
        lo, hi = sorted([v1, v2])
        joint_err = max(v1err, v2err)
        if debug:
            print '      %6.3f +/- %6.3f   %6.3f +/- %6.3f   -->   %6.3f + %3.1f * %6.3f = %6.3f >? %6.3f   %s' % (v1, v1err, v2, v2err, lo, factor, joint_err, lo + factor * joint_err, hi, lo + factor * joint_err > hi)
        return lo + factor * joint_err > hi

    # ----------------------------------------------------------------------------------------
    def consistent_slope_and_y_icpt(self, vals1, vals2, debug=False):
        consistent_slopes = self.consistent(vals1['m'], vals1['m_err'], vals2['m'], vals2['m_err'], debug=debug)
        consistent_y_icpts = self.consistent(vals1['b'], vals1['b_err'], vals2['b'], vals2['b_err'], debug=debug)
        return consistent_slopes and consistent_y_icpts

    # ----------------------------------------------------------------------------------------
    def fit_two_piece_istart(self, gene, istart, positions_to_try_to_fit, fitfo, print_dbg_header=False, debug=False):
        if debug and print_dbg_header:
            print '             position   ratio       (one piece / two pieces)'

        # NOTE I'm including the zero bin here -- do I really want to do that? UPDATE yes, I think so -- we know there will be zero mutations in that bin, but the number of sequences in it still contains information (uh, I think)
        prexyvals = {pos : {k : v[:istart] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}  # arrays up to, but not including, <istart>
        postxyvals = {pos : {k : v[istart : istart + self.max_fit_length] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}  # arrays from <istart> onwards
        bothxyvals = {pos : {k : v[:istart + self.max_fit_length] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}

        candidate_ratios, residfo = {}, {}  # NOTE <residfo> is really just for dbg printing... but we have to sort before we print, so we need to keep the info around
        for pos in positions_to_try_to_fit:
            prevals = prexyvals[pos]
            postvals = postxyvals[pos]
            bothvals = bothxyvals[pos]

            if sum(postvals['obs']) < self.n_muted_min or sum(postvals['total']) < self.n_total_min:
                continue

            # if both slope and intercept are quite close to each other, the fits aren't going to say they're wildly inconsistent
            pre_approx = self.approx_fit_vals(prevals)
            post_approx = self.approx_fit_vals(postvals)
            if pre_approx['m'] > post_approx['m'] or self.consistent_slope_and_y_icpt(pre_approx, post_approx):
                continue

            onefit = self.get_curvefit(bothvals['n_mutelist'], bothvals['freqs'], bothvals['errs'], y_icpt_bounds=(-1., 1.))

            # don't bother with the two-piece fit if the one-piece fit is pretty good
            if onefit['residuals_over_ndof'] < self.min_zero_icpt_residual:
                continue

            prefit = self.get_curvefit(prevals['n_mutelist'], prevals['freqs'], prevals['errs'], y_icpt_bounds=(-1., 1.))
            postfit = self.get_curvefit(postvals['n_mutelist'], postvals['freqs'], postvals['errs'], y_icpt_bounds=(-1., 1.))
            twofit_residuals = prefit['residuals_over_ndof'] * prefit['ndof'] + postfit['residuals_over_ndof'] * postfit['ndof']
            twofit_ndof = prefit['ndof'] + postfit['ndof']
            twofit_residuals_over_ndof = twofit_residuals / twofit_ndof

            # TODO add a requirement that all the candidate slopes are consistent both pre and post

            candidate_ratios[pos] = onefit['residuals_over_ndof'] / twofit_residuals_over_ndof if twofit_residuals_over_ndof > 0. else float('inf')
            residfo[pos] = {'onefo' : onefit, 'prefo' : prefit, 'postfo' : postfit, 'twofo' : {'residuals_over_ndof' : twofit_residuals_over_ndof}}
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
                if candidate_ratios[pos] > self.min_min_candidate_ratio and residfo[pos]['onefo']['residuals_over_ndof'] > self.min_zero_icpt_residual:
                    pos_str = utils.color('yellow', pos_str)
                print_str = ['                 %s    %5s            %5s / %-5s               ' % (pos_str, fstr(candidate_ratios[pos]),
                                                                                                  fstr(residfo[pos]['onefo']['residuals_over_ndof']), fstr(residfo[pos]['twofo']['residuals_over_ndof']))]
                for n_mutes in range(self.n_max_mutations_per_segment + 1):
                    if n_mutes in bothxyvals[pos]['n_mutelist']:
                        inm = bothxyvals[pos]['n_mutelist'].index(n_mutes)
                        print_str.append('%4d / %-4d' % (bothxyvals[pos]['obs'][inm], bothxyvals[pos]['total'][inm]))
                    else:
                        print_str.append('           ')
                print ''.join(print_str)

        if len(candidates) < istart:
            return

        fitfo['min_snp_ratios'][istart] = min(candidate_ratios.values())
        fitfo['mean_snp_ratios'][istart] = numpy.mean(candidate_ratios.values())
        fitfo['candidates'][istart] = candidate_ratios

    # ----------------------------------------------------------------------------------------
    def fit_istart(self, gene, istart, positions_to_try_to_fit, fitfo, print_dbg_header=False, debug=False):
        if debug and print_dbg_header:
            print '             position   ratio    (m=0 / m=big)      big bounds',
            print '%0s %s' % ('', ''.join(['%11d' % nm for nm in range(self.n_max_mutations_per_segment + 1)]))  # NOTE *has* to correspond to line at bottom of fcn below

        subxyvals = {pos : {k : v[istart : istart + self.max_fit_length] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}  # arrays from <istart> onwards

        candidate_ratios, residfo = {}, {}  # NOTE <residfo> is really just for dbg printing... but we have to sort before we print, so we need to keep the info around
        for pos in positions_to_try_to_fit:
            pvals = subxyvals[pos]

            if sum(pvals['obs']) < self.n_muted_min or sum(pvals['total']) < self.n_total_min:
                continue

            big_y_icpt = numpy.average(pvals['freqs'], weights=pvals['weights'])  # corresponds, roughly, to the expression level of the least common allele to which we have sensitivity NOTE <istart> is at index 0
            big_y_icpt_err = numpy.std(pvals['freqs'], ddof=1)  # NOTE this "err" is from the variance over bins, and ignores the sample statistics of each bin. This is a little weird, but good: it captures cases where the points aren't well-fit by a line, either because of multiple alleles with very different prevalences, or because the sequences aren't very independent NOTE the former case is very important)
            big_y_icpt_bounds = (big_y_icpt - 1.5*big_y_icpt_err, big_y_icpt + 1.5*big_y_icpt_err)  # we want the bounds to be lenient enough to accomodate non-zero slopes (in the future, we could do something cleverer like extrapolating with the slope of the line to x=0)

            # if the bounds include zero, there won't be much difference between the two fits
            if big_y_icpt_bounds[0] <= 0.:
                continue

            # if a rough estimate of the x-icpt is less than zero, the zero-icpt fit is probably going to be pretty good
            _, _, approx_y_icpt, _ = self.approx_fit_vals(pvals)
            if approx_y_icpt < 0.:
                continue

            zero_icpt_fit = self.get_curvefit(pvals['n_mutelist'], pvals['freqs'], pvals['errs'], y_icpt_bounds=(0., 0.))

            # don't bother with the big-icpt fit if the zero-icpt fit is pretty good
            if zero_icpt_fit['residuals_over_ndof'] < self.min_zero_icpt_residual:
                continue

            big_icpt_fit = self.get_curvefit(pvals['n_mutelist'], pvals['freqs'], pvals['errs'], y_icpt_bounds=big_y_icpt_bounds)

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
                if candidate_ratios[pos] > self.min_min_candidate_ratio and residfo[pos]['zerofo']['residuals_over_ndof'] > self.min_zero_icpt_residual:
                    pos_str = utils.color('yellow', pos_str)
                print_str = ['                 %s    %5s    %5s / %-5s   [%5.3f - %5.3f] ' % (pos_str, fstr(candidate_ratios[pos]),
                                                                                              fstr(residfo[pos]['zerofo']['residuals_over_ndof']), fstr(residfo[pos]['bigfo']['residuals_over_ndof']),
                                                                                              residfo[pos]['bigfo']['y_icpt_bounds'][0], residfo[pos]['bigfo']['y_icpt_bounds'][1])]
                print_str += '    '
                for n_mutes in range(self.n_max_mutations_per_segment + 1):
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

        start = time.time()
        self.xyvals = {}
        self.positions_to_plot = {gene : set() for gene in self.counts}
        print '%s: looking for new alleles' % utils.color('red', 'try ' + str(self.itry))
        for gene in sorted(self.counts):
            if debug:
                sys.stdout.flush()
                print '  %s observed %d %s (ignoring %d of these that were too highly mutated,' % (utils.color_gene(gene, width=15), self.gene_obs_counts[gene], utils.plural_str('time', self.gene_obs_counts[gene]), self.n_seqs_too_highly_mutated[gene]),
                print '%d that had 5p deletions larger than %d,' % (self.n_big_del_skipped['5p'][gene], self.n_bases_to_exclude['5p'][gene]),
                print 'and %d that had 3p deletions larger than %d)' % (self.n_big_del_skipped['3p'][gene], self.n_bases_to_exclude['3p'][gene])

            if self.gene_obs_counts[gene] < self.n_total_min:
                continue

            positions = sorted(self.counts[gene])
            self.xyvals[gene] = {pos : self.get_allele_finding_xyvals(gene, pos) for pos in positions}
            positions_to_try_to_fit = [pos for pos in positions if sum(self.xyvals[gene][pos]['obs']) > self.n_muted_min or sum(self.xyvals[gene][pos]['total']) > self.n_total_min]  # ignore positions with neither enough mutations nor total observations

            if debug and len(positions) > len(positions_to_try_to_fit):
                self.print_skip_debug(gene, positions, positions_to_try_to_fit)

            if len(positions_to_try_to_fit) < self.n_max_snps:
                if debug:
                    print '          not enough positions with enough observations to fit %s' % utils.color_gene(gene)
                continue

            # loop over each snp hypothesis
            fitfo = {n : {} for n in ('min_snp_ratios', 'mean_snp_ratios', 'candidates')}
            not_enough_candidates = []  # just for dbg printing
            for istart in range(1, self.n_max_snps + 1):
                if istart < self.n_snps_to_switch_to_two_piece_method:
                    self.fit_istart(gene, istart, positions_to_try_to_fit, fitfo, print_dbg_header=(istart==1), debug=debug)
                else:
                    self.fit_two_piece_istart(gene, istart, positions_to_try_to_fit, fitfo, print_dbg_header=(istart==self.n_snps_to_switch_to_two_piece_method), debug=debug)

                if istart not in fitfo['candidates']:  # just for dbg printing
                    not_enough_candidates.append(istart)

            if debug and len(not_enough_candidates) > 0:
                print '      not enough candidates for istarts: %s' % ' '.join([str(i) for i in not_enough_candidates])

            if debug and len(fitfo['candidates']) > 0:
                print '  evaluating each snp hypothesis'
                print '    snps       min ratio'
            istart_candidates = []
            for istart in fitfo['candidates']:  # note that not all <istart>s get added to fitfo
                if debug:
                    print '    %2d     %9s' % (istart, fstr(fitfo['min_snp_ratios'][istart])),
                if self.is_a_candidate(gene, fitfo, istart, debug=debug):
                    if debug and len(istart_candidates) == 0:
                        print '   %s' % utils.color('yellow', '(best)'),
                    istart_candidates.append(istart)
                print ''

            if len(istart_candidates) > 0:
                n_candidate_snps = min(istart_candidates)  # add the candidate with the smallest number of snps to the germline set, and run again
                self.add_new_allele(gene, fitfo, n_candidate_snps, debug=debug)

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
        for gene in self.positions_to_plot:  # we can make plots for the positions we didn't fit, but there's a *lot* of them and they're slow
            for position in self.positions_to_plot[gene]:
                plotting.make_allele_finding_plot(plotdir + '/' + utils.sanitize_name(gene), gene, position, self.xyvals[gene][position])

        check_call(['./bin/permissify-www', plotdir])
        print '(%.1f sec)' % (time.time()-start)
