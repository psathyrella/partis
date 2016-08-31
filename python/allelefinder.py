import time
import math
import sys
import os
import operator
from subprocess import check_call
import scipy
import glob
import numpy

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

        self.n_max_snps = 20  # max number of snps, i.e. try excluding up to this many bins on the left
        self.n_max_mutations_per_segment = 30  # don't look at sequences whose v segments have more than this many mutations
        self.max_fit_length = 99999  # UPDATE nevermind, I no longer think there's a reason not to fit the whole thing OLD: don't fit more than this many bins for each <istart> (the first few positions in the fit are the most important, and if we fit too far to the right these important positions get diluted) UPDATE I'm no longer so sure that I shouldn't fit the whole shebang 

        self.n_snps_to_switch_to_two_piece_method = 5

        self.n_muted_min = 30  # don't fit positions that have fewer total mutations than this (i.e. summed over bins)
        self.n_total_min = 150  # ...or fewer total observations than this

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

        self.n_fits = 0

        self.default_initial_glfo = None
        if self.args.default_initial_germline_dir is not None:  # if this is set, we want to take any new allele names from this directory's glfo if they're in there
            self.default_initial_glfo = glutils.read_glfo(self.args.default_initial_germline_dir, glfo['chain'])

        self.finalized = False

    # ----------------------------------------------------------------------------------------
    def init_gene(self, gene):
        self.counts[gene] = {}
        for igl in range(len(self.glfo['seqs'][utils.get_region(gene)][gene])):
            self.counts[gene][igl] = {}
            for istart in range(self.n_max_mutations_per_segment + 1):  # istart and n_mutes are equivalent
                self.counts[gene][igl][istart] = {n : 0 for n in ['muted', 'total'] + utils.nukes}
        self.gene_obs_counts[gene] = 0
        self.n_seqs_too_highly_mutated[gene] = 0

    # ----------------------------------------------------------------------------------------
    def get_seqs(self, info, region):
        germline_seq = info[region + '_gl_seq']
        assert len(info['seqs']) == 1
        query_seq = info[region + '_qr_seqs'][0]
        assert len(germline_seq) == len(query_seq)

        # don't use the leftmost and rightmost few bases
        germline_seq = germline_seq[self.args.new_allele_excluded_bases[0] : len(germline_seq) - self.args.new_allele_excluded_bases[1]]
        query_seq = query_seq[self.args.new_allele_excluded_bases[0] : len(query_seq) - self.args.new_allele_excluded_bases[1]]
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

            n_mutes, germline_seq, query_seq = self.get_seqs(info, region)

            if n_mutes > self.n_max_mutations_per_segment:
                self.n_seqs_too_highly_mutated[gene] += 1
                continue

            for ipos in range(len(germline_seq)):
                igl = ipos + int(info[region + '_5p_del']) + self.args.new_allele_excluded_bases[0]  # position in original (i.e. complete) germline gene

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
    def get_curvefit(self, n_mutelist, freqs, errs, y_icpt_bounds):
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

        if debug:
            print '    candidate',
        return True

    # ----------------------------------------------------------------------------------------
    def approx_x_icpt(self, pvals, debug=False):
        # NOTE ignores weights... but that's ok
        x, y = pvals['n_mutelist'], pvals['freqs']  # tmp shorthand
        slopes = [(y[i] - y[i-1]) / (x[i] - x[i-1]) for i in range(1, len(x))]  # only uses adjacent points, and double-counts interior points, but we don't care
        mean_slope = numpy.average(slopes)
        mean_xval = numpy.average(pvals['n_mutelist'])
        mean_yval = numpy.average(pvals['freqs'])
        if debug:
            print ' %.3f - %.3f * %.3f = %.3f' % (mean_yval, mean_slope, mean_xval, mean_yval - mean_slope * mean_xval)
        return mean_yval - mean_slope * mean_xval

    # ----------------------------------------------------------------------------------------
    def approx_slope_and_error(self, pvals, debug=False):
        # NOTE ignores weights... but that's ok
        x, y, w = pvals['n_mutelist'], pvals['freqs'], pvals['weights']  # tmp shorthand
        slopes = [(y[i] - y[i-1]) / (x[i] - x[i-1]) for i in range(1, len(x))]  # only uses adjacent points, and double-counts interior points, but we don't care (we don't use steps of two, because then we'd the last one if it's odd-length)
        weights = [numpy.mean([w[i-1], w[i]]) for i in range(1, len(x))]
        mean_slope = numpy.average(slopes, weights=weights)
        err = numpy.std(slopes, ddof=1) / math.sqrt(len(x))  # uh, I think <x> gives us the right number of independent measurements. In any case, this is *very* approximate
        return mean_slope, err

    # ----------------------------------------------------------------------------------------
    def fit_two_piece_istart(self, gene, istart, positions_to_try_to_fit, fitfo, print_dbg_header=False, debug=False):
        if debug and print_dbg_header:
            print '             position   ratio       (one piece / two pieces)'

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

            # if rough estimates of the slopes are roughly compatible (or preslope is bigger than postslope), the fits aren't going to tell a very different story
            preslope, preerr = self.approx_slope_and_error(prevals)
            postslope, posterr = self.approx_slope_and_error(postvals)
            fac = 0.5

            if preslope > postslope or preslope + fac*preerr > postslope - fac*posterr:
                continue

            onefit = self.get_curvefit(bothvals['n_mutelist'], bothvals['freqs'], bothvals['errs'], y_icpt_bounds=(0., 1.))

            # don't bother with the two-piece fit if the one-piece fit is pretty good
            if onefit['residuals_over_ndof'] < self.min_zero_icpt_residual:
                continue

            prefit = self.get_curvefit(prevals['n_mutelist'], prevals['freqs'], prevals['errs'], y_icpt_bounds=(0., 1.))
            postfit = self.get_curvefit(postvals['n_mutelist'], postvals['freqs'], postvals['errs'], y_icpt_bounds=(0., 1.))
            twofit_residuals = prefit['residuals_over_ndof'] * prefit['ndof'] + postfit['residuals_over_ndof'] * postfit['ndof']
            twofit_ndof = prefit['ndof'] + postfit['ndof']
            twofit_residuals_over_ndof = twofit_residuals / twofit_ndof
            # tmpstr = ''
            # if onefit['residuals_over_ndof'] / twofit_residuals_over_ndof > 2.25:
            #     tmpstr = utils.color('yellow', 'yes')
            # print '   %3d  %5s   %5s / %-5s  %3s  (%-5s  %5s)' % (pos, fstr(onefit['residuals_over_ndof'] / twofit_residuals_over_ndof), fstr(onefit['residuals_over_ndof']), fstr(twofit_residuals_over_ndof), tmpstr, fstr(prefit['residuals_over_ndof']), fstr(postfit['residuals_over_ndof']))

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
                for n_mutes in range(1, self.n_max_mutations_per_segment + 1):
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
            print '%0s %s' % ('', ''.join(['%11d' % nm for nm in range(1, self.n_max_mutations_per_segment + 1)]))

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
            if self.approx_x_icpt(pvals) < 0.:
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
                for n_mutes in range(1, self.n_max_mutations_per_segment + 1):
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
        left, right = self.args.new_allele_excluded_bases
        newpos = self.glfo[utils.conserved_codons[chain][region] + '-positions'][template_gene]  # codon position for template gene should be ok
        for oldname_gene, oldname_seq in self.default_initial_glfo['seqs'][region].items():  # NOTE <oldname_{gene,seq}> is the old *name* corresponding to the new (snp'd) allele, whereas <old_seq> is the allele from which we inferred the new (snp'd) allele
            # first see if they match up through the cysteine
            oldpos = self.default_initial_glfo[utils.conserved_codons[chain][region] + '-positions'][oldname_gene]
            if oldname_seq[left : oldpos + 3] != new_seq[left : newpos + 3]:
                continue

            # then require that any bases in common to the right of the cysteine in the new allele match the ones in the old one (where "in common" means either of them can be longer, since this just changes the insertion length)
            bases_to_right_of_cysteine = min(len(oldname_seq) - (oldpos + 3), len(new_seq) - right - (newpos + 3))

            if bases_to_right_of_cysteine > 0 and oldname_seq[oldpos + 3 : oldpos + 3 + bases_to_right_of_cysteine] != new_seq[newpos + 3 : newpos + 3 + bases_to_right_of_cysteine]:
                continue

            print '        using old name %s for new allele %s:' % (utils.color_gene(oldname_gene), utils.color_gene(new_name))
            def print_sequence_chunks(seq, cpos, name):
                print '            %s%s%s%s%s   %s' % (utils.color('red', seq[:left]),
                                                       seq[left : cpos],
                                                       utils.color('reverse_video', seq[cpos : cpos + 3]),
                                                       seq[cpos + 3 : cpos + 3 + bases_to_right_of_cysteine],
                                                       utils.color('red', seq[cpos + 3 + bases_to_right_of_cysteine:]),
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
            new_len_str = utils.color('red', new_seq[len(old_seq):])
            new_seq_for_cf = new_seq[:len(old_seq)]
        elif len(old_seq) > len(new_seq):
            old_len_str = utils.color('red', old_seq[len(new_seq):])
            old_seq_for_cf = old_seq[:len(new_seq)]
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
        too_close_to_ends = range(self.args.new_allele_excluded_bases[0]) + range(len(glseq) - self.args.new_allele_excluded_bases[1], len(glseq))
        not_enough_counts = set(positions) - set(positions_to_try_to_fit) - set(too_close_to_ends)  # well, not enough counts, *and* not too close to the ends

        print '          skipping %d / %d positions:' % (len(positions) - len(positions_to_try_to_fit), len(positions)),
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
                print '  %s observed %d %s, %d too highly mutated' % (utils.color_gene(gene, width=15), self.gene_obs_counts[gene], utils.plural_str('time', self.gene_obs_counts[gene]), self.n_seqs_too_highly_mutated[gene])

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
