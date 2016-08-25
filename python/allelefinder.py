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

        self.n_max_snps = 10  # max number of snps, i.e. try excluding up to this many bins on the left
        self.n_max_mutations_per_segment = 20  # don't look at sequences whose v segments have more than this many mutations
        self.min_fit_length = 5  # don't both fitting an <istart> (i.e. snp hypothesis) if it doesn't have at least this many bins
        self.max_fit_length = 10  # don't fit more than this many bins for each <istart> (the first few positions in the fit are the most important, and if we fit too far to the right these important positions get diluted)

        self.n_muted_min = 10  # don't fit positions that have fewer total mutations than this (i.e. summed over bins)
        self.n_total_min = 50  # ...or fewer total observations than this
        self.n_muted_min_per_bin = 8  # <istart>th bin has to have at least this many mutated sequences (i.e. 2-3 sigma from zero)

        self.min_min_candidate_ratio = 2.25  # every candidate ratio must be greater than this
        self.min_zero_icpt_residual = 2.25

        self.min_min_candidate_ratio_to_plot = 1.25  # don't plot positions that're below this (for all <istart>)

        self.default_slope_bounds = (-0.2, 0.2)  # fitting function needs some reasonable bounds from which to start (I could at some point make slope part of the criteria for candidacy, but it wouldn't add much sensitivity)

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
            'print_str' : '    %9.3f +/- %-9.3f   %7.4f +/- %7.4f    %7.4f' % (y_icpt, y_icpt_err, slope, slope_err, float(residual_sum) / ndof),
            'n_mutelist' : n_mutelist,
            'freqs' : freqs,
            'errs' : errs
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
                print '    min snp ratio %s too small (less than %s)' % (fstr(fitfo['min_snp_ratios'][istart]), fstr(self.min_min_candidate_ratio))
            return False
        for candidate_pos in fitfo['candidates'][istart]:  # return false if any of the candidate positions don't have enough counts with <istart> mutations (probably a homozygous new allele with more than <istart> snps) UPDATE did I mean heterozygous?
            n_istart_muted = self.counts[gene][candidate_pos][istart]['muted']
            if n_istart_muted < self.n_muted_min_per_bin:
                if debug:
                    print '    not enough mutated counts at candidate position %d with %d %s (%s < %s)' % (candidate_pos, istart, utils.plural_str('mutations', n_istart_muted), fstr(n_istart_muted), fstr(self.n_muted_min_per_bin))
                return False

        if debug:
            print '    candidate'
        return True

    # ----------------------------------------------------------------------------------------
    def approx_x_icpt(self, pvals):
        # NOTE ignores weights... but that's ok
        x, y = pvals['n_mutelist'], pvals['freqs']  # tmp shorthand
        slopes = [(y[i] - y[i-1]) / (x[i] - x[i-1]) for i in range(1, len(x))]  # only uses adjacent points, and double-counts interior points, but we don't care
        mean_slope = numpy.average(slopes)
        mean_xval = numpy.average(pvals['n_mutelist'])
        mean_yval = numpy.average(pvals['freqs'])
        return mean_yval - mean_slope * mean_xval

    # ----------------------------------------------------------------------------------------
    def fit_istart(self, gene, istart, positions_to_try_to_fit, fitfo, debug=False):
        subxyvals = {pos : {k : v[istart : istart + self.max_fit_length] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}

        candidate_ratios, residfo = {}, {}  # NOTE <residfo> is really just for dbg printing... but we have to sort before we print, so we need to keep the info around
        for pos in positions_to_try_to_fit:
            pvals = subxyvals[pos]

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

            residfo[pos] = {'zero_icpt_resid' : zero_icpt_fit['residuals_over_ndof'],
                            'big_icpt_resid' :  big_icpt_fit['residuals_over_ndof'],
                            'fixed_y_icpt' : big_y_icpt,
                            'fixed_y_icpt_err' : big_y_icpt_err}

            if residfo[pos]['zero_icpt_resid'] / residfo[pos]['big_icpt_resid'] > self.min_min_candidate_ratio_to_plot:
                self.positions_to_plot[gene].add(pos)  # if we already decided to plot it for another <istart>, it'll already be in there

        candidates = [pos for pos, _ in sorted(candidate_ratios.items(), key=operator.itemgetter(1), reverse=True)]  # sort the candidate positions in decreasing order of residual ratio
        candidates = candidates[:istart]  # remove any extra ones (if we have more than we need)
        candidate_ratios = {pos : ratio for pos, ratio in candidate_ratios.items() if pos in candidates}

        if debug:
            if len(candidates) > 0:
                print '    %d %s' % (istart, utils.plural_str('snp', istart))
            for pos in candidates:
                xtrastrs = (' ', ' ')  #('[', ']') if XXX else (' ', ' ')
                pos_str = '%3s' % str(pos)
                if candidate_ratios[pos] > self.min_min_candidate_ratio and residfo[pos]['zero_icpt_resid'] > self.min_zero_icpt_residual:
                    pos_str = utils.color('yellow', pos_str)
                print_str = ['               %s %s    %5s   (%5s / %-5s)   %5.3f +/- %5.3f %s' % (xtrastrs[0], pos_str, fstr(candidate_ratios[pos]),
                                                                                                  fstr(residfo[pos]['zero_icpt_resid']), fstr(residfo[pos]['big_icpt_resid']),
                                                                                                  residfo[pos]['fixed_y_icpt'], residfo[pos]['fixed_y_icpt_err'],
                                                                                                  xtrastrs[1])]
                print_str += '    '
                for n_mutes in range(1, self.n_max_mutations_per_segment + 1):
                    if n_mutes in subxyvals[pos]['n_mutelist']:
                        inm = subxyvals[pos]['n_mutelist'].index(n_mutes)
                        print_str.append('%4d / %-4d' % (subxyvals[pos]['obs'][inm], subxyvals[pos]['total'][inm]))
                    else:
                        print_str.append('           ')
                print ''.join(print_str)

        if len(candidates) < istart:
            # if debug:
            #     print '      not enough candidates (%d < %d)' % (len(candidates), istart)
            return

        fitfo['min_snp_ratios'][istart] = min(candidate_ratios.values())
        fitfo['candidates'][istart] = candidate_ratios

    # ----------------------------------------------------------------------------------------
    def add_new_allele(self, gene, fitfo, n_candidate_snps, debug=False):
        # figure out what the new nukes are
        old_seq = self.glfo['seqs'][utils.get_region(gene)][gene]
        new_seq = old_seq
        mutfo = {}
        for pos in sorted(fitfo['candidates'][n_candidate_snps]):
            obs_counts = {nuke : self.counts[gene][pos][n_candidate_snps][nuke] for nuke in utils.nukes}  # NOTE it's super important to only use the counts from sequences with <n_candidate_snps> total mutations
            sorted_obs_counts = sorted(obs_counts.items(), key=operator.itemgetter(1), reverse=True)
            original_nuke = self.glfo['seqs'][utils.get_region(gene)][gene][pos]
            new_nuke = None
            for nuke, _ in sorted_obs_counts:  # take the most common one that isn't the existing gl nuke
                if nuke != original_nuke:
                    new_nuke = nuke
                    break
            print '   %3d  (%s --> %s)' % (pos, original_nuke, new_nuke),
            assert old_seq[pos] == original_nuke
            mutfo[pos] = {'original' : original_nuke, 'new' : new_nuke}
            new_seq = new_seq[:pos] + new_nuke + new_seq[pos+1:]
        print ''

        new_name, mutfo = glutils.get_new_allele_name_and_change_mutfo(gene, mutfo)
        if self.default_initial_glfo is not None:  # if this is set, we want to take the names from this directory's glfo
            left, right = self.args.new_allele_excluded_bases
            # see if there's already a name for <new_name>'s sequence
            # NOTE <oldname_{gene,seq}> is the old *name* corresponding to the new (snp'd) allele, whereas <old_seq> is the allele from which we inferred the new (snp'd) allele
            for oldname_gene, oldname_seq in self.default_initial_glfo['seqs'][utils.get_region(new_name)].items():
                if oldname_seq[left : len(oldname_seq) - right] == new_seq[left : len(new_seq) - right]:
                    print '        using old name %s for new allele %s' % (utils.color_gene(oldname_gene), utils.color_gene(new_name))
                    new_name = oldname_gene
                    new_seq = oldname_seq  # in case the excluded positions are different
                    break

        print '          %s   %s' % (old_seq, utils.color_gene(gene))
        print '          %s   %s' % (utils.color_mutants(old_seq, new_seq), utils.color_gene(new_name))

        # and add it to the set of new alleles for this gene
        self.new_allele_info.append({
            'template-gene' : gene,
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
            fitfo = {n : {} for n in ('min_snp_ratios', 'candidates')}
            not_enough_candidates = []  # just for dbg printing
            for istart in range(1, self.n_max_snps + 1):
                if debug and istart == 1:
                    print '                                 resid. / ndof'
                    print '             position   ratio    (m=0 / m=big)          big',
                    print '%5s %s' % ('', ''.join(['%11d' % nm for nm in range(1, self.n_max_mutations_per_segment + 1)]))

                self.fit_istart(gene, istart, positions_to_try_to_fit, fitfo, debug=debug)

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
                    istart_candidates.append(istart)

            if len(istart_candidates) > 0:
                n_candidate_snps = min(istart_candidates)  # add the candidate with the smallest number of snps to the germline set, and run again
                print '\n    found a new allele candidate separated from %s by %d %s at %s:' % (utils.color_gene(gene), n_candidate_snps,
                                                                                                utils.plural_str('snp', n_candidate_snps), utils.plural_str('position', n_candidate_snps)),
                self.add_new_allele(gene, fitfo, n_candidate_snps, debug=debug)

        if debug:
            if len(self.new_allele_info) > 0:
                print '  found %d new %s: %s' % (len(self.new_allele_info), utils.plural_str('allele', len(self.new_allele_info)), ' '.join([utils.color_gene(nfo['gene']) for nfo in self.new_allele_info]))
            else:
                print '    no new alleles'
            print '  allele finding time (%d fits): %.1f' % (self.n_fits, time.time()-start)

        self.finalized = True

    # ----------------------------------------------------------------------------------------
    def plot(self, base_plotdir, only_csv=False):
        if not self.finalized:
            self.finalize(debug=debug)

        plotdir = base_plotdir + '/allele-finding'
        if self.itry is not None:
            plotdir = plotdir + '/try-' + str(self.itry)

        print '    plotting'

        for old_gene_dir in glob.glob(plotdir + '/*'):  # has to be a bit more hackey than elsewhere, since we have no way of knowing what genes might have had their own directories written last time we wrote to this dir
            if not os.path.isdir(old_gene_dir):
                raise Exception('not a directory: %s' % old_gene_dir)
            utils.prep_dir(old_gene_dir, wildlings=('*.csv', '*.svg'))
            os.rmdir(old_gene_dir)
        utils.prep_dir(plotdir, wildlings=('*.csv', '*.svg'))

        if only_csv:  # not implemented
            print '    only_csv not yet implemented in allelefinder'
            return

        start = time.time()
        for gene in self.positions_to_plot:  # we can make plots for the positions we didn't fit, but there's a *lot* of them and they're slow
            for position in self.positions_to_plot[gene]:
                plotting.make_allele_finding_plot(plotdir + '/' + utils.sanitize_name(gene), gene, position, self.xyvals[gene][position])

        check_call(['./bin/permissify-www', plotdir])
        print '      allele finding plot time: %.1f' % (time.time()-start)
