import time
import math
import os
import operator
import scipy
import glob

import fraction_uncertainty
from mutefreqer import MuteFreqer
import plotting
import utils

# ----------------------------------------------------------------------------------------
def fstr(fval):
    if fval < 10:
        return '%.2f' % fval
    elif fval < 1e4:
        return '%.0f.' % fval
    else:
        return '%.0e' % fval

# ----------------------------------------------------------------------------------------
class AlleleFinder(object):
    def __init__(self, glfo, args):
        self.glfo = glfo
        self.args = args

        self.new_allele_info = []

        self.small_number = 1e-5
        self.n_max_mutations_per_segment = 20  # don't look at sequences whose v segments have more than this many mutations
        self.n_max_snps = self.n_max_mutations_per_segment - 9  # try excluding up to this many bins (on the left) when doing the fit (leaves at least 9 points for fit)
        self.n_muted_min = 15  # don't fit positions that have fewer mutations than this
        self.n_total_min = 15  # ...or fewer total observations than this
        self.n_five_prime_positions_to_exclude = 5  # skip positions that are too close to the 5' end of V (misassigned insertions look like snps)
        self.min_non_candidate_positions_to_fit = 10  # always fit at least a few non-candidate positions
        self.min_y_intercept = 0.15  # corresponds, roughly, to the expression level of the least common allele to which we have sensitivity
        self.default_slope_bounds = (-0.2, 0.2)  # fitting function needs some reasonable bounds from which to start
        self.big_y_icpt_bounds = (self.min_y_intercept, 1.5)  # snp-candidate positions should fit well when forced to use these bounds, but non-snp positions should fit like &*@!*
        self.min_score = 2  # (mean ratio over snp candidates) - (first non-candidate ratio) must be greater than this
        self.min_min_candidate_ratio = 2.25  # every candidate ratio must be greater than this
        # self.max_non_candidate_ratio = 2.  # first non-candidate has to be smaller than this
        self.fitted_positions = {}  # positions that, for any <istart>, we have fit info

        self.mfreqer = MuteFreqer(glfo)
        self.gene_obs_counts = {}  # only used for allele-finding
        self.counts = {}
        self.plotvals = {}

        self.finalized = False

    # ----------------------------------------------------------------------------------------
    def increment(self, info):
        self.mfreqer.increment(info)

        for region in utils.regions:
            regional_freq, len_excluding_ambig = utils.get_mutation_rate(self.glfo['seqs'], info, restrict_to_region=region, return_len_excluding_ambig=True)
            n_mutes = regional_freq * len_excluding_ambig  # total number of mutations in the region (for allele finding stuff)
            if abs(n_mutes - int(n_mutes)) > 1e6:
                raise Exception('n mutated %f not an integer' % n_mutes)
            n_mutes = int(n_mutes)

            gene = info[region + '_gene']
            if gene not in self.counts:
                self.counts[gene] = {}
                self.gene_obs_counts[gene] = 0
            self.gene_obs_counts[gene] += 1

            gcts = self.counts[gene]  # shorthand name

            germline_seq = info[region + '_gl_seq']
            query_seq = info[region + '_qr_seq']
            assert len(germline_seq) == len(query_seq)

            for ipos in range(len(germline_seq)):
                igl = ipos + int(info[region + '_5p_del'])  # account for left-side deletions in the indexing

                if germline_seq[ipos] in utils.ambiguous_bases or query_seq[ipos] in utils.ambiguous_bases:  # skip if either germline or query sequence is ambiguous at this position
                    continue

                if igl not in gcts:  # if we have not yet observed this position in a query sequence, initialize it
                    gcts[igl] = {}

                if igl not in gcts:
                    gcts[igl] = {}
                if utils.get_region(gene) == 'v':
                    if n_mutes not in gcts[igl]:
                        gcts[igl][n_mutes] = {n : 0 for n in ['muted', 'total'] + utils.nukes}
                    gcts[igl][n_mutes]['total'] += 1
                    if query_seq[ipos] != germline_seq[ipos]:  # if this position is mutated
                        gcts[igl][n_mutes]['muted'] += 1  # mark that we saw this germline position mutated once in a sequence with <n_mutes> regional mutation frequency
                    gcts[igl][n_mutes][query_seq[ipos]] += 1  # only used to work out what the snp'd base is if there's a new allele

    # ----------------------------------------------------------------------------------------
    def get_residual_sum(self, xvals, yvals, errs, slope, intercept):
        def expected(x):
            return slope * x + intercept
        residual_sum = sum([(y - expected(x))**2 / err**2 for x, y, err in zip(xvals, yvals, errs)])
        return residual_sum

    # ----------------------------------------------------------------------------------------
    def get_curvefit(self, n_mutelist, freqs, errs, y_icpt_bounds=None):
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
    def get_allele_finding_xyvals(self, gene, position):
        gpcounts = self.counts[gene][position]
        iterinfo = gpcounts.items()

        obs = [d['muted'] for nm, d in iterinfo if nm < self.n_max_mutations_per_segment]

        lohis = [fraction_uncertainty.err(d['muted'], d['total'], use_beta=True) for nm, d in iterinfo if nm < self.n_max_mutations_per_segment]
        errs = [(hi - lo) / 2 for lo, hi, _ in lohis]
        weights = [1./(e*e) for e in errs]

        freqs = [float(d['muted']) / d['total'] if d['total'] > 0 else 0. for nm, d in iterinfo if nm < self.n_max_mutations_per_segment]
        total = [d['total'] for nm, d in iterinfo if nm < self.n_max_mutations_per_segment]

        n_mutelist = [nm for nm in gpcounts.keys() if nm < self.n_max_mutations_per_segment]

        return {'obs' : obs, 'total' : total, 'n_mutelist' : n_mutelist, 'freqs' : freqs, 'errs' : errs, 'weights' : weights}

    # ----------------------------------------------------------------------------------------
    def is_a_candidate(self, score, min_snp_ratio, max_non_snp_ratio):
        # if score < self.min_score:  # last snp candidate has to be a lot better than the first non-snp
        #     return False
        if min_snp_ratio < self.min_min_candidate_ratio:  # worst snp candidate has to be pretty good on its own
            return False
        # if max_non_snp_ratio > self.min_min_candidate_ratio:  # first non-snp candidate has to be pretty bad on its own
        #     return False

        return True


    # ----------------------------------------------------------------------------------------
    def get_positions_to_fit(self, gene, gene_results, debug=False):
        self.fitted_positions[gene] = set()

        positions = sorted(self.mfreqer.counts[gene].keys())
        xyvals = {pos : self.get_allele_finding_xyvals(gene, pos) for pos in positions}
        positions_to_try_to_fit = [pos for pos in positions if sum(xyvals[pos]['obs']) > self.n_muted_min or sum(xyvals[pos]['total']) > self.n_total_min]  # ignore positions with neither enough mutations or total observations
        if len(positions_to_try_to_fit) < self.n_max_snps - 1 + self.min_non_candidate_positions_to_fit:
            gene_results['not_enough_obs_to_fit'].add(gene)
            if debug:
                print '          not enough positions with enough observations to fit %s' % utils.color_gene(gene)
                return None, None
        if debug and len(positions) > len(positions_to_try_to_fit):
            print '          skipping %d / %d positions (with fewer than %d mutations and %d observations)' % (len(positions) - len(positions_to_try_to_fit), len(positions), self.n_muted_min, self.n_total_min)

        self.plotvals[gene] = {}
        for pos in positions_to_try_to_fit:
            self.plotvals[gene][pos] = xyvals[pos]

        return positions_to_try_to_fit, xyvals

    # ----------------------------------------------------------------------------------------
    def fit_istart(self, gene, istart, positions_to_try_to_fit, subxyvals, fitfo, debug=False):
        residuals = {}
        for pos in positions_to_try_to_fit:
            # skip positions that are too close to the 5' end of V (misassigned insertions look like snps)
            if pos > len(self.glfo['seqs'][utils.get_region(gene)][gene]) - self.n_five_prime_positions_to_exclude - 1:
                continue

            # as long as we already have a few non-candidate positions, skip positions that have no frequencies greater than the min y intercept (note that they could in principle still have a large y intercept, but we don't really care)
            if len(residuals) > istart + self.min_non_candidate_positions_to_fit and len([f for f in subxyvals[pos]['freqs'] if f > self.min_y_intercept]) == 0:
                continue

            if sum(subxyvals[pos]['total']) < self.n_total_min:
                continue

            # also skip positions that only have a few points to fit (i.e. genes that were very rare, or I guess maybe if they were always eroded past this position)
            if len(subxyvals[pos]['n_mutelist']) < 3:
                continue

            zero_icpt_fit = self.get_curvefit(subxyvals[pos]['n_mutelist'], subxyvals[pos]['freqs'], subxyvals[pos]['errs'], y_icpt_bounds=(0. - self.small_number, 0. + self.small_number))
            big_icpt_fit = self.get_curvefit(subxyvals[pos]['n_mutelist'], subxyvals[pos]['freqs'], subxyvals[pos]['errs'], y_icpt_bounds=self.big_y_icpt_bounds)

            residuals[pos] = {'zero_icpt' : zero_icpt_fit['residuals_over_ndof'], 'big_icpt' : big_icpt_fit['residuals_over_ndof']}

            self.fitted_positions[gene].add(pos)  # if we already did the fit for another <istart>, it'll already be in there

        if len(residuals) <= istart:  # needs to be at least one longer, so we have the first-non-snp
            if debug:
                print '      not enough observations to fit more than %d snps' % (istart - 1)
            return

        residual_ratios = {pos : float('inf') if r['big_icpt'] == 0. else r['zero_icpt'] / r['big_icpt'] for pos, r in residuals.items()}
        sorted_ratios = sorted(residual_ratios.items(), key=operator.itemgetter(1), reverse=True)  # sort the positions in decreasing order of residual ratio
        candidate_snps = [pos for pos, _ in sorted_ratios[:istart]]  # the first <istart> positions are the "candidate snps"
        max_non_snp, max_non_snp_ratio = sorted_ratios[istart]  # position and ratio for largest non-candidate
        min_candidate_ratio = min([residual_ratios[cs] for cs in candidate_snps])

        fitfo['scores'][istart] = (min_candidate_ratio - max_non_snp_ratio) / max(self.small_number, max_non_snp_ratio)
        fitfo['min_snp_ratios'][istart] = min([residual_ratios[cs] for cs in candidate_snps])
        fitfo['max_non_snp_ratios'][istart] = max_non_snp_ratio
        fitfo['candidates'][istart] = {cp : residual_ratios[cp] for cp in candidate_snps}

        if debug:
            for pos in candidate_snps + [max_non_snp, ]:
                xtrastrs = ('[', ']') if pos == max_non_snp else (' ', ' ')
                pos_str = '%3s' % str(pos)
                if residual_ratios[pos] > self.min_min_candidate_ratio:
                    pos_str = utils.color('yellow', pos_str)
                print '               %s %s    %5s   (%5s / %-5s)       %3d / %3d %s' % (xtrastrs[0], pos_str, fstr(residual_ratios[pos]),
                                                                                       fstr(residuals[pos]['zero_icpt']), fstr(residuals[pos]['big_icpt']),
                                                                                       sum(subxyvals[pos]['obs']), sum(subxyvals[pos]['total']), xtrastrs[1])

            print '            %38s score: %-5s = (%-5s - %5s) / %-5s' % ('', fstr(fitfo['scores'][istart]), fstr(min_candidate_ratio), fstr(max_non_snp_ratio), fstr(max_non_snp_ratio))

    # ----------------------------------------------------------------------------------------
    def add_new_allele(self, gene, fitfo, n_candidate_snps, debug=False):
        # figure out what the new nukes are
        old_seq = self.glfo['seqs'][utils.get_region(gene)][gene]
        new_seq = old_seq
        mutfo = {}
        for pos in sorted(fitfo['candidates'][n_candidate_snps]):
            obs_counts = {nuke : self.counts[gene][pos][n_candidate_snps][nuke] for nuke in utils.nukes}  # NOTE it's super important to only use the counts from sequences with <n_candidate_snps> total mutations
            sorted_obs_counts = sorted(obs_counts.items(), key=operator.itemgetter(1), reverse=True)
            original_nuke = self.mfreqer.counts[gene][pos]['gl_nuke']
            new_nuke = None
            for nuke, _ in sorted_obs_counts:  # take the most common one that isn't the existing gl nuke
                if nuke != original_nuke:
                    new_nuke = nuke
                    break
            print '   %3d  (%s --> %s)' % (pos, original_nuke, new_nuke),
            assert old_seq[pos] == original_nuke
            mutfo[pos] = {'original' : original_nuke, 'new' : new_nuke}
            new_seq = new_seq[:pos] + new_nuke + new_seq[pos+1:]

        new_name = utils.get_new_allele_name(gene, mutfo, new_seq)
        print ''
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
    def finalize(self, debug=False):
        assert not self.finalized

        self.mfreqer.finalize()

        start = time.time()
        gene_results = {'not_enough_obs_to_fit' : set(), 'didnt_find_anything_with_fit' : set(), 'new_allele' : set()}
        if debug:
            print '\nlooking for new alleles:'
        for gene in sorted(self.mfreqer.counts):
            if utils.get_region(gene) != 'v':
                continue
            if debug:
                print '\n%s (observed %d %s)' % (utils.color_gene(gene), self.gene_obs_counts[gene], utils.plural_str('time', self.gene_obs_counts[gene]))

            positions_to_try_to_fit, xyvals = self.get_positions_to_fit(gene, gene_results, debug=debug)
            if positions_to_try_to_fit is None:
                continue

            fitfo = {n : {} for n in ('scores', 'min_snp_ratios', 'max_non_snp_ratios', 'candidates')}
            for istart in range(1, self.n_max_snps):
                if debug:
                    if istart == 1:
                        print '                                 resid. / ndof'
                        print '             position   ratio   (m=0 / m>%5.2f)       muted / obs ' % self.big_y_icpt_bounds[0]
                    print '  %d %s' % (istart, utils.plural_str('snp', istart))

                subxyvals = {pos : {k : v[istart:] for k, v in xyvals[pos].items()} for pos in positions_to_try_to_fit}
                self.fit_istart(gene, istart, positions_to_try_to_fit, subxyvals, fitfo, debug=debug)
                if istart not in fitfo['scores']:  # if it didn't get filled, we didn't have enough observations to do the fit
                    break

            istart_candidates = []
            for istart, score in sorted(fitfo['scores'].items(), key=operator.itemgetter(1), reverse=True):
                if self.is_a_candidate(score, fitfo['min_snp_ratios'][istart], fitfo['max_non_snp_ratios'][istart]):
                    istart_candidates.append(istart)

            # if debug:
            #     print '\n  fit results for each snp hypothesis:'
            #     print '    snps      score       min snp     max non-snp'
            #     for istart, score in sorted(fitfo['scores'].items(), key=operator.itemgetter(1), reverse=True):
            #         print_str = '    %2d     %9s   %9s   %9s' % (istart, fstr(score), fstr(fitfo['min_snp_ratios'][istart]), fstr(fitfo['max_non_snp_ratios'][istart]))
            #         if istart == n_candidate_snps:
            #             print_str = utils.color('red', print_str)
            #         print print_str

            if len(istart_candidates) > 0:
                n_candidate_snps = min(istart_candidates)  # add the candidate with the smallest number of snps to the germline set, and run again (if the firs
                gene_results['new_allele'].add(gene)
                print '\n    found a new allele candidate separated from %s by %d %s at %s:' % (utils.color_gene(gene), n_candidate_snps,
                                                                                                utils.plural_str('snp', n_candidate_snps), utils.plural_str('position', n_candidate_snps)),
                self.add_new_allele(gene, fitfo, n_candidate_snps, debug=debug)
            else:
                gene_results['didnt_find_anything_with_fit'].add(gene)
                if debug:
                    print '  no new alleles'

        if debug:
            print 'found new alleles for %d %s (there were also %d without new alleles, and %d without enough observations to fit)' % (len(gene_results['new_allele']), utils.plural_str('gene', len(gene_results['new_allele'])),
                                                                                                                                       len(gene_results['didnt_find_anything_with_fit']), len(gene_results['not_enough_obs_to_fit']))
            print '      allele finding time: %.1f' % (time.time()-start)

        self.finalized = True

    # ----------------------------------------------------------------------------------------
    def write(self, base_outdir, debug=False):
        if not self.finalized:
            self.finalize(debug=debug)

        if self.args.new_allele_fname is not None:
            n_new_alleles = len(self.new_allele_info)
            print '  writing %d new %s to %s' % (n_new_alleles, utils.plural_str('allele', n_new_alleles), self.args.new_allele_fname)
            with open(self.args.new_allele_fname, 'w') as outfile:
                for allele_info in self.new_allele_info:
                    outfile.write('>%s\n' % allele_info['gene'])
                    outfile.write('%s\n' % allele_info['seq'])

    # ----------------------------------------------------------------------------------------
    def plot(self, base_plotdir, only_csv=False):
        if not self.finalized:
            self.finalize(debug=debug)

        plotdir = base_plotdir + '/allele-finding'

        print '\nneed to make sure subdir preping is ok\n'
        for old_gene_dir in glob.glob(plotdir + '/*'):
            if not os.path.isdir(old_gene_dir):
                raise Exception('not a directory: %s' % old_gene_dir)
            utils.prep_dir(old_gene_dir, multilings=('*.csv', '*.svg'))
            os.rmdir(old_gene_dir)

        if only_csv:  # not implemented
            return

        for gene in self.plotvals:
            if utils.get_region(gene) != 'v':
                continue

            start = time.time()
            for position in self.plotvals[gene]:
                if position not in self.fitted_positions[gene]:  # we can make plots for the positions we didn't fit, but there's a *lot* of them and they're slow
                    continue
                # if 'allele-finding' not in self.TMPxyvals[gene][position] or self.TMPxyvals[gene][position]['allele-finding'] is None:
                #     continue
                plotting.make_allele_finding_plot(plotdir + '/' + utils.sanitize_name(gene), gene, position, self.plotvals[gene][position])
            print '      allele finding plot time: %.1f' % (time.time()-start)

