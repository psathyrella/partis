from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import random
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

from .hist import Hist
from . import utils
from . import glutils

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
    def big_discontinuity_factor(self, istart):
        if istart == 1:
            return 2.25  # since the <istart - 1>th bin is the zero bin, in which we sometimes expect a very small number of sequences, this needs to be smaller here
        elif istart < 4:
            return 2.5  # i.e. check everything that's more than <factor> sigma away (where "check" means actually do the fits, as long as it passes all the other prefiltering steps)
        else:
            return 3.5  # i.e. check everything that's more than <factor> sigma away (where "check" means actually do the fits, as long as it passes all the other prefiltering steps)

    def __init__(self, glfo, args, itry=0):
        self.region = 'v'
        self.glfo = glfo
        self.args = args
        self.itry = itry

        self.param_order = ['slope', 'y_icpt']

        self.fraction_of_seqs_to_exclude = 0.01  # exclude the fraction of sequences with largest v_{5,3}p deletions whose counts add up to this fraction of total sequences NOTE you don't want to make this too big, because although you'll be removing all the seqs with large 4p deletions, this number also gets used when you're deciding whether your new allele is in the default glfo
        self.n_bases_to_exclude = {'5p' : {}, '3p' : {}}  # i.e. on all the seqs we keep, we exclude this many bases; and any sequences that have larger deletions than this are not kept
        self.genes_to_exclude = set()  # genes that, with the above restrictions, are entirely excluded

        self.n_max_queries = 100000  # the few times i've run on samples with like a million seqs, it starts kicking a lot of extra novel alleles that are clearly spurious, so for now just subsample down to 100k, which is plenty for decent germline inference

        self.max_fit_length = 10

        self.hard_code_five = 5
        self.hard_code_three = 3
        self.hard_code_two = 2

        self.n_muted_min = 30  # don't fit positions that have fewer total mutations than this (i.e. summed over bins)
        self.n_total_min = 150  # ...or fewer total observations than this
        self.n_muted_min_per_bin = 7  # <istart>th bin has to have at least this many mutated sequences (i.e. 2-3 sigma from zero)
        self.min_fraction_per_bin = 0.005  # require that every bin (i.e. n_muted) from 0 through <self.args.max_n_snps> has at least 1% of the total (unless overridden)

        self.min_min_candidate_ratio = 2.25  # every candidate ratio must be greater than this
        self.min_mean_candidate_ratio = 2.75  # mean of candidate ratios must be greater than this
        self.min_bad_fit_residual = 1.95
        self.max_good_fit_residual = 4.5  # since this is unbounded above (unlike the min bad fit number), it needs to depend on how bad the bad fit/good fit ratio is (although, this starts making it hard to distinguish this from the ratio criterion, but see next parameter below))
        self.large_residual_ratio = 4.25  # if the ratio's bigger than this, we don't apply the max good fit residual criterion (i.e. if the ratio is a total slam dunk, it's ok if the good fit is shitty) UPDATE this is dumb, I should just find away around this bullshit
        self.default_consistency_sigmas = 3.  # default number of sigma for the boundary between consistent and inconsistent fits
        self.discontinuity_consistency_sigma = 4.5
        # self.max_consistent_candidate_fit_sigma = 10.  # this is extremely permissiive, since we don't expect  them to actually the same -- in particular, the slopes are given by the position's mutation rate (among I think maybe other things)
        self.min_discontinuity_slope_ratio = 2.5  # the fractional difference between the slope *at* the discontinuity and that on either side has to be at least this big

        self.min_min_candidate_ratio_to_plot = 1.5  # don't plot positions that're below this (for all <istart>)

        self.n_warn_alleles_per_gene = 2

        self.default_slope_bounds = (-0.1, 1.)  # fitting function needs some reasonable bounds from which to start (I could at some point make slope part of the criteria for candidacy, but it wouldn't add much sensitivity)
        self.unbounded_y_icpt_bounds = (-1., 1.5)

        self.seq_info = {}

        self.counts, self.fitfos = {}, {}
        self.inferred_allele_info = []  # new alleles with respect to the template genes from which we actually inferred them, for usage internal to allelefinder
        self.new_allele_info = []  #  new alleles with respect to original template genes, for external use (distinction is important if we infer a new allele from another previously-inferred new allele)
        self.positions_to_plot = {}
        self.n_seqs_too_highly_mutated = {}  # sequences (per-gene) that had more than <self.args.n_max_mutations_per_segment> mutations
        self.gene_obs_counts = {}  # NOTE same as n_clonal_representatives
        self.overall_mute_counts = Hist(self.args.n_max_mutations_per_segment - 1, 0.5, self.args.n_max_mutations_per_segment - 0.5)  # i.e. 0th (underflow) bin corresponds to zero mutations
        self.per_gene_mute_counts = {}  # crappy name -- this is the denominator for each position in <self.counts>. Which is usually the same for most positions, but it's cleaner to keep it separate than choose one of 'em.
        self.n_big_del_skipped = {s : {} for s in self.n_bases_to_exclude}

        self.n_fits = 0

        self.default_initial_glfo = None
        if self.args.default_initial_germline_dir is not None:  # if this is set, we want to take any new allele names from this directory's glfo if they're in there
            self.default_initial_glfo = glutils.read_glfo(self.args.default_initial_germline_dir, glfo['locus'])

        self.n_excluded_clonal_queries = {}
        self.n_clones = {}
        self.n_clonal_representatives = {}  # NOTE same as gene_obs_counts

        self.finalized = False

        self.reflengths = {}
        # self.alleles_with_evidence = set()

    # ----------------------------------------------------------------------------------------
    def dbgfcn(self, pos, istart, pos_2=None):
        dbg_positions = None  # [238, 226]
        dbg_istart = None if dbg_positions is None else len(dbg_positions) #None # 4

        if dbg_positions is None or dbg_istart is None:
            return False

        is_dbg = istart == dbg_istart
        if pos_2 is None:
            is_dbg &= pos in dbg_positions
        else:
            is_dbg &= pos in dbg_positions and pos_2 in dbg_positions
        return is_dbg

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
        self.n_excluded_clonal_queries[gene] = 0
        self.n_clones[gene] = 0
        self.n_clonal_representatives[gene] = 0

    # ----------------------------------------------------------------------------------------
    def choose_cluster_representatives(self, swfo, cluster, debug=False):  # NOTE there's somewhat similar code in AlleleClusterer
        # I could probably make this a lot faster by just always only taking one sequence from highly mutated clusters (they have so many shared mutations anyway)
        # also, now that I'm clustering only within each gene, I don't need the checks about all the sequences having the same gene
        # note that when there *is* a new allele, this is wildly over-conservative (i.e. it collapses more than we need to), since (almost) *everyone* in the cluster will [appear to] have the mutation. I don't think there's really a way around that, though
        # TODO the per-gene counters in here are deprecated now that I'm clustering per-gene in increment_and_finalize()

        assert len(cluster) > 0
        if len(cluster) == 1:
            self.n_clonal_representatives[swfo[cluster[0]][self.region + '_gene']] += 1
            self.n_clones[swfo[cluster[0]][self.region + '_gene']] += 1
            if debug:
                print('  singleton %s' % cluster[0])
            return cluster

        n_muted = {q : {'n' : self.seq_info[q]['n_mutes'], 'positions' : self.seq_info[q]['positions']} for q in cluster}
        sorted_cluster = sorted(cluster, key=lambda q: n_muted[q]['n'])
        genes = [swfo[q][self.region + '_gene'] for q in sorted_cluster]
        if genes.count(genes[0]) != len(genes):  # should only happen if the part that differentiates the genes was deleted or the read didn't extend that far, in which case I think it's fine
            print('  %s genes in cluster %s not all the same: %s' % (utils.color('yellow', 'warning'), ' '.join(sorted_cluster), ' '.join([utils.color_gene(g) for g in genes])))
        if debug:
            print('%s' % utils.color_gene(genes[0]))
            for q in sorted_cluster:
                print('    %20s  %3d  %s' % (q, n_muted[q]['n'], ' '.join([str(p) for p in sorted(n_muted[q]['positions'])])))

        # ----------------------------------------------------------------------------------------
        def no_shared_mutations(q1, q2):  # true if q1 and q2 don't have mutations at any of the same positions
            if n_muted[q1]['n'] == 0 or n_muted[q2]['n'] == 0:
                return True
            p1, p2 = n_muted[q1]['positions'], n_muted[q2]['positions']
            return len(p1 & p2) == 0

        chosen_queries = [q for q in sorted_cluster if n_muted[q]['n'] == 0]  # just throw all the unmutated ones in to start with
        sorted_cluster = [q for q in sorted_cluster if n_muted[q]['n'] > 0]  # and remove them from <sorted_cluster>
        while len(sorted_cluster) > 0:
            qchosen = sorted_cluster[0]
            chosen_queries.append(qchosen)
            sorted_cluster = [q for q in sorted_cluster if no_shared_mutations(q, qchosen)]  # remove everybody that shares any mutations with sorted_cluster[0] (including the one we just chose)

        self.n_excluded_clonal_queries[genes[0]] += len(cluster) - len(chosen_queries)
        self.n_clones[genes[0]] += 1
        self.n_clonal_representatives[genes[0]] += len(chosen_queries)

        if debug:
            n_chosen_str = str(len(chosen_queries))
            if len(chosen_queries) > 1:
                n_chosen_str = utils.color('blue', n_chosen_str)
            print('  %s  %s' % (n_chosen_str, ' '.join(chosen_queries)))

            print('      %s      naive' % swfo[cluster[0]]['naive_seq'])
            for query in sorted(cluster, key=lambda q: n_muted[q]['n']):  # can't use <sorted_cluster> since it's empty now
                print('    %s %s  %2d  %s' % (utils.color('blue', 'x') if query in chosen_queries else ' ', utils.color_mutants(swfo[query]['naive_seq'], swfo[query]['seqs'][0]), self.seq_info[query]['n_mutes'], utils.color('blue', query) if query in chosen_queries else query))
            print('')

        return chosen_queries

    # ----------------------------------------------------------------------------------------
    def set_excluded_bases(self, swfo, debug=False):
        # debug = 2
        for side in self.n_bases_to_exclude:
            # first, for each gene, count how many times we saw each deletion length
            dcounts = {}
            for query in swfo['queries']:
                gene = swfo[query][self.region + '_gene']
                if gene not in dcounts:
                    dcounts[gene] = {}
                dlen = swfo[query][self.region + '_' + side + '_del']
                if dlen not in dcounts[gene]:
                    dcounts[gene][dlen] = 0
                dcounts[gene][dlen] += 1

            # then, for each gene, find the deletion length such that only a few (self.fraction_of_seqs_to_exclude) of the sequences have longer deletions
            #  - exclude any sequences with deletions longer than this
            #  - in sequences which we keep, exclude the positions which this deletion length deletes
            for gene in dcounts:
                observed_deletions = sorted(dcounts[gene].keys())
                total_obs = sum(dcounts[gene].values())
                running_sum = 0
                if debug > 1:
                    print(gene)
                    print('  observed %s deletions: %s (counts %s)' % (side, ' '.join([str(d) for d in observed_deletions]), ' '.join([str(c) for c in dcounts[gene].values()])))
                    print('     len   fraction')
                for dlen in observed_deletions:
                    self.n_bases_to_exclude[side][gene] = dlen  # setting this before the "if" means that if we fall through (e.g. if there aren't enough sequences to get above the threshold) we'll still have a reasonable default
                    running_sum += dcounts[gene][dlen]
                    if debug > 1:
                        print('    %4d    %5.3f' % (dlen, float(running_sum) / total_obs))
                    if float(running_sum) / total_obs > 1. - self.fraction_of_seqs_to_exclude:  # if we've already added deletion lengths accounting for most of the sequences, ignore the rest of 'em
                        break
                if debug > 1:
                    print('     choose', self.n_bases_to_exclude[side][gene])

        # print choices and check consistency
        if debug > 1:
            print('    exclusions:  5p   3p')
        for gene in sorted(dcounts.keys()):
            total_exclusion_length = self.n_bases_to_exclude['5p'][gene] + self.n_bases_to_exclude['3p'][gene]
            if len(self.glfo['seqs'][self.region][gene]) - total_exclusion_length < self.args.min_allele_finding_gene_length:  # if the non-excluded part of the gene is too short, don't even bother looking for new alleles with/on it
                self.genes_to_exclude.add(gene)
            if debug > 1:
                print('                %3d  %3d  %s' % (self.n_bases_to_exclude['5p'][gene], self.n_bases_to_exclude['3p'][gene], utils.color_gene(gene, width=15)), end=' ')
                if gene in self.genes_to_exclude:
                    print('%s excluding from analysis' % utils.color('red', 'too long:'), end=' ')
                print('')

        if debug and len(self.genes_to_exclude) > 0:
            print('    excluding %d / %d genes whose reads are too short (adjust with --min-allele-finding-gene-length) %s' % (len(self.genes_to_exclude), len(dcounts), '' if len(self.genes_to_exclude) > 10 else ' '.join([utils.color_gene(g) for g in self.genes_to_exclude])))

    # ----------------------------------------------------------------------------------------
    def get_seqs_for_query(self, info, gene):
        germline_seq = info[self.region + '_gl_seq']
        assert len(info['seqs']) == 1
        query_seq = info[self.region + '_qr_seqs'][0]
        assert len(germline_seq) == len(query_seq)

        left_exclusion = self.n_bases_to_exclude['5p'][gene] - info[self.region + '_5p_del']
        right_exclusion = self.n_bases_to_exclude['3p'][gene] - info[self.region + '_3p_del']
        assert left_exclusion >= 0  # internal consistency check -- we should've already removed all the sequences with bigger deletions
        assert right_exclusion >= 0
        if left_exclusion + right_exclusion >= len(germline_seq):
            raise Exception('excluded all bases for %s: %d + %d >= %d' % (info['unique_ids'], left_exclusion, right_exclusion, len(germline_seq)))
        germline_seq = germline_seq[left_exclusion : len(germline_seq) - right_exclusion]
        query_seq = query_seq[left_exclusion : len(query_seq) - right_exclusion]
        # NOTE <germline_seq> and <query_seq> no longer correspond to <info>, but that should be ok
        if gene not in self.reflengths:  # make sure that all query sequences for this gene are precisely the same length
            self.reflengths[gene] = len(query_seq)
        assert self.reflengths[gene] == len(query_seq)  # just an internal consistency check now -- they should all be identical

        n_mutes, muted_positions = utils.hamming_distance(germline_seq, query_seq, return_mutated_positions=True)
        muted_positions = [p + self.n_bases_to_exclude['5p'][gene] for p in muted_positions]
        # double check positions are correct:
        # full_gl_seq = self.glfo['seqs'][self.region][gene]
        # print muted_positions
        # print full_gl_seq
        # print utils.color_mutants(full_gl_seq, self.n_bases_to_exclude['5p'][gene] * 'N' + query_seq + self.n_bases_to_exclude['3p'][gene] * 'N', print_isnps=True)
        # for ipos in range(len(full_gl_seq)):
        #     if ipos < self.n_bases_to_exclude['5p'][gene]:
        #         continue
        #     if ipos >= len(query_seq):
        #         continue
        #     if ipos in muted_positions:
        #         assert query_seq[ipos - self.n_bases_to_exclude['5p'][gene]] != full_gl_seq[ipos]
        #     else:
        #         assert query_seq[ipos - self.n_bases_to_exclude['5p'][gene]] == full_gl_seq[ipos]

        self.seq_info[info['unique_ids'][0]] = {'n_mutes' : n_mutes, 'positions' : set(muted_positions), 'gl_seq' : germline_seq, 'qr_seq' : query_seq}

    # ----------------------------------------------------------------------------------------
    def skip_query(self, uid, info):  # and fill self.seq_info
        gene = info[self.region + '_gene']

        if gene in self.genes_to_exclude:
            return True

        if gene not in self.counts:
            self.init_gene(gene)

        for side in self.n_bases_to_exclude:  # NOTE this is important, because if there is a snp to the right of the cysteine, all the sequences in which it is removed by the v_3p deletion will be shifted one bin leftward, screwing everything up (same goes now that we're doing the same thing on the left side)
            if info[self.region + '_' + side + '_del'] > self.n_bases_to_exclude[side][gene]:
                self.n_big_del_skipped[side][gene] += 1
                return True  # if you continue with the loop, the counters will be off

        self.get_seqs_for_query(info, gene)

        if self.seq_info[uid]['n_mutes'] > self.args.n_max_mutations_per_segment:
            self.n_seqs_too_highly_mutated[gene] += 1
            return True

        return False

    # ----------------------------------------------------------------------------------------
    def increment_query(self, uid, gene):
        self.gene_obs_counts[gene] += 1

        n_mutes, germline_seq, query_seq = [self.seq_info[uid][tn] for tn in ['n_mutes', 'gl_seq', 'qr_seq']]  # historical...
        self.overall_mute_counts.fill(n_mutes)  # NOTE this is almost the same as the hists in mutefreqer.py, except those divide by sequence length, and more importantly, they don't make sure all the sequences are the same length (i.e. the base exclusions stuff)
        self.per_gene_mute_counts[gene].fill(n_mutes)  # NOTE this is almost the same as the hists in mutefreqer.py, except those divide by sequence length, and more importantly, they don't make sure all the sequences are the same length (i.e. the base exclusions stuff)

        assert len(germline_seq) == len(self.glfo['seqs'][self.region][gene]) - self.n_bases_to_exclude['5p'][gene] - self.n_bases_to_exclude['3p'][gene]
        for ipos in range(len(germline_seq)):
            igl = ipos + self.n_bases_to_exclude['5p'][gene]  # position in original (i.e. complete) germline gene

            if germline_seq[ipos] == utils.ambig_base or query_seq[ipos] == utils.ambig_base:  # skip if either germline or query sequence is ambiguous at this position
                continue

            self.counts[gene][igl][n_mutes]['total'] += 1
            if query_seq[ipos] != germline_seq[ipos]:  # if this position is mutated
                self.counts[gene][igl][n_mutes]['muted'] += 1  # mark that we saw this germline position mutated once in a sequence with <n_mutes> regional mutation frequency
            self.counts[gene][igl][n_mutes][query_seq[ipos]] += 1  # if there's a new allele, we need this to work out what the snp'd base is

    # ----------------------------------------------------------------------------------------
    def get_residual_sum(self, xvals, yvals, errs, slope, intercept, debug=False):
        def expected(x):
            return slope * x + intercept
        residuals = [(y - expected(x))**2 / err**2 for x, y, err in zip(xvals, yvals, errs)]
        if debug:
            print('resid: ' + ' '.join(['%5.3f' % r for r in residuals]))
        return sum(residuals)

    # ----------------------------------------------------------------------------------------
    def dbgstr(self, fitfo, extra_str='', pvals=None):
        return_strs = []
        if pvals is not None:
            return_strs.append('  ' + ' '.join(['%5d' % n for n in pvals['n_mutelist']]))
            return_strs.append('  ' + ' '.join(['%5.3f' % f for f in pvals['freqs']]))
            return_strs.append('  ' + ' '.join(['%5.3f' % e for e in pvals['errs']]))
        return_strs.append('        %s  m: %5.3f +/- %5.4f   b: %5.3f +/- %5.4f     %5.3f / %-3d = %5.3f' % (extra_str, fitfo['slope'], fitfo['slope_err'], fitfo['y_icpt'], fitfo['y_icpt_err'],
                                                                                                             fitfo['residual_sum'], fitfo['ndof'], fitfo['residuals_over_ndof']))
        return '\n'.join(return_strs)

    # ----------------------------------------------------------------------------------------
    def default_fitfo(self, xvals, yvals, errs, ndof=0, y_icpt_bounds=None):
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
            n_mutelist.append(pvals['n_mutelist'][im])
            freqs.append(pvals['freqs'][im])
            errs.append(pvals['errs'][im])
        return n_mutelist, freqs, errs

    # ----------------------------------------------------------------------------------------
    def cov_err_ok(self, cov_err, errs):
        non_null_errs = [e for e in errs if e > utils.eps]
        if len(non_null_errs) == 0:
            print('%s all errs very small for:\n%s\n%s\n%s' % (utils.color('red', 'error'), n_mutelist, freqs, errs))
        min_point_err = min(non_null_errs)
        if cov_err / min_point_err < 1e-5:  # if the covariance from scipy is much smaller than the uncertainty on all the points, it was probably a "perfect" fit (i.e. cov_err is basically or exactly zero 'cause the line goes exactly through every point)
            return False
        return True

    # ----------------------------------------------------------------------------------------
    def get_scipy_results(self, fixed_y_icpt, linefunc, n_mutelist, freqs, errs, bounds, fitfo):
        params, cov = scipy.optimize.curve_fit(linefunc, n_mutelist, freqs, sigma=errs, bounds=bounds)
        self.n_fits += 1

        def retrieve(param):
            iparam = self.param_order.index(param)
            cov_err = math.sqrt(cov[iparam][iparam])
            if self.cov_err_ok(cov_err, errs):
                err = cov_err
            else:
                err = self.hack_err(param, errs, n_mutelist)
            return params[iparam], err

        fitfo['slope'], fitfo['slope_err'] = retrieve('slope')
        if not fixed_y_icpt:
            fitfo['y_icpt'], fitfo['y_icpt_err'] = retrieve('y_icpt')

    # ----------------------------------------------------------------------------------------
    def get_curvefit(self, pvals, y_icpt_bounds, debug=False):
        n_mutelist, freqs, errs = self.get_tmp_fitvals(pvals)  # this is probably kind of slow
        if y_icpt_bounds[0] == y_icpt_bounds[1]:  # fixed y-icpt
            fitfo = self.default_fitfo(n_mutelist, freqs, errs, ndof=len(n_mutelist) - 1, y_icpt_bounds=y_icpt_bounds)
            fitfo['y_icpt'] = y_icpt_bounds[0]
            if fitfo['ndof'] > 0:
                def linefunc(x, slope):
                    return slope*x + fitfo['y_icpt']
                bounds = self.default_slope_bounds
                self.get_scipy_results(True, linefunc, n_mutelist, freqs, errs, bounds, fitfo)  # modifies <fitfo>
            elif fitfo['ndof'] == 0:
                fitfo = self.approx_fit_vals(pvals, fixed_y_icpt=fitfo['y_icpt'], debug=debug)
        else:  # floating y-icpt
            fitfo = self.default_fitfo(n_mutelist, freqs, errs, ndof=len(n_mutelist) - 2, y_icpt_bounds=y_icpt_bounds)
            if fitfo['ndof'] > 0:
                def linefunc(x, slope, y_icpt):
                    return slope*x + y_icpt
                bounds = [[s, y] for s, y in zip(self.default_slope_bounds, y_icpt_bounds)]
                self.get_scipy_results(False, linefunc, n_mutelist, freqs, errs, bounds, fitfo)  # modifies <fitfo>
            elif fitfo['ndof'] == 0:
                fitfo = self.approx_fit_vals(pvals, fixed_y_icpt=None, debug=debug)

        if fitfo['ndof'] > 0:
            fitfo['residual_sum'] = self.get_residual_sum(n_mutelist, freqs, errs, fitfo['slope'], fitfo['y_icpt'])
            fitfo['residuals_over_ndof'] = float(fitfo['residual_sum']) / fitfo['ndof']
        else:
            fitfo['residual_sum'] = 1.
            fitfo['residuals_over_ndof'] = 0.

        if debug:
            print(self.dbgstr(fitfo, extra_str='fit', pvals={'n_mutelist' : n_mutelist, 'freqs' : freqs, 'errs': errs}))  # not necessarily the same as <pvals>

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
            print(' ', ' '.join(['%5d' % n for n in gcts]))
            for nuke in utils.nukes:
                print(nuke, ' '.join(['%5d' % gcts[n][nuke] for n in gcts]), overall_nuke_totals[nuke])
            print(' ', ' '.join(['%5.3f' % f for f in freqs]))

        reweighted_freqs = []
        for n_muted in gcts:
            freq, total = 0., 0.
            for nuke in utils.nukes:
                reweight = 0.
                if per_bin_nuke_totals[n_muted][nuke] > 0:
                    reweight = float(overall_nuke_totals[nuke]) / per_bin_nuke_totals[n_muted][nuke]
                total += reweight * gcts[n_muted][nuke]
                if nuke != self.glfo['seqs'][self.region][gene][position]:
                    freq += reweight * gcts[n_muted][nuke]
            reweighted_freqs.append(freq / total if total > 0. else 0.)
        if debug:
            print(' ', ' '.join(['%5.3f' % f for f in reweighted_freqs]))
        return reweighted_freqs

    # ----------------------------------------------------------------------------------------
    def get_allele_finding_xyvals(self, gene, position):
        from . import fraction_uncertainty
        gcts = self.counts[gene][position]  # shorthand name

        obs = [d['muted'] for d in gcts.values()]

        lohis = [fraction_uncertainty.err(d['muted'], d['total']) if d['total'] > 0. else (0., 1.) for d in gcts.values()]  # set uncertainty bounds to (0., 1.) for zero-denominator bins
        errs = [(hi - lo) / 2 for lo, hi in lohis]
        weights = [1./(e*e) for e in errs]

        freqs = [float(d['muted']) / d['total'] if d['total'] > 0 else 0. for d in gcts.values()]
        total = [d['total'] for d in gcts.values()]

        n_mutelist = list(gcts.keys())

        return {'obs' : obs, 'total' : total, 'n_mutelist' : n_mutelist, 'freqs' : freqs, 'errs' : errs, 'weights' : weights}

    # ----------------------------------------------------------------------------------------
    def get_both_pre_post_vals(self, gene, istart, positions_to_try_to_fit):
        # NOTE I'm including the zero bin here -- do I really want to do that? UPDATE yes, I think so -- we know there will be zero mutations in that bin, but the number of sequences in it still contains information (uh, I think)
        min_ibin = max(0, istart - self.max_fit_length)
        max_ibin = min(self.args.n_max_mutations_per_segment, istart + self.max_fit_length)
        bothxyvals = {pos : {k : v[min_ibin : max_ibin] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}
        prexyvals = {pos : {k : v[min_ibin : istart] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}  # arrays up to, but not including, <istart>
        postxyvals = {pos : {k : v[istart : max_ibin] for k, v in self.xyvals[gene][pos].items()} for pos in positions_to_try_to_fit}  # arrays from <istart> onwards
        return bothxyvals, prexyvals, postxyvals

    # ----------------------------------------------------------------------------------------
    def is_a_candidate(self, candidfo, debug=False):
        if candidfo['min_snp_ratio'] < self.min_min_candidate_ratio:  # worst snp candidate has to be pretty good on its own
            if debug:
                print('    min snp ratio %s too small (less than %s)' % (fstr(candidfo['min_snp_ratio']), fstr(self.min_min_candidate_ratio)), end=' ')
            return False
        if candidfo['mean_snp_ratio'] < self.min_mean_candidate_ratio:  # mean of snp candidate ratios has to be even better
            if debug:
                print('    mean snp ratio %s too small (less than %s)' % (fstr(candidfo['mean_snp_ratio']), fstr(self.min_mean_candidate_ratio)), end=' ')
            return False

        # make sure all the snp positions have similar fits (the bin totals for all snp positions should be highly correlated, since they should ~all be present in ~all sequences [that stem from the new allele])
        # NOTE that with multiple multi-snp new alleles that share some, but not all, positions, we don't expect consistency. In particular, at shared positions, the nsnp bin for the other allele will be high, and the prevalence will be off.
        for pos_1, pos_2 in itertools.combinations(candidfo['positions'], 2):
            fitfo_1, fitfo_2 = candidfo['fitfos'][pos_1], candidfo['fitfos'][pos_2]
            if not self.consistent_discontinuities(fitfo_1['onefo'], fitfo_2['onefo'], candidfo['istart'], debug=self.dbgfcn(pos_1, candidfo['istart'], pos_2=pos_2)):
                if debug:
                    print('    positions %d and %d have inconsistent discontinuities' % (pos_1, pos_2))
                return False

            # if not self.consistent_fits(fitfo_1['postfo'], fitfo_2['postfo'], factor=self.max_consistent_candidate_fit_sigma, debug=self.dbgfcn(pos_1, candidfo['istart'], pos_2=pos_2)):
            #     if debug:
            #         print '    positions %d and %d have inconsistent post-istart fits' % (pos_1, pos_2)
            #     return False

            # if this nsnp is less than 3, and there's a second new allele with smaller nsnp, the pre-fit will be super inconsistent, but this is really rare in actual data... so I'm commenting it since I care more about avoiding false positives

            # if not self.consistent_fits(fitfo_1['prefo'], fitfo_2['prefo'], factor=self.max_consistent_candidate_fit_sigma, debug=self.dbgfcn(pos_1, candidfo['istart'], pos_2=pos_2)):
            #     if debug:
            #         print '    positions %d and %d have inconsistent pre-istart fits' % (pos_1, pos_2)
            #         self.consistent_fits(fitfo_1['prefo'], fitfo_2['prefo'], factor=self.max_consistent_candidate_fit_sigma, debug=True)
            #     return False

        if debug:
            print('    candidate', end=' ')
        return True

    # ----------------------------------------------------------------------------------------
    def hack_err(self, param, yerrs, xvals):  # median of the errors divided by the square root of the sample size (should be reasonable...)
        y_err = float(numpy.median(yerrs)) / math.sqrt(len(yerrs))
        if param == 'slope':
            x_length = max(1, (xvals[-1] - xvals[0]) / 2.)  # full x length seems to kinda underestimate the error
            return y_err / x_length  # divide by the x-axis length, to convert to units of slope
        elif param == 'y_icpt':
            return y_err
        else:
            assert False

    # ----------------------------------------------------------------------------------------
    def approx_fit_vals(self, pvals, fixed_y_icpt=None, debug=False):
        # NOTE uncertainties are kinda complicated if you do the weighted mean, so screw it, it works fine with the plain mean

        def getslope(i1, i2, shift=False):  # return two-point slope between indices i1 and i2 (if <shift>, we replace <yv> with a vector in which each value is shifted alternately up or down by its uncertainty)
            if not shift:
                tmp_y = yv
            else:
                tmp_y = [yv[i] + (-1) ** (i%2) * ev[i] for i in range(len(yv))]  # alternately shifted up/down
            return (tmp_y[i2] - tmp_y[i1]) / (xv[i2] - xv[i1])

        xv, yv, ev = pvals['n_mutelist'], pvals['freqs'], pvals['errs']  # tmp shorthand
        fitfo = self.default_fitfo(xv, yv, ev)
        if fixed_y_icpt is not None:  # hurg, can't this be handled by passing as arguments to the fcn? don't want to figure it out just now though
            fitfo['y_icpt'] = fixed_y_icpt

        # if we were only given one point, return the defualt fitfo, possibly setting the slope based on treating the fixed y-icpt as an additional point
        assert len(xv) > 0
        if len(xv) == 1:
            if fixed_y_icpt is not None and xv[0] > 0.:  # <xv[0]> can't be able to be negative, but if it's zero, then the fit values aren't well-defined
                fitfo['slope'] = (yv[0] - fixed_y_icpt) / (xv[0] - 0.)
            return fitfo

        # first set slope and slope error
        slopes = [getslope(i-1, i) for i in range(1, len(xv))]  # only uses adjacent points, and double-counts interior points, but we don't care (we don't use steps of two, because then we'd lose the last one if it's odd-length)
        if len(xv) == 2:  # add a slope for an additional, hypothetical, pair of points, obtained by shifting one point each direction by its uncertainty
            slopes.append(getslope(0, 1, shift=True))
        fitfo['slope'] = numpy.average(slopes)
        var_err = numpy.std(slopes, ddof=1) / math.sqrt(len(xv))  # error based on variance
        if self.cov_err_ok(var_err, ev):
            fitfo['slope_err'] = var_err
        else:
            fitfo['slope_err'] = self.hack_err('slope', ev, xv)

        # then (if necessary) set y-icpt and its error
        if fixed_y_icpt is None:
            y_icpts = [yv[i] - getslope(i-1, i) * xv[i] for i in range(1, len(xv))]
            fitfo['y_icpt'] = numpy.average(y_icpts)
            if len(xv) == 2:
                y_icpts.append(yv[1] - getslope(0, 1, shift=True) * xv[1])
            var_err = numpy.std(y_icpts, ddof=1) / math.sqrt(len(xv))  # error based on variance
            if self.cov_err_ok(var_err, ev):
                fitfo['y_icpt_err'] = var_err
            else:
                fitfo['y_icpt_err'] = self.hack_err('y_icpt', ev, xv)

        if debug:
            print(self.dbgstr(fitfo, extra_str='apr', pvals=pvals))

        return fitfo

    # ----------------------------------------------------------------------------------------
    def consistent(self, v1, v1err, v2, v2err, factor=None, dbgstr='', debug=False):
        # i.e. if both slope and intercept are within <factor> std deviations of each other, don't bother fitting, because the fit isn't going to say they're wildly inconsistent
        if factor is None:
            factor = self.default_consistency_sigmas
        lo, hi = sorted([float(v1), float(v2)])
        joint_err = max(float(v1err), float(v2err))
        if joint_err == float('inf'):
            if debug:
                print('      %6s  (inf err)' % dbgstr)
            return True
        else:
            if debug:
                # print '      %6s  %6.3f +/- %7.4f   %6.3f +/- %7.4f   -->   %6.3f + %3.1f * %7.4f = %7.4f >? %6.3f   %s' % (dbgstr, v1, v1err, v2, v2err, lo, factor, joint_err, lo + factor * joint_err, hi, 'consistent' if (lo + factor * joint_err > hi) else 'nope')
                print('      %6s  %6.3f +/- %7.4f   %6.3f +/- %7.4f   -->  ' % (dbgstr, v1, v1err, v2, v2err), end=' ')
                print('(%-6.3f - %6.3f) / %7.4f = %4.2f <? %3.1f    %s' % (hi, lo, joint_err, (hi - lo) / joint_err, factor, 'consistent' if (lo + factor * joint_err > hi) else 'nope'))
            return lo + factor * joint_err > hi

    # ----------------------------------------------------------------------------------------
    def consistent_fits(self, vals1, vals2, factor=None, debug=False):
        # if I stick with consistent_discontinuities(), I think I can remove a bunch of the fit-consistency infrastructure (although I'm not sure that I really want to)
        if len(vals1['xvals']) < 5 and vals1['xvals'][0] == 0:  # the fit uncertainties are way low in cases where the points have large uncertainties, but line up really well. This only really happens when there's only a few points, though
            return self.consistent_bin_vals(vals1, vals2, factor=factor, debug=debug)
        else:
            consistent = True
            for valname in ['slope', 'y_icpt']:
                consistent &= self.consistent(vals1[valname], vals1[valname + '_err'], vals2[valname], vals2[valname + '_err'], factor=factor, dbgstr=valname, debug=debug)
            return consistent

    # ----------------------------------------------------------------------------------------
    def consistent_discontinuities(self, vals1, vals2, istart, debug=False):  # NOTE this is getting passed the one-piece fitfos, which doesn't necessarily make that much sense, but we just want the x- and y-vals, so it's ok
        # this kind of duplicates consistent_bin_vals(), but oh, well, since I may way to remove one or the other eventually
        def getdiff(vals):
            diff = vals['yvals'][istart] - vals['yvals'][istart - 1]
            err = max(vals['errs'][istart - 1], vals['errs'][istart])
            return diff, err
        diff1, err1 = getdiff(vals1)
        diff2, err2 = getdiff(vals2)
        return self.consistent(diff1, err1, diff2, err2, factor=self.discontinuity_consistency_sigma, dbgstr='discontinuities', debug=debug)

    # ----------------------------------------------------------------------------------------
    def consistent_bin_vals(self, vals1, vals2, factor=None, debug=False):
        # hey, wait, should this be using self.consistent()?
        if factor is None:
            factor = self.default_consistency_sigmas
        net_sigma = 0.
        assert len(vals1['xvals']) == len(vals2['xvals'])
        for ipos in range(len(vals1['xvals'])):
            ydiff = vals2['yvals'][ipos] - vals1['yvals'][ipos]
            joint_err = max(vals1['errs'][ipos], vals2['errs'][ipos])  # at some point i should do something slightly more sensible for my joint errors (maybe geometric mean?, quadrature [but they're not independent]?)
            net_sigma += ydiff / joint_err
            if debug:
                print('    (%6.3f - %6.3f) / %7.4f = %5.2f' % (vals2['yvals'][ipos], vals1['yvals'][ipos], joint_err, ydiff / joint_err))

        if debug:
            print('    net sigma from %d bins: %4.2f ?> %4.2f  %s' % (len(vals1['xvals']), net_sigma, factor, 'consistent' if  factor > net_sigma else 'nope'))
        return factor > net_sigma

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
    def get_big_y_icpt_bounds(self, big_y_icpt, big_y_icpt_err):
        # we want the bounds to be lenient enough to accomodate non-zero slopes (in the future, we could do something cleverer like extrapolating with the slope of the line to x=0)
        return (big_y_icpt - 1.5*big_y_icpt_err, big_y_icpt + 1.5*big_y_icpt_err)

    # ----------------------------------------------------------------------------------------
    def very_different_bin_totals(self, pvals, istart, debug=False):
        # i.e. if there's a homozygous new allele at <istart> + 1  UPDATE wait, why +1?
        last_total = pvals['total'][istart - 1]
        istart_total = pvals['total'][istart]
        joint_total_err = max(math.sqrt(last_total), math.sqrt(istart_total))
        if debug:
            print('    different bin totals:  diff / err = (%.0f - %.0f) / %5.1f = %.1f ?> %.1f'  % (istart_total, last_total, joint_total_err, (istart_total - last_total) / joint_total_err, self.big_discontinuity_factor(istart)))
        return istart_total - last_total > self.big_discontinuity_factor(istart) * joint_total_err  # it the total (denominator) is very different between the two bins

    # ----------------------------------------------------------------------------------------
    def big_discontinuity(self, pvals, istart, debug=False):  # NOTE same as very_different_bin_totals(), except for freqs rather than totals
        # (note that the size of the discontinuity tells us about the allele prevalence -- [ie].[eg]. we can be quite confident of a new allele with a small absolute discontinuity)
        # also note that for large sample sizes, the uncertainties will be small enough that this fcn will in general call it a discontinuity even for non-snpd positions

        if pvals['total'][istart] < 4:  # if there's nothing in this bin, there's certainly not a new allele (although, note, this should have already been checked for)
            return False

        if pvals['total'][istart - 1] < 10:  # if there's hardly any entries in the previous bin (i.e. presumably a homozygous new allele) then just use bin totals UPDATE you know, i'm not really sure what but that it wouldn't make more sense to always use some combination of the two
            return self.very_different_bin_totals(pvals, istart, debug=debug)

        joint_freq_err = max(pvals['errs'][istart - 1], pvals['errs'][istart])
        last_freq = pvals['freqs'][istart - 1]
        istart_freq = pvals['freqs'][istart]
        if debug:
            print('    discontinuity:  diff / err = (%5.3f - %5.3f) / %5.3f = %.1f ?> %.1f'  % (istart_freq, last_freq, joint_freq_err, (istart_freq - last_freq) / joint_freq_err, self.big_discontinuity_factor(istart)))
        return istart_freq - last_freq > self.big_discontinuity_factor(istart) * joint_freq_err

    # ----------------------------------------------------------------------------------------
    def fit_position(self, gene, istart, pos, prevals, postvals, bothvals, candidate_ratios, residfo, debug=False):
        dbg = self.dbgfcn(pos, istart)
        if dbg:
            print('pos %d' % pos)
        big_y_icpt, big_y_icpt_err = self.get_big_y(postvals)
        big_y_icpt_bounds = self.get_big_y_icpt_bounds(big_y_icpt, big_y_icpt_err)  # (big_y_icpt - 1.5*big_y_icpt_err, big_y_icpt + 1.5*big_y_icpt_err)  # we want the bounds to be lenient enough to accomodate non-zero slopes (in the future, we could do something cleverer like extrapolating with the slope of the line to x=0)

        def returnfcn(label):
            if dbg:
                print(label)
                return False  # keep going (don't skip it) if it's a dbg pos/istart
            else:
                return True

        # need to have enough mutated counts in the <istart>th bin (this is particularly important (partly) because it's the handle that tells us it's *this* <istart> that's correct, rather than <istart> + 1)
        if self.counts[gene][pos][istart]['muted'] < self.n_muted_min_per_bin:
            if returnfcn('only %d muted in <istart>th bin' % self.counts[gene][pos][istart]['muted']):
                return

        if sum(postvals['obs']) < self.n_muted_min or sum(postvals['total']) < self.n_total_min:
            if returnfcn('too few overall post-counts'):
                return

        # skip if the discontinuity is less than <factor> sigma, or if hardly any entries at i-1 and the bin totals are closer than <factor> sigma (not actualy OR, but basically)
        if not self.big_discontinuity(bothvals, istart, debug=dbg):
            if returnfcn('no big dicontinuity'):
                return

        if istart <= self.hard_code_three:
            # if the bounds include zero, there won't be much difference between the two fits
            if big_y_icpt_bounds[0] <= 0.:
                if returnfcn('big-y-icpt lower bound %f <= 0.' % big_y_icpt_bounds[0]):
                    return

            # if a rough estimate of the y-icpt is less than zero, the zero-icpt fit is probably going to be pretty good
            approx_fitfo = self.approx_fit_vals(postvals)
            if approx_fitfo['y_icpt'] < 0.:
                if returnfcn('approx post fit y-icpt %f < 0' % approx_fitfo['y_icpt']):
                    return

        # if there's only two points in <prevals>, we can't use the bad fit there to tell us this isn't a candidate, so we check and skip if the <istart - 1>th freq isn't really low
        if istart == 2 and bothvals['freqs'][istart - 1] > big_y_icpt - 1.5 * big_y_icpt_err:  # TODO wait isn't this the same as lower bound?
            if returnfcn('complicated istart = 2 special case'):
                return

        # approximate pre-slope should be smaller than approximate post-slope (for smaller <istart>s, post-slope tends to be flat, so you can't require this)
        if istart >= self.hard_code_five:
            pre_approx = self.approx_fit_vals(prevals)
            post_approx = self.approx_fit_vals(postvals)
            if not self.consistent(pre_approx['slope'], pre_approx['slope_err'], post_approx['slope'], post_approx['slope_err'], dbgstr='slope', debug=dbg) and pre_approx['slope'] > post_approx['slope']:
                if returnfcn('pre approx slope bigger than post approx slope'):
                    return

        onefit = self.get_curvefit(bothvals, y_icpt_bounds=(0., 0.), debug=dbg)

        # don't bother with the two-piece fit if the one-piece fit is pretty good
        if onefit['residuals_over_ndof'] < self.min_bad_fit_residual:
            if returnfcn('one-piece fit is pretty good %f' % onefit['residuals_over_ndof']):
                return

        prefit = self.get_curvefit(prevals, y_icpt_bounds=(0., 0.), debug=dbg)
        postfit = self.get_curvefit(postvals, y_icpt_bounds=big_y_icpt_bounds, debug=dbg)
        twofit_residuals = prefit['residuals_over_ndof'] * prefit['ndof'] + postfit['residuals_over_ndof'] * postfit['ndof']
        twofit_ndof = prefit['ndof'] + postfit['ndof']
        twofit_residuals_over_ndof = twofit_residuals / twofit_ndof

        # pre-slope should be smaller than post-slope also for the full fits (for smaller <istart>s, post-slope tends to be flat, so you can't require this)
        if istart >= self.hard_code_five:
            if not self.consistent(prefit['slope'], prefit['slope_err'], postfit['slope'], postfit['slope_err'], dbgstr='slope', debug=dbg) and prefit['slope'] > postfit['slope']:
                if returnfcn('pre slope %f bigger than post slope %f' % (prefit['slope'], postfit['slope'])):
                    return

        ratio = onefit['residuals_over_ndof'] / twofit_residuals_over_ndof if twofit_residuals_over_ndof > 0. else float('inf')
        if dbg:
            print('  %5.3f / %5.3f = %5.3f' % (onefit['residuals_over_ndof'], twofit_residuals_over_ndof, ratio))

        # make sure two-piece fit is at least ok (unless the residual ratio is incredibly convincing)
        if ratio < self.large_residual_ratio and twofit_residuals_over_ndof > self.max_good_fit_residual:
            if returnfcn('two-piece fit not good enough %f' % twofit_residuals_over_ndof):
                return

        # the slope at the discontinuity should be much larger than on either side
        discontinuity_slope = (bothvals['freqs'][istart] - bothvals['freqs'][istart - 1]) / 1.  # oh,  pedantry
        for side, sideslope, sideslopeerr in [['pre', prefit['slope'], prefit['slope_err']], ['post', postfit['slope'], postfit['slope_err']]]:
            if sideslope == 0. or sideslopeerr == float('inf'):  # the slope error should be set to inf if there's only one point
                continue
            if discontinuity_slope < sideslope:  # shouldn't really happen, but if I do this I don't have to worry about signs in the bit below
                if returnfcn('disc. slope less that %-4s slope %5.3f < %5.3f' % (side, discontinuity_slope, sideslope)):
                    return
            frac_diff = abs((discontinuity_slope - sideslope) / sideslope)
            if frac_diff < self.min_discontinuity_slope_ratio:
                if returnfcn('disc. slope not enough bigger than %-4s slope: abs((%5.3f - %5.3f) / %5.3f) = %4.2f < %3.1f' % (side, discontinuity_slope, sideslope, sideslope, frac_diff, self.min_discontinuity_slope_ratio)):
                    return

        # add it as a candidate
        candidate_ratios[pos] = ratio
        residfo[pos] = {'onefo' : onefit, 'prefo' : prefit, 'postfo' : postfit, 'twofo' : {'residuals_over_ndof' : twofit_residuals_over_ndof}}
        if dbg or candidate_ratios[pos] > self.min_min_candidate_ratio_to_plot:
            self.positions_to_plot[gene].add(pos)  # if we already decided to plot it for another <istart>, it'll already be in there

    # ----------------------------------------------------------------------------------------
    def print_candidates(self, istart, candidates, ratios, residfo, bothxyvals):  # NOTE <ratios> in general contains positions not in <candidates> (i.e. <candidates> is the ones we're supposed to be printing right now)
        if len(candidates) == 0:
            return

        if not self.already_printed_dbg_header:
            print('             position   ratio       (one piece / two pieces)  ', end=' ')
            print('%0s %s' % ('', ''.join(['%11d' % nm for nm in range(self.args.n_max_mutations_per_segment + 1)])))  # NOTE *has* to correspond to line at bottom of fcn below
            self.already_printed_dbg_header = True
        print('    %d %s' % (istart, utils.plural_str('snp', istart)))

        for pos in candidates:
            pos_str = '%3s' % str(pos)
            if ratios[pos] > self.min_min_candidate_ratio and residfo[pos]['onefo']['residuals_over_ndof'] > self.min_bad_fit_residual:
                pos_str = utils.color('yellow', pos_str)
            print_str = ['                 %s    %5s            %5s / %-5s               ' % (pos_str, fstr(ratios[pos]),
                                                                                              fstr(residfo[pos]['onefo']['residuals_over_ndof']), fstr(residfo[pos]['twofo']['residuals_over_ndof']))]
            for n_mutes in range(self.args.n_max_mutations_per_segment + 1):
                if n_mutes in bothxyvals[pos]['n_mutelist']:
                    inm = bothxyvals[pos]['n_mutelist'].index(n_mutes)
                    print_str.append('%4d / %-4d' % (bothxyvals[pos]['obs'][inm], bothxyvals[pos]['total'][inm]))
                else:
                    print_str.append('           ')
            print(''.join(print_str))

    # ----------------------------------------------------------------------------------------
    def fit_istart(self, gene, istart, positions_to_try_to_fit, debug=False):
        bothxyvals, prexyvals, postxyvals = self.get_both_pre_post_vals(gene, istart, positions_to_try_to_fit)
        ratios, residfo = {}, {}
        for pos in positions_to_try_to_fit:
            self.fit_position(gene, istart, pos, prexyvals[pos], postxyvals[pos], bothxyvals[pos], ratios, residfo, debug=debug)
        sorted_positions = sorted(ratios, key=lambda p: ratios[p], reverse=True)  # sort the candidate positions in decreasing order of residual ratio
        sorted_positions = sorted_positions[ : len(sorted_positions) - len(sorted_positions) % istart]  # remove any extra positions
        if len(sorted_positions) >= 2 * istart:  # if there's more than one candidate allele, sorted such that similar positions are together, and maybe we'll get the combinations right
            if istart < self.hard_code_three:  # for small <istart>, y-icpt is a good proxy for allele prevalence, and two new alleles are unlikely to have exactly the same prevalence (also, the within correctly-sorted position groups shit is all highly correlated)
                sorted_positions = sorted(sorted_positions, key=lambda p: residfo[p]['postfo']['y_icpt'], reverse=True)
            else:  # whereas for bigger <istart>, two-fit fit quality works well
                sorted_positions = sorted(sorted_positions, key=lambda p: residfo[p]['twofo']['residuals_over_ndof'], reverse=True)  # 
        position_lists = [sorted_positions[i : i + istart] for i in range(0, len(sorted_positions), istart)]  # and divide them into groups of length <istart> (note that we don't really have a good way of knowing which positions should go together if there's more than one group of candidates (and if <istart> is greater than 1), but since they're sorted by ratio, similar ones are together, which does an ok job)
        for plist in position_lists:
            assert len(plist) == istart  # shouldn't happen any more
            if debug:
                self.print_candidates(istart, plist, ratios, residfo, bothxyvals)
            self.fitfos[gene].append({
                'istart' : istart,
                'positions' : plist,
                'min_snp_ratio' : min([ratios[p] for p in plist]),
                'mean_snp_ratio' : numpy.mean([ratios[p] for p in plist]),
                'fitfos' : {p : residfo[p] for p in plist},
            })

    # ----------------------------------------------------------------------------------------
    def add_allele_to_new_allele_info(self, template_gene, candidfo, debug=False):
        n_candidate_snps = candidfo['istart']

        # figure out what the new nukes are
        old_seq = self.glfo['seqs'][self.region][template_gene]
        new_seq = old_seq
        mutfo = {}
        for pos in sorted(candidfo['positions']):
            obs_counts = {nuke : self.counts[template_gene][pos][n_candidate_snps][nuke] for nuke in utils.nukes}  # NOTE it's super important to only use the counts from sequences with <n_candidate_snps> total mutations
            sorted_obs_counts = sorted(list(obs_counts.items()), key=operator.itemgetter(1), reverse=True)
            original_nuke = self.glfo['seqs'][self.region][template_gene][pos]
            new_nuke = None
            for nuke, _ in sorted_obs_counts:  # take the most common one that isn't the existing gl nuke
                if nuke != original_nuke:
                    new_nuke = nuke
                    break
            assert old_seq[pos] == original_nuke
            mutfo[pos] = {'original' : original_nuke, 'new' : new_nuke}
            new_seq = new_seq[:pos] + new_nuke + new_seq[pos+1:]

        final_name, final_mutfo = glutils.choose_new_allele_name(template_gene, new_seq, snpfo=mutfo)  # final as in destined for <self.new_allele_info>, not for <self.inferred_allele_info>

        # reminder: number of mutations in <final_mutfo> is not necessarily equal to <n_candidate_snps>
        assert len(mutfo) == n_candidate_snps  # ...but this should be the same. Can remove this once I finish fixing the bug

        if self.default_initial_glfo is not None:  # if this is set, we want to take the names from this directory's glfo (i.e. see if there's already a name for <final_name>'s sequence)
            if final_name in self.default_initial_glfo['seqs'][self.region]:
                assert False  # uh... I don't think this can happen and i'm not sure what to do if it does
            else:
                equiv_name, equiv_seq = glutils.find_equivalent_gene_in_glfo(self.default_initial_glfo, new_seq, utils.cdn_pos(self.glfo, self.region, template_gene), new_name=final_name,
                                                                             exclusion_5p=self.n_bases_to_exclude['5p'][template_gene], exclusion_3p=self.n_bases_to_exclude['3p'][template_gene], debug=debug)
                if equiv_name is not None:
                    final_name = equiv_name
                    new_seq = equiv_seq

        if final_name in self.glfo['seqs'][self.region]:
            if debug:
                print('    new gene %s already in glfo (probably 3p end length issues), so skipping it' % utils.color_gene(final_name))
            return

        # we actually expect the slope to be somewhat negative (since as the mutation rate increases a higher fraction of them revert to germline)
        # this is heuristically parameterized by the non-zero values
        # NOTE if there's two new alleles separated from a known allele that isn't in the sample, this doesn't work (this is, of course, quite rare)
        remove_template = True
        homozygous_line = {'slope' : -0.01, 'slope_err' : 0.015, 'y_icpt' : 1.1, 'y_icpt_err' : 0.12}
        for pos in candidfo['fitfos']:  # if every position is consistent with slope = 0, y_icpt = 1, remove the template gene
            if not self.consistent_fits(candidfo['fitfos'][pos]['postfo'], homozygous_line, factor=1.):
                remove_template = False

        if debug:
            print('  %s %s separated from %s by %d snp%s at:  ' % (utils.color('red', 'new'), utils.color_gene(final_name), utils.color_gene(template_gene), n_candidate_snps, utils.plural(n_candidate_snps)), end=' ')
            print('  '.join([('%d (%s --> %s)' % (pos, mutfo[pos]['original'], mutfo[pos]['new'])) for pos in sorted(mutfo)]))
            if mutfo != final_mutfo:
                print('      note: final snp positions (%s) differ from inferred snp positions (%s)' % (' '.join([str(p) for p in sorted(final_mutfo)]), ' '.join([str(p) for p in sorted(mutfo)])))
            # old_len_str, new_len_str = '', ''
            # old_seq_for_cf, new_seq_for_cf = old_seq, new_seq
            # if len(new_seq) > len(old_seq):  # i.e if <old_seq> (the template gene) is shorter than the sequence corresponding to the original name for the new allele that we found from it
            #     new_len_str = utils.color('blue', new_seq[len(old_seq):])
            #     new_seq_for_cf = new_seq[:len(old_seq)]
            #     print '         %d extra (blue) bases in new seq were not considered' % (len(new_seq) - len(old_seq))
            # elif len(old_seq) > len(new_seq):
            #     old_len_str = utils.color('blue', old_seq[len(new_seq):])
            #     old_seq_for_cf = old_seq[:len(new_seq)]
            #     print '         %d extra (blue) bases in old seq were not considered' % (len(old_seq) - len(new_seq))
            # print '          %s%s   %s' % (old_seq_for_cf, old_len_str, utils.color_gene(template_gene))
            # print '          %s%s   %s' % (utils.color_mutants(old_seq_for_cf, new_seq_for_cf), new_len_str, utils.color_gene(final_name))

        # and add it to the list of new alleles for this gene
        self.inferred_allele_info.append({
            'template-gene' : template_gene,  # the "immediate" template, not the ancestral one
            'gene' : final_name,  # reminder: <final_name> doesn't necessarily correspond to 'snp-positions'
            'seq' : new_seq,
            'cpos' : utils.cdn_pos(self.glfo, self.region, template_gene),
            'snp-positions' : list(mutfo.keys()),  # reminder: *not* necessarily the same as <final_mutfo>
            'aligned-seq' : None,
            'plot-paths' : []  # not filled until we plot, since there may not be a plotdir defined
        })
        self.new_allele_info.append({
            'gene' : final_name,
            'seq' : new_seq,
            'cpos' : utils.cdn_pos(self.glfo, self.region, template_gene),
            'template-gene' : template_gene,  # the "immediate" template, not the ancestral one
            'remove-template-gene' : remove_template,
        })

    # ----------------------------------------------------------------------------------------
    def print_skip_debug(self, gene, positions, positions_to_try_to_fit):

        def get_skip_str(skip_positions):
            skip_str = ''
            if len(skip_positions) > 0 and len(skip_positions) < 20:
                skip_str = ' (' + ' '.join([str(p) for p in skip_positions]) + ')'
            return skip_str

        glseq = self.glfo['seqs'][self.region][gene]
        too_close_to_ends = list(range(self.n_bases_to_exclude['5p'][gene])) + list(range(len(glseq) - self.n_bases_to_exclude['3p'][gene], len(glseq)))
        not_enough_counts = set(positions) - set(positions_to_try_to_fit) - set(too_close_to_ends)  # well, not enough counts, *and* not too close to the ends

        print('          skipping', end=' ')
        print('%d / %d positions:' % (len(positions) - len(positions_to_try_to_fit), len(positions)), end=' ')
        print('%d were too close to the ends%s' % (len(too_close_to_ends), get_skip_str(too_close_to_ends)), end=' ')
        print('and %d had fewer than %d mutations and fewer than %d observations%s' % (len(not_enough_counts), self.n_muted_min, self.n_total_min, get_skip_str(not_enough_counts)))

    # ----------------------------------------------------------------------------------------
    def print_summary(self, genes_to_use):
        if len(genes_to_use) == 0:
            print('  no genes with enough counts')
            return

        binline, contents_line = self.overall_mute_counts.horizontal_print(bin_centers=True, bin_decimals=0, contents_decimals=0)
        print('   %s mutations:' % self.region)
        print('              %s' % binline)
        # print '                    %s  overall' % contents_line
        for gene in genes_to_use:
            _, contents_line = self.per_gene_mute_counts[gene].horizontal_print(bin_centers=True, bin_decimals=0, contents_decimals=0)
            print('              %s     %s' % (contents_line, utils.color_gene(gene)))

        print('   sequence counts:')
        print('                 excluded             excluded                 excluded             included                       actually')
        print('              >%2d mutations       5p del (>N bases)        3p del (>N bases)       total seqs        clones          used' % self.args.n_max_mutations_per_segment)
        for gene in genes_to_use:
            print('                 %5d            %5d  %5s             %5d  %4s            %7d          %7d        %7d      %s' % (self.n_seqs_too_highly_mutated[gene],
                                                                                                                                    self.n_big_del_skipped['5p'][gene], '(%d)' % self.n_bases_to_exclude['5p'][gene],
                                                                                                                                    self.n_big_del_skipped['3p'][gene], '(%d)' % self.n_bases_to_exclude['3p'][gene],
                                                                                                                                    self.n_clonal_representatives[gene] + self.n_excluded_clonal_queries[gene],
                                                                                                                                    self.n_clones[gene], self.n_clonal_representatives[gene],
                                                                                                                                    utils.color_gene(gene)))
        print('%s: looking for new alleles among %d gene%s (%d genes didn\'t have enough counts)' % (utils.color('blue', 'try ' + str(self.itry)), len(genes_to_use), utils.plural(len(genes_to_use)), len(self.counts) - len(genes_to_use)))

    # ----------------------------------------------------------------------------------------
    def increment_and_finalize(self, swfo, debug=False):
        assert not self.finalized
        start = time.time()
        if debug:
            print('allele finding')

        # first prepare some things, and increment for each chosen query
        self.set_excluded_bases(swfo, debug=debug)
        queries_to_use = [q for q in swfo['queries'] if not self.skip_query(q, swfo[q])]  # skip_query() also fills self.seq_info if we're not skipping the query (and sometimes also if we do skip it)
        if len(queries_to_use) > self.n_max_queries:
            print('  note: subsampling %d queries down to %d before finding alleles (yes this is a hack, it would be better to fix the issues with false positives on super large samples, but this is also fine for now)' % (len(queries_to_use), self.n_max_queries))
            queries_to_use = numpy.random.choice(queries_to_use, size=self.n_max_queries, replace=False)
        if len(queries_to_use) == 0:
            print('  no queries for allele finding')  # NOTE don't return here -- there's some stuff below that should happen

        if debug:
            print('                        total   clones    representatives')
        n_total_clusters = 0
        def keyfunc(q):
            return swfo[q][self.region + '_gene']
        for gene, gene_queries in itertools.groupby(sorted(queries_to_use, key=keyfunc), key=keyfunc):
            gene_queries = list(gene_queries)  # otherwise i can't print the length down there...
            clusters = utils.collapse_naive_seqs(swfo, queries=gene_queries)
            n_representatives = 0
            for cluster in clusters:
                cluster_representatives = self.choose_cluster_representatives(swfo, cluster)
                for qchosen in cluster_representatives:
                    self.increment_query(qchosen, swfo[qchosen][self.region + '_gene'])
                n_representatives += len(cluster_representatives)
            n_total_clusters += len(clusters)
            if debug:
                print('      %s %6d  %6d    %6d' % (utils.color_gene(gene, width=15), len(gene_queries), len(clusters), n_representatives))
        if debug:
            print('    %d seqs chosen to represent %d clones with %d total seqs' % (sum(self.n_clonal_representatives.values()), n_total_clusters, len(queries_to_use)))

        # then finalize
        genes_to_use = [g for g in sorted(self.counts) if self.gene_obs_counts[g] >= self.n_total_min]
        if debug:
            self.print_summary(genes_to_use)
        # if not self.args.always_find_new_alleles:  # NOTE this is (on purpose) summed over all genes -- genes with homozygous unknown alleles would always fail this criterion
        #     total = int(self.overall_mute_counts.integral(include_overflows=True))  # underflow bin is zero mutations, and we want overflow, too NOTE why the fuck isn't this quite equal to sum(self.gene_obs_counts.values())?
        #     for n_mutes in range(self.args.n_max_snps):
        #         if self.overall_mute_counts.bin_contents[n_mutes] / float(total) < self.min_fraction_per_bin:
        #             print '        not looking for new alleles: not enough counts (%d / %d = %.3f < %.3f)' % (self.overall_mute_counts.bin_contents[n_mutes], total, self.overall_mute_counts.bin_contents[n_mutes] / float(total), self.min_fraction_per_bin),
        #             print 'with %d %s mutations (override this with --always-find-new-alleles, or by reducing --n-max-snps)' % (n_mutes, self.region)
        #             self.finalized = True
        #             return

        self.xyvals = {}
        self.positions_to_plot = {gene : set() for gene in self.counts}
        for gene in genes_to_use:
            if debug:
                print(' %s %3d count%s' % (utils.color_gene(gene, width=21), self.gene_obs_counts[gene], utils.plural(self.gene_obs_counts[gene])))
            positions = sorted(self.counts[gene])
            self.xyvals[gene] = {pos : self.get_allele_finding_xyvals(gene, pos) for pos in positions}
            positions_to_try_to_fit = [pos for pos in positions if sum(self.xyvals[gene][pos]['obs']) > self.n_muted_min or sum(self.xyvals[gene][pos]['total']) > self.n_total_min]  # ignore positions with neither enough mutations nor total observations

            # if debug and len(positions) > len(positions_to_try_to_fit):
            #     self.print_skip_debug(gene, positions, positions_to_try_to_fit)

            if len(positions_to_try_to_fit) < self.args.n_max_snps:
                if debug:
                    print('          not enough positions with enough observations to fit %s' % utils.color_gene(gene))
                continue

            # if self.per_gene_mute_counts[gene].bin_contents[0] > self.n_total_min:  # UPDATE this will have to change to accomodate repertoires with very few unmutated sequences
            #     self.alleles_with_evidence.add(gene)

            # loop over each snp hypothesis
            self.already_printed_dbg_header = False
            self.fitfos[gene] = []
            for istart in range(1, self.args.n_max_snps + 1):
                self.fit_istart(gene, istart, positions_to_try_to_fit, debug=debug)

            candidates = []
            for icand in range(len(self.fitfos[gene])):
                candidfo = self.fitfos[gene][icand]
                if debug:
                    if icand == 0:
                        print('   snps       min ratio')
                    print('   %2d     %9s' % (candidfo['istart'], fstr(candidfo['min_snp_ratio'])), end=' ')
                if self.is_a_candidate(candidfo, debug=debug):
                    candidates.append(candidfo)
                if debug:
                    print('')

            # first take the biggest one, then if there's any others that have entirely non-overlapping positions, we don't need to re-run
            already_used_positions = set()
            for candidfo in reversed(candidates):  # biggest <istart> first
                these_positions = set(candidfo['positions'])
                if len(these_positions & already_used_positions) > 0:
                    continue
                alleles_this_gene = [g for g in (list(self.counts.keys()) + [gfo['gene'] for gfo in self.new_allele_info]) if utils.are_alleles(g, gene)]  # <gene> is the template gene
                if len(alleles_this_gene) >= self.n_warn_alleles_per_gene:
                    print('  %s inferred %d alleles for gene %s' % (utils.color('yellow', 'note:'), len(alleles_this_gene), gene))
                if self.args.n_max_alleles_per_gene is not None:  # NOTE we're not really looping over the most likely new alleles first here (and there isn't really a good way to do that), so it's pretty arbitrary which ones get skipped
                    if len(alleles_this_gene) >= self.args.n_max_alleles_per_gene:
                        print('    --n-max-alleles-per-gene: already have %d allele%s for %s (%s), so skipping new inferred allele' % (len(alleles_this_gene), utils.plural(len(alleles_this_gene)), utils.color_gene(gene), ' '.join(utils.color_gene(g) for g in alleles_this_gene)))
                        continue
                already_used_positions |= these_positions
                self.add_allele_to_new_allele_info(gene, candidfo, debug=debug)

        if debug:
            if len(self.inferred_allele_info) > 0:
                print('  found %d new %s: %s' % (len(self.inferred_allele_info), utils.plural_str('allele', len(self.inferred_allele_info)), ' '.join([utils.color_gene(nfo['gene']) for nfo in self.inferred_allele_info])))
            else:
                print('    no new alleles')
            print('  allele finding: %d fits in %.1f sec' % (self.n_fits, time.time()-start))

        self.finalized = True
        return self.new_allele_info

    # ----------------------------------------------------------------------------------------
    def get_fitfo_for_plotting(self, gene, pos):
        # loop over all the new alleles to figure out which one we're supposed to take the fitfo from
        for newfo in self.inferred_allele_info:
            if newfo['template-gene'] != gene:
                continue
            if pos not in newfo['snp-positions']:
                continue
            newfo['plot-paths'].append(utils.sanitize_name(gene) + '/' + str(pos) + '.svg')
            cfos = [candidfo for candidfo in self.fitfos[gene] if candidfo['istart'] == len(newfo['snp-positions']) and pos in candidfo['positions']]
            if len(cfos) == 0:
                print('  shouldn\'t be able to get here if there\'s no candidfos with the right <istart>')
                continue
            return newfo['gene'], cfos[0]['fitfos'][pos]  # just arbitrarily take the first one (I don't think you can really get two -- that would mean a position was shared by more than one new allele)
        return None, None

    # ----------------------------------------------------------------------------------------
    def plot(self, base_plotdir, only_csv=False):
        from . import plotting
        if not self.finalized:
            self.finalize(debug=debug)

        if only_csv:  # not implemented
            print('    <only_csv> not yet implemented in allelefinder')
            return

        start = time.time()
        print('    plotting allele finding', end=' ')
        sys.stdout.flush()

        plotdir = base_plotdir + '/allele-finding/try-' + str(self.itry)

        for old_gene_dir in glob.glob(plotdir + '/*'):  # has to be a bit more hackey than elsewhere, since we have no way of knowing what genes might have had their own directories written last time we wrote to this dir
            if not os.path.isdir(old_gene_dir):
                raise Exception('not a directory: %s' % old_gene_dir)
            utils.prep_dir(old_gene_dir, wildlings=('*.csv', '*.svg'))
            os.rmdir(old_gene_dir)
        utils.prep_dir(plotdir, wildlings=('*.csv', '*.svg'))


        if self.args.plot_and_fit_absolutely_everything is None:
            for gene in self.positions_to_plot:  # we can make plots for the positions we didn't fit, but there's a *lot* of them and they're slow
                for position in self.positions_to_plot[gene]:
                    new_gene, fitfos = self.get_fitfo_for_plotting(gene, position)
                    plotting.make_allele_finding_plot(plotdir + '/' + utils.sanitize_name(gene), gene, position, self.xyvals[gene][position], xmax=self.args.n_max_mutations_per_segment, fitfos=fitfos, new_gene=new_gene)
        else:
            for gene in self.counts:  # we can make plots for the positions we didn't fit, but there's a *lot* of them and they're slow
                for position in sorted(self.counts[gene]):
                    both, pre, post = self.get_both_pre_post_vals(gene, istart=self.args.plot_and_fit_absolutely_everything, positions_to_try_to_fit=[position])
                    big_y_icpt, big_y_icpt_err = self.get_big_y(post[position])
                    big_y_icpt_bounds = self.get_big_y_icpt_bounds(big_y_icpt, big_y_icpt_err)
                    fitfos = {
                        'onefo' : self.get_curvefit(both[position], y_icpt_bounds=(0., 0.), debug=True),
                        'prefo' : self.get_curvefit(pre[position], y_icpt_bounds=(0., 0.), debug=True),
                        'postfo' : self.get_curvefit(post[position], y_icpt_bounds=big_y_icpt_bounds, debug=True),
                    }
                    plotting.make_allele_finding_plot(plotdir + '/' + utils.sanitize_name(gene), gene, position, self.xyvals[gene][position], xmax=self.args.n_max_mutations_per_segment, fitfos=fitfos)

        plotting.make_html(plotdir, fnames=[newfo['plot-paths'] for newfo in self.inferred_allele_info], title=('inferred alleles'))

        print('(%.1f sec)' % (time.time()-start))
