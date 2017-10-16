import time
import sys
import utils
import numpy
import re
from subprocess import check_call
import copy
import math

import plotconfig
from hist import Hist
import indelutils

class PerformancePlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, name):
        self.name = name
        self.values, self.hists = {}, {}  # the dictionary-based approach in <self.values> is nice because you can decide your hist bounds after filling everything

        for column in plotconfig.gene_usage_columns:
            self.values[column] = {'right' : 0, 'wrong' : 0}
        for column in plotconfig.int_columns:  # it might be nicer to eventually switch these to hists (I think the ony reason they're separte is that they predate the existence of the hist class)
            self.values[column] = {}
        for rstr in plotconfig.rstrings:
            self.values[rstr + 'hamming_to_true_naive'] = {}
            self.values[rstr + 'muted_bases'] = {}
        self.values['shm_indel_length'] = {}

        self.hists['mute_freqs'] = Hist(25, -0.04, 0.04)  # only do mutation frequency for the whole sequence
        # NOTE this hist bounds here are intended to be super inclusive, whereas in compare-plotdirs.py we apply the more-restrictive ones from plotconfig.py (we still shift overflows here, where appropriate, though)
        for region in utils.regions:
            self.hists[region + '_gene_right_vs_mute_freq'] = Hist(25, 0., 0.4)  # correct *up* to allele (i.e. you can get the allele wrong)
            self.hists[region + '_gene_wrong_vs_mute_freq'] = Hist(25, 0., 0.4)
            self.hists[region + '_allele_right_vs_per_gene_support'] = Hist(25, 0., 1.)  # whereas these require the *correct* allele
            self.hists[region + '_allele_wrong_vs_per_gene_support'] = Hist(25, 0., 1.)

        self.subplotdirs = ['gene-call', 'mutation', 'boundaries']

    # ----------------------------------------------------------------------------------------
    def hamming_to_true_naive(self, true_line, line, restrict_to_region=''):
        true_naive_seq = true_line['naive_seq']
        inferred_naive_seq = line['naive_seq']
        if len(true_naive_seq) != len(inferred_naive_seq):
            raise Exception('different length true and inferred naive seqs for %s\n  %s\n  %s' % (' '.join(line['unique_ids']), true_line['naive_seq'], line['naive_seq']))
        if restrict_to_region != '':  # NOTE very similar to utils.get_n_muted(), except, we want to use the true bounds for both true and naive sequences
            if restrict_to_region in utils.regions:
                bounds = true_line['regional_bounds'][restrict_to_region]
            elif restrict_to_region == 'cdr3':
                bounds = (true_line['codon_positions']['v'], true_line['codon_positions']['j'] + 3)
            else:
                print 'invalid regional restriction %s' % restrict_to_region
            if restrict_to_region == 'v':  # NOTE this is kind of hackey, especially treating v differently to d and j, but it kind of makes sense -- v is fundamentally different in that germline v is a real source of diversity, so it makes sense to isolate the v germline accuracy from the boundary-call accuracy like this
                bounds = (bounds[0], bounds[1] - 3)  # most of the boundary uncertainty is in the last three bases
            true_naive_seq = true_naive_seq[bounds[0] : bounds[1]]
            inferred_naive_seq = inferred_naive_seq[bounds[0] : bounds[1]]
        return utils.hamming_distance(true_naive_seq, inferred_naive_seq)

    # ----------------------------------------------------------------------------------------
    def add_fail(self):
        for column in plotconfig.gene_usage_columns:
            self.values[column]['wrong'] += 1

    # ----------------------------------------------------------------------------------------
    def set_bool_column(self, true_line, inf_line, column, overall_mute_freq):
        if utils.are_alleles(true_line[column], inf_line[column]):  # NOTE this doesn't require allele to be correct, but set_per_gene_support() does
            self.values[column]['right'] += 1
            self.hists[column + '_right_vs_mute_freq'].fill(overall_mute_freq)  # NOTE this'll toss a KeyError if you add bool column that aren't [vdj]_gene
        else:
            self.values[column]['wrong'] += 1
            self.hists[column + '_wrong_vs_mute_freq'].fill(overall_mute_freq)

    # ----------------------------------------------------------------------------------------
    def set_per_gene_support(self, true_line, inf_line, region):
        # if inf_line[region + '_per_gene_support'].keys()[0] != inf_line[region + '_gene']:
        #     print '   WARNING best-supported gene %s not same as viterbi gene %s' % (utils.color_gene(inf_line[region + '_per_gene_support'].keys()[0]), utils.color_gene(inf_line[region + '_gene']))
        support = inf_line[region + '_per_gene_support'].values()[0]  # sorted, ordered dict with gene : logprob key-val pairs
        if true_line[region + '_gene'] == inf_line[region + '_gene']:  # NOTE this requires allele to be correct, but set_bool_column() does not
            self.hists[region + '_allele_right_vs_per_gene_support'].fill(support)
        else:
            self.hists[region + '_allele_wrong_vs_per_gene_support'].fill(support)

    # ----------------------------------------------------------------------------------------
    def evaluate(self, true_line, inf_line, simglfo=None):
        if len(inf_line['unique_ids']) > 1:
            raise Exception('mutli-seq lines not yet handled')
        iseq = 0

        def addval(col, simval, infval):
            if col[2:] == '_insertion':  # stored as the actual inserted bases
                simval = len(simval)
                infval = len(infval)
            diff = infval - simval
            if diff not in self.values[col]:
                self.values[col][diff] = 0
            self.values[col][diff] += 1

        if indelutils.has_indels(true_line['indelfos'][iseq]) or indelutils.has_indels(inf_line['indelfos'][iseq]):
            simlen = indelutils.net_length(true_line['indelfos'][iseq])
            inflen = indelutils.net_length(inf_line['indelfos'][iseq])
            addval('shm_indel_length', simlen, inflen)
            if simlen != inflen:  # this is probably because the simulated shm indel was within the cdr3, so we attempt to fix it by switching the sim line to non-reversed
                print '    %s true and inferred shm net indel lengths different, so skipping rest of performance evaluation' % ' '.join(inf_line['unique_ids'])
                # note that you can't really evaluate the rest of the performance vars in a particularly meaningful when the indel info is different (like I tried to do below) since you have to decide how to assign the indel'd bases (like, is it correct to assign the indel'd bases to a deletion? or to an insertion? or to the j?)
                return
                # true_line = copy.deepcopy(true_line)
                # utils.remove_all_implicit_info(true_line)
                # true_line['indelfos'][iseq] = indelutils.get_empty_indel()
                # true_line['seqs'][iseq] = true_line['input_seqs'][iseq]
                # utils.add_implicit_info(simglfo, true_line)

        mutfo = {lt : {mt : {} for mt in ['freq', 'total']} for lt in ['sim', 'inf']}
        for rstr in plotconfig.rstrings:
            if rstr == '':  # these are already in the <line>s, so may as well not recalculate
                mutfo['sim']['freq'][rstr], mutfo['sim']['total'][rstr] = true_line['mut_freqs'][iseq], true_line['n_mutations'][iseq]
                mutfo['inf']['freq'][rstr], mutfo['inf']['total'][rstr] = inf_line['mut_freqs'][iseq], inf_line['n_mutations'][iseq]
            else:
                mutfo['sim']['freq'][rstr], mutfo['sim']['total'][rstr] = utils.get_mutation_rate_and_n_muted(true_line, iseq=iseq, restrict_to_region=rstr.rstrip('_'))
                mutfo['inf']['freq'][rstr], mutfo['inf']['total'][rstr] = utils.get_mutation_rate_and_n_muted(inf_line, iseq=iseq, restrict_to_region=rstr.rstrip('_'))

        for col in plotconfig.gene_usage_columns:
            self.set_bool_column(true_line, inf_line, col, mutfo['sim']['freq'][''])  # this also sets the fraction-correct-vs-mute-freq hists

        for column in plotconfig.int_columns:
            addval(column, true_line[column], inf_line[column])

        for rstr in plotconfig.rstrings:
            addval(rstr + 'hamming_to_true_naive', 0, self.hamming_to_true_naive(true_line, inf_line, restrict_to_region=rstr.rstrip('_')))
            addval(rstr + 'muted_bases', mutfo['sim']['total'][rstr], mutfo['inf']['total'][rstr])

        # if self.hamming_to_true_naive(true_line, inf_line, restrict_to_region='v') > 5:
        #     print '%20s %2d  %s %s' % (inf_line['unique_ids'][0], self.hamming_to_true_naive(true_line, inf_line, restrict_to_region='v'),
        #                                utils.color_gene(true_line['v_gene'], width=20), 'ok' if inf_line['v_gene'] == true_line['v_gene'] else utils.color_gene(inf_line['v_gene'], width=20))

        for region in utils.regions:
            if region + '_per_gene_support' in inf_line:
                self.set_per_gene_support(true_line, inf_line, region)

        self.hists['mute_freqs'].fill(mutfo['inf']['freq'][''] - mutfo['sim']['freq'][''])  # when we're evaluating on multi-seq hmm output, we synthesize single-sequence lines for each sequence

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, only_csv=False):
        print '  plotting performance',
        import fraction_uncertainty
        import plotting
        start = time.time()
        for substr in self.subplotdirs:
            utils.prep_dir(plotdir + '/' + substr, wildlings=('*.csv', '*.svg'))

        for column in self.values:
            if column in plotconfig.gene_usage_columns:
                right = self.values[column]['right']
                wrong = self.values[column]['wrong']
                lo, hi = fraction_uncertainty.err(right, right + wrong)
                hist = plotting.make_bool_hist(right, wrong, self.name + '-' + column)
                plotting.draw_no_root(hist, plotname=column, plotdir=plotdir + '/gene-call', write_csv=True, stats='0-bin', only_csv=only_csv)
            else:
                hist = plotting.make_hist_from_dict_of_counts(self.values[column], 'int', self.name + '-' + column, normalize=False)
                if 'hamming_to_true_naive' in column:
                    xtitle = 'hamming distance'
                    tmpplotdir = plotdir + '/mutation'
                else:
                    xtitle = 'inferred - true'
                    if 'muted' in column:
                        tmpplotdir = plotdir + '/mutation'
                    else:
                        tmpplotdir = plotdir + '/boundaries'
                plotting.draw_no_root(hist, plotname=column, plotdir=tmpplotdir, write_csv=True, only_csv=only_csv, xtitle=xtitle, shift_overflows=True)

        for column in self.hists:
            if '_vs_mute_freq' in column or '_vs_per_gene_support' in column:  # only really care about the fraction, which we plot below
                continue
            plotting.draw_no_root(self.hists[column], plotname=column, plotdir=plotdir + '/mutation', write_csv=True, only_csv=only_csv, ytitle='counts', xtitle='inferred - true', shift_overflows=True)

        # fraction correct vs mute freq
        for region in utils.regions:
            hright = self.hists[region + '_gene_right_vs_mute_freq']
            hwrong = self.hists[region + '_gene_wrong_vs_mute_freq']
            if hright.integral(include_overflows=True) == 0:
                continue
            plotting.make_fraction_plot(hright, hwrong, plotdir + '/gene-call', region + '_fraction_correct_vs_mute_freq', xlabel='mut freq', ylabel='fraction correct up to allele', xbounds=(0., 0.5), only_csv=only_csv, write_csv=True)

        # per-gene support stuff
        for region in utils.regions:
            if self.hists[region + '_allele_right_vs_per_gene_support'].integral(include_overflows=True) == 0:
                continue
            hright = self.hists[region + '_allele_right_vs_per_gene_support']
            hwrong = self.hists[region + '_allele_wrong_vs_per_gene_support']
            plotting.make_fraction_plot(hright, hwrong, plotdir + '/gene-call', region + '_allele_fraction_correct_vs_per_gene_support', xlabel='support', ylabel='fraction with correct allele', xbounds=(-0.1, 1.1), only_csv=only_csv, write_csv=True)

        if not only_csv:  # write html file and fix permissiions
            for substr in self.subplotdirs:
                plotting.make_html(plotdir + '/' + substr, n_columns=4)

        print '(%.1f sec)' % (time.time()-start)
