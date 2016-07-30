import sys
import utils
import numpy
import plotting
import re
from hist import Hist
from subprocess import check_call
import fraction_uncertainty
import copy
import math

# Columns for which we just want to know, Did we guess the right value? (for other columns, we store guess - true)
bool_columns = ('v_gene', 'd_gene', 'j_gene')

class PerformancePlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, name, only_correct_gene_fractions=False):
        self.name = name
        self.values, self.hists = {}, {}
        self.only_correct_gene_fractions = only_correct_gene_fractions
        for column in tuple(list(utils.index_columns) + ['cdr3_length', ]):
            if column == 'cdr3_length':  # kind of finicky to figure out what this is, so I don't always set it
                continue
            self.values[column] = {}
            if column in bool_columns:
                self.values[column]['right'] = 0
                self.values[column]['wrong'] = 0
        self.values['hamming_to_true_naive'] = {}
        self.hists['hamming_to_true_naive_normed'] = Hist(25, 0., 0.5)
        for region in utils.regions:
            self.values[region + '_hamming_to_true_naive'] = {}
            self.hists[region + '_hamming_to_true_naive_normed'] = Hist(25, 0., 0.5)
        self.hists['mute_freqs'] = Hist(30, -0.05, 0.05, xtitle='inferred - true')
        for region in utils.regions:
            self.hists[region + '_mute_freqs'] = Hist(30, -0.05, 0.05, xtitle='inferred - true')
        for region in utils.regions:  # plots of correct gene calls vs mute freq
            self.hists[region + '_gene_right_vs_mute_freq'] = Hist(25, 0., 0.4)
            self.hists[region + '_gene_wrong_vs_mute_freq'] = Hist(25, 0., 0.4)
            self.hists[region + '_allele_right_vs_per_gene_support'] = Hist(25, 0., 1.)  # NOTE these require the correct *allele*, while the other ones don't
            self.hists[region + '_allele_wrong_vs_per_gene_support'] = Hist(25, 0., 1.)

    # ----------------------------------------------------------------------------------------
    def hamming_distance_to_true_naive(self, true_line, line, restrict_to_region='', normalize=False, padfo=None, debug=False):
        """
        Hamming distance between the inferred naive sequence and the tue naive sequence.
        <restrict_to_region> if set, restrict the comparison to the section of the *true* sequence assigned to the given region.
        NOTE this will not in general correspond to the similarly-assigned region in the inferred naive sequence.
        if <normalize> divide by sequence length
        """

        true_naive_seq = true_line['naive_seq']
        inferred_naive_seq = line['naive_seq']
        if len(true_naive_seq) != len(inferred_naive_seq):
            print '%20s    true      inf' % ''
            for k in true_line:
                print '%20s   %s' % (k, true_line[k]),
                if k in line:
                    print '   %s' % line[k]
                else:
                    print '    NOPE'
            for k in line:
                if k not in true_line:
                    print '  not in true line   %20s    %s' % (k, line[k])
            raise Exception('%s true and inferred sequences not the same length\n   %s\n   %s\n' % (line['unique_id'], true_naive_seq, inferred_naive_seq))

        # assert False # read through this whole damn thing and make sure it's ok

        left_hack_add_on = ''
        right_hack_add_on = ''
        # if len(true_line['seq']) > len(utils.remove_ambiguous_ends(line['seq'], line['fv_insertion'], line['jf_insertion'])):  # ihhhmmm doesn't report the bits of the sequence it erodes off the ends, so we have to add them back on
        # # if len(true_naive_seq) > len(inferred_naive_seq):  # hm, now why did I use line['seq'] stuff before?
        #     assert False
        #     start = true_line['seq'].find(line['seq'])
        #     assert start >= 0
        #     end = len(line['seq']) + start
        #     left_hack_add_on = true_line['seq'][: start]
        #     right_hack_add_on = true_line['seq'][ end :]
        #     # extra_penalty = len(left_hack_add_on) + len(right_hack_add_on)
        #     inferred_naive_seq = 'N'*len(left_hack_add_on) + inferred_naive_seq + 'N'*len(right_hack_add_on)
        #     if debug:
        #         print '  adding to inferred naive seq'


        if padfo is not None:  # remove N padding from the inferred sequence
            if debug:
                print 'removing padfo'
                print inferred_naive_seq
            if inferred_naive_seq[padfo['padleft'] : ].count('N') == padfo['padleft']:  # this fails to happen if reset_effective_erosions_and_effective_insertions already removed the Ns
                inferred_naive_seq = inferred_naive_seq[padfo['padleft'] : ]
            elif debug:  # NOTE if no debug, we just fall through, which isok
                print 'tried to remove non Ns!\n   %s\n   padleft %d\n' % (inferred_naive_seq, padfo['padleft'])
            if padfo['padright'] > 0:
                if inferred_naive_seq[ : padfo['padright']].count('N') == padfo['padright']:  # this fails to happen if reset_effective_erosions_and_effective_insertions already removed the Ns
                    inferred_naive_seq = inferred_naive_seq[ : -padfo['padright']]
                elif debug:  # NOTE if no debug, we just fall through, which isok
                    print 'tried to remove non Ns!\n   %s\n   padright %d\n' % (inferred_naive_seq, padfo['padright'])
            if debug:
                print padfo['padleft'] * ' ' + inferred_naive_seq + padfo['padleft'] * ' '

        bounds = None
        if restrict_to_region != '':
            bounds = true_line['regional_bounds'][restrict_to_region]
            if debug:
                print 'restrict to %s' % restrict_to_region
                utils.color_mutants(true_naive_seq, inferred_naive_seq, print_result=True, extra_str='      ')
                utils.color_mutants(true_naive_seq[bounds[0] : bounds[1]], inferred_naive_seq[bounds[0] : bounds[1]], print_result=True, extra_str='      ' + bounds[0]*' ')
            true_naive_seq = true_naive_seq[bounds[0] : bounds[1]]
            inferred_naive_seq = inferred_naive_seq[bounds[0] : bounds[1]]

        if len(true_naive_seq) != len(inferred_naive_seq):
            raise Exception('still not the same lengths for %s\n  %s\n  %s' % (line['unique_ids'][0], true_naive_seq, inferred_naive_seq))
        fraction, len_excluding_ambig = utils.hamming_fraction(true_naive_seq, inferred_naive_seq, return_len_excluding_ambig=True)
        total_distance = int(fraction * len_excluding_ambig)
        if len(true_naive_seq) == 0:
            if not (restrict_to_region == 'd' and utils.get_chain(true_line['v_gene']) != 'h'):
                print 'WARNING zero length sequence in hamming_distance_to_true_naive'
            return 0
        if normalize:
            return float(total_distance) / len(true_naive_seq)
        else:
            return total_distance

    # ----------------------------------------------------------------------------------------
    def add_fail(self):
        for column in self.values:
            if column in bool_columns:
                self.values[column]['wrong'] += 1
            else:
                pass

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
        if inf_line[region + '_per_gene_support'].keys()[0] != inf_line[region + '_gene']:
            print '   WARNING best-supported gene %s not same as viterbi gene %s' % (utils.color_gene(inf_line[region + '_per_gene_support'].keys()[0]), utils.color_gene(inf_line[region + '_gene']))
        support = inf_line[region + '_per_gene_support'].values()[0]  # sorted, ordered dict with gene : logprob key-val pairs
        if true_line[region + '_gene'] == inf_line[region + '_gene']:  # NOTE this requires allele to be correct, but set_bool_column() does not
            self.hists[region + '_allele_right_vs_per_gene_support'].fill(support)
        else:
            self.hists[region + '_allele_wrong_vs_per_gene_support'].fill(support)

    # ----------------------------------------------------------------------------------------
    def add_partial_fail(self, true_line, line):
        # NOTE does not fill all the hists ('cause it kind of can't, right?)

        overall_mute_freq = utils.get_mutation_rate(true_line, iseq=0)  # true value

        for column in self.values:
            if column in bool_columns:
                if column in line:
                    self.set_bool_column(true_line, line, column, overall_mute_freq)
            else:
                pass

        for region in utils.regions:
            if region + '_per_gene_support' in inf_line:
                self.set_per_gene_support(true_line, inf_line, region)

    # ----------------------------------------------------------------------------------------
    def evaluate(self, true_line, inf_line, padfo=None):

        overall_mute_freq = utils.get_mutation_rate(true_line, iseq=0)  # true value

        for column in self.values:
            if self.only_correct_gene_fractions and column not in bool_columns:
                continue
            if column in bool_columns:
                self.set_bool_column(true_line, inf_line, column, overall_mute_freq)
            else:
                trueval, guessval = 0, 0
                if column[2:] == '_insertion':  # insertion length
                    trueval = len(true_line[column])
                    guessval = len(inf_line[column])
                # elif '_content' in column:
                #     seq_to_use = inf_line[column[ : column.find('_', 3)]]  # NOTE has to work for seq_content *and* vd_insertion_content, hence the 3
                #         for nuke in seq_to_use:
                #             self.counts[col][nuke] += 1
                elif 'hamming_to_true_naive' in column:
                    trueval = 0  # NOTE this is a kind of weird way to do it, since diff ends up as really just the guessval, but it does the job
                    restrict_to_region = column[0].replace('h', '')  # if fist char in <column> is not an 'h', restrict to that region
                    assert '_norm' not in column  # moved these to <self.hists>
                    guessval = self.hamming_distance_to_true_naive(true_line, inf_line, normalize=False, restrict_to_region=restrict_to_region, padfo=padfo)
                else:
                    trueval = int(true_line[column])
                    guessval = int(inf_line[column])

                diff = guessval - trueval
                if diff not in self.values[column]:
                    self.values[column][diff] = 0
                self.values[column][diff] += 1

        for region in utils.regions:
            if region + '_per_gene_support' in inf_line:
                self.set_per_gene_support(true_line, inf_line, region)

        for column in [c for c in self.hists if 'hamming_to_true' in c]:
            restrict_to_region = column[0] if column[0] in utils.regions else ''
            hfrac = self.hamming_distance_to_true_naive(true_line, inf_line, normalize=True, restrict_to_region=restrict_to_region, padfo=padfo)
            self.hists[column].fill(hfrac)

        for column in [c for c in self.hists if 'mute_freqs' in c]:  # NOTE not '_vs_mute_freq' (plural is important)
            if '_vs_mute_freq' in column or '_per_gene_support' in column:  # fill these above
                continue
            if len(re.findall('[vdj]_', column)) == 1:
                region = re.findall('[vdj]_', column)[0][0]
            else:
                region = ''
            trueval = utils.get_mutation_rate(true_line, iseq=0, restrict_to_region=region)
            guessval = utils.get_mutation_rate(inf_line, iseq=0, restrict_to_region=region)
            self.hists[column].fill(guessval - trueval)

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, only_csv=False):
        utils.prep_dir(plotdir, wildlings=('*.csv', '*.svg'))
        for column in self.values:
            if self.only_correct_gene_fractions and column not in bool_columns:
                continue
            if column in bool_columns:
                right = self.values[column]['right']
                wrong = self.values[column]['wrong']
                errs = fraction_uncertainty.err(right, right+wrong)
                print '  %s\n    correct up to allele: %4d / %-4d = %4.4f (-%.3f, +%.3f)' % (column, right, right+wrong, float(right) / (right + wrong), errs[0], errs[1])
                hist = plotting.make_bool_hist(right, wrong, self.name + '-' + column)
                plotting.draw_no_root(hist, plotname=column, plotdir=plotdir, write_csv=True, stats='0-bin', only_csv=only_csv)
            else:
                # TODO this is dumb... I should make the integer-valued ones histograms as well
                hist = plotting.make_hist_from_dict_of_counts(self.values[column], 'int', self.name + '-' + column, normalize=True)
                log = ''
                xtitle = 'hamming distance' if 'hamming_to_true_naive' in column else 'inferred - true'
                plotting.draw_no_root(hist, plotname=column, plotdir=plotdir, write_csv=True, log=log, only_csv=only_csv, xtitle=xtitle)

        for column in self.hists:
            if '_vs_mute_freq' in column or '_vs_per_gene_support' in column:  # only really care about the fraction, which we plot below
                continue
            plotting.draw_no_root(self.hists[column], plotname=column, plotdir=plotdir, write_csv=True, log=log, only_csv=only_csv)

        # fraction correct vs mute freq
        for region in utils.regions:
            hright = self.hists[region + '_gene_right_vs_mute_freq']
            hwrong = self.hists[region + '_gene_wrong_vs_mute_freq']
            if hright.integral(include_overflows=True) == 0:
                continue
            plotting.make_fraction_plot(hright, hwrong, plotdir, region + '_fraction_correct_vs_mute_freq', xlabel='mut freq', ylabel='fraction correct up to allele', xbounds=(0., 0.5), write_csv=True)

        # per-gene support crap
        for region in utils.regions:
            if self.hists[region + '_allele_right_vs_per_gene_support'].integral(include_overflows=True) == 0:
                continue
            hright = self.hists[region + '_allele_right_vs_per_gene_support']
            hwrong = self.hists[region + '_allele_wrong_vs_per_gene_support']
            plotting.make_fraction_plot(hright, hwrong, plotdir, region + '_allele_fraction_correct_vs_per_gene_support', xlabel='support', ylabel='fraction with correct allele', xbounds=(-0.1, 1.1), write_csv=True)

        if not only_csv:
            plotting.make_html(plotdir)
