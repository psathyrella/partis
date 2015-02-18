import sys
import utils
import plotting
from hist import Hist
from subprocess import check_call

assert plotting.check_root()

# Columns for which we just want to know, Did we guess the right value? (for other columns, we store guess - true)
bool_columns = ('v_gene', 'd_gene', 'j_gene')

class PerformancePlotter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, germlines, plotdir, name):
        self.germlines = germlines
        self.plotdir = plotdir
        self.name = name
        utils.prep_dir(self.plotdir + '/plots', wildling=None, multilings=['*.csv', '*.svg', '*.root'])
        self.values = {}
        for column in utils.index_columns:
            if column == 'cdr3_length':  # kind of finicky to figure out what this is, so I don't always set it
                continue
            self.values[column] = {}
            if column in bool_columns:
                self.values[column]['right'] = 0
                self.values[column]['wrong'] = 0
        self.values['hamming_to_true_naive'] = {}
        for region in utils.regions:
            self.values[region + '_hamming_to_true_naive'] = {}
            self.values[region + '_hamming_to_true_naive_normed'] = {}
        # for bound in utils.boundaries:
        #     self.counts[bound + '_insertion_content'] = {'A':0, 'C':0, 'G':0, 'T':0}  # base content of each insertion
        # self.counts['seq_content'] = {'A':0, 'C':0, 'G':0, 'T':0}
        # n_bins, xmin, xmax = 100, 0.0, 1.0
        self.hists = {}
        self.hists['mute_freqs'] = Hist(30, -0.05, 0.05)

    # ----------------------------------------------------------------------------------------
    def hamming_distance_to_true_naive(self, true_line, line, query_name, restrict_to_region='', normalize=False, debug=False):
        """
        Hamming distance between the inferred naive sequence and the tue naive sequence.
        <restrict_to_region> if set, restrict the comparison to the section of the *true* sequence assigned to the given region.
        NOTE this will not in general correspond to the similarly-assigned region in the inferred naive sequence.
        if <normalize> divide by sequence length
        """

        true_naive_seq = utils.get_full_naive_seq(self.germlines, true_line)
        inferred_naive_seq = utils.get_full_naive_seq(self.germlines, line)

        left_hack_add_on = ''
        right_hack_add_on = ''
        if len(true_line['seq']) > len(line['seq']):  # ihhhmmm doesn't report the bits of the sequence it erodes off the ends, so we have to add them back on
        # if len(true_naive_seq) > len(inferred_naive_seq):  # hm, now why I did use line['seq'] stuff before?
            start = true_line['seq'].find(line['seq'])
            assert start >= 0
            end = len(line['seq']) + start
            left_hack_add_on = true_line['seq'][: start]
            right_hack_add_on = true_line['seq'][ end :]
            # extra_penalty = len(left_hack_add_on) + len(right_hack_add_on)
            inferred_naive_seq = 'x'*len(left_hack_add_on) + inferred_naive_seq + 'x'*len(right_hack_add_on)
            if debug:
                print '  adding to inferred naive seq'

        bounds = None
        if restrict_to_region != '':
            bounds = utils.get_regional_naive_seq_bounds(restrict_to_region, self.germlines, true_line)  # get the bounds of this *true* region
            true_naive_seq = true_naive_seq[bounds[0] : bounds[1]]
            inferred_naive_seq = inferred_naive_seq[bounds[0] : bounds[1]]


        # if len(true_naive_seq) > len(inferred_naive_seq):
            

        if debug:
            print restrict_to_region, 'region, bounds', bounds
            print '  true ', true_naive_seq
            print '  infer', inferred_naive_seq

        if len(true_naive_seq) != len(inferred_naive_seq):
            print 'ERROR still not the same lengths for %s' % query_name
            print '  true ', true_naive_seq
            print '  infer', inferred_naive_seq
            sys.exit()
        total_distance = utils.hamming(true_naive_seq, inferred_naive_seq)
        if len(true_naive_seq) == 0:
            print 'WARNING zero length sequence in hamming_distance_to_true_naive'
            return 0
        if normalize:
            return int(100 * (float(total_distance) / len(true_naive_seq)))
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
    def add_partial_fail(self, true_line, line):
        for column in self.values:
            if column in bool_columns:
                if column in line and utils.are_alleles(true_line[column], line[column]):  # NOTE you have to change this below as well!
                    self.values[column]['right'] += 1
                    # if column == 'v_gene':
                    #     print '  partial right ', true_line[column], line[column]
                else:
                    self.values[column]['wrong'] += 1
                    # if column == 'v_gene':
                    #     print '  partial wrong ', true_line[column], line[column] if column in line else 'FOO'
            else:
                pass

    # ----------------------------------------------------------------------------------------
    def evaluate(self, true_line, inf_line):
        for column in self.values:
            if column in bool_columns:
                if utils.are_alleles(true_line[column], inf_line[column]):  # NOTE you have to change this above as well!
                # if true_line[column] == inf_line[column]:
                    self.values[column]['right'] += 1
                else:
                    self.values[column]['wrong'] += 1
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
                    normalize = '_norm' in column
                    guessval = self.hamming_distance_to_true_naive(true_line, inf_line, inf_line['unique_id'], restrict_to_region=restrict_to_region, normalize=normalize)
                else:
                    trueval = int(true_line[column])
                    guessval = int(inf_line[column])

                diff = guessval - trueval
                if diff not in self.values[column]:
                    self.values[column][diff] = 0
                self.values[column][diff] += 1

        for column in self.hists:
            trueval = utils.get_mutation_rate(self.germlines, true_line)
            guessval = utils.get_mutation_rate(self.germlines, inf_line)
            self.hists[column].fill(guessval - trueval)

    # ----------------------------------------------------------------------------------------
    def plot(self):
        for column in self.values:
            if column in bool_columns:
                right = self.values[column]['right']
                wrong = self.values[column]['wrong']
                print '  %s\n    correct up to allele: %4d / %-4d = %4.2f' % (column, right, right+wrong, float(right) / (right + wrong))
                hist = plotting.make_bool_hist(right, wrong, self.name + '-' + column)
                plotting.draw(hist, 'bool', plotname=column, plotdir=self.plotdir, write_csv=True)
            else:
                hist = plotting.make_hist_from_dict_of_counts(self.values[column], 'int', self.name + '-' + column, normalize=True)
                log = ''
                if column.find('hamming_to_true_naive') >= 0:
                    hist.GetXaxis().SetTitle('hamming distance')
                else:
                    hist.GetXaxis().SetTitle('inferred - true')
                plotting.draw(hist, 'int', plotname=column, plotdir=self.plotdir, write_csv=True, log=log)
        for column in self.hists:
            hist = plotting.make_hist_from_my_hist_class(self.hists[column], 'mute_freqs')
            plotting.draw(hist, 'float', plotname=column, plotdir=self.plotdir, write_csv=True, log=log)
        
        check_call(['./bin/makeHtml', self.plotdir, '3', 'null', 'svg'])
        check_call(['./bin/permissify-www', self.plotdir])  # NOTE this should really permissify starting a few directories higher up
