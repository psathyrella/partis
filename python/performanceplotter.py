import sys
import utils
import plotting
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
        self.values['mute_freqs'] = {}
    # ----------------------------------------------------------------------------------------
    def hamming_distance_to_true_naive(self, true_line, line, query_name):
        """ hamming distance between the inferred naive sequence and the tue naive sequence """
        true_naive_seq = utils.get_full_naive_seq(self.germlines, true_line)
        inferred_naive_seq = utils.get_full_naive_seq(self.germlines, line)
        if line['fv_insertion'] != '':
            assert true_line['fv_insertion'] == ''
        # if len(true_naive_seq) != len(inferred_naive_seq):  # probably because we inferred a v left or j right deletion that wasn't there
        #     assert False
        #     # TODO should really penalize if the length is wrong, otherwise you can improve your hamming distance by making your match really short
        #     print 'v_5p true %d inferred %d' % (int(true_line['v_5p_del']), int(line['v_5p_del']))
        #     if int(line['v_5p_del']) > 0 and int(true_line['v_5p_del']) == 0:
        #         true_naive_seq = true_naive_seq[int(line['v_5p_del']) : ]

        #     print 'j_3p true %d inferred %d' % (int(true_line['j_3p_del']), int(line['j_3p_del']))
        #     if int(line['j_3p_del']) > 0 and int(true_line['j_3p_del']) == 0:
        #         true_naive_seq = true_naive_seq[ : -int(line['j_3p_del'])]

        #     if len(true_naive_seq) != len(inferred_naive_seq):
        #         print 'ERROR lengths still not equal in %s' % query_name
        #         print true_naive_seq
        #         print inferred_naive_seq

        if utils.hamming(true_naive_seq, inferred_naive_seq) > 200:
            print 'huge hamming %d for %s' % (utils.hamming(true_naive_seq, inferred_naive_seq), query_name)
            utils.color_mutants(true_naive_seq, inferred_naive_seq, True)

        return utils.hamming(true_naive_seq, inferred_naive_seq)

    # ----------------------------------------------------------------------------------------
    def mutation_rate(self, line):
        naive_seq = utils.get_full_naive_seq(self.germlines, line)
        muted_seq = line['seq']
        # print ''
        # utils.color_mutants(naive_seq, muted_seq, True)
        n_mutes = utils.hamming(naive_seq, muted_seq)
        return int(100 * float(n_mutes) / len(naive_seq))  # utils.hamming() asserts they're the same length

    # ----------------------------------------------------------------------------------------
    def add_fail(self):
        for column in self.values:
            if column in bool_columns:
                self.values[column]['wrong'] += 1
            else:
                print 'perfplotter: not sure what to do with a fail'

    # ----------------------------------------------------------------------------------------
    def evaluate(self, true_line, line, query_name):
        for column in self.values:
            if column in bool_columns:
                if true_line[column] == line[column]:
                    self.values[column]['right'] += 1
                else:
                    self.values[column]['wrong'] += 1
                    # if column == 'v_gene':
                    #     print 'wrong', query_name
            else:
                trueval, guessval = 0, 0
                if 'insertion' in column:
                    trueval = len(true_line[column])
                    guessval = len(line[column])
                elif column == 'hamming_to_true_naive':
                    trueval = 0  # NOTE this is a kind of weird way to do it, since diff ends up as really just the guessval, but it's ok for now
                    guessval = self.hamming_distance_to_true_naive(true_line, line, query_name)
                elif column == 'mute_freqs':
                    trueval = self.mutation_rate(true_line)
                    guessval = self.mutation_rate(line)
                else:
                    trueval = int(true_line[column])
                    guessval = int(line[column])
                diff = guessval - trueval
                if diff not in self.values[column]:
                    self.values[column][diff] = 0
                self.values[column][diff] += 1

    # ----------------------------------------------------------------------------------------
    def plot(self):
        for column in self.values:
            if column in bool_columns:
                right = self.values[column]['right']
                wrong = self.values[column]['wrong']
                print '  %s\n    correct: %4d / %-4d = %4.2f' % (column, right, right+wrong, float(right) / (right + wrong))
                hist = plotting.make_bool_hist(right, wrong, self.name + '-' + column)
                plotting.draw(hist, 'bool', plotname=column, plotdir=self.plotdir, write_csv=True)
            else:
                hist = plotting.make_hist(self.values[column], 'int', self.name + '-' + column, normalize=True)
                if column == 'hamming_to_true_naive':
                    hist.GetXaxis().SetTitle('hamming distance')
                else:
                    hist.GetXaxis().SetTitle('inferred - true')
                plotting.draw(hist, 'int', plotname=column, plotdir=self.plotdir, write_csv=True)
        check_call(['./permissify-www', self.plotdir])  # NOTE this should really permissify starting a few directories higher up
        check_call(['makeHtml', self.plotdir, '3', 'null', 'svg'])
