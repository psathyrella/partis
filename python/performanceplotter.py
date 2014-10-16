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
        utils.prep_dir(self.plotdir + '/plots', wildling=None, multilings=['*.csv','*.svg'])
        self.values = {}
        for column in utils.index_columns:
            if column == 'cdr3_length':  # kind of finicky to figure out what this is, so I don't always set it
                continue
            self.values[column] = {}
            if column in bool_columns:
                self.values[column]['right'] = 0
                self.values[column]['wrong'] = 0
            # else:
            #     self.values[column] = {}
        self.values['hamming_to_true_naive'] = {}
    # ----------------------------------------------------------------------------------------
    def hamming_distance_to_true_naive(self, true_line, line, query_name):
        """ hamming distance between the inferred naive sequence and the tue naive sequence """
        true_naive_seq = utils.get_full_naive_seq(self.germlines, true_line)
        inferred_naive_seq = utils.get_full_naive_seq(self.germlines, line)
        # utils.color_mutants(true_naive_seq, inferred_naive_seq, True)
        if len(true_naive_seq) != len(inferred_naive_seq):
            print 'ERROR %s different naive lengths' % query_name
            print '    ', true_naive_seq
            print '    ', inferred_naive_seq
            print '   cropping the longer one (it\'s probably because the j match doesn\'t go all the way to the right side of the sequence)'
            # TODO should really penalize if it's wrong
            min_length = min(len(true_naive_seq), len(inferred_naive_seq))
            true_naive_seq = true_naive_seq[:min_length]
            inferred_naive_seq = inferred_naive_seq[:min_length]
        return utils.hamming(true_naive_seq, inferred_naive_seq)

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
            else:
                trueval, guessval = 0, 0
                if 'insertion' in column:
                    trueval = len(true_line[column])
                    guessval = len(line[column])
                elif column == 'hamming_to_true_naive':
                    trueval = 0  # NOTE this is a kind of weird way to do it, since diff ends up as really just the guessval, but it's ok for now
                    guessval = self.hamming_distance_to_true_naive(true_line, line, query_name)
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
            else:
                hist = plotting.make_hist(self.values[column], 'int', self.name + '-' + column, normalize=True)
                if column == 'hamming_to_true_naive':
                    hist.GetXaxis().SetTitle('hamming distance')
                else:
                    hist.GetXaxis().SetTitle('inferred - true')
                plotting.draw(hist, 'int', plotname=column, plotdir=self.plotdir)
        check_call(['./permissify-www', self.plotdir])  # NOTE this should really permissify starting a few directories higher up
        check_call(['makeHtml', self.plotdir, '3', 'null', 'svg'])
