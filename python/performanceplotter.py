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
            else:
                self.values[column] = {}
    # ----------------------------------------------------------------------------------------
    def hamming_distance_to_naive(self, true_line, line):
        original_seqs = {}  # original (non-eroded) germline seqs
        lengths = {}  # length of each match (including erosion)
        eroded_seqs = {}  # eroded germline seqs
        utils.get_reco_event_seqs(self.germlines, line, original_seqs, lengths, eroded_seqs)
        for seq in eroded_seqs.items():
            print seq

    # ----------------------------------------------------------------------------------------
    def evaluate(self, true_line, line):
        self.hamming_distance_to_naive(true_line, line)
        for column in self.values:
            trueval = true_line[column]
            guessval = line[column]
            if column in bool_columns:
                # if utils.are_alleles(guessval, trueval):
                if guessval == trueval:
                    self.values[column]['right'] += 1
                else:
                    self.values[column]['wrong'] += 1
            else:
                if 'insertion' in column:
                    trueval = len(trueval)
                    guessval = len(guessval)
                else:
                    trueval = int(trueval)
                    guessval = int(guessval)
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
                hist.GetXaxis().SetTitle('inferred - true')
                plotting.draw(hist, 'int', plotname=column, plotdir=self.plotdir)
        check_call(['./permissify-www', self.plotdir])  # NOTE this should really permissify starting a few directories higher up
        check_call(['makeHtml', self.plotdir, '3', 'null', 'svg'])
