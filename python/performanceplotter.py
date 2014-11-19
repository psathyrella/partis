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
        for region in utils.regions:
            self.values[region + '_hamming_to_true_naive'] = {}
        # for bound in utils.boundaries:
        #     self.counts[bound + '_insertion_content'] = {'A':0, 'C':0, 'G':0, 'T':0}  # base content of each insertion TODO add correlation to previous base
        # self.counts['seq_content'] = {'A':0, 'C':0, 'G':0, 'T':0}
        self.values['mute_freqs'] = {}
    # ----------------------------------------------------------------------------------------
    def hamming_distance_to_true_naive(self, true_line, line, query_name, restrict_to_region='', normalize=False):
        """
        Hamming distance between the inferred naive sequence and the tue naive sequence.
        <restrict_to_region> if set, restrict the comparison to the section of the *true* sequence assigned to the given region.
        NOTE this will not in general correspond to the similarly-assigned region in the inferred naive sequence.
        if <normalize> divide by sequence length
        """
        true_naive_seq = utils.get_full_naive_seq(self.germlines, true_line)
        inferred_naive_seq = utils.get_full_naive_seq(self.germlines, line)

        extra_penalty = 0
        if len(true_naive_seq) > len(inferred_naive_seq):  # probably the method added v left or j right erosions that aren't there
            if line['v_5p_del'] > 0:
                true_naive_seq = true_naive_seq[line['v_5p_del'] : ]
                extra_penalty = line['v_5p_del']  # treat it as if the bases were wrong
            if line['j_3p_del'] > 0:
                true_naive_seq = true_naive_seq[ : -line['j_3p_del']]
                extra_penalty = line['j_3p_del']  # treat it as if the bases were wrong

        if restrict_to_region != '':
            bounds = utils.get_regional_naive_seq_bounds(restrict_to_region, self.germlines, true_line)  # get the bounds of this *true* region
            true_naive_seq = true_naive_seq[bounds[0] : bounds[1]]
            inferred_naive_seq = inferred_naive_seq[bounds[0] : bounds[1]]

        total_distance = utils.hamming(true_naive_seq, inferred_naive_seq)
        total_distance += extra_penalty
        if normalize:
            return int(100 * (float(total_distance) / len(true_naive_seq)))
        else:
            return total_distance

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
    def evaluate(self, true_line, line):
        for column in self.values:
            if column in bool_columns:
                if utils.are_alleles(true_line[column], line[column]):  # NOTE you have to change this above as well!
                # if true_line[column] == line[column]:
                    self.values[column]['right'] += 1
                else:
                    self.values[column]['wrong'] += 1
            else:
                trueval, guessval = 0, 0
                if column[2:] == '_insertion':  # insertion length
                    trueval = len(true_line[column])
                    guessval = len(line[column])
                # elif '_content' in column:  # TODO ARG this is harder than I thought. I need to modify to allow distributions rather than just [guess - true]
                #     seq_to_use = line[column[ : column.find('_', 3)]]  # NOTE has to work for seq_content *and* vd_insertion_content, hence the 3
                #         for nuke in seq_to_use:
                #             self.counts[col][nuke] += 1
                elif column == 'hamming_to_true_naive':
                    trueval = 0  # NOTE this is a kind of weird way to do it, since diff ends up as really just the guessval, but it's ok for now
                    guessval = self.hamming_distance_to_true_naive(true_line, line, line['unique_id'])
                elif column.find('hamming_to_true_naive') == 2:  # i.e. it's '[vdj]_hamming_to_true_naive'
                    trueval = 0  # NOTE this is a kind of weird way to do it, since diff ends up as really just the guessval, but it's ok for now
                    guessval = self.hamming_distance_to_true_naive(true_line, line, line['unique_id'], restrict_to_region=column[0], normalize=True)
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
                print '  %s\n    correct up to allele: %4d / %-4d = %4.2f' % (column, right, right+wrong, float(right) / (right + wrong))
                hist = plotting.make_bool_hist(right, wrong, self.name + '-' + column)
                plotting.draw(hist, 'bool', plotname=column, plotdir=self.plotdir, write_csv=True)
            else:
                hist = plotting.make_hist(self.values[column], 'int', self.name + '-' + column, normalize=True)
                log = ''
                if column.find('hamming_to_true_naive') >= 0:
                    hist.GetXaxis().SetTitle('hamming distance')
                    if column.find('hamming_to_true_naive') > 0:
                        log = 'y'
                else:
                    hist.GetXaxis().SetTitle('inferred - true')
                plotting.draw(hist, 'int', plotname=column, plotdir=self.plotdir, write_csv=True, log=log)
        check_call(['./permissify-www', self.plotdir])  # NOTE this should really permissify starting a few directories higher up
        check_call(['makeHtml', self.plotdir, '3', 'null', 'svg'])
