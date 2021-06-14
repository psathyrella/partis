import sys
import itertools
import csv
import os

from hist import Hist
import utils
import glutils

# ----------------------------------------------------------------------------------------
class CorrCounter(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self):
        combo_headers = ['v_gene', 'd_gene', 'j_gene']  # we want all pairwise combos of these with each other
        pair_headers = [('v_3p_del', 'v_gene'),  # whereas for these we specify each individual pair
                        ('d_5p_del', 'd_gene'), ('d_3p_del', 'd_gene'),
                        ('j_5p_del', 'j_gene'),
                        ('v_gene', 'cdr3_length'), ('d_gene', 'cdr3_length'), ('j_gene', 'cdr3_length'),
                        ('vd_insertion', 'v_gene'), ('vd_insertion', 'd_gene'), ('dj_insertion', 'd_gene'), ('dj_insertion', 'j_gene')]
        self.hpairs = list(itertools.combinations(combo_headers, 2)) + pair_headers
        self.cvecs = {(a, b) : [] for a, b in self.hpairs}

    # ----------------------------------------------------------------------------------------
    def increment(self, info):
        for h_a, h_b in self.hpairs:
            def vfcn(h, v): return len(v) if '_insertion' in h else v
            self.cvecs[(h_a, h_b)].append(tuple([vfcn(h, info[h]) for h in (h_a, h_b)]))

    # ----------------------------------------------------------------------------------------
    def clean_plots(self, plotdir):
        utils.prep_dir(plotdir, wildlings=(('*.csv', '*.svg')))

    # ----------------------------------------------------------------------------------------
    def plot(self, plotdir, only_csv=False, debug=False):
        if only_csv:
            return

        import plotting
        import plotconfig

        for hk in self.cvecs:
            xvals, yvals = zip(*self.cvecs[hk])
            lfcn, y_lfcn = [utils.shorten_gene_name if '_gene' in h else str for h in hk]
            plotting.plot_smatrix(plotdir, '-vs-'.join(hk), xylists=(xvals, yvals), xlabel=plotconfig.plot_titles.get(hk[0], hk[0]), ylabel=plotconfig.plot_titles.get(hk[1], hk[1]), lfcn=lfcn, y_lfcn=y_lfcn, tdbg=debug)

        n_per_row = 3
        fnames = ['%s.svg' % '-vs-'.join(hk) for hk in self.hpairs]
        fnames = [fnames[i:i+n_per_row] for i in range(0, len(fnames), n_per_row)]
        plotting.make_html(plotdir, fnames=fnames)
