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
        self.all_headers = combo_headers + list(sorted(set([h for hp in pair_headers for h in hp if h not in combo_headers])))
        self.hpairs = list(itertools.combinations(combo_headers, 2)) + pair_headers  # pairs with some chance of being correlated, i.e. for which it's worth making a pairwise plot
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
        # ----------------------------------------------------------------------------------------
        def get_corr(hk, xvals, yvals):
            xbins, ybins = [sorted(set(tl)) for tl in (xvals, yvals)]  # duplicates code in plotting.plot_smatrix()
            xlabels, ylabels = [xbins.index(v) for v in xvals], [ybins.index(v) for v in yvals]
            # pcorr = numpy.corrcoef(xlabels, ylabels)[0, 1]
            from sklearn import metrics
            # plain_mi = metrics.cluster.mutual_info_score(xlabels, ylabels)
            norm_mi = metrics.cluster.normalized_mutual_info_score(xlabels, ylabels)
            # adj_mi = metrics.cluster.adjusted_mutual_info_score(xlabels, ylabels)
            # print '     %6.2f  %6.2f  %6.2f  %s' % (plain_mi, norm_mi, adj_mi, ' '.join(hk))
            return norm_mi
        # ----------------------------------------------------------------------------------------
        if not only_csv:
            import plotting
            import plotconfig

        corr_vals = [[None for _ in self.all_headers] for _ in self.all_headers]
        for hk in self.cvecs:
            xvals, yvals = zip(*self.cvecs[hk])
            if len(set(xvals)) > 1 and len(set(yvals)) > 1:
                iheads = [self.all_headers.index(h) for h in hk]
                corr_vals[iheads[0]][iheads[1]] = get_corr(hk, xvals, yvals)
                corr_vals[iheads[1]][iheads[0]] = corr_vals[iheads[0]][iheads[1]]  # just to make it symmetric, which it isn't because of the hand-ordering in the headers above (it might be nicer to only fill in one half, but whatever)

            if not only_csv:
                lfcn, y_lfcn = [utils.shorten_gene_name if '_gene' in h else str for h in hk]
                plotting.plot_smatrix(plotdir, '-vs-'.join(hk), xylists=(xvals, yvals), xlabel=plotconfig.plot_titles.get(hk[0], hk[0]), ylabel=plotconfig.plot_titles.get(hk[1], hk[1]), lfcn=lfcn, y_lfcn=y_lfcn, tdbg=debug)

        if not only_csv:
            plotting.plot_smatrix(plotdir, 'mutual-infos', smatrix=corr_vals, xybins=[[plotconfig.plot_titles.get(h, h) for h in self.all_headers] for _ in range(2)], float_vals=True, blabel='norm. mutual info.', tdbg=debug)

            n_per_row = 3
            fnames = ['%s.svg' % '-vs-'.join(hk) for hk in self.hpairs]
            fnames = [['mutual-infos.svg']] + [fnames[i:i+n_per_row] for i in range(0, len(fnames), n_per_row)]
            plotting.make_html(plotdir, fnames=fnames)
