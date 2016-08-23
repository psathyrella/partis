import sys
import glob
import time
from subprocess import check_call
import os

import utils
import plotting
from hist import Hist

# ----------------------------------------------------------------------------------------
class AlleleRemover(object):
    def __init__(self, glfo, args):
        self.glfo = glfo
        self.args = args

        self.n_five_prime_positions_to_exclude = 2  # skip positions that are too close to the 5' end of V (depending on sequence method, these can be unreliable. uh, I think?)
        self.n_three_prime_positions_to_exclude = 4  # skip positions that are too close to the 3' end of V (misassigned insertions look like snps)

        self.n_bins = 30
        self.xmin, self.xmax = 0., 1.

        self.counts = {}

        self.finalized = False

    # ----------------------------------------------------------------------------------------
    def increment(self, qinfo, line_info, debug=False):
        for region in ['v', ]:
            best_hfrac = None
            second_closest_gene, second_smallest_hfrac = None, None  # NOTE I'm not actually positive that the best s-w match has the smallest hfrac, but I don't really care, either
            for score, gene in qinfo['matches'][region]:
                germline_seq = self.glfo['seqs'][region][gene][qinfo['glbounds'][gene][0] : qinfo['glbounds'][gene][1]]
                assert len(line_info['seqs']) == 1
                query_seq = line_info['seqs'][0][qinfo['qrbounds'][gene][0] : qinfo['qrbounds'][gene][1]]  # NOTE this means we're in general using different length matches for each gene, but it seems the least bad thing to do, at least a.t.m.
                assert len(germline_seq) == len(query_seq)

                if gene == line_info[region + '_gene']:  # if this is the best match
                    assert germline_seq == line_info[region + '_gl_seq']  # TODO remove this
                    assert query_seq == line_info[region + '_qr_seqs'][0]  # TODO remove this

                germline_seq = germline_seq[self.n_five_prime_positions_to_exclude : -self.n_three_prime_positions_to_exclude]
                query_seq = query_seq[self.n_five_prime_positions_to_exclude : -self.n_three_prime_positions_to_exclude]
                hfrac = utils.hamming_fraction(germline_seq, query_seq)

                # print '    %5d  %3d / %3d = %5.3f %s' % (score, utils.hamming_distance(germline_seq, query_seq), len(germline_seq), hfrac, gene)

                if gene == line_info[region + '_gene']:  # if this is the best match
                    best_hfrac = hfrac
                    continue

                if second_smallest_hfrac is None or hfrac < second_smallest_hfrac:
                    second_smallest_hfrac = hfrac
                    second_closest_gene = gene

            if second_smallest_hfrac is None:
                print '    no other matches for ', line_info['unique_ids'][0]
                if len(qinfo['matches'][region]) > 1:
                       raise Exception(line_info['unique_ids'][0])
                return

            if debug:
                print 'best %.3f   second %.3f' % (best_hfrac, second_smallest_hfrac)

            if gene not in self.counts:
                self.counts[gene] = {'best' : Hist(self.n_bins, self.xmin, self.xmax, title='best'), 'second' : Hist(self.n_bins, self.xmin, self.xmax, title='second')}
            self.counts[gene]['best'].fill(best_hfrac)
            self.counts[gene]['second'].fill(second_smallest_hfrac)

    # ----------------------------------------------------------------------------------------
    def finalize(self, debug=False):
        assert not self.finalized

        self.finalized = True

    # ----------------------------------------------------------------------------------------
    def plot(self, base_plotdir, only_csv=False):
        if not self.finalized:
            self.finalize(debug=debug)

        print '    plotting'
        plotdir = base_plotdir + '/allele-removing'

        utils.prep_dir(plotdir, wildlings=('*.csv', '*.svg'))

        if only_csv:  # not implemented
            print '    only_csv not yet implemented in allelefinder'
            return

        start = time.time()
        for gene, hists in self.counts.items():
            plotting.draw_no_root(hists['best'], plotname=utils.sanitize_name(gene), plotdir=plotdir, more_hists=[hists['second'], ], errors=True, plottitle=gene, xtitle='counts', linewidths=[5, 3])

        plotting.make_html(plotdir)
        check_call(['./bin/permissify-www', plotdir])
        print '      allele removing plot time: %.1f' % (time.time()-start)
