import sys

import utils

# ----------------------------------------------------------------------------------------
class AlleleRemover(object):
    def __init__(self, glfo, args):
        self.glfo = glfo
        self.args = args

        self.n_five_prime_positions_to_exclude = 2  # skip positions that are too close to the 5' end of V (depending on sequence method, these can be unreliable. uh, I think?)
        self.n_three_prime_positions_to_exclude = 4  # skip positions that are too close to the 3' end of V (misassigned insertions look like snps)

        self.counts = {}

        self.finalized = False

    # ----------------------------------------------------------------------------------------
    def increment(self, qinfo, line_info):
        for region in ['v', ]:
            best_hfrac = None
            second_closest_gene, second_smallest_hfrac = None, None  # NOTE I'm not actually positive that the best s-w match has the smallest hfrac, but I don't really care, either
            for _, gene in qinfo['matches'][region]:
                assert len(line_info['seqs']) == 1
                query_seq = line_info[region + '_qr_seqs'][0]  # NOTE keep it inside the loop since we shorten it below

                germline_seq = self.glfo['seqs'][region][gene][qinfo['glbounds'][gene][0] : qinfo['glbounds'][gene][1]]
                assert len(germline_seq) == len(query_seq)

                if gene == line_info[region + '_gene']:  # if this is the best match
                    assert germline_seq == line_info[region + '_gl_seq']  # TODO remove this

                germline_seq = germline_seq[self.n_five_prime_positions_to_exclude : -self.n_three_prime_positions_to_exclude]
                query_seq = query_seq[self.n_five_prime_positions_to_exclude : -self.n_three_prime_positions_to_exclude]
                hfrac = utils.hamming_fraction(germline_seq, query_seq)

                print '   ', gene, hfrac

                if gene == line_info[region + '_gene']:  # if this is the best match
                    best_hfrac = hfrac
                    continue

                if second_smallest_hfrac is None or hfrac < second_smallest_hfrac:
                    second_smallest_hfrac = hfrac
                    second_closest_gene = gene

            print 'best %.3f   second %.3f' % (best_hfrac, second_smallest_hfrac)

            sys.exit()

        #     gene = line_info[region + '_gene']
        #     if gene not in self.counts:
        #         self.counts[gene] = {}

        #     germline_seq = line_info[region + '_gl_seq']
        #     assert len(line_info['seqs']) == 1
        #     query_seq = line_info[region + '_qr_seqs'][0]
        #     assert len(germline_seq) == len(query_seq)
        #     germline_seq = germline_seq[self.n_five_prime_positions_to_exclude : -self.n_three_prime_positions_to_exclude]
        #     query_seq = query_seq[self.n_five_prime_positions_to_exclude : -self.n_three_prime_positions_to_exclude]
        #     n_mutes = utils.hamming_distance(germline_seq, query_seq)

    # ----------------------------------------------------------------------------------------
    def finalize(self, debug=False):
        assert not self.finalized


        self.finalized = True

    # ----------------------------------------------------------------------------------------
    def plot(self, base_plotdir, only_csv=False):
        if not self.finalized:
            self.finalize(debug=debug)
