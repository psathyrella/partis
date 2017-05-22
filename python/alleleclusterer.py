import time
import sys
import os

import utils

# ----------------------------------------------------------------------------------------
class AlleleClusterer(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, args):
        self.args = args
        self.region = 'v'
        self.min_seqs = 15

    # ----------------------------------------------------------------------------------------
    def get_alleles(self, swfo, glfo, debug=True):
        qr_seqs = {query : swfo[query][self.region + '_qr_seqs'][0] for query in swfo['queries']}
        consensus_fname = self.args.workdir + '/cluster-consensus-seqs.fa'
        threshold = swfo['mute-freqs'][self.region]
        print '  running vsearch with threshold %.2f (*300 = %d)' % (threshold, int(threshold * 300))
        start = time.time()
        partition = utils.run_vsearch(qr_seqs, self.args.workdir + '/vsearch', threshold=threshold, n_procs=self.args.n_procs, consensus_fname=consensus_fname)
        print '  vsearch: %.1f sec' % (time.time()-start)

        # read vsearch output (and convert consensus file to a more useful format)
        consensus_seqs = utils.read_fastx(consensus_fname)
        os.remove(consensus_fname)
        consensus_info = {}
        for seqfo in consensus_seqs:
            namefo = {k : v for k, v in [kvstr.split('=') for kvstr in seqfo['name'].rstrip(';').split(';')]}
            consensus_info[namefo['centroid']] = {'size' : int(namefo['seqs']), 'cons_seq' : seqfo['seq']}

        #
        final_alleles = {}
        for centroid, cinfo in consensus_info.items():
            if cinfo['size'] < self.min_seqs:
                continue

            clusters = [cl for cl in partition if centroid in cl]
            if len(clusters) != 1:
                raise Exception('didn\'t find one cluster for %s (found %d)' % (centroid, len(clusters)))
            cluster = clusters[0]

            geneset = set([swfo[q][self.region + '_gene'] for q in cluster])
            glseqs = {g : glfo['seqs'][self.region][g] for g in geneset}
            if debug:
                print '%15s: %d seqs' % (centroid, cinfo['size'])
                print '  %-12s   %s' % ('consensus', cinfo['cons_seq'])
                for gene, glseq in glseqs.items():
                    print '  %-12s   %s' % (utils.color_gene(gene, width=12), utils.color_mutants(cinfo['cons_seq'], glseq, print_isnps=True, align=True))

            newname = hash(cinfo['cons_seq'])
            final_alleles[newname] = cinfo['cons_seq']

        return final_alleles
