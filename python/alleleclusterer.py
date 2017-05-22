import sys
import os

import utils

# ----------------------------------------------------------------------------------------
class AlleleClusterer(object):
    def __init__(self, args):
        self.region = 'v'
        self.args = args

    def get_alleles(self, swfo, glfo):
        naive_seqs = {query : swfo[query]['naive_seq'] for query in swfo['queries']}
        threshold = 3 * swfo['mute-freqs'][self.region]
        consensus_fname = self.args.workdir + '/cluster-consensus-seqs.fa'
        utils.run_vsearch(naive_seqs, self.args.workdir + '/vsearch', threshold, n_procs=self.args.n_procs, consensus_fname=consensus_fname)
        consensus_seqs = utils.read_fastx(consensus_fname)
        consensus_info = {}
        for seqfo in consensus_seqs:
            namefo = {k : v for k, v in [kvstr.split('=') for kvstr in seqfo['name'].rstrip(';').split(';')]}
            consensus_info[namefo['centroid']] = seqfo['seq']
        for centroid, cons_seq in consensus_info.items():
            glseq = glfo['seqs'][self.region][swfo[centroid][self.region + '_gene']]
            min_len = min(len(glseq), len(cons_seq))
            print '%15s  %s' % (centroid, utils.color_gene(swfo[centroid][self.region + '_gene']))
            utils.color_mutants(glseq[:min_len], cons_seq[:min_len], print_result=True, print_isnps=True, extra_str='    ')
        os.remove(consensus_fname)
