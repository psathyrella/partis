import operator
import time
import sys
import os

import utils
import glutils

# ----------------------------------------------------------------------------------------
class AlleleClusterer(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, args):
        self.args = args
        self.region = 'v'
        self.min_seqs = 8

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
        icluster = 0
        for seqfo in consensus_seqs:
            namefo = {k : v for k, v in [kvstr.split('=') for kvstr in seqfo['name'].rstrip(';').split(';')]}
            if namefo['centroid'] not in partition[icluster]:  # <partition> should be in the same order
                raise Exception('ack!')
            consensus_info[namefo['centroid']] = {'size' : int(namefo['seqs']), 'cons_seq' : seqfo['seq'], 'icluster' : icluster}
            icluster += 1

        final_alleles = {}
        for centroid, cinfo in consensus_info.items():
            if cinfo['size'] < self.min_seqs:
                continue

            cluster = partition[cinfo['icluster']]
            glcounts = {}
            for query in cluster:
                gene = swfo[query][self.region + '_gene']
                if gene not in glcounts:
                    glcounts[gene] = 0
                glcounts[gene] += 1
            most_common_gene, _ = sorted(glcounts.items(), key=operator.itemgetter(1), reverse=True)[0]  # not sure what but the most similar gene wouldn't be a better choice, but I think it doesn't matter since it's just getting used for find_new_allele_in_existing_glfo()
            if debug:
                print '%15s: %d seqs' % (centroid, cinfo['size'])
                print '  %-12s  %4s   %s' % ('consensus', '', cinfo['cons_seq'])
                for gene in glcounts:
                    print '  %-12s  %4d   %s' % (utils.color_gene(gene, width=12), glcounts[gene], utils.color_mutants(cinfo['cons_seq'], glfo['seqs'][self.region][gene], print_isnps=True, align=True))

            new_name = most_common_gene + '+XXX'
            new_seq = cinfo['cons_seq']
            template_cpos = glfo[utils.conserved_codons[glfo['locus']][self.region] + '-positions'][most_common_gene]
            new_name, new_seq = glutils.find_new_allele_in_existing_glfo(glfo, self.region, new_name, new_seq, template_cpos)
            if new_name not in glfo['seqs'][self.region]:
                print '  actual new allele %s' % utils.color_gene(new_name)
                final_alleles[new_name] = new_seq

        return final_alleles
