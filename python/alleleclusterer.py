import numpy
import itertools
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
        self.absolute_min_seqs = 15

    # ----------------------------------------------------------------------------------------
    def get_alleles(self, swfo, glfo, debug=True):
        qr_seqs = {query : swfo[query][self.region + '_qr_seqs'][0] for query in swfo['queries']}
        consensus_fname = self.args.workdir + '/cluster-consensus-seqs.fa'
        msa_fname = self.args.workdir + '/msa.fa'
        threshold = swfo['mute-freqs'][self.region]
        print '  running vsearch with threshold %.2f (*300 = %d)' % (threshold, int(threshold * 300))
        start = time.time()
        _ = utils.run_vsearch(qr_seqs, self.args.workdir + '/vsearch', threshold=threshold, n_procs=self.args.n_procs, msa_fname=msa_fname)
        print '  vsearch: %.1f sec' % (time.time()-start)

        # consensus_seqs = utils.read_fastx(consensus_fname)
        # os.remove(consensus_fname)
        # consensus_info = {}
        # icluster = 0
        # for seqfo in consensus_seqs:
        #     namefo = {k : v for k, v in [kvstr.split('=') for kvstr in seqfo['name'].rstrip(';').split(';')]}
        #     if namefo['centroid'] not in partition[icluster]:  # <partition> should be in the same order
        #         raise Exception('ack!')
        #     consensus_info[namefo['centroid']] = {'size' : int(namefo['seqs']), 'cons_seq' : seqfo['seq'], 'icluster' : icluster}
        #     icluster += 1

        # could go back to using just the consensus file if I"m not going to calculate the dot product below
        msa_seqs = utils.read_fastx(msa_fname)
        msa_info = []
        for seqfo in msa_seqs:
            if seqfo['name'][0] == '*':  # start of new cluster (centroid is first, and is marked with a '*')
                centroid = seqfo['name'].lstrip('*')
                msa_info.append({'centroid' : centroid, 'seqfos' : [{'name' : centroid, 'seq' : seqfo['seq']}]})
            elif seqfo['name'] == 'consensus':
                msa_info[-1]['cons_seq'] = seqfo['seq']
            else:
                msa_info[-1]['seqfos'].append(seqfo)
        os.remove(msa_fname)

        final_alleles = {}
        min_seqs = max(self.absolute_min_seqs, self.args.min_allele_prevalence_fraction * len(swfo['queries']))
        for clusterfo in msa_info:
            if len(clusterfo['seqfos']) < min_seqs:
                continue

            glcounts = {}
            for seqfo in clusterfo['seqfos']:
                gene = swfo[seqfo['name']][self.region + '_gene']
                if gene not in glcounts:
                    glcounts[gene] = 0
                glcounts[gene] += 1
            # dot_products = [utils.dot_product(clusterfo['cons_seq'], seq1, seq2) for seq1, seq2 in itertools.combinations([seqfo['seq'] for seqfo in clusterfo['seqfos']], 2)]
            # mean_dot_product = numpy.average(dot_products)
            sorted_glcounts = sorted(glcounts.items(), key=operator.itemgetter(1), reverse=True)
            if debug:
                print '%d seqs' % len(clusterfo['seqfos'])
                print '    %-12s  %4s   %s' % ('consensus', '', clusterfo['cons_seq'])
                for gene, counts in sorted_glcounts:
                    print '    %-12s  %4d   %s' % (utils.color_gene(gene, width=12), counts, utils.color_mutants(clusterfo['cons_seq'], glfo['seqs'][self.region][gene], print_isnps=True, align=True))

            most_common_gene, _ = sorted_glcounts[0]  # not sure what but the most similar gene wouldn't be a better choice, but I think it doesn't matter since it's just getting used for find_new_allele_in_existing_glfo()
            new_name = most_common_gene + '+' + str(random.randint(0, 999999))  # er, not sure what to use here, but at least this probably won't result in collisions
            new_seq = clusterfo['cons_seq']
            template_cpos = glfo[utils.conserved_codons[glfo['locus']][self.region] + '-positions'][most_common_gene]
            new_name, new_seq = glutils.find_new_allele_in_existing_glfo(glfo, self.region, new_name, new_seq, template_cpos)
            if new_name in glfo['seqs'][self.region]:
                print '  existing gene %s' % utils.color_gene(new_name)
                continue
            most_common_glseq = glfo['seqs'][self.region][most_common_gene]  # not necessarily the most similar (but it probably is)
            if len(new_seq[:template_cpos]) == len(most_common_glseq[:template_cpos]):
                n_snps = utils.hamming_distance(new_seq[:template_cpos], most_common_glseq[:template_cpos])
                if n_snps < self.args.n_max_snps:
                    print '  too close (%d snps) to existing gene %s' % (n_snps, utils.color_gene(most_common_gene))
                    continue
            print '  actual new allele %s' % utils.color_gene(new_name)
            final_alleles[new_name] = new_seq

        return final_alleles
