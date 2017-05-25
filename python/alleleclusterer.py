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
        self.absolute_n_seqs_min = 15

    # ----------------------------------------------------------------------------------------
    def get_alleles(self, queryfo, threshold, glfo, swfo=None, debug=True):
        # NOTE do *not* modify <glfo> (in the future it would be nice to just modify <glfo>, but for now we need it to be super clear in partitiondriver what is happening to <glfo>)
        if swfo is None:
            qr_seqs = {q : qfo['qr_seq'] for q, qfo in queryfo.items()}
            gene_info = {q : qfo['gene'] for q, qfo in queryfo.items()}
        else:
            assert queryfo is None
            qr_seqs = {query : utils.get_qr_seqs_with_indels_reinstated(swfo[query], iseq=0)[self.region] for query in swfo['queries']}
            gene_info = {q : swfo[q][self.region + '_gene'] for q in swfo['queries']}

        msa_fname = self.args.workdir + '/msa.fa'
        print '  running vsearch with threshold %.2f (*300 = %d)' % (threshold, int(threshold * 300))
        _ = utils.run_vsearch('cluster', qr_seqs, self.args.workdir + '/vsearch', threshold=threshold, n_procs=self.args.n_procs, msa_fname=msa_fname)

        msa_info = []
        msa_seqs = utils.read_fastx(msa_fname)
        for seqfo in msa_seqs:
            if seqfo['name'][0] == '*':  # start of new cluster (centroid is first, and is marked with a '*')
                centroid = seqfo['name'].lstrip('*')
                msa_info.append({'centroid' : centroid, 'seqfos' : [{'name' : centroid, 'seq' : seqfo['seq']}]})
            elif seqfo['name'] == 'consensus':
                msa_info[-1]['cons_seq'] = seqfo['seq'].replace('+', '')  # gaaaaah not sure what the +s mean
            else:
                msa_info[-1]['seqfos'].append(seqfo)
        os.remove(msa_fname)

        new_alleles = {}
        n_seqs_min = max(self.absolute_n_seqs_min, self.args.min_allele_prevalence_fraction * len(qr_seqs))
        for clusterfo in sorted(msa_info, key=lambda cfo: len(cfo['seqfos']), reverse=True):
            if len(clusterfo['seqfos']) < n_seqs_min:
                break

            # dot_products = [utils.dot_product(clusterfo['cons_seq'], seq1, seq2) for seq1, seq2 in itertools.combinations([seqfo['seq'] for seqfo in clusterfo['seqfos']], 2)]
            # mean_dot_product = numpy.average(dot_products)

            glcounts = {}
            for seqfo in clusterfo['seqfos']:
                gene = gene_info[seqfo['name']]
                if gene not in glcounts:
                    glcounts[gene] = 0
                glcounts[gene] += 1
            sorted_glcounts = sorted(glcounts.items(), key=operator.itemgetter(1), reverse=True)
            if debug:
                print '%d seqs' % len(clusterfo['seqfos'])
                print '    %-12s  %4s   %s' % ('consensus', '', clusterfo['cons_seq'])
                for gene, counts in sorted_glcounts:
                    print '    %-12s  %4d   %s' % (utils.color_gene(gene, width=12), counts, utils.color_mutants(clusterfo['cons_seq'], glfo['seqs'][self.region][gene], print_isnps=True, align=True))

            # choose some existing gene to use as a template (the most similar gene might be a better choice, but deciding on "most similar" would involve adjudicating between snps and indels, and it shouldn't really matter)
            template_gene, _ = sorted_glcounts[0]
            template_seq = glfo['seqs'][self.region][template_gene]  # not necessarily the most similar (but it probably is)
            template_cpos = utils.cpos(glfo, self.region, template_gene)

            new_name = template_gene + '+' + ''.join(numpy.random.choice(list('xyz'), size=8))  # er, not sure what to use here, but at least this probably won't result in collisions
            new_seq = clusterfo['cons_seq'].replace('-', '')
            new_name, new_seq = glutils.find_new_allele_in_existing_glfo(glfo, self.region, new_name, new_seq, template_cpos)

            if new_name in glfo['seqs'][self.region]:
                print '    existing gene %s' % utils.color_gene(new_name)
                continue

            if len(new_seq[:template_cpos]) == len(template_seq[:template_cpos]):
                n_snps = utils.hamming_distance(new_seq[:template_cpos], template_seq[:template_cpos])
                if n_snps < self.args.n_max_snps:
                    print '    too close (%d snp%s) to existing gene %s' % (n_snps, utils.plural(n_snps), utils.color_gene(template_gene))
                    continue

            if self.too_close_to_already_added_gene(new_seq, new_alleles):
                continue

            print '  %s new allele %s' % (utils.color('bold', utils.color('blue', '-->')), utils.color_gene(new_name))
            new_alleles[new_name] = {'template-gene' : template_gene, 'gene' : new_name, 'seq' : new_seq}

        return new_alleles

    # ----------------------------------------------------------------------------------------
    def too_close_to_already_added_gene(self, new_seq, new_alleles):
        for added_name, added_info in new_alleles.items():
            _, isnps = utils.color_mutants(added_info['seq'], new_seq, return_isnps=True, align=True)  # oh man that could be cleaner
            if len(isnps) < self.args.n_max_snps:
                print '    too close (%d snp%s) to gene we just added %s' % (len(isnps), utils.plural(len(isnps)), utils.color_gene(added_name))
                return True
        return False
