import numpy
import itertools
import operator
import time
import sys
import os

import utils
import glutils
import indelutils

# ----------------------------------------------------------------------------------------
class AlleleClusterer(object):
    # ----------------------------------------------------------------------------------------
    def __init__(self, args):
        self.args = args
        self.region = 'v'
        self.other_region = 'j'
        self.absolute_n_seqs_min = 15
        self.min_cluster_fraction = 0.005
        self.max_j_mutations = 8
        self.small_number_of_j_mutations = 3
        self.all_j_mutations = None

    # ----------------------------------------------------------------------------------------
    def get_alleles(self, queryfo, threshold, glfo, swfo=None, reco_info=None, simglfo=None, debug=False):
        debug = True
        if debug:
            print 'clustering for new alleles'

        default_initial_glfo = glfo
        if self.args.default_initial_germline_dir is not None:  # if this is set, we want to take any new allele names from this directory's glfo if they're in there
            default_initial_glfo = glutils.read_glfo(self.args.default_initial_germline_dir, glfo['locus'])
        else:
            print '  %s --default-initial-germline-dir isn\'t set, so new allele names won\'t correspond to existing names' % utils.color('yellow', 'warning')

        # NOTE do *not* modify <glfo> (in the future it would be nice to just modify <glfo>, but for now we need it to be super clear in partitiondriver what is happening to <glfo>)
        if swfo is None:
            print '  note: not collapsing clones, since we\'re working from vsearch v-only info'
            qr_seqs = {q : qfo['qr_seq'] for q, qfo in queryfo.items()}
            gene_info = {q : qfo['gene'] for q, qfo in queryfo.items()}
        else:
            assert queryfo is None
            clusters = utils.collapse_naive_seqs(swfo)
            qr_seqs = {}
            self.all_j_mutations = {}
            for cluster in clusters:  # take the sequence with the lowest j mutation for each cluster, if it doesn't have too many j mutations NOTE choose_cluster_representatives() in allelefinder is somewhat similar
                clusterstr = ':'.join(cluster)
                j_mutations = {q : utils.get_n_muted(swfo[q], iseq=0, restrict_to_region=self.other_region) for q in cluster}
                best_query, smallest_j_mutations = sorted(j_mutations.items(), key=operator.itemgetter(1))[0]
                if smallest_j_mutations < self.max_j_mutations:
                    qr_seqs[best_query] = indelutils.get_qr_seqs_with_indels_reinstated(swfo[best_query], iseq=0)[self.region]
                for query in cluster:
                    self.all_j_mutations[query] = j_mutations[query]  # I don't think I can key by the cluster str, since here things correspond to the naive-seq-collapsed clusters, then we remove some of the clusters, and then cluster with vsearch
            print '   collapsed %d input sequences into %d representatives from %d clones (removed %d clones with >= %d j mutations)' % (len(swfo['queries']), len(qr_seqs), len(clusters), len(clusters) - len(qr_seqs), self.max_j_mutations)
            gene_info = {q : swfo[q][self.region + '_gene'] for q in qr_seqs}

        msa_fname = self.args.workdir + '/msa.fa'
        print '   vsearch clustering %d %s segments with threshold %.2f (*300 = %d)' % (len(qr_seqs), self.region, threshold, int(threshold * 300))
        assert self.region == 'v'  # would need to change the 300
        _ = utils.run_vsearch('cluster', qr_seqs, self.args.workdir + '/vsearch', threshold=threshold, msa_fname=msa_fname)

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

        n_initial_clusters = len(msa_info)
        print '     read %d vsearch clusters (%d sequences))' % (n_initial_clusters, sum([len(cfo['seqfos']) for cfo in msa_info]))

        n_seqs_min = max(self.absolute_n_seqs_min, self.min_cluster_fraction * len(qr_seqs))
        msa_info = [cfo for cfo in msa_info if len(cfo['seqfos']) >= n_seqs_min]
        print '     removed %d clusters with fewer than %d sequences' % (n_initial_clusters - len(msa_info), n_seqs_min)
        msa_info = sorted(msa_info, key=lambda cfo: len(cfo['seqfos']), reverse=True)

        n_existing_gene_clusters = 0
        # n_clusters_for_new_alleles = 0
        new_alleles = {}
        if debug:
            print '  looping over %d clusters with %d sequences' % (len(msa_info), sum([len(cfo['seqfos']) for cfo in msa_info]))
            print '   seqs   %s mutations (mean)' % self.other_region
        for clusterfo in msa_info:
            assert len(clusterfo['seqfos']) >= n_seqs_min

            # dot_products = [utils.dot_product(clusterfo['cons_seq'], seq1, seq2) for seq1, seq2 in itertools.combinations([seqfo['seq'] for seqfo in clusterfo['seqfos']], 2)]
            # mean_dot_product = numpy.average(dot_products)

            glcounts, true_glcounts = {}, {}
            for seqfo in clusterfo['seqfos']:
                gene = gene_info[seqfo['name']]
                if gene not in glcounts:
                    glcounts[gene] = 0
                glcounts[gene] += 1
                if reco_info is not None:
                    gene = reco_info[seqfo['name']][self.region + '_gene']
                    if gene not in true_glcounts:
                        true_glcounts[gene] = 0
                    true_glcounts[gene] += 1
            sorted_glcounts = sorted(glcounts.items(), key=operator.itemgetter(1), reverse=True)
            if reco_info is not None:
                true_sorted_glcounts = sorted(true_glcounts.items(), key=operator.itemgetter(1), reverse=True)

            mean_j_mutations = numpy.mean([self.all_j_mutations[seqfo['name']] for seqfo in clusterfo['seqfos']])

            # choose the most common existing gene to use as a template (the most similar gene might be a better choice, but deciding on "most similar" would involve adjudicating between snps and indels, and it shouldn't really matter)
            template_gene, _ = sorted_glcounts[0]
            template_seq = glfo['seqs'][self.region][template_gene]
            template_cpos = utils.cdn_pos(glfo, self.region, template_gene)

            new_seq = clusterfo['cons_seq'].replace('-', '')  # I don't really completely understand the dashes in this sequence, but it seems to be right to just remove 'em

            equiv_name, equiv_seq = glutils.find_equivalent_gene_in_glfo(default_initial_glfo, new_seq, template_cpos)
            if equiv_name is not None:
                new_name = equiv_name
                new_seq = equiv_seq
            else:
                new_name, _ = glutils.choose_new_allele_name(template_gene, new_seq)

            if new_name in glfo['seqs'][self.region]:  # note that this only looks in <glfo>, not in <new_alleles>
                n_existing_gene_clusters += 1
                continue

            if debug:
                print '    %-4d' % len(clusterfo['seqfos']),
                if self.all_j_mutations is not None:
                    print '   %5.1f' % mean_j_mutations,
                print ''
                print '           %-20s  %4s   %s' % ('consensus', '', new_seq)
                for gene, counts in sorted_glcounts:
                    print '           %-20s  %4d   %s' % (utils.color_gene(gene, width=20), counts, utils.color_mutants(new_seq, glfo['seqs'][self.region][gene], print_isnps=True, align=True))
                if reco_info is not None:
                    print '       %s   %.1f' % (utils.color('green', 'true'), numpy.mean([utils.get_n_muted(reco_info[seqfo['name']], iseq=0, restrict_to_region=self.other_region) for seqfo in clusterfo['seqfos']]))
                    print '           %-20s  %4s   %s' % ('consensus', '', new_seq)
                    for gene, counts in true_sorted_glcounts:
                        print '           %-12s  %4d   %s' % (utils.color_gene(gene, width=20), counts, utils.color_mutants(new_seq, simglfo['seqs'][self.region][gene], print_isnps=True, align=True))

            if new_name in new_alleles:  # already added it
                # n_clusters_for_new_alleles += 1
                continue
            assert new_seq not in new_alleles.values()  # if it's the same seq, it should've got the same damn name

            if len(new_seq[:template_cpos]) == len(template_seq[:template_cpos]):
                n_snps = utils.hamming_distance(new_seq[:template_cpos], template_seq[:template_cpos])
                # if n_snps < self.args.n_max_snps and mean_j_mutations > self.small_number_of_j_mutations:
                factor = 1.75
                if n_snps < factor * mean_j_mutations:  # i.e. we keep if it's *further* than factor * <number of j mutations> from the closest existing allele (should presumably rescale by some factor to go from j --> v, but it seems like the factor's near to 1.)
                    if debug:
                        print '      too close (%d snp%s < %.2f = %.2f * %d j mutation%s) to existing glfo gene %s' % (n_snps, utils.plural(n_snps), factor * mean_j_mutations, factor, mean_j_mutations, utils.plural(mean_j_mutations), utils.color_gene(template_gene))
                    continue

            if self.too_close_to_already_added_gene(new_seq, new_alleles, debug=debug):
                continue

            print '  %s allele %s%s' % (utils.color('red', 'new'), utils.color_gene(new_name), ' (exists in default germline dir)' if new_name in default_initial_glfo['seqs'][self.region] else '')
            new_alleles[new_name] = {'template-gene' : template_gene, 'gene' : new_name, 'seq' : new_seq}

        if debug:
            print '  %d / %d clusters consensed to existing genes' % (n_existing_gene_clusters, len(msa_info))

        return new_alleles

    # ----------------------------------------------------------------------------------------
    def too_close_to_already_added_gene(self, new_seq, new_alleles, debug=False):
        for added_name, added_info in new_alleles.items():
            _, isnps = utils.color_mutants(added_info['seq'], new_seq, return_isnps=True, align=True)  # oh man that could be cleaner
            if len(isnps) < self.args.n_max_snps:
                if debug:
                    print '      too close (%d snp%s) to gene we just added %s' % (len(isnps), utils.plural(len(isnps)), utils.color_gene(added_name))
                return True
        return False
