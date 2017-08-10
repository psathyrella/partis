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
    def __init__(self, args, glfo, reco_info=None, simglfo=None):
        self.args = args
        self.glfo = glfo
        self.reco_info = reco_info
        self.simglfo = simglfo

        self.region = 'v'
        self.other_region = 'j'
        self.absolute_n_seqs_min = 15
        # self.min_cluster_fraction = 0.004
        self.max_number_of_clusters = 50  # of order the number of expected "V gene clusters", i.e. counting groups of very nearby genes as one cluster
        self.max_j_mutations = 8
        self.small_number_of_j_mutations = 3
        self.min_n_snps = 5

        self.all_j_mutations = None


    # ----------------------------------------------------------------------------------------
    def get_glcounts(self, clusterfo, gene_info):  # return map from gene : counts for this cluster
        glcounts, true_glcounts = {}, {}
        for seqfo in clusterfo['seqfos']:
            gene = gene_info[seqfo['name']]
            if gene not in glcounts:
                glcounts[gene] = 0
            glcounts[gene] += 1
            if self.reco_info is not None:
                gene = self.reco_info[seqfo['name']][self.region + '_gene']
                if gene not in true_glcounts:
                    true_glcounts[gene] = 0
                true_glcounts[gene] += 1
        sorted_glcounts = sorted(glcounts.items(), key=operator.itemgetter(1), reverse=True)
        true_sorted_glcounts = None
        if self.reco_info is not None:
            true_sorted_glcounts = sorted(true_glcounts.items(), key=operator.itemgetter(1), reverse=True)
        return sorted_glcounts, true_sorted_glcounts

    # ----------------------------------------------------------------------------------------
    def print_cluster(self, clusterfo, sorted_glcounts, new_seq, mean_j_mutations, true_sorted_glcounts):
        if self.all_j_mutations is not None:
            print '   %5.1f' % mean_j_mutations,
        print ''
        print '           %-20s  %4s   %s' % ('consensus', '', new_seq)
        for gene, counts in sorted_glcounts:
            print '           %-20s  %4d   %s' % (utils.color_gene(gene, width=20), counts, utils.color_mutants(new_seq, self.glfo['seqs'][self.region][gene], print_isnps=True, align=True))
        if self.reco_info is not None:
            print '       %s   %.1f' % (utils.color('green', 'true'), numpy.mean([utils.get_n_muted(self.reco_info[seqfo['name']], iseq=0, restrict_to_region=self.other_region) for seqfo in clusterfo['seqfos']]))
            print '           %-20s  %4s   %s' % ('consensus', '', new_seq)
            for gene, counts in true_sorted_glcounts:
                print '           %-12s  %4d   %s' % (utils.color_gene(gene, width=20), counts, utils.color_mutants(new_seq, self.simglfo['seqs'][self.region][gene], print_isnps=True, align=True))

    # ----------------------------------------------------------------------------------------
    def decide_whether_to_remove_template_genes(self, msa_info, gene_info, new_alleles, debug=True):
        if len(new_alleles) == 0:
            return

        if debug:
            print '              template  new'
            print '      size      snps    snps    assigned',
            if self.reco_info is not None:
                print '         true',
            print ''

        dbg_print = debug  # don't print all the tiny clusters
        templates = {newfo['template-gene'] : newfo['gene'] for newfo in new_alleles.values()}
        all_glcounts = {}
        for clusterfo in sorted(msa_info, key=lambda cfo: len(cfo['seqfos']), reverse=True):
            sorted_glcounts, true_sorted_glcounts = self.get_glcounts(clusterfo, gene_info)  # it would be nice to not re-call this for the clusters we already called it on above
            for gene, counts in sorted_glcounts:  # <gene> is the one assigned by sw before allele clustering
                if debug and len(clusterfo['seqfos']) < 5:
                    if dbg_print:
                        print '     not printing clusters smaller than 5'
                    dbg_print = False

                if gene not in all_glcounts:  # add it before we decide whether to switch it, so a template gene with zero counts will be in there with zero counts
                    all_glcounts[gene] = 0
                if gene in templates:  # if this was a template for a new allele, we have to decide whether to apportion some or all of the sequences in this cluster to that new allele
                    template_gene = gene
                    template_cpos = utils.cdn_pos(self.glfo, self.region, template_gene)
                    cons_seq = clusterfo['cons_seq']
                    template_seq = self.glfo['seqs'][self.region][template_gene]
                    new_allele_seq = new_alleles[templates[template_gene]]['seq']

                    compare_len = min([template_cpos, len(cons_seq), len(template_seq), len(new_allele_seq)])  # NOTE this doesn't account for indels, i.e. the template and consensus sequences are in general different lengths, but that's ok, it'll just inflate the hamming distance for sequences that differ from consensus by indels, and all we care is finding the one that doesn't have any indels
                    n_template_snps = utils.hamming_distance(cons_seq[:compare_len], template_seq[:compare_len])
                    n_new_snps = utils.hamming_distance(cons_seq[:compare_len], new_allele_seq[:compare_len])

                    if debug and dbg_print:
                        print '    %5d       %2d      %2d' % (len(clusterfo['seqfos']), n_template_snps, n_new_snps),

                    if n_new_snps < n_template_snps:  # reassign to the new allele
                        gene = templates[template_gene]
                        if gene not in all_glcounts:  # add it before we decide whether to switch it, so a template gene with zero counts will be in there with zero counts
                            all_glcounts[gene] = 0

                    if debug and dbg_print:
                        print '    %s' % utils.color_gene(gene, width=15),
                        if self.reco_info is not None:
                            true_gene = true_sorted_glcounts[0][0]  # NOTE this is the most *common* simulated gene in the cluster, not necessarily the one corresponding to these particular sequences... but clusters with new alleles should generally be dominated by one gene, so oh, well
                            if true_gene == gene:
                                print '    %s' % utils.color('green', 'ok'),
                            else:
                                print '    %s' % utils.color_gene(true_gene, width=15),
                        print ''

                all_glcounts[gene] += counts

        if debug:
            print '  final counts:'
            for gene, counts in sorted(all_glcounts.items(), key=operator.itemgetter(1), reverse=True):
                print '    %4d  %s' % (counts, utils.color_gene(gene))

        for new_name, newfo in new_alleles.items():
            # print '%s  %s  %.1f / %.1f = %.4f' % (new_name, newfo['template-gene'], all_glcounts[newfo['template-gene']], float(sum(all_glcounts.values())), all_glcounts[newfo['template-gene']] / float(sum(all_glcounts.values())))
            if all_glcounts[newfo['template-gene']] / float(sum(all_glcounts.values())) < self.args.min_allele_prevalence_fraction:  # NOTE all_glcounts only includes large clusters, and the constituents of those clusters are clonal representatives, so this isn't quite the same as in alleleremover
                newfo['remove-template-gene'] = True

    # ----------------------------------------------------------------------------------------
    def get_alleles(self, queryfo, swfo=None, debug=False):
        print 'clustering for new alleles'

        default_initial_glfo = self.glfo
        if self.args.default_initial_germline_dir is not None:  # if this is set, we want to take any new allele names from this directory's glfo if they're in there
            default_initial_glfo = glutils.read_glfo(self.args.default_initial_germline_dir, self.glfo['locus'])
        else:
            print '  %s --default-initial-germline-dir isn\'t set, so new allele names won\'t correspond to existing names' % utils.color('yellow', 'warning')

        # NOTE do *not* modify <self.glfo> (in the future it would be nice to just modify <self.glfo>, but for now we need it to be super clear in partitiondriver what is happening to <self.glfo>)
        if swfo is None:
            assert False  # need to figure out a threshold
            # threshold=self.sw_info['mute-freqs']['v'],
            print '  note: not collapsing clones, and not removing non-full-length sequences, since we\'re working from vsearch v-only info'
            qr_seqs = {q : qfo['qr_seq'] for q, qfo in queryfo.items()}
            gene_info = {q : qfo['gene'] for q, qfo in queryfo.items()}
        else:
            assert queryfo is None

            # first remove non-full-length sequences
            full_length_queries = [q for q in swfo['queries'] if swfo[q]['v_5p_del'] == 0 and swfo[q]['j_3p_del'] == 0]
            print '   removing %d/%d sequences with v_5p or j_3p deletions' % (len(swfo['queries']) - len(full_length_queries), len(swfo['queries']))
            if len(full_length_queries) == 0:
                return {}

            # then cluster by full-length (v+d+j) naive sequence
            clusters = utils.collapse_naive_seqs(swfo, queries=full_length_queries)

            # then build <qr_seqs> from the v sequences corresponding to the least-j-mutated sequence in each of these clusters (skipping clusterings that are too mutated)
            qr_seqs = {}
            self.all_j_mutations = {}
            for cluster in clusters:
                clusterstr = ':'.join(cluster)
                j_mutations = {q : utils.get_n_muted(swfo[q], iseq=0, restrict_to_region=self.other_region) for q in cluster}
                best_query, smallest_j_mutations = sorted(j_mutations.items(), key=operator.itemgetter(1))[0]  # take the sequence with the lowest j mutation for each cluster, if it doesn't have too many j mutations NOTE choose_cluster_representatives() in allelefinder is somewhat similar
                if smallest_j_mutations < self.max_j_mutations:
                    qr_seqs[best_query] = indelutils.get_qr_seqs_with_indels_reinstated(swfo[best_query], iseq=0)[self.region]
                for query in cluster:
                    self.all_j_mutations[query] = j_mutations[query]  # I don't think I can key by the cluster str, since here things correspond to the naive-seq-collapsed clusters, then we remove some of the clusters, and then cluster with vsearch
            print '   collapsed %d input sequences into %d representatives from %d clones (removed %d clones with >= %d j mutations)' % (len(full_length_queries), len(qr_seqs), len(clusters), len(clusters) - len(qr_seqs), self.max_j_mutations)
            gene_info = {q : swfo[q][self.region + '_gene'] for q in qr_seqs}  # assigned gene for the clonal representative from each cluster that we used (i.e. *not* from every sequence in the sample)

            assert self.region == 'v'  # need to think about whether this should always be j, or if it should be self.other_region
            j_mfreqs = [utils.get_mutation_rate(swfo[q], iseq=0, restrict_to_region='j') for q in qr_seqs]
            threshold = numpy.mean(j_mfreqs) / 1.5  # v mut freq will be way off for any very different new alleles

        # then vsearch cluster the v-sequences in <qr_seqs> using a heuristic j-mutation-based threshold
        msa_fname = self.args.workdir + '/msa.fa'
        print '   vsearch clustering %d %s segments with threshold %.2f (*300 = %d)' % (len(qr_seqs), self.region, threshold, int(threshold * 300))
        assert self.region == 'v'  # would need to change the 300
        _ = utils.run_vsearch('cluster', qr_seqs, self.args.workdir + '/vsearch', threshold=threshold, msa_fname=msa_fname, vsearch_binary=self.args.vsearch_binary)
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

        # then throw out smaller clusters
        # n_seqs_min = max(self.absolute_n_seqs_min, self.min_cluster_fraction * len(msa_info))
        n_seqs_min = self.absolute_n_seqs_min
        clusterfos = [cfo for cfo in msa_info if len(cfo['seqfos']) >= n_seqs_min]
        print '     removed %d clusters with fewer than %d sequences' % (n_initial_clusters - len(clusterfos), n_seqs_min)
        clusterfos = sorted(clusterfos, key=lambda cfo: len(cfo['seqfos']), reverse=True)
        if len(clusterfos) > self.max_number_of_clusters:
            print '     taking the %d largest clusters (removing %d)' % (self.max_number_of_clusters, len(clusterfos) - self.max_number_of_clusters)
            clusterfos = clusterfos[:self.max_number_of_clusters]

        # and finally loop over eah cluster, deciding if it corresponds to a new allele
        if debug:
            print '  looping over %d clusters with %d sequences' % (len(clusterfos), sum([len(cfo['seqfos']) for cfo in clusterfos]))
            print '   rank  seqs   %s mutations (mean)' % self.other_region
        new_alleles = {}
        n_existing_gene_clusters = 0
        for iclust in range(len(clusterfos)):
            clusterfo = clusterfos[iclust]

            # dot_products = [utils.dot_product(clusterfo['cons_seq'], seq1, seq2) for seq1, seq2 in itertools.combinations([seqfo['seq'] for seqfo in clusterfo['seqfos']], 2)]
            # mean_dot_product = numpy.average(dot_products)

            # choose the most common existing gene to use as a template (the most similar gene might be a better choice, but deciding on "most similar" would involve adjudicating between snps and indels, and it shouldn't really matter)
            sorted_glcounts, true_sorted_glcounts = self.get_glcounts(clusterfo, gene_info)
            template_gene, template_counts = sorted_glcounts[0]
            template_seq = self.glfo['seqs'][self.region][template_gene]
            template_cpos = utils.cdn_pos(self.glfo, self.region, template_gene)
            mean_j_mutations = numpy.mean([self.all_j_mutations[seqfo['name']] for seqfo in clusterfo['seqfos']])

            new_seq = clusterfo['cons_seq'].replace('-', '')  # I'm not sure that I completely understand the dashes in this sequence, but it seems to be right to just remove 'em

            equiv_name, equiv_seq = glutils.find_equivalent_gene_in_glfo(default_initial_glfo, new_seq, template_cpos)
            if equiv_name is not None:
                new_name = equiv_name
                new_seq = equiv_seq
            else:
                new_name, _ = glutils.choose_new_allele_name(template_gene, new_seq)  # TODO it would be nice to pass in indel info, so the name matches the simulation name when there's indels

            if debug:
                print '    %-3d  %4d' % (iclust, len(clusterfo['seqfos'])),

            if new_name in self.glfo['seqs'][self.region]:  # note that this only looks in <self.glfo>, not in <new_alleles>
                n_existing_gene_clusters += 1
                if debug:
                    print '    %s' % utils.color_gene(new_name)
                continue

            if new_name in new_alleles:  # already added it
                if debug:
                    print '    %s (%s)' % (utils.color_gene(new_name), utils.color('red', 'new'))
                continue
            assert new_seq not in new_alleles.values()  # if it's the same seq, it should've got the same damn name

            if debug:
                self.print_cluster(clusterfo, sorted_glcounts, new_seq, mean_j_mutations, true_sorted_glcounts)

            if len(new_seq[:template_cpos]) == len(template_seq[:template_cpos]):
                n_snps = utils.hamming_distance(new_seq[:template_cpos], template_seq[:template_cpos])  # TODO should probably update this to do the same thing (with min([])) as up in decide_whether_to_remove_template_genes()
                # if n_snps < self.args.n_max_snps and mean_j_mutations > self.small_number_of_j_mutations:
                factor = 1.75
                if n_snps < self.min_n_snps or n_snps < factor * mean_j_mutations:  # i.e. we keep if it's *further* than factor * <number of j mutations> from the closest existing allele (should presumably rescale by some factor to go from j --> v, but it seems like the factor's near to 1.)
                    if debug:
                        print '      too close (%d snp%s < %.2f = %.2f * %.1f mean j mutation%s) to existing glfo gene %s' % (n_snps, utils.plural(n_snps), factor * mean_j_mutations, factor, mean_j_mutations, utils.plural(mean_j_mutations), utils.color_gene(template_gene))
                    continue

            if self.too_close_to_already_added_gene(new_seq, new_alleles, debug=debug):
                continue

            print '       %s %s%s' % (utils.color('red', 'new'), utils.color_gene(new_name), ' (exists in default germline dir)' if new_name in default_initial_glfo['seqs'][self.region] else '')
            new_alleles[new_name] = {'template-gene' : template_gene, 'gene' : new_name, 'seq' : new_seq}

        if debug:
            print '  %d / %d clusters consensed to existing genes' % (n_existing_gene_clusters, len(msa_info))

        self.decide_whether_to_remove_template_genes(msa_info, gene_info, new_alleles)

        return new_alleles

    # ----------------------------------------------------------------------------------------
    def too_close_to_already_added_gene(self, new_seq, new_alleles, debug=False):
        for added_name, added_info in new_alleles.items():
            _, isnps = utils.color_mutants(added_info['seq'], new_seq, return_isnps=True, align=True)  # oh man that could be cleaner
            if len(isnps) < self.min_n_snps or len(isnps) < self.args.n_max_snps:
                if debug:
                    print '      too close (%d snp%s) to gene we just added %s' % (len(isnps), utils.plural(len(isnps)), utils.color_gene(added_name))
                return True
        return False
