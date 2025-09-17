from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import copy
import numpy
import itertools
import operator
import time
import sys
import os

from . import utils
from . import glutils
from . import indelutils
from .hist import Hist
from . import mds

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
        self.max_mutations = {
            'j': 8,
            # 'v' : 25,  # only used for mfreq hists at the moment
        }
        self.min_n_snps = 5
        self.mfreq_ratio_threshold = 1.25
# ----------------------------------------------------------------------------------------
        self.n_mds_components = 3
        self.XXX_n_kmeans_clusters = 3
# ----------------------------------------------------------------------------------------

        # it's not a super sensible distinction, but for now <qr_seqs> is passed around (i.e. not a member variable) since it's what's actually looped over for clusters, and everything else is a member variable
        self.all_j_mutations = None
        self.gene_info = None
        self.mfreqs = None
        self.mean_mfreqs = None
        # self.mfreq_hists = None  # TODO remove this if you're not going to use it (or maybe just print the results of the donut test) hopefully don't need it after switch to kmeans
        self.adjusted_glcounts = None  # gene assignment counts after reassignment for new + template genes

        assert len(utils.gap_chars) == 2  # i just want to hard code it below, dammit
        assert '-' in utils.gap_chars
        assert '.' in utils.gap_chars

    # ----------------------------------------------------------------------------------------
    def get_glcounts(self, clusterfo):  # return map from gene : counts for this cluster
        glcounts, true_glcounts = {}, {}
        for seqfo in clusterfo['seqfos']:
            gene = self.gene_info[seqfo['name']]
            if gene not in glcounts:
                glcounts[gene] = 0
            glcounts[gene] += 1
            if self.reco_info is not None:
                gene = self.reco_info[seqfo['name']][self.region + '_gene']
                if gene not in true_glcounts:
                    true_glcounts[gene] = 0
                true_glcounts[gene] += 1
        sorted_glcounts = sorted(list(glcounts.items()), key=operator.itemgetter(1), reverse=True)
        true_sorted_glcounts = None
        if self.reco_info is not None:
            true_sorted_glcounts = sorted(list(true_glcounts.items()), key=operator.itemgetter(1), reverse=True)
        return sorted_glcounts, true_sorted_glcounts

    # ----------------------------------------------------------------------------------------
    def choose_clonal_representatives(self, swfo, debug=False):
        # NOTE do *not* modify <self.glfo> (in the future it would be nice to just modify <self.glfo>, but for now we need it to be super clear in partitiondriver what is happening to <self.glfo>)

        # first remove non-full-length sequences
        full_length_queries = [q for q in swfo['queries'] if swfo[q]['v_5p_del'] == 0 and swfo[q]['j_3p_del'] == 0]
        print('   removing %d/%d sequences with v_5p or j_3p deletions' % (len(swfo['queries']) - len(full_length_queries), len(swfo['queries'])))
        if len(full_length_queries) == 0:
            return None, None

        # then cluster by full-length (v+d+j) naive sequence
        clusters = utils.collapse_naive_seqs(swfo, queries=full_length_queries)

        # then build <qr_seqs> from the v sequences corresponding to the least-j-mutated sequence in each of these clusters (skipping clusterings that are too mutated)
        qr_seqs = {}
        self.all_j_mutations = {}
        for cluster in clusters:
            clusterstr = ':'.join(cluster)
            j_mutations = {q : utils.get_n_muted(swfo[q], iseq=0, restrict_to_region=self.other_region) for q in cluster}
            best_query, smallest_j_mutations = sorted(list(j_mutations.items()), key=operator.itemgetter(1))[0]  # take the sequence with the lowest j mutation for each cluster, if it doesn't have too many j mutations NOTE choose_cluster_representatives() in allelefinder is somewhat similar
            if smallest_j_mutations < self.max_mutations['j']:
                qr_seqs[best_query] = indelutils.get_qr_seqs_with_indels_reinstated(swfo[best_query], iseq=0)[self.region]
            for query in cluster:
                self.all_j_mutations[query] = j_mutations[query]  # I don't think I can key by the cluster str, since here things correspond to the naive-seq-collapsed clusters, then we remove some of the clusters, and then cluster with vsearch
        print('   collapsed %d input sequences into %d representatives from %d clones (removed %d clones with >= %d j mutations)' % (len(full_length_queries), len(qr_seqs), len(clusters), len(clusters) - len(qr_seqs), self.max_mutations['j']))
        if len(qr_seqs) == 0:
            return None, None

        self.gene_info = {q : swfo[q][self.region + '_gene'] for q in qr_seqs}  # assigned gene for the clonal representative from each cluster that we used (i.e. *not* from every sequence in the sample)
        self.mfreqs = {  # NOTE only includes cluster representatives, i.e. it's biased towards sequences with low overall mutation, and low j mutation
            'v' : {q : utils.get_mutation_rate(swfo[q], iseq=0, restrict_to_region='v') for q in qr_seqs},
            'j' : {q : utils.get_mutation_rate(swfo[q], iseq=0, restrict_to_region='j') for q in qr_seqs},
        }
        self.mean_mfreqs = {r : numpy.mean(list(self.mfreqs[r].values())) for r in self.mfreqs}
        # assert self.region == 'v'  # this won't work if our region is j, since it's too short; there's always/often a dip/gap between 0 mutations and the rest of the distribution
        # self.mfreq_hists = {self.region : Hist(30, 0., 0.3)}  # not reall sure whether it's better to use n_mutes or mfreq, but I already have mfreq
        # for query in qr_seqs:
        #     for region in self.mfreq_hists:
        #         self.mfreq_hists[region].fill(self.mfreqs[region][query])
        print('    mutation among all cluster representatives:   v / j = %6.3f / %6.3f = %6.3f' % (self.mean_mfreqs['v'], self.mean_mfreqs['j'], self.mean_mfreqs['v'] / self.mean_mfreqs['j']))

        assert self.region == 'v'  # need to think about whether this should always be j, or if it should be self.other_region
        j_mfreqs = [utils.get_mutation_rate(swfo[q], iseq=0, restrict_to_region='j') for q in qr_seqs]
        threshold = numpy.mean(j_mfreqs) / 1.5  # v mut freq will be way off for any very different new alleles

        return qr_seqs, threshold

    # ----------------------------------------------------------------------------------------
    def vsearch_cluster_v_seqs(self, qr_seqs, threshold, debug=False):
        # then vsearch cluster the v-sequences in <qr_seqs> using a heuristic j-mutation-based threshold
        msa_fname = self.args.workdir + '/msa.fa'
        print('   vsearch clustering %d %s segments with threshold %.2f (*300 = %d)' % (len(qr_seqs), self.region, threshold, int(threshold * 300)))
        assert self.region == 'v'  # would need to change the 300
        _ = utils.run_vsearch('cluster', qr_seqs, self.args.workdir + '/vsearch', threshold=threshold, msa_fname=msa_fname, vsearch_binary=self.args.vsearch_binary)
        msa_info = []
        msa_seqs = utils.read_fastx(msa_fname)
        for seqfo in msa_seqs:
            if seqfo['name'][0] == '*':  # start of new cluster (centroid is first, and is marked with a '*')
                centroid = seqfo['name'].lstrip('*')
                msa_info.append({'centroid' : centroid, 'seqfos' : [{'name' : centroid, 'seq' : seqfo['seq']}]})  # I don't seem to actually be using the identity of the centroid sequence for anything
            elif seqfo['name'] == 'consensus':
                msa_info[-1]['cons_seq'] = seqfo['seq'].replace('+', '')  # gaaaaah not sure what the +s mean
            else:
                msa_info[-1]['seqfos'].append(seqfo)
        os.remove(msa_fname)
        n_initial_clusters = len(msa_info)
        print('     read %d vsearch clusters (%d sequences))' % (n_initial_clusters, sum([len(cfo['seqfos']) for cfo in msa_info])))

        # then throw out smaller clusters
        # n_seqs_min = max(self.absolute_n_seqs_min, self.min_cluster_fraction * len(msa_info))
        n_seqs_min = self.absolute_n_seqs_min
        clusterfos = [cfo for cfo in msa_info if len(cfo['seqfos']) >= n_seqs_min]
        print('     removed %d clusters with fewer than %d sequences' % (n_initial_clusters - len(clusterfos), n_seqs_min))
        clusterfos = sorted(clusterfos, key=lambda cfo: len(cfo['seqfos']), reverse=True)
        if len(clusterfos) > self.max_number_of_clusters:
            print('     taking the %d largest clusters (removing %d)' % (self.max_number_of_clusters, len(clusterfos) - self.max_number_of_clusters))
            clusterfos = clusterfos[:self.max_number_of_clusters]

        return clusterfos, msa_info

    # ----------------------------------------------------------------------------------------
    def get_family_groups(self, qr_seqs, swfo):
        family_groups = {}  # this is kinda wasteful to copy all the sequences as well
        for name, seq in qr_seqs.items():
            # fam = utils.gene_family(adjusted_XXX[name][self.region + '_gene'])
            fam = utils.gene_family(swfo[name][self.region + '_gene'])
            if fam not in family_groups:
                family_groups[fam] = []
            family_groups[fam].append({'name' : name, 'seq' : seq})
        return family_groups

    # ----------------------------------------------------------------------------------------
    def get_clusterfos_from_partition(self, partition, all_seqs):
        clusterfos = []
        for cluster in partition:
            cfo = {'seqfos' : [{'name' : uid, 'seq' : all_seqs[uid]} for uid in cluster]}  # note that vsearch clustering also adds 'centroid', but I think it isn't subsequently used
            cfo['cons_seq'] = utils.cons_seq(unaligned_seqfos=cfo['seqfos'])  # NOTE switching to the new fcn from old_bio_cons_seq() without testing
            clusterfos.append(cfo)
        return clusterfos

    # ----------------------------------------------------------------------------------------
    def kmeans_cluster_v_seqs(self, qr_seqs, swfo, plotdir=None, debug=False):
        if plotdir is not None:
            utils.prep_dir(plotdir, wildlings=['*.svg'], subdirs=[d for d in os.listdir(plotdir) if os.path.isdir(plotdir + '/' + d)], rm_subdirs=True)

        clusterfos = []
        if debug:
            print('kmeans clustering')
            print('  seqs    family')
        for family, seqfos in self.get_family_groups(qr_seqs, swfo).items():
            if debug:
                print('  %5d     %s' % (len(seqfos), family))
            partition = mds.run_bios2mds(self.n_mds_components, self.XXX_n_kmeans_clusters, seqfos, self.args.workdir + '/mds', self.args.random_seed, reco_info=self.reco_info, region=self.region, plotdir=plotdir + '/' + family if plotdir is not None else None)
            # partition = mds.run_sklearn_mds(self.n_mds_components, self.XXX_n_kmeans_clusters, seqfos, self.args.random_seed, reco_info=self.reco_info, region=self.region, plotdir=plotdir + '/' + family if plotdir is not None else None)
            clusterfos += self.get_clusterfos_from_partition(partition, qr_seqs)

        clusterfos = sorted(clusterfos, key=lambda c: len(c['seqfos']), reverse=True)
        return clusterfos

    # ----------------------------------------------------------------------------------------
    def print_cluster(self, iclust, clusterfo, sorted_glcounts, new_seq, true_sorted_glcounts, mean_cluster_mfreqs, has_indels):
        if iclust > 0:
            print('')
        print('    %-3d  %4d   %6.3f' % (iclust, len(clusterfo['seqfos']), mean_cluster_mfreqs['v'] / mean_cluster_mfreqs['j'] if mean_cluster_mfreqs['j'] > 0. else 0.), end=' ')
        for igene in range(len(sorted_glcounts)):
            if igene > 0:
                print('%22s' % '', end=' ')
            gene, counts = sorted_glcounts[igene]
            print('   %-s %4d      %2d%s' % (utils.color_gene(gene, width=20), counts, utils.hamming_distance(new_seq, self.glfo['seqs'][self.region][gene], align=True), ' (%s)' % utils.color('blue', 'x') if has_indels else '   '), end=' ')
            if igene < len(sorted_glcounts) - 1 or self.reco_info is not None:
                print('')
        if self.reco_info is not None:
            for igene in range(len(true_sorted_glcounts)):
                gene, counts = true_sorted_glcounts[igene]
                print('%17s       %s %-s %4d %s    %2d   ' % ('', utils.color('green', '['), utils.color_gene(gene[:23], width=20), counts, utils.color('green', ']'), utils.hamming_distance(new_seq, self.simglfo['seqs'][self.region][gene], align=True)), end=' ')
                if igene < len(true_sorted_glcounts) - 1:
                    print('')

        # print '           %-20s  %4s   %s' % ('consensus', '', new_seq)
        # for gene, counts in sorted_glcounts:
        #     print '           %-20s  %4d   %s' % (utils.color_gene(gene, width=20), counts, utils.color_mutants(new_seq, self.glfo['seqs'][self.region][gene], print_isnps=True, align=True))
        # if self.reco_info is not None:
        #     print '       %s   %.1f' % (utils.color('green', 'true'), numpy.mean([utils.get_n_muted(self.reco_info[seqfo['name']], iseq=0, restrict_to_region=self.other_region) for seqfo in clusterfo['seqfos']]))
        #     print '           %-20s  %4s   %s' % ('consensus', '', new_seq)
        #     for gene, counts in true_sorted_glcounts:
        #         print '           %-12s  %4d   %s' % (utils.color_gene(gene, width=20), counts, utils.color_mutants(new_seq, self.simglfo['seqs'][self.region][gene], print_isnps=True, align=True))


    # ----------------------------------------------------------------------------------------
    def reassign_template_counts(self, msa_info, new_alleles, debug=False):
        # XXX need to update family_groups here
        if len(new_alleles) == 0:
            return

        if debug:
            print('              template  new')
            print('      size      snps    snps    assigned', end=' ')
            if self.reco_info is not None:
                print('         true', end=' ')
            print('')

        dbg_print = debug  # don't print all the tiny clusters
        templates = {newfo['template-gene'] : newfo['gene'] for newfo in new_alleles.values()}
        self.adjusted_glcounts = {}
        for clusterfo in sorted(msa_info, key=lambda cfo: len(cfo['seqfos']), reverse=True):
            sorted_glcounts, true_sorted_glcounts = self.get_glcounts(clusterfo)  # it would be nice to not re-call this for the clusters we already called it on above
            for gene, counts in sorted_glcounts:  # <gene> is the one assigned by sw before allele clustering
                if debug and len(clusterfo['seqfos']) < 5:
                    if dbg_print:
                        print('     not printing clusters smaller than 5')
                    dbg_print = False

                if gene not in self.adjusted_glcounts:  # add it before we decide whether to switch it, so a template gene with zero counts will be in there with zero counts
                    self.adjusted_glcounts[gene] = 0
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
                        print('    %5d      %3d     %3d' % (len(clusterfo['seqfos']), n_template_snps, n_new_snps), end=' ')

                    if n_new_snps < n_template_snps:  # reassign to the new allele
                        gene = templates[template_gene]
                        if gene not in self.adjusted_glcounts:  # add it before we decide whether to switch it, so a template gene with zero counts will be in there with zero counts
                            self.adjusted_glcounts[gene] = 0

                    if debug and dbg_print:
                        print('    %s' % utils.color_gene(gene, width=15), end=' ')
                        if self.reco_info is not None:
                            true_gene = true_sorted_glcounts[0][0]  # NOTE this is the most *common* simulated gene in the cluster, not necessarily the one corresponding to these particular sequences... but clusters with new alleles should generally be dominated by one gene, so oh, well
                            if true_gene == gene:
                                print('    %s' % utils.color('green', 'ok'), end=' ')
                            else:
                                print('    %s' % utils.color_gene(true_gene, width=15), end=' ')
                        print('')

                self.adjusted_glcounts[gene] += counts

        if debug:
            print('  final counts:')
            for gene, counts in sorted(list(self.adjusted_glcounts.items()), key=operator.itemgetter(1), reverse=True):
                print('    %4d  %s' % (counts, utils.color_gene(gene)))

    # # ----------------------------------------------------------------------------------------
    # def check_for_donuts(self, debug=False):
    #     hist = self.mfreq_hists['v']
    #     fac = 1.5  # one-half the number of sigmas
    #     print hist

    #     min_freq, min_freq_err, min_ibin = None, None, None
    #     for ibin in range(1, hist.n_bins):  # don't care about underflow or overflow
    #         if min_freq is None or hist.bin_contents[ibin] < min_freq:
    #             min_freq = hist.bin_contents[ibin]
    #             min_freq_err = hist.errors[ibin]
    #             min_ibin = ibin

    #         print '  %2d  %9.3f - %3.1f * %7.2f ?> %9.3f - %3.1f * %7.2f' % (ibin, hist.bin_contents[ibin], fac, hist.errors[ibin], min_freq, fac, min_freq_err),
    #         print '--> %9.3f ?> %9.3f' % (hist.bin_contents[ibin] - fac * hist.errors[ibin], min_freq + fac * min_freq_err),
    #         if hist.bin_contents[ibin] - fac * hist.errors[ibin] > min_freq + fac * min_freq_err:
    #             print 'ack!',
    #         print ''

    # ----------------------------------------------------------------------------------------
    def get_alleles(self, swfo, plotdir=None, debug=False):
        print('clustering for new alleles')

        # NOTE do *not* modify <self.glfo> (in the future it would be nice to just modify <self.glfo>, but for now we need it to be super clear in partitiondriver what is happening to <self.glfo>)
        default_initial_glfo = self.glfo
        if self.args.default_initial_germline_dir is not None:  # if this is set (and it essentially always is, since it has a default value), we want to take any new allele names from this directory's glfo if they're in there
            default_initial_glfo = glutils.read_glfo(self.args.default_initial_germline_dir, self.glfo['locus'])
        else:
            print('  %s --default-initial-germline-dir isn\'t set, so new allele names won\'t correspond to existing names' % utils.color('yellow', 'warning'))
        glfo_to_modify = copy.deepcopy(default_initial_glfo)  # so we can add new genes to it, so we can check for equivalency more easily TODO fix that shit, obviously

        qr_seqs, threshold = self.choose_clonal_representatives(swfo, debug=debug)
        if qr_seqs is None:
            return {}

        # self.check_for_donuts(debug=debug)

        if not self.args.kmeans_allele_cluster:
            clusterfos, msa_info = self.vsearch_cluster_v_seqs(qr_seqs, threshold, debug=debug)
        else:
            clusterfos = self.kmeans_cluster_v_seqs(qr_seqs, swfo, plotdir=plotdir, debug=debug)
            msa_info = clusterfos

        # and finally loop over each cluster, deciding if it corresponds to a new allele
        if debug:
            print('  looping over %d clusters with %d sequences' % (len(clusterfos), sum([len(cfo['seqfos']) for cfo in clusterfos])))
            print('   rank  seqs   v/j mfreq                 seqs      snps (%s)' % utils.color('blue', 'indels'))
        new_alleles = {}
        n_existing_gene_clusters = 0
        for iclust in range(len(clusterfos)):
            clusterfo = clusterfos[iclust]

            # dot_products = [utils.dot_product(clusterfo['cons_seq'], seq1, seq2) for seq1, seq2 in itertools.combinations([seqfo['seq'] for seqfo in clusterfo['seqfos']], 2)]
            # mean_dot_product = numpy.average(dot_products)

            # choose the most common existing gene to use as a template (the most similar gene might be a better choice, but deciding on "most similar" would involve adjudicating between snps and indels, and it shouldn't really matter)
            sorted_glcounts, true_sorted_glcounts = self.get_glcounts(clusterfo)
            template_gene, template_counts = sorted_glcounts[0]
            template_seq = self.glfo['seqs'][self.region][template_gene]
            template_cpos = utils.cdn_pos(self.glfo, self.region, template_gene)

            assert '.' not in clusterfo['cons_seq']  # make sure you haven't switched to something that doesn't use '-' for gap chars
            new_seq = clusterfo['cons_seq'].replace('-', '')  # I'm not sure that I completely understand the dashes in this sequence, but it seems to be right to just remove 'em

            aligned_template_seq, aligned_new_seq = utils.align_seqs(template_seq, clusterfo['cons_seq'])
            has_indels = '-' in aligned_template_seq.strip('-') or '-' in aligned_new_seq.strip('-')  # only counts internal indels
            cluster_mfreqs = {r : [self.mfreqs[r][seqfo['name']] for seqfo in clusterfo['seqfos']] for r in self.mfreqs}  # regional mfreqs for each sequence in the cluster corresponding to the initially-assigned existing gene
            mean_cluster_mfreqs = {r : numpy.mean(cluster_mfreqs[r]) for r in cluster_mfreqs}

            equiv_name, equiv_seq = glutils.find_equivalent_gene_in_glfo(glfo_to_modify, new_seq, template_cpos)
            if equiv_name is not None:
                new_name = equiv_name
                new_seq = equiv_seq
            else:
                new_name, _ = glutils.choose_new_allele_name(template_gene, new_seq, indelfo={'indels' : ['xxx', 'xxx', 'xxx']} if has_indels else None)  # the fcn just checks to see if it's non-None and of length greater than zero...TODO it would be nice to figure out actual snp and indel info

            if debug:
                self.print_cluster(iclust, clusterfo, sorted_glcounts, new_seq, true_sorted_glcounts, mean_cluster_mfreqs, has_indels)

            if new_name in self.glfo['seqs'][self.region]:  # note that this only looks in <self.glfo>, not in <new_alleles>
                n_existing_gene_clusters += 1
                if debug:
                    print('existing %s' % utils.color_gene(new_name))
                continue

            if new_name in new_alleles:  # already added it NOTE might make more sense to use <glfo_to_modify> here instead of <new_alleles> (or just not have freaking both of them)
                if debug:
                    print('%s (%s)' % (utils.color_gene(new_name), utils.color('red', 'new')))
                continue
            assert new_seq not in list(new_alleles.values())  # if it's the same seq, it should've got the same damn name

            if not has_indels:  # we assume that the presence of indels somewhat precludes false positives, which is equivalent to an assumption about the rarity of shm indels
                if self.too_close_to_existing_glfo_gene(clusterfo, new_seq, template_seq, template_cpos, template_gene, debug=debug):  # presumably if it were really close to another (non-template) existing glfo gene, that one would've been the template
                    continue

                if mean_cluster_mfreqs['j'] > 0. and self.mean_mfreqs['j'] > 0.:
                    this_cluster_ratio = mean_cluster_mfreqs['v'] / mean_cluster_mfreqs['j']
                    overall_ratio = self.mean_mfreqs['v'] / self.mean_mfreqs['j']
                    if this_cluster_ratio / overall_ratio < self.mfreq_ratio_threshold:
                        if debug:
                            print('v / j cluster mfreqs too small %6.3f / %6.3f = %6.3f < %6.3f' % (this_cluster_ratio, overall_ratio, this_cluster_ratio / overall_ratio, self.mfreq_ratio_threshold))
                        continue

            if self.too_close_to_already_added_gene(new_seq, new_alleles, debug=debug):  # this needs to be applied even if there are indels, since the indels are with respect to the (existing glfo) template gene, not to the [potentially] previously-added gene
                continue

            print('%s %s%s' % (utils.color('red', 'new'), utils.color_gene(new_name), ' (exists in default germline dir)' if new_name in default_initial_glfo['seqs'][self.region] else ''))
            new_alleles[new_name] = {'template-gene' : template_gene, 'gene' : new_name, 'seq' : new_seq, 'cpos' : template_cpos}
            if new_alleles[new_name]['gene'] not in glfo_to_modify['seqs'][self.region]:  # if it's in <default_initial_glfo> it'll already be in there
                glutils.add_new_allele(glfo_to_modify, new_alleles[new_name], use_template_for_codon_info=False)  # just so we can check for equivalency next time through (we don't really actually do anything else with <glfo_to_modify>)

        if debug:
            print('  %d / %d clusters consensed to existing genes' % (n_existing_gene_clusters, len(msa_info)))

        self.reassign_template_counts(msa_info, new_alleles, debug=False)
        for new_name, newfo in new_alleles.items():
            # print '%s  %s  %.1f / %.1f = %.4f' % (new_name, newfo['template-gene'], self.adjusted_glcounts[newfo['template-gene']], float(sum(self.adjusted_glcounts.values())), self.adjusted_glcounts[newfo['template-gene']] / float(sum(self.adjusted_glcounts.values())))
            if self.adjusted_glcounts[newfo['template-gene']] / float(sum(self.adjusted_glcounts.values())) < self.args.min_allele_prevalence_fractions[self.region]:  # NOTE self.adjusted_glcounts only includes large clusters, and the constituents of those clusters are clonal representatives, so this isn't quite the same as in alleleremover
                newfo['remove-template-gene'] = True

        return new_alleles

    # ----------------------------------------------------------------------------------------
    def too_close_to_existing_glfo_gene(self, clusterfo, new_seq, template_seq, template_cpos, template_gene, debug=False):
        if len(new_seq[:template_cpos]) != len(template_seq[:template_cpos]):  # TODO update this to use the new n_snps from the aligned template/new seqs
            return False

        mean_j_mutations = numpy.mean([self.all_j_mutations[seqfo['name']] for seqfo in clusterfo['seqfos']])  # TODO <self.all_j_mutations> uses everybody in the cluster, rather than just the representative. It'd be nice to synchronize this with other things
        # TODO should probably update this to do the same thing (with min([])) as up in decide_whether_to_remove_template_genes(), or just use the new <align> option to utils.hamming_distance (although that would be slower, and this seems to work ok)
        pre_cpos_snps = utils.hamming_distance(new_seq[:template_cpos], template_seq[:template_cpos])
        factor = 1.75
        if pre_cpos_snps < self.min_n_snps or pre_cpos_snps < factor * mean_j_mutations:  # i.e. we keep if it's *further* than factor * <number of j mutations> from the closest existing allele (should presumably rescale by some factor to go from j --> v, but it seems like the factor's near to 1.)
            if debug:
                print('too close to existing glfo gene %s (%d snp%s < %.2f = %.2f * %.1f mean j mutation%s)' % (utils.color_gene(template_gene), pre_cpos_snps, utils.plural(pre_cpos_snps), factor * mean_j_mutations, factor, mean_j_mutations, utils.plural(mean_j_mutations)))
            return True

        return False

    # ----------------------------------------------------------------------------------------
    def too_close_to_already_added_gene(self, new_seq, new_alleles, debug=False):
        for added_name, added_info in new_alleles.items():
            n_snps = utils.hamming_distance(added_info['seq'], new_seq, align=True)
            if n_snps < self.min_n_snps or n_snps < self.args.n_max_snps:
                if debug:
                    print('too close (%d snp%s) to gene we just added %s' % (n_snps, utils.plural(n_snps), utils.color_gene(added_name)))
                return True
        return False
