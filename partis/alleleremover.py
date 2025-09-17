from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import sys
import glob
import time
import operator
from subprocess import check_call
import copy
import os

from . import utils
from . import glutils
from .hist import Hist

# ----------------------------------------------------------------------------------------
class AlleleRemover(object):
    def __init__(self, glfo, args, simglfo=None, reco_info=None):
        self.glfo = glfo
        self.args = args
        self.simglfo = simglfo
        self.reco_info = reco_info
        if self.reco_info is not None:
            self.simcounts = {r : {g : 0 for g in self.simglfo['seqs'][r]} for r in utils.regions}
            for uid in self.reco_info:
                for tmpreg in utils.regions:
                    self.simcounts[tmpreg][self.reco_info[uid][tmpreg + '_gene']] += 1

        self.genes_to_keep = set()  # NOTE if we're doing several regions, these include genes from all of 'em
        self.genes_to_remove = set()
        self.dbg_strings = {}

        self.finalized = False

        self.alleles_per_group = 2 if self.args.dont_find_new_alleles else 1
        self.n_max_snps = {
            'v' : self.args.n_max_snps - 2,  # not sure why I originally put the 2 there, maybe just to adjust the space between where allelefinder stops looking and where alleleclusterer starts?
            'd' : 5,  # if you make it simply proportional to the v one, it ends up less than 1 (also, you really want it to depend both on the sequence length and the density of the imgt set for the region [and maybe relative shm rates in the region])
            'j' : 5,
        }

    # ----------------------------------------------------------------------------------------
    def finalize(self, gene_counts, annotations=None, regions=None, debug=False):
        # NOTE a.t.m. this is frequently using vsearch V alignments only, so we can't collapse clones beforehand. It would be more accurate, though, if we could collapse clones first (if using vsearch alignments, it also includes a lot of crap sequences that later get removed)
        # NOTE <gene_counts> is in general floats, not integers
        assert not self.finalized

        if regions is None:
            regions = utils.regions if gene_counts is None else [r for r in utils.regions if r in gene_counts]

        # if we pass in <annotations> instead of <gene_counts> we have to transfer the info (if we pass in both, we assume gene counts was already done correctly, and just use the annotations for simulation debug printing)
        if gene_counts is None:
            assert annotations is not None
            gene_counts = utils.get_gene_counts_from_annotations(annotations)

        if debug:
            total_counts = sum([sum(gene_counts[r].values()) for r in gene_counts]) / float(len(gene_counts))  # ok it's weird to take the average, but they should all be the same, and it's cleaner than choosing one of 'em
            print('  removing least likely genes with total counts %.1f%s' % (total_counts, '' if self.simglfo is None else (' (%d in simulation)' % sum(self.simcounts[regions[0]].values()))))  # should really take the average for the sim one as well

        for region in regions:
            sorted_gene_counts = sorted(list(gene_counts[region].items()), key=lambda x: (x[1], x[0]), reverse=True)  # sort first by count, break ties alphabetically by gene name (the latter to get a repeatable order)
            self.finalize_region(region, sorted_gene_counts, annotations=annotations, debug=debug)

    # ----------------------------------------------------------------------------------------
    def separate_into_classes(self, region, sorted_gene_counts, easycounts):  # where each class contains all alleles with the same distance from start to cyst, and within a hamming distance of <self.args.n_max_snps>
        snp_groups = utils.separate_into_snp_groups(self.glfo, region, self.n_max_snps[region], genelist=[g for g, _ in sorted_gene_counts])
        class_counts = []  # this could stand to be cleaned up... it's kind of a holdover from before I moved the separating fcn to utils
        for sgroup in snp_groups:
            class_counts.append([{'gene' : gfo['gene'], 'counts' : easycounts[gfo['gene']], 'seq' : gfo['seq'], 'hdist' : gfo['hdist']} for gfo in sgroup])
        return class_counts

    # ----------------------------------------------------------------------------------------
    def finalize_region(self, region, sorted_gene_counts, annotations=None, debug=False):
        easycounts = {gene : counts for gene, counts in sorted_gene_counts}
        total_counts = sum([counts for counts in easycounts.values()])
        class_counts = self.separate_into_classes(region, sorted_gene_counts, easycounts)

        genes_to_keep = set()

        if debug:
            print('   %s groups separated by %d snps  (-: same group as previous kept gene)' % (utils.color('blue', region), self.n_max_snps[region]))
            print('     %-20s       %5s %s        removed genes (snps counts%s)%s%s' % ('genes to keep', 'counts',
                                                                                        '' if self.simglfo is None else utils.color('blue', 'sim'),
                                                                                        '' if self.simglfo is None else utils.color('blue', ' sim counts'),
                                                                                        '' if self.simglfo is None else ('  ' + utils.color('red', 'x:') + ' not in simulation'),
                                                                                        '' if (annotations is None or self.reco_info is None) else ('               %s sim counts/genes for the queries assigned to this kept gene %s' % (utils.color('blue', '['), utils.color('blue', ']'))),
            ), end=' ')
            def count_str(cnt):
                if cnt < 10.:
                    return '%.1f' % cnt
                else:
                    return '%.0f' % cnt
            def simcountstr(gene, ws):  # counts in simulation for <gene> (note that this is _not_ the same as sim_gene_count_str(), since this takes no account of _which_ queries these counts occur in [plus it's coming from the opposite point of view])
                if self.simglfo is None:
                    rstr = ''
                elif gene in self.simglfo['seqs'][utils.get_region(gene)]:
                    rstr = utils.color('blue', (' %' + ws + 'd') % self.simcounts[utils.get_region(gene)][gene])
                else:
                    rstr = utils.color('red', (' %' + ws + 's') % 'x')
                return rstr
            def sim_gene_count_str(kgene):  # figure out simulation genes and counts for the uids assigned to <kgene>
                if annotations is None or self.reco_info is None:
                    return ''
                uids_this_gene = [uid for uid, line in annotations.items() if line[region + '_gene'] == kgene]
                sim_genes = {}  # simulation genes for the uids that we assigned to <kgene> (note that self.simcounts doesn't have this per-uid information)
                for uid in uids_this_gene:
                    sgene = self.reco_info[uid][region + '_gene']
                    if sgene not in sim_genes:
                        sim_genes[sgene] = 0
                    sim_genes[sgene] += 1
                sorted_sim_gene_counts = sorted(list(sim_genes.items()), key=operator.itemgetter(1), reverse=True)
                count_str = ' '.join([utils.color('blue' if sg == kgene else 'red', str(c)) for sg, c in sorted_sim_gene_counts])
                sgene_str = ' '.join([utils.color_gene(sg) for sg, _ in sorted_sim_gene_counts])
                return '%s   %s' % (count_str, sgene_str)

        for iclass in range(len(class_counts)):
            gclass = class_counts[iclass]
            kept_this_class = []
            for ig in range(len(gclass)):
                gfo = gclass[ig]

                if float(gfo['counts']) / total_counts < self.args.min_allele_prevalence_fractions[region]:  # always skip everybody that's super uncommon
                    pass  # don't keep it
                elif ig == 0:  # keep the first one from this class
                    genes_to_keep.add(gfo['gene'])
                    kept_this_class.append(gfo['gene'])
                elif utils.hamming_distance(gclass[0]['seq'], gclass[ig]['seq']) == 0:  # don't keep it if it's indistinguishable from the most common one (the matches are probably mostly really the best one)
                    pass  # don't keep it
                elif len(kept_this_class) < self.alleles_per_group:  # always keep the most common <self.alleles_per_group> in each class [note: defaults to 1 if looking for new alleles, otherwise 2]
                    genes_to_keep.add(gfo['gene'])
                    kept_this_class.append(gfo['gene'])
                else:
                    pass  # don't keep it

                if debug and gfo['gene'] in genes_to_keep:
                    snpstr = ' ' if ig == 0 else '(%d)' % utils.hamming_distance(gclass[0]['seq'], gfo['seq'])  # only happens if we keep more than one from this class
                    print('\n      %s%-s  %7s%s  %-3s' % ('- ' if ig > 0 else '  ', utils.color_gene(gfo['gene'], width=20), count_str(gfo['counts']), simcountstr(gfo['gene'], '4'), snpstr), end=' ')
            if debug:
                if len(kept_this_class) == 0:
                    print('\n      %s%-s  %7s%4s  %-3s' % ('  ', utils.color('blue', 'none', width=20, padside='right'), '-', '', ''), end=' ')
                removedfo = [gfo for gfo in gclass if gfo['gene'] not in genes_to_keep]
                removed_str = ''
                if len(removedfo) > 0:
                    number_strs = ['(%d %3s%s)' % (gfo['hdist'], count_str(gfo['counts']), simcountstr(gfo['gene'], '1')) for gfo in removedfo]
                    name_strs = ['%s' % (utils.color_gene(gfo['gene'])) for gfo in removedfo]
                    removed_str = '%s  %s' % (' '.join(number_strs), ' '.join(name_strs))
                annotation_str = ''
                if (annotations is not None and self.reco_info is not None) and len(kept_this_class) > 0:
                    annotation_str = '%s %s %s' % (utils.color('blue', '['), sim_gene_count_str(kept_this_class[-1]), utils.color('blue', ']'))
                print('     %s  %s  %s' % (removed_str, (70 - utils.len_excluding_colors(removed_str)) * ' ', annotation_str), end=' ')
        if debug:
            print('')

        genes_to_remove = set(self.glfo['seqs'][region]) - genes_to_keep

        print('    keeping %d / %d %s gene%s' % (len(genes_to_keep), len(self.glfo['seqs'][region]), region, utils.plural(len(genes_to_keep))))
        if len(genes_to_keep) == 0:
            print('   would\'ve kept zero genes, instead keeping all of them')
            genes_to_keep = copy.deepcopy(genes_to_remove)
            genes_to_remove.clear()

        if self.simglfo is not None:
            missing_genes = set(self.simglfo['seqs'][region]) - genes_to_keep
            if len(missing_genes) > 0:
                print('    %s %d simulation genes (counts): %s' % (utils.color('red', 'missing'), len(missing_genes), '  '.join([('%s %d' % (utils.color_gene(g), self.simcounts[region][g])) for g in sorted(missing_genes)])))
            completely_absent_genes = missing_genes - genes_to_remove
            if len(completely_absent_genes) > 0:
                print('%s %d simulation genes completely absent: %s' % (utils.color('red', 'warning'), len(completely_absent_genes),  '  '.join([('%s %d' % (utils.color_gene(g), self.simcounts[region][g])) for g in sorted(completely_absent_genes)])))

        self.genes_to_keep |= genes_to_keep  # add the ones from _this_ region (rhs) to the ones from all regions (lhs)
        self.genes_to_remove |= genes_to_remove

        self.finalized = True
