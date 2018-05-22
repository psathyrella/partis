import sys
import glob
import time
import operator
from subprocess import check_call
import copy
import os

import utils
from hist import Hist

# ----------------------------------------------------------------------------------------
class AlleleRemover(object):
    def __init__(self, glfo, args):
        self.glfo = glfo
        self.args = args

        self.genes_to_keep = None
        self.genes_to_remove = None
        self.dbg_strings = {}

        self.finalized = False

        self.n_max_snps = {
            'v' : self.args.n_max_snps - 2,  # not sure why I originally put the 2 there, maybe just to adjust the space between where allelefinder stops looking and where alleleclusterer starts?
            'd' : 5,  # if you make it simply proportional to the v one, it ends up less than 1 (also, you really want it to depend both on the sequence length and the density of the imgt set for the region [and maybe relative shm rates in the region])
            'j' : 5,
        }

    # ----------------------------------------------------------------------------------------
    def separate_into_classes(self, region, sorted_gene_counts, easycounts):  # where each class contains all alleles with the same distance from start to cyst, and within a hamming distance of <self.args.n_max_snps>
        snp_groups = utils.separate_into_snp_groups(self.glfo, region, self.n_max_snps[region], genelist=[g for g, _ in sorted_gene_counts])
        class_counts = []  # this could stand to be cleaned up... it's kind of a holdover from before I moved the separating fcn to utils
        for sgroup in snp_groups:
            class_counts.append([{'gene' : gfo['gene'], 'counts' : easycounts[gfo['gene']], 'seq' : gfo['seq'], 'hdist' : gfo['hdist']} for gfo in sgroup])
        return class_counts

    # ----------------------------------------------------------------------------------------
    def finalize(self, gene_counts, debug=False):
        # NOTE a.t.m. this is frequently using vsearch V alignments only, so we can't collapse clones beforehand. It would be more accurate, though, if we could collapse clones first (if using vsearch alignments, it also includes a lot of crap sequences that later get removed)
        # NOTE <gene_counts> is in general floats, not integers
        assert not self.finalized
        if debug:
            print '  removing least likely genes (%.1f total counts)' % (sum([sum(gene_counts[r].values()) for r in gene_counts]) / float(len(gene_counts)))  # ok it's weird to take they average, but they should all be the same, and it's cleaner than choosing one of 'em

        for region in [r for r in utils.regions if r in gene_counts]:
            self.finalize_region(region, sorted(gene_counts[region].items(), key=operator.itemgetter(1), reverse=True), debug=debug)

    # ----------------------------------------------------------------------------------------
    def finalize_region(self, region, sorted_gene_counts, debug=False):
        easycounts = {gene : counts for gene, counts in sorted_gene_counts}
        total_counts = sum([counts for counts in easycounts.values()])
        class_counts = self.separate_into_classes(region, sorted_gene_counts, easycounts)

        self.genes_to_keep = set()

        if debug:
            print '  %s (groups separated by %d snps)' % (utils.color('blue', region), self.n_max_snps[region])
            print '     %-20s    %5s (%s)      removed genes (snps counts)  names' % ('genes to keep', 'counts', 'snps'),
            def count_str(cnt):
                if cnt < 10.:
                    return '%.1f' % cnt
                else:
                    return '%.0f' % cnt

        for iclass in range(len(class_counts)):
            gclass = class_counts[iclass]
            n_from_this_class = 0
            for ig in range(len(gclass)):
                gfo = gclass[ig]

                if float(gfo['counts']) / total_counts < self.args.min_allele_prevalence_fraction:  # always skip everybody that's super uncommon
                    pass  # don't keep it
                elif ig == 0:  # keep the first one from this class
                    self.genes_to_keep.add(gfo['gene'])
                    n_from_this_class += 1
                elif utils.hamming_distance(gclass[0]['seq'], gclass[ig]['seq']) == 0:  # don't keep it if it's indistinguishable from the most common one (the matches are probably mostly really the best one)
                    pass  # don't keep it
                elif n_from_this_class < self.args.n_alleles_per_gene:  # always keep the most common <self.args.n_alleles_per_gene> in each class [note: defaults to 1 if looking for new alleles, otherwise 2]
                    self.genes_to_keep.add(gfo['gene'])
                    n_from_this_class += 1
                else:
                    pass  # don't keep it

                if debug and gfo['gene'] in self.genes_to_keep:
                    snpstr = ' ' if ig == 0 else '(%d)' % utils.hamming_distance(gclass[0]['seq'], gfo['seq'])
                    print '\n       %-s  %7s  %-3s' % (utils.color_gene(gfo['gene'], width=20), count_str(gfo['counts']), snpstr),
            if debug:
                if n_from_this_class == 0:
                    print '\n       %-s  %7s  %-3s' % (utils.color('blue', 'none', width=20, padside='right'), '-', ''),
                removedfo = [gfo for gfo in gclass if gfo['gene'] not in self.genes_to_keep]
                if len(removedfo) > 0:
                    number_strs = ['(%d %s)' % (gfo['hdist'], count_str(gfo['counts'])) for gfo in removedfo]
                    name_strs = ['%s' % utils.color_gene(gfo['gene']) for gfo in removedfo]
                    print '        %s  %s' % (' '.join(number_strs), ' '.join(name_strs)),
        if debug:
            print ''

        self.genes_to_remove = set(self.glfo['seqs'][region]) - self.genes_to_keep

        print '    keeping %d / %d %s gene%s' % (len(self.genes_to_keep), len(self.glfo['seqs'][region]), region, utils.plural(len(self.genes_to_keep)))
        # print '    removing %d %s genes: %d with no matches, %d with unconvincing matches' % (len(self.genes_to_remove), region, len(set(self.glfo['seqs'][region]) - set(easycounts)), len(set(easycounts) - self.genes_to_keep))
        if len(self.genes_to_keep) == 0:
            print '   would\'ve kept zero genes, instead keeping all of them'
            self.genes_to_keep = copy.deepcopy(self.genes_to_remove)
            self.genes_to_remove.clear()

        self.finalized = True
