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
    def __init__(self, glfo, args, alfinder):
        self.glfo = glfo
        self.args = args
        self.alfinder = alfinder

        self.genes_to_keep = None
        self.genes_to_remove = None
        self.dbg_strings = {}
        self.region = 'v'

        self.finalized = False

    # ----------------------------------------------------------------------------------------
    def separate_into_classes(self, sorted_gene_counts, easycounts):  # where each class contains all alleles with the same distance from start to cyst, and within a hamming distance of <self.args.n_max_snps>
        snp_groups = utils.separate_into_snp_groups(self.glfo, self.region, self.args.n_max_snps, genelist=[g for g, _ in sorted_gene_counts])
        class_counts = []  # this could stand to be cleaned up... it's kind of a holdover from before I moved the separating fcn to utils
        for sgroup in snp_groups:
            class_counts.append([{'gene' : gfo['gene'], 'counts' : easycounts[gfo['gene']], 'seq' : gfo['seq']} for gfo in sgroup])
        return class_counts

    # ----------------------------------------------------------------------------------------
    def finalize(self, sorted_gene_counts, debug=False):
        # NOTE a.t.m. this is using vsearch V alignments only, so we can't collapse clones beforehand. It would be more accurate, though, if we could collapse clones first
        # NOTE <sorted_gene_counts> is usually/always floats instead of integers
        assert not self.finalized
        easycounts = {gene : counts for gene, counts in sorted_gene_counts}
        total_counts = sum([counts for counts in easycounts.values()])

        self.genes_to_keep = set()

        if debug:
            print '  removing least likely genes (%.1f total counts)' % total_counts
            print '     %-20s    %5s (%s)      removed genes (counts)' % ('genes to keep', 'counts', 'snps'),
            def count_str(cnt):
                if cnt < 10.:
                    return '%.1f' % cnt
                else:
                    return '%.0f' % cnt

        class_counts = self.separate_into_classes(sorted_gene_counts, easycounts)
        for iclass in range(len(class_counts)):
            gclass = class_counts[iclass]
            n_from_this_class = 0
            for ig in range(len(gclass)):
                gfo = gclass[ig]
                if self.args.n_max_total_alleles is not None and len(self.genes_to_keep) >= self.args.n_max_total_alleles:  # command line can specify the total number of alleles
                    break

                if float(gfo['counts']) / total_counts < self.args.min_allele_prevalence_fraction:  # always skip everybody that's super uncommon
                    pass
                elif ig == 0:  # keep the first one from this class
                    self.genes_to_keep.add(gfo['gene'])
                    n_from_this_class += 1
                elif utils.hamming_distance(gclass[0]['seq'], gclass[ig]['seq']) == 0:  # don't keep it if it's indistinguishable from the most common one (the matches are probably mostly really the best one)
                    pass  # don't keep it
                elif n_from_this_class < self.args.n_alleles_per_gene:  # always keep the most common <self.args.n_alleles_per_gene> in each class
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
                    removal_strs = ['%s (%s)' % (utils.color_gene(gfo['gene']), count_str(gfo['counts'])) for gfo in removedfo]
                    print '        %s' % '  '.join(removal_strs),
        if debug:
            print ''

        self.genes_to_remove = set(self.glfo['seqs'][self.region]) - self.genes_to_keep

        print '    keeping %d / %d %s gene%s' % (len(self.genes_to_keep), len(self.glfo['seqs'][self.region]), self.region, utils.plural(len(self.genes_to_keep)))
        # print '    removing %d %s genes: %d with no matches, %d with unconvincing matches' % (len(self.genes_to_remove), self.region, len(set(self.glfo['seqs'][self.region]) - set(easycounts)), len(set(easycounts) - self.genes_to_keep))

        self.finalized = True
