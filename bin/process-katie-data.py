#!/usr/bin/env python
import numpy
import csv
import operator
import sys
import os

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils

# ----------------------------------------------------------------------------------------
def build_v_gene_set(glfo, introns):
    total_d_counts = {}
    refseqs = {}
    for d_gene, counts in introns.items():
        total_d_counts[d_gene] = sum(counts.values())
    for d_gene, _ in sorted(total_d_counts.items(), key=operator.itemgetter(1), reverse=True):
        counts = introns[d_gene]

        # first decide on the reference sequences
        refseq, column_counts = None, None
        for seq in sorted(counts, key=len, reverse=True):
            if refseq is None:  # first one, i.e. the longest
                refseq = seq
                column_counts = [{n : 0 for n in utils.nukes} for i in range(len(refseq))]
            ioffset = len(refseq) - len(seq)
            partial_refseq = refseq[ioffset:]
            assert len(partial_refseq) == len(seq)
            for ibase in range(ioffset, len(refseq)):
                column_counts[ibase][seq[ibase - ioffset]] += counts[seq]

        refseqs[d_gene] = []
        for basecounts in column_counts:
            most_common_base = sorted(basecounts.items(), key=operator.itemgetter(1), reverse=True)[0][0]
            refseqs[d_gene].append(most_common_base)
        refseqs[d_gene] = ''.join(refseqs[d_gene])

        n_ok = 0
        mutecounts = {}
        for seq in sorted(counts, key=len, reverse=True):
            # print '    %3d   %150s' % (count, seq)
            partial_refseq = refseqs[d_gene][len(refseqs[d_gene]) - len(seq):]
            if seq == partial_refseq:
                n_ok += counts[seq]
            else:
                # utils.color_mutants(partial_refseq, seq, print_result=True, extra_str='                ')
                n_mutes = utils.hamming_distance(partial_refseq, seq)
                if n_mutes not in mutecounts:
                    mutecounts[n_mutes] = 0
                mutecounts[n_mutes] += counts[seq]
        print '  %s   %4d / %-4d ok' % (utils.color_gene(d_gene, width=10), n_ok, n_ok + sum(mutecounts.values())),
        if len(mutecounts) > 0:
            print '(mean of %.1f mutations among the other %d' % (numpy.average(mutecounts.keys(), weights=mutecounts.values()), sum(mutecounts.values())),
        print ''

    # add the intronic v genes to glfo
    for d_gene, refseq in refseqs.items():
        glfo['seqs']['v'][utils.generate_dummy_v(d_gene)] = refseq
        glfo['cyst-positions'][utils.generate_dummy_v(d_gene)] = len(refseq) - 3

    # write a glfo dir with everything
    glutils.write_glfo(outdir + '/germlines/imgt-and-intronic', glfo, debug=True)

    # remove the original v genes, and write a glfo dir with just the intronic ones
    glutils.remove_genes(glfo, [g for g in glfo['seqs']['v'] if 'xDx' not in g], debug=True)
    glutils.write_glfo(outdir + '/germlines/intronic', glfo, debug=True)


# tmpglfo = glutils.read_glfo('tmp-germlines', 'h')
glfo = glutils.read_glfo('data/germlines/human', 'h')

infname = '/fh/fast/matsen_e/data/2016-06-02-katie/VMO_Memory-3/VMO_Memory-3.tsv'
outdir = '/fh/fast/matsen_e/processed-data/partis/2016-06-02-katie'

n_failed, n_v_ok, n_total = 0, 0, 0
introns = {}

with open(infname) as infile:
    reader = csv.DictReader(infile, delimiter='\t')
    iline = 0
    introns = {}
    for line in reader:
        iline += 1
        # if iline > 1000:
        #     break
        n_total += 1
        if line[utils.adaptive_headers['v_gene']] != 'unresolved':  # for the moment, skip seqs that have v genes
            n_v_ok += 1
            continue
        line = utils.convert_from_adaptive_headers(glfo, line, uid='%09d' % iline, only_dj_rearrangements=True)
        if line['failed']:
            n_failed += 1
            continue

        d_gene = line['d_gene']
        if d_gene not in introns:
            introns[d_gene] = {}
        iseq = glfo['seqs']['v'][line['v_gene']]
        if iseq not in introns[d_gene]:
            introns[d_gene][iseq] = 0
        introns[d_gene][iseq] += 1

build_v_gene_set(glfo, introns)

print 'processed %d (of these:   %d failed   %d v ok)' % (n_total, n_failed, n_v_ok)
