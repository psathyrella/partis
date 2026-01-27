#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import csv
import os
import sys
import argparse
import colored_traceback.always
import glob
import subprocess
import operator
from io import open

import partis.utils as utils
import partis.glutils as glutils
from partis.clusterpath import ClusterPath

parser = argparse.ArgumentParser()
# NOTE to compare multiple output runs see  datascripts/meta/qa013-synth/read-lh-cf.py
parser.add_argument('--partis-file', required=True, help='partis yaml partition output file that includes alternative annotation information (i.e. --calculate-alternative-annotations was set while partitioning)')
parser.add_argument('--linearham-dir', required=True, help='linearham output dir (the main/parent dir))')
parser.add_argument('--prob-to-ignore', default=0.15, help='don\'t print sequences with probabilities smaller than this')
parser.add_argument('--outdir', help='if set, write csv info for the printed naive sequences to here')
args = parser.parse_args()

# ----------------------------------------------------------------------------------------
def read_linearham_output():
    lh_info = {}
    glbfns = [f for gstr in ['cluster', 'iclust-'] for f in glob.glob('%s/%s*' % (args.linearham_dir, gstr))]
    clusterdirs = [d for d in (glbfns) if os.path.isdir(d)]  # works for old style 'clusterN' and new-style 'cluster-N'
    if len(clusterdirs) == 0:
        raise Exception('no linearham cluster subdirs (of form clusterN/ or cluster-N/ found in %s' % args.linearham_dir)
    print('  reading linearham info for %d cluster%s: %s' % (len(clusterdirs), utils.plural(len(clusterdirs)),' '.join(os.path.basename(cd) for cd in clusterdirs)))
    for cdir in clusterdirs:
        sfnames = [f for gstr in ['', 'lineage*/'] for f in glob.glob('%s/%scluster_seqs.fasta' % (cdir, gstr))]
        if len(sfnames) == 0:
            raise Exception('no sequence files \'*_seqs.fasta\' found in %s' % cdir)
        sfn = utils.get_single_entry(sfnames)
        input_seqfos = utils.read_fastx(sfn)
        input_uids = [sfo['name'] for sfo in input_seqfos if sfo['name'] != 'naive']
        # aa_naive_seqs.fasta:   prob of each aa naive seq
        # aa_naive_seqs.dnamap:  prob of each nuc naive seq contributing to each of those aa naive seqs
        aasfn = '%s/aa_naive_seqs.fasta'%os.path.dirname(sfn)
        aa_seq_infos = utils.read_fastx(aasfn)
        for iseq, sfo in enumerate(aa_seq_infos):
            tlist = sfo['name'].split('_')
            assert len(tlist) == 3
            assert int(tlist[1]) == iseq
            sfo['prob'] = float(tlist[2])
        with open(aasfn.replace('.fasta', '.dnamap')) as outfile:  # this is some weird bastardization of a fasta file
            iseq = -1
            for line in outfile:
                if line[0] == '>':
                    iseq += 1
                    assert line.strip().lstrip('>') == aa_seq_infos[iseq]['name']
                    aa_seq_infos[iseq]['nuc_seqs_probs'] = []
                    continue
                prob, naive_nuc_seq = line.strip().split(',')
                aa_seq_infos[iseq]['nuc_seqs_probs'].append((naive_nuc_seq, float(prob)))
        lh_info[':'.join(input_uids)] = aa_seq_infos
    return lh_info

# ----------------------------------------------------------------------------------------
def print_naive_seq_lines(nseq_info, namestr, namecolor, ref_seq=None, amino_acid=False, writefo=None):
    def i_aa_color(i_aa):
        tmpcolors = ['purple', 'yellow', 'red', 'blue', 'green']
        return tmpcolors[i_aa % len(tmpcolors)]
    total_prob = 0.
    breaking = False
    for naive_seq, prob, i_aa_seq in sorted(nseq_info, key=operator.itemgetter(1), reverse=True):
        if ref_seq is None:
            ref_seq = naive_seq
        breakstr = ''
        if 1. - total_prob < args.prob_to_ignore:
            breaking = True
            breakstr = 'total: %5.2f (breaking after %.2f)' % (prob+total_prob, 1. - args.prob_to_ignore)
        print('     %s %s    %5.2f    %s   %s' % (utils.color_mutants(ref_seq, naive_seq, amino_acid=amino_acid, align_if_necessary=True),
                                                  utils.color(i_aa_color(i_aa_seq), str(i_aa_seq), width=2), prob, utils.color(namecolor, namestr, width=9, padside='right'), breakstr))
        if writefo is not None:
            writefo.append({'method' : namestr, 'prob' : prob, 'seq' : naive_seq})
        if breaking:
            break
        total_prob += prob
    print('')
    return ref_seq

# ----------------------------------------------------------------------------------------
def print_all_lines(lh_aa_seq_infos, pline, amino_acid=False):
    seq_len = len(lh_aa_seq_infos[0]['seq'] if amino_acid else pline['naive_seq'])
    anstr = '%s %s naive seqs' % (headstr('1.' if amino_acid else '3.'), 'amino acid' if amino_acid else 'nucleotide')
    print('  %s:%s aa seq' % (anstr, (seq_len - utils.len_excluding_colors(anstr)) * ' '))
    if amino_acid:
        codon_str = utils.color('reverse_video', 'X')
        vpos, jpos = [pline['codon_positions'][r] // 3 for r in ['v', 'j']]
        cdstr = '%s%s%s%s%s' % (' '*vpos, codon_str, '-'*(jpos - vpos - 1), codon_str, ' '*(seq_len - jpos - 1))
    else:
        cdstr = seq_len*' '
    print('    %s index  prob' % cdstr)
    writefo=[]
    ref_seq = print_naive_seq_lines(get_lh_nsinfo(lh_aa_seq_infos, amino_acid=amino_acid), 'linearham', 'green', amino_acid=amino_acid, writefo=writefo)
    _ = print_naive_seq_lines(get_partis_nsinfo(pline, amino_acid=amino_acid), 'partis', 'blue', ref_seq=ref_seq, amino_acid=amino_acid, writefo=writefo)  # use the linearham naive seq as ref_seq also for partis
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    if args.outdir is not None:
        print('  writing %s seqs to %s' % ('aa' if amino_acid else 'nuc', args.outdir))
        with open('%s/%s-seqs.csv'%(args.outdir, 'aa' if amino_acid else 'nuc'), 'w') as jfile:
            writer = csv.DictWriter(jfile, writefo[0].keys())
            writer.writeheader()
            for wfo in writefo:
                writer.writerow(wfo)

# ----------------------------------------------------------------------------------------
def print_gene_calls(pline):
# TODO would be nice to read the linearham annotation so that we can highlight the v/d/j genes that linearham assumed were correct
    print('  %s partis gene calls (linearham only considers one gene combo):' % headstr('4.'))
    print('        prob   gene')
    for region in utils.regions:
        print('      %s' % utils.color('blue', region))
        for gene, prob in pline['alternative-annotations']['gene-calls'][region]:
            print('        %4.2f  %s' % (prob, utils.color_gene(gene, width=15)))

# ----------------------------------------------------------------------------------------
def get_lh_nsinfo(lh_aa_seq_infos, amino_acid=False):
    if amino_acid:
        return [(s['seq'], s['prob'], i) for i, s in enumerate(lh_aa_seq_infos)]
    else:
        return [(ns, np, i) for i, s in enumerate(lh_aa_seq_infos) for ns, np in s['nuc_seqs_probs']]

# ----------------------------------------------------------------------------------------
def get_partis_nsinfo(pline, amino_acid=False):
    nuc_naive_seqs = pline['alternative-annotations']['naive-seqs'] if 'alternative-annotations' in pline else [(pline['naive_seq'], 1.), ]
    pdict, sdict = {}, {}
    for aseq, nseq, prob in [(utils.ltranslate(nseq, trim=True), nseq, prob) for nseq, prob in nuc_naive_seqs]:  # add up the probs for any nuc seqs that code for the same aa seq
        if aseq not in pdict:
            pdict[aseq] = 0.
            sdict[aseq] = []
        pdict[aseq] += prob
        sdict[aseq].append(nseq)
    aa_naive_seqs = sorted(list(pdict.items()), key=operator.itemgetter(1), reverse=True)
    if amino_acid:
        return [(s, p, i) for i, (s, p) in enumerate(aa_naive_seqs)]
    else:
        nuc_seq_aa_indices = {}
        for iseq, (aseq, _) in enumerate(aa_naive_seqs):
            for nseq in sdict[aseq]:
                nuc_seq_aa_indices[nseq] = iseq
        return [(s, p, nuc_seq_aa_indices[s]) for s, p in nuc_naive_seqs]

# ----------------------------------------------------------------------------------------
def headstr(hstr):
    return utils.color('green', hstr)

# ----------------------------------------------------------------------------------------
glfo, annotation_list, cpath = utils.read_output(args.partis_file)
lh_info = read_linearham_output()

annotations = {':'.join(adict['unique_ids']) : adict for adict in annotation_list}  # collect the annotations in a dictionary so they're easier to access
most_likely_partition = cpath.partitions[cpath.i_best]  # a partition is represented as a list of lists of strings, with each string a sequence id
print('  %d (of %d) clusters from partis file share uids with %d linearham cluster%s' % (len([c for c in most_likely_partition if any(len(set(c) & set(lc.split(':'))) > 0 for lc in lh_info)]), len(most_likely_partition), len(lh_info), utils.plural(len(lh_info))))
sorted_clusters = sorted(most_likely_partition, key=len, reverse=True)
for cluster in sorted_clusters:
    pline = annotations[':'.join(cluster)]
    utils.trim_fwk_insertions(glfo, pline, modify_alternative_annotations=True)  # linearham probably won't have them, so we need to remove them so things line up
    if 'alternative-annotations' not in pline:
        print('  note: no alternative annotations in %s, so can\'t print partis alternative naive sequences' % args.partis_file)

    p_uids = set(pline['unique_ids'])
    lh_clusters = [(uidstr, cfo) for uidstr, cfo in lh_info.items() if set(uidstr.split(':')) & p_uids]  # lh clusters with any uids in common iwth this partis <cluster> (there should only be 1)
    lh_aa_seq_infos = []
    if len(lh_clusters) == 0:
        # print '  no linearham clusters with any of these uids' % utils.color('red', 'error')
        continue
    elif len(lh_clusters) != 1:
        raise Exception('should have only one linearham cluster with uids in common with this cluster, but found %d' % len(lh_clusters))
    else:
        lh_uidstr, lh_aa_seq_infos = lh_clusters[0]
        lh_uids = set(lh_uidstr.split(':'))
        print('%s with sizes: partis %d   linearham %d (%d in common)' % (utils.color('blue', 'starting clusters'), len(pline['unique_ids']), len(lh_uids), len(lh_uids & p_uids)))
        if len(lh_uids - p_uids) > 0:
            print('  %s %d extra uids in linearham cluster' % (utils.color('yellow', 'warning'), len(lh_uids - p_uids)))
        if len(p_uids - lh_uids) > 0:
            print('  %s %d extra uids in partis cluster' % (utils.color('yellow', 'warning'), len(p_uids - lh_uids)))
    if len(lh_aa_seq_infos) == 0:
        print('  %s no prob/naive_seq pairs in linearham for this cluster' % utils.color('red', 'error'))

    print_all_lines(lh_aa_seq_infos, pline, amino_acid=True)
    utils.print_reco_event(utils.synthesize_single_seq_line(pline, iseq=0), extra_str='    ', label='%s annotation for a single (arbitrary) sequence from the cluster:'%headstr('2.'))
    print_all_lines(lh_aa_seq_infos, pline, amino_acid=False)

    if 'alternative-annotations' in pline:
        print_gene_calls(pline)

    print('')
