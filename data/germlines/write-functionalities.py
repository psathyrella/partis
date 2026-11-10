#!/usr/bin/env python3
"""Write <gldir>/functionalities.csv: the functionality ('F', 'ORF', 'P', or '' for unknown) of every gene in a germline set, with where it came from.

Labels come only from sequence, never from names: a name says nothing about a sequence (the same name can have different sequences in different
databases, and imgt has renamed some genes). For each gene, in order:
  - imgt-seq: exact sequence match to IMGT GENE-DB entries (all functionalities, from imgt-download/genedb/get-genedb.sh). If the sequence occurs
      at several loci with different functionality (e.g. IGHD4-4*01 F and IGHD4-11*01 ORF), we take the best (F > ORF > P), since a read can't tell them apart.
  - ogrdb (macaque only): exact sequence match to an OGRDB allele (ogrdb-download/macaque/get-functional.py); OGRDB's 'functional' is a boolean ("false if
      it is a pseudogene"), so we use F or P.
  - imgt-nearest: otherwise the label of the nearest IMGT GENE-DB entry of the same locus and region, if it's within --max-nearest-distance (mismatches
      plus gap positions; bases of our gene that the imgt entry doesn't cover count as differences, but not vice versa). Ties take the best label. For a v gene, an in-frame stop codon before
      the cysteine makes it P regardless.
  - otherwise '' (unknown).
We don't otherwise judge genes here (see data/germlines/README.md and issue #421).

Run from the partis main dir, e.g.:
  ./data/germlines/write-functionalities.py --species human --write
"""

import argparse
import csv
import os
import sys
from collections import defaultdict, Counter

from Bio import Align

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/data/germlines', '')
if partis_dir not in sys.path:
    sys.path.insert(1, partis_dir)
from partis import utils
from partis import glutils

func_rank = {'F' : 0, 'ORF' : 1, 'P' : 2}

# ----------------------------------------------------------------------------------------
def best_label(labels):
    return min(labels, key=lambda l: func_rank[l])

# ----------------------------------------------------------------------------------------
def read_imgt(fname):
    """ return {locus : {region : [(name, functionality, seq), ...]}} for imgt GENE-DB fasta <fname> """
    entries = defaultdict(lambda: defaultdict(list))
    name, info = None, None
    seqs = {}
    with open(fname) as ffile:
        for line in ffile:
            line = line.strip()
            if line[:1] == '>':
                info = line[1:].split('|')
                name = info[1]
                seqs[name] = [glutils.strip_functionality(info[3]), '']
            elif name is not None:
                seqs[name][1] += line.upper()
    for name, (functionality, seq) in seqs.items():
        entries[name[:3].lower()][name[3].lower()].append((name, functionality, seq))
    return entries

# ----------------------------------------------------------------------------------------
def read_ogrdb(fname):
    """ return {seq : [(label, functionality), ...]} """
    ogrdb = defaultdict(list)
    with open(fname) as ffile:
        for line in csv.DictReader(ffile):
            ogrdb[line['sequence']].append((line['label'], 'F' if line['functional'] == 'True' else 'P'))
    return ogrdb

# ----------------------------------------------------------------------------------------
def has_stop_before_cysteine(seq, cpos):
    return any(seq[i : i + 3] in utils.codon_table['stop'] for i in range(0, cpos, 3))

# ----------------------------------------------------------------------------------------
def get_label(glfo, region, gene, seq, imgt_entries, imgt_by_seq, ogrdb, aligner, args):
    if seq in imgt_by_seq:
        names, labels = zip(*imgt_by_seq[seq])
        return best_label(labels), 'imgt-seq:%s' % '/'.join(names)
    if ogrdb is not None and seq in ogrdb:
        names, labels = zip(*ogrdb[seq])
        return best_label(labels), 'ogrdb:%s' % '/'.join(names)
    dists = [(int(-aligner.score(seq, iseq)), iname, ilabel) for iname, ilabel, iseq in imgt_entries]
    if len(dists) > 0:
        min_dist = min(d for d, _, _ in dists)
        if min_dist <= args.max_nearest_distance:
            nearest = [(n, l) for d, n, l in dists if d == min_dist]
            label, source = best_label([l for _, l in nearest]), 'imgt-nearest:%s:%d' % ('/'.join(n for n, _ in nearest), min_dist)
            if region == 'v' and has_stop_before_cysteine(seq, glfo['cyst-positions'][gene]):
                label, source = 'P', source + ':stop-codon'
            return label, source
    return '', 'none'

# ----------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--species', required=True, choices=['human', 'macaque'])
    parser.add_argument('--gldir', help='germline dir to label (default data/germlines/<species>)')
    parser.add_argument('--imgt-fname', help='imgt GENE-DB fasta (default data/germlines/imgt-download/genedb/<species>.fasta)')
    parser.add_argument('--ogrdb-fname', help='OGRDB functional csv (default, for macaque only, data/germlines/ogrdb-download/macaque/ogrdb-functional.csv)')
    parser.add_argument('--loci', default='igh:igk:igl')
    parser.add_argument('--max-nearest-distance', type=int, default=8, help='see docstring (default is the same as partis\'s --n-max-snps for new alleles)')
    parser.add_argument('--write', action='store_true', help='write <gldir>/functionalities.csv (default: dry run, print summary only)')
    args = parser.parse_args()
    args.gldir = utils.non_none([args.gldir, 'data/germlines/%s' % args.species])
    args.imgt_fname = utils.non_none([args.imgt_fname, 'data/germlines/imgt-download/genedb/%s.fasta' % args.species])
    if args.ogrdb_fname is None and args.species == 'macaque':
        args.ogrdb_fname = 'data/germlines/ogrdb-download/macaque/ogrdb-functional.csv'
    args.loci = utils.get_arg_list(args.loci)

    imgt = read_imgt(args.imgt_fname)
    ogrdb = read_ogrdb(args.ogrdb_fname) if args.ogrdb_fname is not None else None
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score, aligner.mismatch_score = 0, -1
    aligner.open_gap_score, aligner.extend_gap_score = -1, -1
    aligner.target_end_gap_score = 0  # our gene (the target) may be shorter than the imgt entry at either end without penalty...
    aligner.query_end_gap_score = -1  # ...but parts of our gene that the imgt entry doesn't cover count as differences (otherwise a short partial imgt entry would be 'near' anything that contains it)

    rows = []
    for locus in args.loci:
        glfo = glutils.read_glfo(args.gldir, locus)
        for region in utils.getregions(locus):
            imgt_by_seq = defaultdict(list)
            for iname, ilabel, iseq in imgt[locus][region]:
                imgt_by_seq[iseq].append((iname, ilabel))
            counts = Counter()
            for gene, seq in glfo['seqs'][region].items():
                if gene == glutils.dummy_d_genes.get(locus):
                    continue
                label, source = get_label(glfo, region, gene, seq, imgt[locus][region], imgt_by_seq, ogrdb, aligner, args)
                rows.append({'gene' : gene, 'functionality' : label, 'source' : source, 'seq' : seq})
                counts[(label if label != '' else 'unknown', source.split(':')[0])] += 1
            print('  %s%s %4d: %s' % (locus, region, sum(counts.values()), '  '.join('%s %s %d' % (l, s, n) for (l, s), n in sorted(counts.items()))))

    print('  total: %s' % '  '.join('%s %d' % (l if l != '' else 'unknown', n) for l, n in sorted(Counter(r['functionality'] for r in rows).items())))
    if args.write:
        outfname = '%s/functionalities.csv' % args.gldir
        with open(outfname, 'w') as ofile:
            writer = csv.DictWriter(ofile, glutils.functionality_headers)
            writer.writeheader()
            writer.writerows(rows)
        print('  wrote %d genes to %s' % (len(rows), outfname))

# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
