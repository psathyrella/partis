#!/usr/bin/env python3
"""Remove genes from a germline set whose conserved codon (v cysteine, j tryptophan/phenylalanine) can't be used: the sequence ends before the codon
('truncated', e.g. imgt "partial in 3'" entries) or, for v, the codon is out of frame with the start of the gene ('out-of-frame', e.g. imgt "partial in 5'"
entries, or a frameshift). For a truncated v, extras.csv puts the cysteine at the end of the sequence, i.e. short of the real cysteine, so reads assigned to it
get the wrong cdr3 length and are called nonproductive (issue #421).

Genes whose codon is intact but mutated (e.g. IGHV1-38-4*01, cysteine TAT) are real germline sequences that partis annotates correctly, so they're not
removed here; the default sets drop them anyway, since they contain only functional genes (see remove-nonfunctional-genes.py).

Genes in <extra_genes_to_remove> are removed as well, with the reason given there.

Run from the partis main dir, e.g.:
  ./data/germlines/remove-unusable-genes.py --species human --write
"""

import argparse
import os
import sys

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/data/germlines', '')
if partis_dir not in sys.path:
    sys.path.insert(1, partis_dir)
from partis import utils
from partis import glutils

extra_genes_to_remove = {
    'human' : {
        'IGLJ2*01_8-------------------------------38' : 'vdjbase deletion allele with only 7 bases left (TGTGGTA)',
    },
}

# ----------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--species', required=True)
    parser.add_argument('--gldir', help='default data/germlines/<species>')
    parser.add_argument('--loci', default='igh:igk:igl')
    parser.add_argument('--write', action='store_true', help='write the germline set back to --gldir (default: dry run)')
    args = parser.parse_args()
    args.gldir = utils.non_none([args.gldir, 'data/germlines/%s' % args.species])
    args.loci = utils.get_arg_list(args.loci)

    for locus in args.loci:
        glfo = glutils.read_glfo(args.gldir, locus)
        removed = glutils.remove_genes_with_bad_codons(glfo, statuses=['truncated', 'out-of-frame'])
        extra_genes = [g for g in extra_genes_to_remove.get(args.species, {}) if g in glfo['seqs'][utils.get_region(g)]]
        glutils.remove_genes(glfo, extra_genes)
        removed.update((g, extra_genes_to_remove[args.species][g]) for g in extra_genes)
        print('%s: removing %d gene%s' % (utils.color('blue', locus), len(removed), utils.plural(len(removed))))
        for gene, reason in removed.items():
            print('    %-50s %s' % (gene, reason))
        if args.write and len(removed) > 0:
            glutils.write_glfo(args.gldir, glfo)
            print('  wrote %s to %s' % (locus, args.gldir))
    if not args.write:
        print('\n%s dry run (no --write): nothing written' % utils.color('yellow', 'note'))

# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
