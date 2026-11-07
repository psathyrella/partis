#!/usr/bin/env python3
"""Remove from a default germline set every gene that <gldir>/functionalities.csv doesn't label F, i.e. ORF, P, and unknown.

The default sets contain only functional genes (issue #421, https://github.com/psathyrella/partis/issues/421#issuecomment-6007296922):
  - anyone who cares about non-functional germline gene assignments is doing something uncommon, so is probably paying attention to the
    germline set, and can easily use whatever set they want (with --initial-germline-dir)
  - whereas non-functional genes in the default set will for sure give some spurious assignments of non-functional genes to functional BCRs
  - nonproductive BCRs come mostly (and maybe only) from out of frame rearrangements and stop codons, not from non-functional germline genes
Genes with unknown functionality are removed too: they're alleles that matched nothing in imgt or ogrdb, and non-ogrdb sources have lots of
poorly-supported alleles.

functionalities.csv has to be up to date (run write-functionalities.py first); afterwards we rewrite it for the remaining genes.

Run from the partis main dir, e.g.:
  ./data/germlines/remove-nonfunctional-genes.py --species human --write
"""

import argparse
import collections
import os
import sys

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/data/germlines', '')
if partis_dir not in sys.path:
    sys.path.insert(1, partis_dir)
from partis import utils
from partis import glutils

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
        functionalities = glutils.read_functionalities(args.gldir, glfo)
        genes = [g for r in utils.regions for g in glfo['seqs'][r] if g != glutils.dummy_d_genes.get(locus)]
        missing_genes = [g for g in genes if g not in functionalities]
        if len(missing_genes) > 0:
            raise Exception('%d %s gene%s missing from (or with a different sequence in) %s/functionalities.csv, so it\'s out of date (run write-functionalities.py first): %s' % (len(missing_genes), locus, utils.plural(len(missing_genes)), args.gldir, utils.color_genes(missing_genes)))
        removed = collections.OrderedDict((g, functionalities[g] if functionalities[g] != '' else 'unknown') for g in genes if functionalities[g] != 'F')
        glutils.remove_genes(glfo, list(removed))
        counts = collections.Counter(removed.values())
        print('%s: removing %d gene%s (%s)' % (utils.color('blue', locus), len(removed), utils.plural(len(removed)), ', '.join('%d %s' % (n, f) for f, n in sorted(counts.items()))))
        for gene, functionality in removed.items():
            print('    %-50s %s' % (gene, functionality))
        if args.write and len(removed) > 0:
            glutils.write_glfo(args.gldir, glfo)
            print('  wrote %s to %s' % (locus, args.gldir))
    if args.write:
        print('\n%s now rewrite the functionality file: ./data/germlines/write-functionalities.py --species %s --gldir %s --write' % (utils.color('yellow', 'note'), args.species, args.gldir))
    else:
        print('\n%s dry run (no --write): nothing written' % utils.color('yellow', 'note'))

# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
