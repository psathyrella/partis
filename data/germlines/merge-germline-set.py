#!/usr/bin/env python3
"""Merge a freshly-downloaded germline set into an existing partis germline dir.

Policy (see data/germlines/README.md and issue #396):
  - We do NOT judge the quality of individual genes here -- that is the job of
    whoever curated the set we downloaded. We just merge in what we downloaded.
  - The merge is always a UNION (glutils.get_merged_glfo): we never silently drop
    a gene we already had. Honoring an upstream removal is a separate, deliberate,
    documented decision -- not something this script does automatically.
  - It prints exactly what changed: per region, the counts plus the names of the genes
    added and the genes only in our old set; get_merged_glfo prints the two conflict
    types (same name / different seq, and same seq / different name), so we don't repeat
    those. Redirect stdout to a file for a full record (see data/germlines/*/changelog/).
  - The one warning this script adds of its own:
      * genes present in old but missing from new, when --expected replacement
        (a re-download of a same-source set legitimately may drop genes, but a big
        missing list usually means a botched download, so a human must look)
  - the two modes differ by where the new set comes from, which sets name/seq precedence
    on a conflict: --expected augmentation is genes from a *different* source (e.g.
    engelbrecht/rodriguez on top of imgt), so we keep our established (old) names;
    --expected replacement is a newer version of the *same* source (e.g. a fresh ogrdb
    pull), so it's authoritative and its (new) names/seqs win.

This supersedes the old mouse-only, python2 parse-ogrdb.py. Run with the venv active
(so `import partis...` resolves); from the partis main dir, e.g.:

  ./data/germlines/merge-germline-set.py --species macaque \\
      --new-dir data/germlines/ogrdb-download/macaque --sanitize \\
      --expected replacement --write

  ./data/germlines/merge-germline-set.py --species human \\
      --new-dir /fh/fast/matsen_e/data/vanwinkle-170/germlines/processed \\
      --expected augmentation --write
"""

import argparse
import os
import sys

# resolve imports whether or not partis is pip-installed: fall back to the repo root
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/data/germlines', '')
if partis_dir not in sys.path:
    sys.path.insert(1, partis_dir)
from partis import utils
from partis import glutils


# ----------------------------------------------------------------------------------------
def classify_region(old_glfo, new_glfo, region):
    """Split the genes in <region> into added / removed / seq_changed / name_changed / unchanged.

    A germline "gene" is a (name, sequence) pair, so a difference can be in either
    field. We key added/removed on *both* name and sequence so a gene that merely got
    renamed (or re-sequenced) upstream is reported as a change rather than as a
    simultaneous add + remove.
    """
    old_by_name = old_glfo['seqs'][region]          # {name: seq}
    new_by_name = new_glfo['seqs'][region]
    old_by_seq = {s: n for n, s in old_by_name.items()}  # {seq: name}
    new_by_seq = {s: n for n, s in new_by_name.items()}

    seq_changed = [(n, old_by_name[n], new_by_name[n])   # same name, different seq
                   for n in new_by_name if n in old_by_name and old_by_name[n] != new_by_name[n]]
    name_changed = [(old_by_seq[s], new_by_seq[s])       # same seq, different name
                    for s, nn in new_by_seq.items() if s in old_by_seq and old_by_seq[s] != nn]
    added = sorted(n for n, s in new_by_name.items()     # genuinely new: neither name nor seq seen before
                   if n not in old_by_name and s not in old_by_seq)
    removed = sorted(n for n, s in old_by_name.items()   # ours, and the new set has neither this name nor this seq
                     if n not in new_by_name and s not in new_by_seq)
    unchanged = [n for n in new_by_name if n in old_by_name and old_by_name[n] == new_by_name[n]]
    return {'added': added, 'removed': removed, 'seq_changed': sorted(seq_changed),
            'name_changed': sorted(name_changed), 'unchanged': unchanged}


# ----------------------------------------------------------------------------------------
def read_locus_glfo(gldir, locus, template_glfo=None, sanitize_names=False, debug=False):
    if not os.path.exists('%s/%s' % (gldir, locus)):
        return None
    # <template_glfo> (the old set) lets read_glfo infer codon positions for any new genes
    # that lack them (get_missing_codon_info); it's harmless when extras.csv is present.
    return glutils.read_glfo(gldir, locus, template_glfo=template_glfo,
                             add_dummy_name_components=sanitize_names, debug=debug)


# ----------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--species', required=True, help='used for default --old-dir/--out-dir and for messages')
    parser.add_argument('--new-dir', required=True, help='partis-layout germline dir of the freshly-downloaded set (<dir>/<locus>/<locus>{v,d,j}.fasta [+ extras.csv])')
    parser.add_argument('--old-dir', help='existing set to merge into (default data/germlines/<species>)')
    parser.add_argument('--out-dir', help='where to write the merged set (default = --old-dir, i.e. update in place). Only written with --write')
    parser.add_argument('--loci', default='igh:igk:igl', help='colon-separated loci to process')
    parser.add_argument('--expected', required=True, choices=['augmentation', 'replacement'],
                        help="where the new set comes from (sets which name wins on a conflict): 'augmentation' = genes from a DIFFERENT source (e.g. engelbrecht/rodriguez on top of imgt), keep our old names; 'replacement' = a newer version of the SAME source (e.g. a fresh ogrdb pull), take the new names, and warn about old genes missing from it")
    parser.add_argument('--sanitize-names', action='store_true', help="the new set's fasta ids aren't valid gene names (e.g. arbitrary ids from a raw download): rename them with add_dummy_name_components. Not needed if the ids are already valid gene names. Codon positions are always inferred from --old-dir when the new set lacks an extras.csv, regardless of this flag.")
    parser.add_argument('--write', action='store_true', help='actually write the merged glfo (default: dry run -- print the summary + warnings only)')
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()

    args.old_dir = utils.non_none([args.old_dir, 'data/germlines/%s' % args.species])
    args.out_dir = utils.non_none([args.out_dir, args.old_dir])
    args.loci = utils.get_arg_list(args.loci)

    any_warning = False

    for locus in args.loci:
        print('%s' % utils.color('blue', locus, width=6))
        old_glfo = read_locus_glfo(args.old_dir, locus, debug=args.debug)
        if old_glfo is None:
            print('  no old %s dir under %s, skipping' % (locus, args.old_dir))
            continue
        new_glfo = read_locus_glfo(args.new_dir, locus, template_glfo=old_glfo,
                                   sanitize_names=args.sanitize_names, debug=args.debug)
        if new_glfo is None:
            print('  %s no new %s dir under %s (old set has %d %s genes that would be untouched)'
                  % (utils.wrnstr(), locus, args.new_dir, sum(len(old_glfo['seqs'][r]) for r in utils.regions), locus))
            continue

        for region in utils.regions:
            if region == 'd' and not utils.has_d_gene(locus):  # dummy d gene injected for loci without d -- nothing to report
                continue
            cats = classify_region(old_glfo, new_glfo, region)
            nold, nnew = len(old_glfo['seqs'][region]), len(new_glfo['seqs'][region])
            n_merged = nold + len(cats['added'])
            print('  %s: old %d, new %d  ->  merged %d   (+%d added, %d old-only, %d seq-changed, %d name-changed)'
                  % (region, nold, nnew, n_merged, len(cats['added']), len(cats['removed']),
                     len(cats['seq_changed']), len(cats['name_changed'])))
            if cats['added']:  # (the seq-changed and name-changed genes are listed by get_merged_glfo below)
                print('      added (%d): %s' % (len(cats['added']), utils.color_genes(cats['added'])))
            if cats['removed']:
                if args.expected == 'replacement':  # a re-download shouldn't be missing much, so warn loudly
                    any_warning = True
                    print('      %s old-only (%d): expected a replacement, but these are missing from the new set -- confirm real upstream removals, not a bad download (kept by the union merge): %s'
                          % (utils.wrnstr(), len(cats['removed']), utils.color_genes(cats['removed'])))
                else:  # augmentation from a different source: missing genes are normal, just record them
                    print('      old-only (%d, in the old set but not the new one, kept by the union merge): %s' % (len(cats['removed']), utils.color_genes(cats['removed'])))

        # name/seq precedence: for an augmentation (new genes from a different source) we keep our
        # established names; for a replacement (re-download of a newer version of the same source)
        # the new set is authoritative, so its names/seqs win. get_merged_glfo keeps glfo_a's, so we
        # pass whichever should win as the first arg.
        if args.expected == 'augmentation':
            merged_glfo, _ = glutils.get_merged_glfo(old_glfo, new_glfo, debug=args.debug, glfo_labels=('old', 'new'))
        else:
            merged_glfo, _ = glutils.get_merged_glfo(new_glfo, old_glfo, debug=args.debug, glfo_labels=('new', 'old'))
        if args.write:
            glutils.write_glfo(args.out_dir, merged_glfo, debug=args.debug)
            print('  wrote merged %s to %s' % (locus, args.out_dir))

    if not args.write:
        print('\n%s dry run (no --write): merged glfo not written' % utils.color('yellow', 'note'))
    if any_warning:
        print('%s see the warnings above before committing this update' % utils.wrnstr())


# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
