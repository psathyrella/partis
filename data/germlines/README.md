## Updating a germline set (general recipe)

Use `merge-germline-set.py` to merge a freshly-downloaded set into an existing one. It's
always a union (we never drop a gene we already had). The two modes differ in which name
wins when both sets have the same sequence under different names:
  - augmentation: the new genes are from a *different* source (e.g. engelbrecht/rodriguez
    on top of imgt), so we keep our established (old) names.
  - replacement: the new set is a newer version of the *same* source (e.g. a fresh ogrdb
    pull), so it's authoritative -- its (new) names win, and any old genes missing from it
    get a warning (maybe real upstream removals, maybe a botched download).
  # augmentation (new genes from a different source), new set already in partis layout:
  ./data/germlines/merge-germline-set.py --species human --new-dir <partis-layout-dir> --expected augmentation --write
  # replacement (newer version of the same source):
  ./data/germlines/merge-germline-set.py --species macaque --new-dir data/germlines/ogrdb-download/macaque --expected replacement --write

`--new-dir` must be in partis layout (`<dir>/<locus>/<locus>{v,d,j}.fasta` [+ extras.csv]).
For an ogrdb download that is a single combined fasta, first split it per region/locus
(see the ogrdb download + split notes below). Codon positions for new genes are inferred
by aligning to the old set, so an extras.csv in the new set is optional; add
`--sanitize-names` only if the new set's fasta ids aren't already valid gene names.
Default is a dry run; add `--write` to write the merged set. It prints the added and
old-only gene names plus the conflicts (redirect stdout to a file for a record -- see the
per-species changelog/ dirs). Review the output and every warning before committing.

---

human: download from https://ogrdb.airr-community.org/germline_sets/Human
  - use the "AIRR-seq reference set" (not the "source set" -- this has duplicates, among other things)
  - make a dir, separate v/d/j files
  - run (separately for each locus): ./bin/partis cache-parameters --locus igh --infname test/mishmash.fa --n-max-queries 1 --sanitize-input-germlines --initial-germline-dir <ogrdb-dir>
    - input file doesn't matter/isn't used
    - add something like: `glutils.write_glfo('ogrdb-ref/ogrdb-sanitized', glfo, debug=True)` to `bin/partis` where `args.sanitize_input_germlines` gets used
  - also added genes from genomic sequencing in engelbrecht/rodriguez papers, downloaded from vdjbase

macaque: from ogrdb https://ogrdb.airr-community.org/germline_sets/Macaca%20mulatta
  - sets G00091 (IGH), G00092 (IGK), G00093 (IGL); pulled via the api, e.g.
    curl https://ogrdb.airr-community.org/api/germline/set/G00091/published/ungapped
  - staged (split per region) in ogrdb-download/macaque/, merged in with
    merge-germline-set.py --species macaque --expected replacement (issue #396)

mouse/ is merged from imgt, ogrdb c57bl, and ogrdb balbc

c57bl/ and balbc/ are just from ogrdb

to update mouse germlines, download from: https://ogrdb.airr-community.org/germline_sets/Mouse
  1) cp unaligned fastas to appropriate dir in ogrdb-downloads/
  2) split apart igh files, cp igk/igl files from germlines/mouse/
  3) merge in with merge-germline-set.py (see "Updating a germline set" above), once per mouse type
  4) after checking that everything went well, replace old mouse/ dirs with new merged dirs

leader and constant genes from imgt, with extra c genes from the flairr-seq authors
  - Any allele with a name of <gene>_FL_<number> or <gene>*N<number> is novel and not found in the IMGT database.
  - Any allele with the name <gene>-FL or *<number>_ext<number> is an extension of a documented IMGT allele. Some alleles with the first syntax have a letter signifying it is an extension of an IMGT allele but contains unique sequence from other extensions of that allele. 
  - i.e. there's tons of duplicate stuff in here and i'm not sorting through it, whatever
