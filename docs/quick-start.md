### Quick start

Once you have partis installed, to annotate a file `/path/to/yourseqs.fa` with BCR sequences you'd run

```./bin/partis annotate --infname /path/to/yourseqs.fa --outfname /path/to/yourseqs-annotate.csv```.

To separate them into clonal families, replace `annotate` with `partition`.
Note that all input must be plus strand sequences.
If you're using Docker, and you mounted your host filesystem as described above, you should replace each `/path/to` with `/host/path/to`.
To parallelize on your local machine, add `--n-procs N`.
The slurm and sge batch systems are also supported (see [below](#parallelization)).
There are also several [faster](#faster-methods) partitioning methods.

Typically, you can expect to annotate 10 thousand sequences on an 8-core desktop in about five minutes, and partition in 25 minutes.
Note, though, that partitioning run times depend heavily on the lineage structure of your sample (for instance a sample with fewer, larger clonal families is usually slower than a more diverse sample).
Memory is not usually a limiting factor, although very large samples (millions of sequences) with very large clonal families can provide exceptions.
In such cases, running on subsamples using `--n-random-queries N` or `--n-max-queries N` is frequently the best option.

By default, it assumes the sequences are human igh.
To change this, use the `--species {human,mouse}` and `--locus {tra,trb,trd,trg,igl,igk,igh}` options.

There are also some example sequences you can run on in `test/example.fa`.

In addition to any output files, partis writes to two directories on your file system.
Temporary working files go in `--workdir`, which is entirely removed upon successful completion.
This defaults into a directory under `/tmp`, and shouldn't need to be modified unless you're running on multiple nodes (see below), in which case it needs to be on a network mount they can all see.
Permament parameter files, on the other hand, are written to `--parameter-dir`, which defaults to a subdir of the current directory (see below).

###### Text output and logging

  - to properly display the ansi color codes in stdout:
    - less -RS (S disables line wrapping -- use left/right arrows to move side-to-side)
    - M-x display-ansi-colors (emacs)
  - it's typically best to view the ascii annotation outputs by making your text size as small as necessary to fit it all on your screen (typically ctrl-<minus>). You don't usually need to view individual bases, and you can zoom back in for the rest of the stdout.
