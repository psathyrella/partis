### Quick start

Once you have partis installed, to annotate a file `/path/to/yourseqs.fa` with BCR sequences you'd run

```/path/to/<partis_dir>/bin/partis annotate --infname /path/to/yourseqs.fa --outfname /path/to/yourseqs-annotate.csv```.

To group them into clonal families, replace `annotate` with `partition`.
Note that all input must be plus strand sequences.
If you're using Docker, and you mounted your host filesystem as described [here](install.md), you should replace `/path/to` with the appropriate host mount point within Docker.
To parallelize on your local machine, add `--n-procs N`.
The slurm and sge batch systems are also supported (details [here](parallele.md)).

Typically, you can expect to annotate 10 thousand sequences on an 8-core desktop in about five minutes, and partition in 25 minutes.

By default, it assumes the sequences are human igh.
To change this, use the `--species {human,mouse,macaque}` and `--locus {tra,trb,trd,trg,igl,igk,igh}` options.

In addition to any output files specified with `--oufname`, partis writes to two directories on your file system.
Temporary working files go in `--workdir`, which is entirely removed upon successful completion.
The workdir defaults to a subdirectory of `/tmp` (`/tmp/$USER/hmms/<random.randint>`), and shouldn't need to be modified unless you're running on multiple nodes (see below), in which case it needs to be on a network mount they can all see.

Permament parameter files, on the other hand, are written to `--parameter-dir`, which defaults to a subdirectory of the current directory (see [here](subcommands.md#cache-parameters)).

###### Text output and logging
  - If you're looking at partis ascii art annotation output, you'll want to make sure the ansi color codes in stdout are displayed properly:
    - less -RS (S disables line wrapping -- use left/right arrows to move side-to-side)
    - M-x display-ansi-colors (emacs)
  - it's typically best to view the ascii annotation outputs by making your text size as small as necessary to fit it all on your screen (typically ctrl-<minus>). You don't usually need to view individual bases, and you can zoom back in for the rest of the stdout.
