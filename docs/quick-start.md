### Quick start

There are a lot of different partis [subcommands](subcommands.md), but probably the first thing you want to do is run `partition` on a fasta input file `/path/to/yourseqs.fa`.
This will group the sequences into clonal families, and then annotate each family with V gene, naive sequence, etc.
Assuming you don't pass it a parameter dir, it will also first infer parameters, including germline inference.

```/path/to/<partis_dir>/bin/partis partition --infname /path/to/yourseqs.fa --outfname /path/to/yourseqs-partition.yaml```.

Note that all input must be plus strand sequences, and if this isn't human igh you must set the `--species {human,mouse,macaque}` and/or `--locus {tra,trb,trd,trg,igl,igk,igh}` options.
If you're using Docker, and you mounted your host filesystem as described [here](install.md#installation-with-docker), you should replace `/path/to` with the appropriate host mount point within Docker.
To parallelize on your local machine, add `--n-procs N`; to paralellize over many machines, the slurm and sge batch systems are currently supported (details [here](parallel.md)).

Typically, you can expect to annotate 10,000 sequences on an 8-core desktop in about five minutes, and partition in 25 minutes.

In addition to any output files specified with `--oufname` (described [here](output-formats.md)), partis writes to two directories on your file system.
Temporary working files go in `--workdir`, which is entirely removed upon successful completion.
The workdir defaults to a subdirectory of `/tmp` (`/tmp/$USER/hmms/<random.randint>`), and this default shouldn't need to be changed unless you're using a batch system to run on multiple machines, in which case it needs to be on a network mount that they can all see.
Permament parameter files are written to `--parameter-dir`, which defaults to a subdirectory of the current directory (see [here](subcommands.md#cache-parameters)).

You can add `--debug 1` (or `--debug 2`) to print a lot of additional information about what's going on.
The output of these should usually be viewed with `less -RS` either directly by piping `|` or after redirecting `>` to a log file (S disables line wrapping -- use left/right arrows to move side-to-side).

For details on the large number of available partition options, run `partis partition --help`.

A variety of overview plots will be written to disk if you set `--plotdir <plotdir>`. Details on their content can be found [here](plotting.md).

If it's taking too long, or using too much memory, try the suggestions [here](subcommands.md#partition), [here](parallel.md) and [here](https://groups.google.com/forum/#!topic/partis/1IEfLapbStw).

After you've partitioned your sample, you might want to view an ascii-art representation of the resulting clusters and annotations with [view-output](subcommands.md#view-output), or calculate selection metrics to predict affinity with [get-selection-metrics](subcommands.md#get-selection-metrics).
You might also want to use the [linearham](https://github.com/matsengrp/linearham/) package for accurate Bayesian infererence of trees and naive sequences.
And for rich, browser-based visualization of families, trees, and annotations we recommend [Olmsted](https://github.com/matsengrp/olmsted/).
