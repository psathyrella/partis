[Up to table of contents](contents.md)

### Quick start

While there are many different partis [subcommands](subcommands.md), likely the first thing you want to do is run `partition` on a fasta input file `/path/to/yourseqs.fa`.
The following command will first infer a set of parameters (including germline genes) on the sample, then group sequences into clonal families and annotate each family with V gene, naive sequence, etc:

```partis partition --infname /path/to/yourseqs.fa --outfname /path/to/yourseqs-partition.yaml```.

If partis is not installed in your path, you will need to use the full path in this command.
This by default assumes input of single-chain positive sense human igh.
If you have another species or locus, set the `--species {human,mouse,macaque,c57bl,balbc}` (details [here](../data/germlines/README.md)) and/or `--locus {tra,trb,trd,trg,igl,igk,igh}` options.
If you have different loci mixed together, you'll need to either set `--paired-loci`, or first run `split-loci.py`, both of which split different loci into separate files (the former runs the latter in the course of handling data with pairing information).
Both of these take the argument `--reverse-negative-strands` if you have a mix of positive and negative sense sequences.
If you have heavy/light pairing information, for instance from 10x single cell data, it is important to incorporate it as described [here](paired-loci.md) in order to take advantage of (among other things) dramatically improved partitioning accuracy.
If you're using Docker, and you mounted your host filesystem as described [here](install.md#installation-with-docker), you should replace `/path/to` with the appropriate host mount point within Docker.

If you prefer, you can import partis as a module rather than running command line scripts, for instance:
```
import partis.utils as utils
glfo, annotation_list, cpath = utils.read_output(yaml_output_file_path)  # read germline info, annotations, and list of partitions from partis yaml output file
antn_dict = utils.get_annotation_dict(annotation_list)  # convert list to dict for convenient access
for cluster in cpath.best():  # loop over clusters in best partition
    annotation = antn_dict[':'.join(cluster)]  # retrieve annotation for cluster
    utils.print_reco_event(annotation)  # print ascii representation of annotation
```
While we do not have automatically generated documentation pages for the functions thereby exposed, we have had very good results asking LLMs what functions to use for a task and how to use them.

If you want it to run faster (by a factor of 5 to 10, and at the cost of some accuracy) add `--fast`, which uses [naive vsearch clustering](subcommands#--fast-synonym---naive-vsearch).
The number of processes on your local machine defaults to the number of cpus, so shouldn't need to be adjusted (with `--n-procs N`) unless you're running several things at once.
To paralellize over many machines, the slurm and sge batch systems are currently supported (details [here](parallel.md)).

Typically, you can expect to annotate 10,000 sequences on an 8-core desktop in about five minutes, and partition in 25 minutes.
If it's taking too long, or using too much memory, there are many ways to improve these by orders of magnitude by trading some accuracy, or by focusing only on specific aspects of the repertoire, described [here](subcommands.md#partition), [here](parallel.md) and [here](https://groups.google.com/forum/#!topic/partis/1IEfLapbStw).

In addition to any output files specified with `--outfname` (described [here](output-formats.md)), partis writes to two directories on your file system.
Temporary working files go in `--workdir`, which is entirely removed upon successful completion.
The workdir defaults to a subdirectory of `/tmp` (`/tmp/$USER/hmms/<random.randint>`), and this default shouldn't need to be changed unless you're using a batch system to run on multiple machines, in which case it needs to be on a network mount that they can all see.
Permament parameter files are written to `--parameter-dir`, which defaults to a subdirectory of the current directory (see [here](subcommands.md#cache-parameters)).

You can add `--debug 1` (or `--debug 2`) to print a lot of additional information about what's going on.
The output of these should usually be viewed with `less -RS` either directly by piping `|` or after redirecting `>` to a log file (`-S` disables line wrapping -- use left/right arrows to move side-to-side).

For details on the large number of available partition options, run `partis partition --help`.

A variety of overview plots will be written to disk if you set `--plotdir <plotdir>`. Details on their content can be found [here](plotting.md).

After you've partitioned your sample, you might want to view an ascii-art representation of the resulting clusters and annotations with [view-output](subcommands.md#view-output), infer [trees](subcommands.md#tree-inference), or calculate selection metrics to predict affinity with [get-selection-metrics](subcommands.md#get-selection-metrics).
You might also want to use the [linearham](https://github.com/matsengrp/linearham/) package for accurate Bayesian infererence of trees and naive sequences.
And for rich, browser-based visualization of families, trees, and annotations we recommend [Olmsted](https://github.com/matsengrp/olmsted/).
