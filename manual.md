### Introduction

Partis is an HMM-based framework for B- and T-cell receptor annotation, simulation, and partitioning.
It is built on top of the [ham](https://github.com/psathyrella/ham) HMM compiler, and also uses the [ig-sw](https://github.com/matsengrp/ig-sw) set of Smith-Waterman annotation tools.
Partis is free software under the GPL v3.

The following two papers describe the annotation and clonal family inference functionality of partis, respectively:

* Ralph, DK, & Matsen IV, FA (2016). [Consistency of VDJ Rearrangement and Substitution Parameters Enables Accurate B Cell Receptor Sequence Annotation.](http://doi.org/10.1371/journal.pcbi.1004409) *PLOS Computational Biology*, 12(1), e1004409.
* Ralph, DK, & Matsen IV, FA (2016). [Likelihood-based Inference of B-cell Clonal Families.](http://arxiv.org/abs/1603.08127) In press, *PLOS Computational Biology*.

This manual is organized into the following sections:

  * [Installation with Docker](#installation-with-docker)
  * [Installation from scratch](#installation-from-scratch)
  * [Quick start](#quick-start)
  * [Subcommands](#subcommands) how to navigate the various `partis` actions
    - [annotate](#annotate) find most likely annotations
	- [partition](#partition) cluster sequences into clonally-related families
	  - [faster methods](#faster-methods)
	  - [cluster annotations](#cluster-annotations)
	- [view-annotations](#view-annotations) Print (to std out) the annotations from an existing annotation output csv.
	- [view-partitions](#view-partitions)  Print (to std out) the partitions from an existing partition output csv.
	- [cache-parameters](#cache-parameters) write parameter values and HMM model files for a given data set
	  - [germline sets](#germline-sets)
	- [simulate](#simulate) make simulated sequences
	- [run-forward](#run-forward) find total probability of sequences
  * [Parallelization](#parallelization)
  * [Text output and logging](#text-output-and-logging)

There are also many flags and optional parameters; unless mentioned below these are beyond the scope of this manual.
Details concerning their purpose, however, may be gleaned by means of the following incantation: `./bin/partis --help`.
In general, we assume that the reader is familiar with the papers describing [annotation](http://arxiv.org/abs/1503.04224) and [clustering](https://arxiv.org/abs/1603.08127) with partis.

To ask questions, or search through past discussions, please use the [google group](https://groups.google.com/forum/#!forum/partis).
For specific issues with the software, e.g. bug reports or feature requests, on the other hand, [submit an issue](https://github.com/psathyrella/partis/issues?utf8=%E2%9C%93&q=) on github.

### Installation with Docker

The easiest way to install partis is with the [Docker image](https://registry.hub.docker.com/u/psathyrella/partis/).
Docker images are kind of like lightweight virtual machines, and as such all the dependencies are taken care of automatically.
If, however, you'll be mucking about under the hood, or you just don't want to deal with Docker, plain installation might be preferable (see [Installation from scratch](#installation-from-scratch) below).

You'll first want install Docker using their [installation instructions](https://docs.docker.com) for your particular system.
Once Docker's installed, pull the partis image from dockerhub, start up a container from this image and attach yourself to it interactively, and compile:

```
sudo docker pull psathyrella/partis
sudo docker run -it -v /:/host psathyrella/partis /bin/bash
./bin/build.sh
```
Depending on your system, the `sudo` may be unnecessary.
Note the `-v`, which mounts the root of the host filesystem to `/host` inside the container.

To detach from the docker container without stopping it (and you don't want to stop it!), hit `ctrl-p ctrl-q`.

###### Docker tips
Docker containers and images are kinda-sorta like virtual machines, only different, so a few things:
  - We use `docker run` above: this creates a new container from (i.e. a new instance of) the partis image
  - If you exit (ctrl-d or `exit`) and then do `docker run` again, that'll create a new container. But most of the time you want to reattach to the one you made before.
    - to reattach to the same container (after detaching with `ctrl-p ctrl-q`):
      - run `docker ps -a` (lists all running and stopped containers) to get the right container ID
      - run `docker attach <ID>`
    - Hence the `-it` and `/bin/bash` options we used above for `docker run`: these allocate a pseudo-tty, keep STDIN open, and run bash instead of the default command, without all of which you can't reattach
    - the Docker docs are good, but googling on stackoverflow is frequently better

### Installation from scratch

Given the wide variety of hardware and operating systems, installing without Docker comes with no guarantees.
That said, you should be able to get partis running with just a few `apt-get` and `pip` commands, at least on recent versions of ubuntu.
We haven't tried on os-x, but if you do, let us know -- especially if you find any tricks.
This docker-free approach is frequently a lot faster (if you already have numpy and scipy installed, for instance, you don't have to wait for docker to compile and install them from scratch).
You also don't have to deal with the additional complications of being inside docker, perhaps most importantly docker's likely incompatibility with your batch queueing system.

Basically, you just need to install a few extra packages in addition to the dependencies listed in the [Dockerfile](https://github.com/psathyrella/partis/blob/master/Dockerfile):

```
sudo apt-get install python-pip cmake scons libgsl0-dev libncurses5-dev libxml2-dev libxslt1-dev mafft r-base
pip install --user numpy scipy matplotlib pandas biopython dendropy==3.12.3 pysam pyyaml seaborn colored_traceback
R --vanilla --slave -e 'install.packages("TreeSim", repos="http://cran.rstudio.com/")'  # optional -- only used for simulation
```

If you need to sort out versions (unlikely, especially if you start from 14.04/16.04, since that's what we're running), follow the Dockerfile chain beginning [here](https://registry.hub.docker.com/u/psathyrella/partis/dockerfile/) and [here](https://github.com/matsengrp/dockerfiles/blob/master/cpp/Dockerfile).

Once you've got all the necessary things on your system, you can proceed to clone the repository and compile:

```
git clone git@github.com:psathyrella/partis
cd partis
./bin/build.sh
```

### Quick start

Once you have partis installed, to annotate a file `/path/to/yourseqs.fa` with BCR sequences you'd run

```./bin/partis annotate --infname /path/to/yourseqs.fa --outfname /path/to/yourseqs-annotate.csv```.

To separate them into clonal families, replace `annotate` with `partition`.

(If you're using Docker, and you mounted your host filesystem as described above, you should replace each `/path/to` with `/host/path/to`.)

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

Besides any output files, partis writes to two directories on your file system.
Temporary working files go in `--workdir`, which is entirely removed upon successful completetion.
This defaults into a random subdir under `/tmp`, and shouldn't need to be modified unless you're running on multiple nodes (see below), in which case it needs to be a network mount they can all see.
Permament parameter files are written to `--parameter-dir`, which defaults to a subdir of the current directory (see below).

### Subcommands

The main script has a number of actions:

```./bin/partis annotate | partition | view-annotations | view-partitions | cache-parameters | simulate | run-forward```,

each of which is described in more detail below.
For more information you can also type `./bin/partis --help` and `./bin/partis <subcommand> --help`.

The fist step in all cases is to infer a set of parameters particular to the input sequence file.
These are written to `--parameter-dir`, and then used for all subsequent runs.
If you don't specify `--parameter-dir`, it defaults to a location in the current directory that amounts to a slight bastardization of your input file path (parameters for `path/to/seqs.fa` will go in `_output/path_to_seqs/`).
This default is designed such that with typical workflows, if your input files have different paths, their parameters will go in different places.

That said, the consequences of using the wrong parameter directory for a set of sequences are potentially dire.
So if you're doing any monkey business, you need to be aware of where partis is telling you that it's putting parameters (it's printed to stdout).
For instance, if you run with one set of sequences in an input file, and then toss some **other** sequences into the same file, partis won't know anything about it, and will use the same (now-inappropriate) parameters.

If `--parameter-dir` (whether explicitly set or left as default) doesn't exist, partis assumes that it needs to cache parameters, and does that before running the requested action.

Whether caching parameters or running on pre-existing parameters, the hmm needs smith-waterman annotations as input.
While this preliminary smith-waterman step is fairly fast, it's also easy to cache the results so you only have to do it once.
By default these smith-waterman annotations are written to a csv file in `--parameter-dir` during parameter caching.
The default filename is a hash of the concatenated input sequence id strings
(Because all sequences need to be aligned and padded to the same length before partititioning, the smith-waterman annotation information for each sequence depends slightly on all the other sequences in the file, hence the hash.)
These defaults should ensure that with typical workflows, smith-waterman only runs once.
If however, you're doing less typical things (running on a subset of sequences in the file), if you want smith-waterman results to be cached you'll need to specify `--sw-cachefname` explicitly, and it'll write it if it doesn't exist, and read from it if it does.

#### annotate

Finds the Viterbi path (i.e., the most likely annotation/alignment) for each sequence, for example:

```./bin/partis annotate --infname test/example.fa --outfname _output/example.csv```

The output csv headers are listed in the table below, and you can view a colored ascii representation of the rearrangement events with the `view-annotations` action.
An example of how to parse this output csv (say, if you want to further process the results) is in `bin/example-output-processing.py`.

All columns listed as "colon-separated lists" are trivial/length one for single hmm, and contain actual colons for multi-hmm.

|   column header        |  description
|------------------------|----------------------------------------------------------------------------
| unique_ids             |  colon-separated list of sequence identification strings
| v_gene         |  V gene in most likely annotation
| d_gene         |  D gene in most likely annotation
| j_gene         |  J gene in most likely annotation
| cdr3_length    |  CDR3 length of most likely annotation (IMGT scheme, i.e. including both codons in their entirety)
| mut_freqs      |  colon-separated list of sequence mutation frequencies
| input_seqs     |  colon-separated list of input sequences, with constant regions (fv/jf insertions) removed
| naive_seq      |  naive (unmutated ancestor) sequence corresponding to most likely annotation
| v_3p_del       |  length of V 3' deletion in most likely annotation
| d_5p_del       |  length of D 5' deletion in most likely annotation
| d_3p_del       |  length of D 3' deletion in most likely annotation
| j_5p_del       |  length of J 5' deletion in most likely annotation
| v_5p_del       |  length of an "effective" V 5' deletion in the most likely annotation, corresponding to a read which does not extend through the entire V segment
| j_3p_del       |  length of an "effective" J 3' deletion in the most likely annotation, corresponding to a read which does not extend through the entire J segment
| vd_insertion       |  sequence of nucleotides corresponding to the non-templated insertion between the V and D segments
| dj_insertion       |  sequence of nucleotides corresponding to the non-templated insertion between the D and J segments
| fv_insertion       |  constant region on the 5' side of the V (accounts for reads which extend beyond the 5' end of V)
| jf_insertion       |  constant region on the 3' side of the J (accounts for reads which extend beyond the 3' end of J)
| mutated_invariants |  true if the conserved codons corresponding to the start and end of the CDR3 code for the same amino acid as in their original germline (cyst and tryp/phen, in IMGT numbering)
| in_frames          |  true if the net effect of VDJ rearrangement and SHM indels leaves both the start and end of the CDR3 (IMGT cyst and tryp/phen) in frame with respect to the start of the germline V sequence
| stops              |  true if there's a stop codon in frame with respect to the start of the germline V sequence
| v_per_gene_support |  approximate probability supporting the top V gene matches, as a semicolon-separated list of colon-separated gene:probability pairs (approximate: monotonically related to the actual probability, but not exactly one-to-one)
| d_per_gene_support |  approximate probability supporting the top D gene matches, as a semicolon-separated list of colon-separated gene:probability pairs (approximate: monotonically related to the actual probability, but not exactly one-to-one)
| j_per_gene_support |  approximate probability supporting the top J gene matches, as a semicolon-separated list of colon-separated gene:probability pairs (approximate: monotonically related to the actual probability, but not exactly one-to-one)
| indelfos       |  colon-separated list of information on any SHM indels that were inferred in the Smith-Waterman step. Written as a literal python dict; can be read in python with `ast.literal_eval(line['indelfo'])`
| indel_reversed_seqs  |  colon-separated list of input sequences with indels "reversed", and with constant regions (fv/jf insertions) removed. Empty string if there are no indels, i.e. if it's the same as 'input_seqs'
| duplicates     |  colon-separated list of "duplicate" sequences for each sequence, i.e. sequences which, after trimming fv/jf insertions, were identical and were thus collapsed.

Note that `utils.process_input_line()` and `utils.get_line_for_output()` can be used to automate input/output (see for example `bin/example-output-processing.py`).

Annotation is the algorithm of choice for annotating sequences where the clonal relationship is different i.e. no sequence in the dataset are from the same germinal center, and therefore are not related by having the same naive sequence. Examples of such datasets could be pooled datasets with BCR sequences from many individuals, where clonal relationship cannot be present.

However for many applications sequence data is created unspecifically for a large amount of BCRs and will contain many sequences being from the same germinal center, hence also sharing the same naive sequence. Using this prior knowledge can greatly improve inference of VDJ gene combination and reconstruction of the naive sequence, and therefore when datasets allow for partitioning, the annotations from the partitioning algorithm should be preferred over the single-sequence annotation results.

#### partition

Example invocation:

```./bin/partis partition --infname test/example.fa --outfname _output/example.csv```

The output csv file headers are listed in the table below, and you can view a colored ascii representation of the clusters with the `view-partitions` action.
We write one line for the most likely partition (with the lowest logprob), as well as a number of lines for the surrounding less-likely partitions (set with `--n-partitions-to-write`)
An example of how to parse this output csv (say, if you want to further process the results) is in `bin/example-output-processing.py`.

|   column header  |  description
|------------------|----------------------------------------------------------------------------
| logprob          |  Total log probability of this partition
| n_clusters       |  Number of clusters (clonal families in this partition)
| partition        |  String representing the clusters, where clusters are separated by ; and sequences within clusters by :, e.g. 'a:b;c:d:e'
| n_procs          |  Number of processes which were simultaneously running for this clusterpath. In practice, final output is usually only written for n_procs = 1

By default, this uses the most accurate and slowest method: hierarchical agglomeration with, first, hamming distance between naive sequences for distant clsuters, and full likelihood calculation for more similar clusters.
Like most clustering algorithms, this scales rather closer than you'd like to quadratically than it does to linearly.
We thus also have two faster and more approximate methods.

##### faster methods

*--seed-unique-id*

Only looks for sequences which are clonally related to the sequence id specified by `--seed-unique-id` (you can also specify `--seed-seq`).
Because this completely avoids the all-against-all comparisons which render full data set partitioning depressingly close to quadratic in sample size, it is vastly faster
It is an excellent option if you have, ahead of time, sequences about which you are particularly interested.
You can also combine it with, for instance `--naive-vsearch`, by seeding with sequences from particularly interesting (probably large) clusters from that approximate method.

*--naive-vsearch*

First calculates naive (unmutated ancestor) for each sequence, then passes these into vsearch for very fast, very heuristic clustering.
The initial naive sequence calculation is easy to parallelize, so is fast if you have access to a fair number of processors.
Vsearch is also very fast, because it makes a large number of heuristic approximations to avoid all-against-all comparison, and thus scales significantly better than quadratically.
With `--n-procs` around 10 for the vsearch step, this should take only of order minutes for a million sequences.
If your sample is too large to comfortably partition with the default method using your computational resources, one option is to run naive vsearch clustering, and use clusters from this step in which you're particularly interested to seed `--seed-unique-id`.

*ignore smaller clusters*

Since in many cases we are most interested in the largest clonal families, you can also specify that, after a few partitition steps, smaller families are discarded.
The idea is that any large clusters will accumulate appreciable size within the first few partition steps, and since most real repertoires are dominated by small clusters (especially singletons), we can learn what we need from these larger clusters without spending time trying to accurately infer the details of the very small clusters.
What counts as a "small" cluster is controlled by `--small-clusters-to-ignore`, which is either a colon-separated list of clusters sizes (e.g. `1:2:3`) or an inclusive range (e.g. `1-10`).
The number of steps after which these small clusters are removed is set with `--n-steps-after-which-to-ignore-small-clusters` (default 3).

##### cluster annotations

If --outfname is set, in addition to the clusters in that file, the most likely annotation for each final cluster are written to --cluster-annotation-fname (default `<--outfname>.replace('.csv', '-cluster-annotations.csv')`).

To annotate an arbitrary collection of sequences using simultaneous multi-HMM inference (which is much, much more accurate than annotating the sequences individually), you can combine the `--queries` and `--n-simultaneous-seqs` arguments.
For instance, if you knew from partitioning that three sequences `a`, `b`, and `c` were clonal, you could run:

``` ./bin/partis annotate --infname in.fa --queries a:b:c --n-simultaneous-seqs 3 --outfname abc-annotation.csv```

In order to get an idea of the uncertainty on a given cluster's naive sequence, you can specify `--calculate-alternative-naive-seqs` during the partition step.
This will save all the naive sequences for intermediate sub-clusters to a cache file so that, afterwards, you can view the alternative naive sequences for the sub-clusters.
For instance:

```
./bin/partis partition --infname test/example.fa --outfname _output/example.csv --calculate-alternative-naive-seqs
./bin/partis view-alternative-naive-seqs --outfname _output/example.csv --queries <queries of interest>  # try piping this to less by adding "| less -RS"
```

if you don't specify `--queries` in the second step, it's ok, since it will print the partitions in the output file before exiting, so you can copy and paste a cluster into the `--queries` argument.

##### cpu and memory usage
Because, at least to a first approximation, accurate clustering entails all-against-all comparison, partitioning is in a fundamentally different computational regime than are single-sequence problems such as annotation.
In order to arrive at a method that can be useful in practice, we have tried to combat this inherent difficulty with a number of different levels of both approximations and caching, several of which are described in the papers.
As in many such cases, this frequently amounts to an attempt to make judicious compromises between cpu and memory usage which are appropriate in each individual circumstance.
Because the computational difficulty of the clustering problem is entirely dependent on the detailed structure of each repertoire, however, it is not always possible to make the optimal choice ahead of time.
For instance, five hundred thousand sequences from a repertoire that consists almost entirely of a single, highly-mutated clone will have very different computational requirements (more!) to one which is more evenly spread among different naive rearrangements and has more typical mutation levels.

Which is a roundabout way of saying: if you find that the sample which you'd like to run on is taking forever with the number of cores you have available, or if it's exhausting your available memory, you have a few options.
First, consider whether it's possible to run separately on independent subsamples of your full sample -- if you divide into ten subsets, each one will typically run a hundred times faster.
You can then, for instance, pull from these ten sub-partitions only the clusters you're interested in (e.g. the biggest ones) and run only on those with the `--queries` option.
Another option, as mentioned above, is to run the very fast but somewhat less accurate `--naive-vsearch` method on the full sample, and then either pass only the most interesting clusters into the full algorithm (using `--queries`), or use `--seed-unique-id` with some sequences from the most interesting clusters.

We're working to include these two strategies as command line options, but want to do some more validation before we settle on the particulars which are best for a general use case.

#### view-annotations

To, e.g. run on the output csv from the `annotate` action:

``` ./bin/partis view-annotations --outfname annotate-output.csv```

#### view-partitions

To, e.g. run on the output csv from the `partition` action:

``` ./bin/partis view-partitions --outfname partition-output.csv```

#### cache-parameters

This is run automatically if `--parameter-dir` doesn't exist (whether this directory is specified explicitly, or is left as default).
So you do not, generally, need to run it on its own.

When presented with a new data set, the first thing we need to do is infer a set of parameters, a task for which we need a preliminary annotation.
As such, partis first runs ig-sw's Smith-Waterman algorithm on the data set.
The smith-waterman annotations are used to build and write out a parameter set, which is in turn used to make a set of HMM model files for each observed allele.
These files are then passed as input to a second, HMM-based, annotation step, which again outputs (more accurate) parameter values and HMM model files.

For example:

``` ./bin/partis cache-parameters --infname test/example.fa --parameter-dir _output/example```

When caching parameters, the parameter csvs from Smith-Waterman and the HMM are put into `/sw` and `/hmm` subdirectories of `--parameter-dir`.
Within each of these, there are a bunch of csv files with (hopefully) self-explanatory names, e.g. `j_gene-j_5p_del-probs.csv` has counts for J 5' deletions subset by J gene.
The hmm model files go in the `hmms` subdirectory, which contains yaml HMM model files for each observed allele.

##### germline sets

By default partis uses the set of germline V, D, and J genes (the "germline set") in `data/germlines`, which is from imgt.
If you have another set you'd like to use, you can do so by setting the `--initial-germline-dir` option to cache-parameters.
A number of studies have shown that, at least for V genes, this set (and all others) is both missing many alleles which occur in real populations and contains many alleles which do not occur in any real sample.
In addition, any individual sample will contain only a small fraction of the genes and alleles in the default germline set.

By default partis handles these problems when caching parameters with a two step process during the preliminary smith-waterman annotation.
The V alleles in the initial germline set which were the best match for at least one sequence, and which accounted for more than `--min-allele-prevalence-fraction` (default 0.0005) of the repertoire, are then divided into classes, such that all alleles within each class are the same length, and separated by less than `--n-max-snps` (default 8) SNPs.
To a first approximation, these classes are the same as imgt gene designations, i.e. the bit before the `*`.
Within each class, we then keep only the most common `--n-alleles-per-gene` alleles (default 2) which are distinguishable from each other in expressed samples (e.g. which don't differ by one base at codon 108).
This amounts to the application of a strong diploid prior.
Which, of course, is not particularly well-justified given the extensive gene deletion and duplication in the BCR locus, and the fact that typical SHM levels are much larger than the typical inter-V-gene distance in the imgt repertoire.
But, it is a reasonable first approximation.
And, more importantly, it is vastly more accurate than simply keeping every single allele which is a best match.

Second, if `--find-new-alleles` is specified, we re-run smith-waterman using this reduced V germline set, and apply a new method for finding any previously unknown V alleles, or V alleles which were over-aggressively pruned in the previous step.
This method is too complex to describe here in detail, but can perhaps best be described as a principled, hypothesis-testing-based generalization of the [tigger](https://www.ncbi.nlm.nih.gov/pubmed/25675496) method.
The entire resulting germline set is written to a subdirectory of `--parameter-dir`, and is used in subsequent runs instead of the initial germline set.

We could probably do the same thing for J, but there doesn't seem to be much polymorphism, so it's probably not worthwhile.
Doing it for D is probably not realistic.

You can disable the allele-removal step with `--dont-remove-unlikely-alleles`.

#### simulate

For testing purposes, we can use the parameters in `--parameter-dir` to create simulated sequences that mimic the data as closely as possible.
This allows us to test how well our algorithms work on a set of sequences for which we know the correct annotations and clonal relationships.
Note that while the parameters describe a very detailed picture of the rearrangement and mutation patterns, we do not (yet!) parametrize either the emipical clonal size distributions or phylogenetic relationships.
These can, however, be specified from a number of realistic options (see `--help` for the simulate subcommand).

For example:

```./bin/partis simulate --outfname _output/examle-simu.csv --parameter-dir test/reference-results/test/parameters/data --n-sim-events 50```.

This will spit out simulated sequences to `--outfname` using some test parameters in `test/reference-results`.
We also specify that we want it to simulate 50 rearrangement events.
To get the actual number of sequences, we multiply this by the mean number of leaves per tree (set with `--n-leaves`, at the moment with a default of five).
At the start of a simulation run, TreeSim generates a set of `--n-trees` trees (default 500 at the moment), and each tree has a number of leaves drawn from an exponential (or zipf, or box, or...) with mean/exponent `--n-leaves`.
Throughout the run, we sample a tree at random from this set for each rearrangement event.

### Parallelization

###### In general

The number of processes on your local machine is set with `--n-procs N`.

In order to parallelize over more processes than the local machine can handle, we currently support slurm and sge: specify one or the other with `--batch-system`.
The default options for each should work, but if you need to add extras (for instance to reserve particular memory requirements) use `--batch-options`, e.g. `--batch-options="--foo bar"`.
Note that you should not specify stdout/stderr locations (`-e` or `-o`) for sge -- we add these options automatically because we need to be able to parse them for each job.
By default, partis writes its temporary files to a working directory of the form `/tmp/$USER/hmms/$RANDOM`.
If you're running on a batch system, though, you need the working directory to be a network mount that every node can see; set this with `--workdir`, e.g. `--workdir /path/to/nfs/$USER/hmms/$RANDOM`.

###### annotate

Sequence annotation lends itself quite readily to independent parallelization.
It should take 0.1-1 second per sequence; if you want it to go faster, just increase `--n-procs`.

###### partition

Clustering, however, really doesn't lend itself at all to independent parallelization -- we need to, at least approximately, compare each sequence to every other one.
For the full and point partis methods, we get around this by starting with `--n-procs` processes.
The input sequences are split evenly among these, and each process does all-against-all comparison (with many optimizations to avoid the full likelihood calculation) of all of its allotted sequences.
The results of this first round are collected and merged together, and then reapportioned among a new, smaller, number of processes.
This is continued until we arrive at one final process which is comparing all sequences.
Since at each stage we cache every calculated log probability, while the later steps have more sequences to compare, they also have more cached numbers at their disposal, and so it's possible to make each step take about the same amount of time.
We currently reduce the number of processes by about 1.6 at each step, as long as the previous step didn't have to calculate too many numbers.

With typical mutation levels, lineage structures, and cluster size distributions (all of which strongly affect clustering time), it's currently best to start with `--n-procs` set so you have about 300 sequences per process.

#### Text output and logging

  - to properly display the ansi color codes in stdoutput:
    - less -RS (S disables line wrapping -- use left/right arrows to move side-to-side)
    - M-x display-ansi-colors (emacs)
  - it's typically best to view the ascii annotation outputs by making your text size as small as necessary to fit it all on your screen (typically ctrl-<minus>). You don't usually need to view individual bases, and you can zoom back in for the rest of the stdout.
