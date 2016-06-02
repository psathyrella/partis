### Introduction

Partis is an HMM-based framework for B-cell receptor annotation, simulation, and partitioning.
It is built on top of the [ham](https://github.com/psathyrella/ham) HMM compiler, and also uses the [ighutil](https://github.com/cmccoy/ighutil) set of Smith-Waterman annotation tools.
Partis is free software under the GPL v3.

This manual is organized into the following sections:

  * [Quick Start](#quick-start) install/run with Docker
  * [Slow Start](#slow-start) install from scratch
  * [Details](#details) how to navigate the various `partis` subcommands

There are also many flags and optional parameters; unless mentioned below these are beyond the scope of this manual.
Details concerning their purpose, however, may be gleaned by means of the following incantation: `./bin/partis --help`.
In general, we assume that the reader is familiar with the papers describing [annotation](http://arxiv.org/abs/1503.04224) and [clustering](https://arxiv.org/abs/1603.08127) with partis.

### Quick Start

Because partis has a lot of dependencies, you'll likely have an easier time of it using the [Docker image](https://registry.hub.docker.com/u/psathyrella/partis/) rather than installing from scratch.
Docker images are kind of like lightweight virtual machines, and as such all the dependencies are taken care of automatically.
If, however, you'll be doing a lot of mucking about under the hood, bare and dockerless installation might be preferable.

You'll first want install Docker using their [installation instructions](https://docs.docker.com) for your particular system.
Once Docker's installed, pull the partis image from dockerhub, start up a container from this image and attach yourself to it interactively, and compile:

```
sudo docker pull psathyrella/partis
sudo docker run -it -v /:/host psathyrella/partis /bin/bash
./bin/build.sh
```
Depending on your system, the `sudo` may be unnecessary.
Note the `-v`, which mounts the root of the host filesystem to `/host` inside the container.

Now you can run individual partis commands (described [below](#details)), or poke around in the code.
If you just want to annotate a small file with BCR sequences, say on your machine at `/path/to/yourseqs.fa`, run

```./bin/partis run-viterbi --infname /host/path/to/yourseqs.fa --outfname /host/path/to/yourseqs-run-viterbi.csv```

Whereas if you'd like to separate them into clonal families, run

```./bin/partis partition --infname /host/path/to/yourseqs.fa --outfname /host/path/to/yourseqs-partition.csv```

Note that now we're inside the container, we access the fasta file at the original path on your host system, but with `/host` tacked on the front (as we specified in `docker run` above).
There's also some example sequences you can run on in `test/example.fa`.

Also note that partis by default infers its HMM parameters on the fly for each data set; while on larger samples this is substantially more accurate than using population-wide averages, if you have fewer than, say, 50 sequences, it is not a particularly meaningful exercise.
Until there exists enough public data such that it is possible to build good population-wide parameter priors, the best thing to do in such cases is to find a larger data set that you think is similar (e.g. same patient, so it has the same germline genes) to the one you're interested in, and infer parameters using that larger set.

Depending on your system, in ten minutes a single process can annotate perhaps 5000 sequences or partition a few hundred.
To parallelize on your local machine, just add `--n-procs N`.
You can also use the approximate clustering methods: point/naive (`--naive-hamming`) or vsearch (`--naive-vsearch`).
The naive-hamming method is perhaps twice as fast as the full method, while the vsearch method is much faster -- it typicallly takes longer to do the pre-annotation to get the naive sequences than it does for vsearch to run.

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

### Slow Start

As noted above, for most use cases you'll likely be happier using the [Docker image](https://registry.hub.docker.com/u/psathyrella/partis/).
But if you're doing lots of development it's nice to be able to work outside of Docker, and installing from scratch isn't really so bad.
You'll need to have recent versions of the following (to see which versions, follow the Dockerfile chain beginning [here](https://registry.hub.docker.com/u/psathyrella/partis/dockerfile/) and [here](https://github.com/matsengrp/dockerfiles/blob/master/cpp/Dockerfile)):

##### required
  - scons
  - pip
  - java 7 (for ighutil)
  - python packages
    - pysam
    - pyyaml
    - cython
    - pandas
    - seaborn

##### required only for simulator
  - dendropy
  - TreeSim

The following packages are also used by partis, but they're included as `git subtree`s in the source code, so you don't need to do anything:
  - tclap
  - yaml-cpp
  - ham
  - samtools
  - bppseqgen
  - vsearch

Once you've got all the necessary things on your system, you can proceed to clone the repository:

```
git clone git@github.com:psathyrella/partis
cd partis
```

And then build:

```
./bin/build.sh
```

### Details

#### Subcommands

`./bin/partis` runs everything, and has a number of potential actions:

```./bin/partis run-viterbi|partition|cache-parameters|simulate|run-forward```.

Each of these subcommands is described in detail below.
For more information you can also type `./bin/partis --help` or `./bin/partis <subcommand> --help`.

The fist step in all cases is to infer a set of parameters particular to the input sequence file.
These are cached as csv and yaml files in `--parameter-dir`, which defaults to a location in the current directory.
You only need to cache these once, then can use them for all subsequent runs.
If `--parameter-dir` (whether explicitly set or left as default) doesn't exist, partis assumes that it needs to cache parameters, and does that before running the requested action.

##### `run-viterbi`: find most likely annotations/alignments

Finds the Viterbi path (i.e., the most likely annotation/alignment) for each sequence, for example:

```./bin/partis run-viterbi --infname test/example.fa --outfname _output/example.csv```

The output csv headers are listed in the table below, and you can view a colored ascii representation of the rearrangement events with the `view-annotations` action.

|   column header        |  description
|------------------------|----------------------------------------------------------------------------
| unique_ids             |  colon-separated list of sequence identification strings (of length 1 if multi-hmm isn't used)
| v_gene         |  V gene in most likely annotation
| d_gene         |  D gene in most likely annotation
| j_gene         |  J gene in most likely annotation
| cdr3_length        |  CDR3 length of most likely annotation (IMGT scheme, i.e. including both codons in their entirety)
| seqs           |  colon-separated list of input sequences (of length 1 if multi-hmm isn't used)
| naive_seq      |  naive (unmutated ancestor) sequence corresponding to most likely annotation
| v_3p_del       |  length of V 3' deletion in most likely annotation
| d_5p_del       |  length of D 5' deletion in most likely annotation
| d_3p_del       |  length of D 3' deletion in most likely annotation
| j_5p_del       |  length of J 5' deletion in most likely annotation
| v_5p_del       |  length of an "effective" V 5' deletion in the most likely annotation, corresponding to a read which does not extend through the entire V segment
| j_3p_del       |  length of an "effective" J 3' deletion in the most likely annotation, corresponding to a read which does not extend through the entire J segment
| vd_insertion       |  sequence of nucleotides corresponding to the non-templated insertion between the V and D segments
| dj_insertion       |  sequence of nucleotides corresponding to the non-templated insertion between the D and J segments
| fv_insertion       |  sequence of nucleotides corresponding to any "effective" non-templated insertion on the 5' side of the V (accounts for reads which extend beyond the 5' end of V)
| jf_insertion       |  sequence of nucleotides corresponding to any "effective" non-templated insertion on the 3' side of the J (accounts for reads which extend beyond the 3' end of J)
| mutated_invariants     |  true if the conserved cysteine or tryptophan (IMGT numbering) were mutated (colon-separated list if multi-hmm)
| in_frames      |  true if conserved cysteine and tryptophan (IMGT numbering) are in the same frame (colon-separated list if multi-hmm)
| stops                  |  true if stop codon was found in the query sequence (colon-separated list if multi-hmm)
| v_per_gene_support     |  approximate probability supporting the top V gene matches (semicolon-separated list of colon-separated gene:probability pairs)
| d_per_gene_support     |  approximate probability supporting the top D gene matches (semicolon-separated list of colon-separated gene:probability pairs)
| j_per_gene_support     |  approximate probability supporting the top J gene matches (semicolon-separated list of colon-separated gene:probability pairs)
| indelfos       |  information on any SHM indels that were inferred in the Smith-Waterman step. Written as a literal python dict; can be read in python with `ast.literal_eval(line['indelfo'])` (colon-separated list if multi-hmm)
| aligned_v_seqs     |  do not use. will soon be removed (see issue #179)
| aligned_d_seqs     |  do not use. will soon be removed (see issue #179)
| aligned_j_seqs     |  do not use. will soon be removed (see issue #179)

Note that `utils.process_input_line()` and `utils.get_line_for_output()` can be used to automate input/output.

##### `partition`: cluster sequences into clonally-related families

Example invocation:

```./bin/partis partition --infname test/example.fa --outfname _output/example.csv```

The output csv file headers are listed in the table below, and you can view a colored ascii representation of the clusters with the `view-partitions` action.
We write one line for the most likely partition (with the lowest logprob), as well as a number of lines for the surrounding less-likely partitions (set with `--n-partitions-to-write`)

|   column header  |  description
|------------------|----------------------------------------------------------------------------
| logprob          |  Total log probability of this partition
| n_clusters       |  Number of clusters (clonal families in this partition)
| partition        |  String representing the clusters, where clusters are separated by ; and sequences within clusters by :, e.g. 'a:b;c:d:e'
| n_procs          |  Number of processes which were simultaneously running for this clusterpath. In practice, final output is usually only written for n_procs = 1

To help visualize the clusters, you can tell it to print the most likely annotations for the final clusters with `--print-cluster-annotations`.
If you specify both `--print-cluster-annotations` and `--outfname`, the annotations will be written to a file name generated from `--outfname` (which can be viewed as other annotations, with `view-annotations`).

By default, this uses the most accurate and slowest method: hierarchical agglomeration with, first, hamming distance between naive sequences for distant clsuters, and full likelihood calculation for more similar clusters.
Like most clustering algorithms, this scales rather closer than you'd like to quadratically than it does to linearly.
We thus also have two faster and more approximate methods.

###### `--naive-hamming`

Use hard boundaries on hamming distance between naive sequences alone, with no likelihood calculation, to cluster with hierarchical agglomeration.
This is perhaps as much as twice as fast as the full method, but sacrifices significant accuracy because it doesn't know (for example) about the differing probabilities of mutations at different points in the sequence, and because it uses a single, fixed annotation rather than integrating over all likely annotations.

###### `--naive-vsearch`

First calculate naive (unmutated ancestor) for each sequence, then pass these into vsearch for very fast, very heuristic clustering.
The naive sequence calculation is easy to parallelize, so is fast if you have access to a fair number of processors.
Vsearch is also very fast, because it makes a large number of heuristic approximations to avoid all-against-all comparison, and thus scales significantly better than quadratically.
With `--n-procs` around 10 for the vsearch step, this should take only of order minutes for a million sequences.
Since it's entirely unprincipled, this of course sacrifices significant accuracy; but since we're using inferred naive sequences it's still much, much more accurate than clustering on SHM'd sequences.

##### `view-annotations`: Print (to std out) the annotations from an existing annotation output csv.

To, e.g. run on the output csv from the `run-viterbi` action:

``` ./bin/partis view-annotations --outfname run-viterbi-output.csv```

##### `view-partitions`: Print (to std out) the partitions from an existing partition output csv.

To, e.g. run on the output csv from the `partition` action:

``` ./bin/partis view-partitions --outfname partition-output.csv```

##### `cache-parameters`: write out parameter values and HMM model files for a given data set

The parameter-caching step is run automatically by the other actions if `--parameter-dir` doesn't exist (whether this directory is specified explicitly, or is left as default).
So you do not, generally, need to run it on its own except for testing.

When presented with a new data set, the first thing we need to do is infer a set of parameters, a task for which we need a preliminary annotation.
As such, partis first runs ighutil's Smith-Waterman algorithm on the data set.
The ighutil annotations are used to build and write out a parameter set, which is in turn used to make a set of HMM model files for each observed allele.
These files are then passed as input to a second, HMM-based, annotation step, which again outputs (more accurate) parameter values and HMM model files.

For example:

``` ./bin/partis cache-parameters --infname test/example.fa --parameter-dir _output/example```

When caching parameters, the parameter csvs from Smith-Waterman and the HMM are put into `/sw` and `/hmm` subdirectories of `--parameter-dir`.
Within each of these, there are a bunch of csv files with (hopefully) self-explanatory names, e.g. `j_gene-j_5p_del-probs.csv` has counts for J 5' deletions subset by J gene.
The hmm model files go in the `hmms` subdirectory, which contains yaml HMM model files for each observed allele.

When reading parameters (e.g. with run-viterbi or partition), the parameters are read directly from `--parameter-dir`, so you need to tack on either `/hmm` or `/sw` to what was passed during parameter caching.

###### Finding New Alleles

By default partis uses the set of germline V, D, and J genes ("germline set") in `data/imgt`.
If you have another set you'd like to use, you can do so by setting `--datadir`.
The default set from imgt, as is well known, is missing many real alleles, and as such, annotations for any individual who has these missing alleles will be wrong.
So we've implemented a method of finding new alleles on the fly when running on a new input data set.
We basically use the idea from [tigger](http://tigger.readthedocs.io/en/latest/), but generalize and robustify it a bit: in particular, we do something more akin to a simultaneous fit over all positions at once.

To try this, use the `--find-new-alleles` option to the `cache-parameters` action.
At the moment, since we don't want to modify users' germline sets without them realizing it, this simply prints out information on any new alleles that it finds.
If you want to incorporate these alleles into your germline set, just choose a name and add it to your germline fasta file.

To see how this works, you can run `./bin/find-new-alleles.py`.
This generates a hypothetical new allele separated by a few snps from an existing allele, simulates a data set with this allele, then attempts to identify this new allele in the simulated sample.

##### `simulate`: make simulated sequences

For testing purposes, we can use the parameters in `--parameter-dir` to create simulated sequences that mimic the data as closely as possible.
This allows us to test how well our algorithms work on a set of sequences for which we know the correct annotations and clonal relationships.
Note that while the parameters describe a very detailed picture of the rearrangement and mutation patterns, we do not (yet!) parametrize either the emipical clonal size distributions or phylogenetic relationships.
These can, however, be specified from a number of realistic options (see `--help` for the simulate subcommand).

For example:

```./bin/partis simulate --outfname _output/example/simu.csv --parameter-dir _output/example/data/hmm --n-sim-events 50```.

This will spit out simulated sequences to `--outfname` using the parameters we just made in `_output/example/data/hmm`.
We also specify that we want it to simulate 50 rearrangement events.
To get the actual number of sequences, we multiply this by the mean number of leaves per tree (5).
At the start of a simulation run, TreeSim generates a set of `--n-trees` trees (default 500 at the moment), and each tree has a number of leaves drawn from an exponential (or zipf, or box, or...) with mean/exponent `--n-leaves`.
Throughout the run, we sample a tree at random from this set for each rearrangement event.

##### `run-forward`: find total probability of sequences

Same as `run-viterbi`, except with the forward algorithm, i.e. it sums over all possible rearrangement events to get the total log probability of the sequence.
Probably mostly useful for testing.

### Parallelization

###### In general

The number of processes on your local machine is set with `--n-procs N`.
Since Smith-Waterman is so much faster than the hmm stuff, depending on your system it can make sense to use fewer process for that step -- you can thus specify a second, smaller number of processes `M` as `--n-procs N:M`.

In order to parallelize over more processes than the local machine can handle, we currently only support slurm.
This is specified by the flag `--slurm`, which runs all subsidiary processes by tacking `srun` onto the front of the command line.
Now, partis writes a lot of temporary files to a working directory, which is by default a random name under `/tmp/$USER`.
If you're running with slurm, though, you need the working directory to be a network mount that everybody can see, so if you're slurming you'll need to set `--workdir` to something visible by your batch nodes.
A suitable choice on our system is `_tmp/$RANDOM`; note, however, that for large data sets (of order, say a million sequences) i/o operations can become expensive, so one is behooved to use as fast a filesystem as possible.

###### run-viterbi

Sequence annotation (run-viterbi) lends itself quite readily to independent parallelization.
It should take 0.1-1 second per sequence; if you want it to go faster, just increase `--n-procs`.

###### partition

Clustering, however, really doesn't lend itself at all to independent parallelization -- we need to, at least approximately, compare each sequence to every other one.
For the full and point partis methods, we get around this by starting with `--n-procs` processes.
The input sequences are split evenly among these, and each process does all-against-all comparison (with many optimizations to avoid the full likelihood calculation) of all of its allotted sequences.
The results of this first round are collected and merged together, and then reapportioned among a new, smaller, number of processes.
This is continued until we arrive at one final process which is comparing all sequences.
Since at each stage we cache every calculated log probability, while the later steps have more sequences to compare, they also have more cached numbers at their disposal, and so it's possible to make each step take about the same amount of time.
We currently reduce the number of processes by about 1.6 at each step, as long as the previous step didn't have to calculate too many numbers.

With typical mutation levels, lineage structures, and cluster size distributions (all of which strongly affect clustering time), it's currently best to start with ``--n-procs` set so you have about 300 sequences per process.

#### Testing

The script `test/test.py` runs quite a few things and compares their outputs to reference values.
At the moment normal installation (e.g. described above) just runs a quick annotation test to make sure things are installed correctly (`--quick` option).
You can, however, run all the others if you like, although they're mainly designed to make sure code modifications haven't broken anything.

