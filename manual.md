### Introduction

Partis is an HMM-based framework for B-cell receptor annotation, simulation, and partitioning.
It is built on top of the [ham](https://github.com/psathyrella/ham) HMM compiler, and also uses the [ighutil](https://github.com/cmccoy/ighutil) set of Smith-Waterman annotation tools.
Partis is free software under the GPL v3.

This manual is organized into the following sections:

  * [Quick Start](#quick-start) install/run with Docker
  * [Slow Start](#slow-start) install from scratch
  * [Details](#details) how to navigate the various `partis.py` subcommands

There are also many flags and optional parameters; unless mentioned below these are tautologically beyond the scope of this manual.
Details concerning their purpose, however, may be gleaned by means of the following incantation: `./bin/partis.py --help`.
In general, we will assume that the reader is familiar with the [paper](http://arxiv.org/abs/1503.04224) describing partis.

### Quick Start

Because partis has a lot of dependencies, you'll likely have an easier time of it using the [Docker image](https://registry.hub.docker.com/u/psathyrella/partis/) rather than installing from scratch.
Docker images are kind of like lightweight virtual machines, and as such all the dependencies are taken care of automatically.
If, however, you'll be doing a lot of mucking about under the hood, bare installation might be preferable.

You'll first want install Docker using their [installation instructions](https://docs.docker.com) for your particular system.
Once Docker's installed, pull the partis image from dockerhub, start up a container from this image and attach yourself to it interactively, and compile:

```
sudo docker pull psathyrella/partis
sudo docker run -it -v /:/host psathyrella/partis /bin/bash
./bin/build.sh
export PATH=$PATH:$PWD/packages/samtools
```
Depending on your setup, the `sudo` may be unnecessary.
Note the `-v`, which mounts the root of the host filesystem to `/host` inside the container.

Now you can run individual partis commands (described below), poke around in the code, or run the scons targets `test` or `validate`.
If you just want to annotate a set of BCR sequences, say at `/path/to/yourseqs.fa`, run

``` ./bin/annotate --infname /host/path/to/yourseqs.fa```

Note that now we're inside the container, we access the fasta file at the original path on your host system, but with `/host` tacked on the front (as we specified in `docker run` above).
This command by default writes the output csv to the directory that `--infname` came from.
Depending on your system, 5000 sequences will take perhaps ten minutes -- if your ratio of patience to sequences is quite different to this, you should look through the parallelization options below.

To detach from the docker container without stopping it, hit `ctrl-p ctrl-q`.

###### Docker tips
Docker containers and images are kinda-sorta like virtual machines, only different, so a few things:
  - We use `docker run` above: this creates a new container from (instance of ) the partis image
  - If you exit (ctrl-d or `exit`) and then do `docker run` again, that'll create a new container. But most of the time you want to reattach to the one you made before.
    - to reattach to the same container:
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
  - ROOT
  - python packages
    - pysam
    - pyyaml
    - cython
    - networkx
    - decorator
    - lxml
    - bs4 (beautiful soup)
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
export PATH=$PATH:$PWD/packages/samtools
```

### Details

#### Subcommands

`partis.py`, in `bin/`, is the script that drives everything.
Every time you invoke it, you choose from a list of actions you want it to perform

```./bin/partis.py --action cache-parameters|simulate|run-viterbi|run-forward|partition|build-hmms|generate-trees ```.

Each of these subcommands is described in detail below.

##### `cache-parameters`: write out parameter values and HMM model files for a given data set

When presented with a new data set, the first thing we need to do is infer a set of parameters, a task for which we need a preliminary annotation.
As such, partis first runs ighutil's Smith-Waterman algorithm on the data set.
The ighutil annotations are used to build and write out a parameter set, which is in turn used to make a set of HMM model files for each observed allele.
These files are then passed as input to a second, HMM-based, annotation step, which again outputs (more accurate) parameter values and HMM model files.

The full command you'd need to cache parameters would look like this:

``` ./bin/partis.py --action cache-parameters --seqfile test/example.fa --is-data --parameter-dir _output/example/data```

`--seqfile` location of input sequences (.[ct]sv, .f[aq]{sta,stq})
`--is-data` lets it know not to expect simulation information (which would be used for validation purposes).
`--skip-unproductive` Smith-Waterman step skips sequences that are annotated to be unproductive rearrangements.
When caching parameters, the parameter csvs from Smith-Waterman and the HMM are put into `/sw` and `/hmm` subdirectories of `--parameter-dir`.
Within each of these, there are a bunch of csv files with (hopefully) self-explanatory names, e.g. `j_gene-j_5p_del-probs.csv` has counts for J 5' deletions subset by J gene.
There is also a `mute-freqs` directory, which has per-position mutation frequencies for each observed allele; and an `hmms` directory, which contains yaml HMM model files for each observed allele.

When reading parameters, the parameters are read from `--parameter-dir`, so tack on either `/hmm` or `/sw` to what was passed during parameter caching.

##### `simulate`: make simulated sequences

Now that we've got a set of parameters for this data set, we can use it to create simulated sequences that mimic the data as closely as possible.
This will allow us to test how well our algorithms work on a set of sequences for which we know the correct annotations.
The basic command to run simulation is

```./bin/partis.py --action simulate --outfname _output/example/simu.csv --parameter-dir _output/example/data/hmm --n-sim-events 50 --n-leaves 5```.

This will spit out simulated sequences to `--outfname` using the parameters we just made in `_output/example/data/hmm`.
We also specify that we want it to simulate 50 rearrangement events.
To get the actual number of sequences, we multiply this by the mean number of leaves per tree (5).
At the start of a simulation run, TreeSim generates a set of `--n-trees` trees (default 500 at the moment), and each tree has a number of leaves drawn from an exponential with mean `--n-leaves`.
Throughout the run, we sample a tree at random from this set for each rearrangement event.

##### `run-viterbi`: find most likely annotations

If you already have parameters and HMM files cached from a previous run, you can just run the Viterbi algorithm by itself.
As an example invocation, we could find the 5 most likely annotations for the first sequence in the previous data set

```./bin/partis.py --action run-viterbi --seqfile test/example.fa --is-data --parameter-dir _output/example/data/hmm --n-best-events 5 --n-max-queries 1 --debug 1```

Here I've kicked the debug level up to 1, so it prints out a colored summary of the candidate rearrangement events.
If you want to save csv output for later, add `--outfname <outfname>`.

##### `run-forward`: find total probability of sequences

Exactly the same as `run-viterbi`, except with the forward algorithm, i.e. it sums over all possible rearrangement events to get the total log probability of the sequence.

### Parallelization

Happily enough, sequence annotation lends itself quite readily to independent parallelization.
You specify the number of processes on your local machine with `--n-procs`.
If you also throw in the `--slurm` flag, subsidiary processes will be run with slurm (under the hood this just adds 'srun' to the front of the commands -- so if you're inside docker you'll also need to install `srun`).
Now, partis writes a lot of temporary files to a working directory, which is by default under `/tmp/$USER`.
If you're running with slurm, though, you need the working directory to be a network mount everybody can see, so if you're slurming you'll need to set `--workdir` to something visible by your batch nodes.
A suitable choice on our system is `_tmp/$RANDOM`.

#### Higher Abstractions

The script `bin/run-driver.py` can help to automate some of these steps.
The command

```./bin/run-driver.py --label example --datafname test/A-every-100-subset-0.tsv.bz2```

will cache data parameters, run simulation, cache simulation parameters, and then run annotation a final time in order to plot performance.
