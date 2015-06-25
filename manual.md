### Introduction

Partis is an HMM-based framework for B-cell receptor annotation, simulation, and partitioning.
It is built on top of the [ham](https://github.com/psathyrella/ham) HMM compiler, and also uses the [ighutil](https://github.com/cmccoy/ighutil) set of Smith-Waterman annotation tools.
Partis is free software under the GPL v3.

This manual is organized into the following sections: Quick Start (to start doing things without understanding what they do); Installation; Subcommands (a rundown on each of the actions `partis.py` can run); Parallelization; and Higher Abstractions (for a description of scripts that automate a number of `partis.py` actions).
There are also many flags and optional parameters; unless mentioned below these are beyond the scope of this manual.
Details concerning their purpose, however, may be gleaned by means of the following incantation: `./bin/partis.py --help`.
In general, we will assume that the reader is familiar with the [paper](http://arxiv.org/abs/1503.04224) describing partis.

### Quick Start

Because partis has a lot of dependencies, you'll likely have an easier time of it using the [Docker image](https://registry.hub.docker.com/u/psathyrella/partis/) rather than installing from scratch.
Docker images are kind of like lightweight virtual machines, and as such all the dependencies are taken care of automatically.
If, however, you'll be doing a lot of mucking about under the hood, bare installation might be preferable.

You'll first want install Docker using their [installation instructions](https://docs.docker.com) for your particular system.
Once Docker's installed, pull the partis image from dockerhub (depending on your setup, the `sudo`s may be unnecessary), start up a container from this image and attach yourself to it interactively, and compile:

```
sudo docker pull psathyrella/partis
sudo docker run -i -t psathyrella/partis /bin/bash
./bin/build.sh && export PATH=$PATH:$PWD/packages/samtools
```

Then, you can run individual partis commands (described below), poke around in the code, or run the scons targets `test` or `validate`.
If you just want to annotate a set of BCR sequences, say that you've got in `yourseqs.fa`
``` ./bin/annotation yourseqs.fa ```
To detach from the docker container without stopping it, hit `ctrl-p ctrl-q`.

###### Docker tips
Docker is really awesome.
But it's confusing at first, mostly because it's kinda-sorta like virtual machines, only different.
The thing to remember is that everything that isn't working for you has in the past also confused many other people on stackoverflow, so it's easy to google what you're doing wrong.

A few tips:
  - We use `docker run` above: this creates a new "container" from the current image (kindasorta: "creates a new instance of the partis image").
  - If you exit and then do `docker run` again, that'll create a new container. But most of the time you want to reattach to the one you made before.
    - to reattach to the same container:
      - run `docker ps -a` (lists all running and stopped containers) to get the right container ID
      - run `docker attach <ID>`
    - Hence the `-i -t` and `/bin/bash` options: this allocates a pseudo-tty keeps STDIN open, and runs bash instead of the default command, without all of which you can't reattach
    - [this](http://stackoverflow.com/questions/26153686/how-to-run-a-command-on-an-already-existing-docker-container) is a helpful discussion. Also, look up the difference between `attach` and `run` in the docker docs.

### Installation

As noted above, for most use cases you'll likely be happier using the [Docker image](https://registry.hub.docker.com/u/psathyrella/partis/).
But installing from scratch isn't really so bad, you just need to make sure you have the following dependencies:

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

The following packages are also used by partis, but they're included as sub trees in the source code, so you don't need to do anything:
  - tclap
  - yaml-cpp
  - ham
  - samtools
  - bppseqgen

Once you've got all the necessary things on your system, you can proceed to clone the repository:

```
git clone git@github.com:psathyrella/partis
cd partis
```

And then build:

```
source ./bin/handbuild.sh
```

### Subcommands

`partis.py`, in `bin/`, is the script that drives everything.
Every time you invoke it, you choose from a list of actions you want it to perform

```./bin/partis.py --action cache-parameters|simulate|run-viterbi|run-forward|partition|build-hmms|generate-trees ```.

Each of these subcommands is described in detail below.

##### `cache-parameters`: write out parameter values and HMM model files for a given data set

When presented with a new data set, the first thing we need to do is infer a set of parameters, a task for which we need a preliminary annotation.
As such, partis first runs ighutil's Smith-Waterman algorithm on the data set.
The ighutil annotations are used to build and write out a parameter set, which is in turn used to make a set of HMM model files for each observed allele.
These files are then passed as input to a second, HMM-based, annotation step, which again outputs (more accurate) parameter values and HMM model files.

The full command you'd need to cache parameters would look like this: TODO use a different example data file here

``` ./bin/partis.py --action cache-parameters --seqfile test/A-every-100-subset-0.tsv.bz2 --is-data --skip-unproductive --parameter-dir _output/example/data --plotdir _plots/example```

`--seqfile` tells it where we've put the input sequences (this example is a subset of the Adaptive data set).
`--is-data` lets it know not to expect simulation information (which would be used for validation purposes).
`--skip-unproductive` ignores sequences that are annotated to be unproductive rearrangements.
The final two arguments tell it where to put output: parameter values and HMM model files will go in `--parameter-dir`, and plots in `--plotdir`.

If you poke around in the output directory after running, you'll find `sw` (for ighutil parameters) and `hmm` subdirectories.
Within each of these, there are a bunch of csv files with (hopefully) self-explanatory names, e.g. `j_gene-j_5p_del-probs.csv` has counts for J 5' deletions subset by J gene.
There is also a `mute-freqs` directory, which has per-position mutation frequencies for each observed allele; and an `hmms` directory, which contains yaml HMM model files for each observed allele.

##### `simulate`: make simulated sequences

Now that we've got a set of parameters for this data set, we can use it to create simulated sequences that mimic the data as closely as possible.
This will allow us to test how well our algorithms work on set of sequences for which we know the correct annotations.
The basic command to run simulation is

```./bin/partis.py --action simulate --outfname _output/example/simu.csv --parameter-dir _output/example/data/hmm --n-max-queries 10```.

This will spit out simulated sequences to `<outfname>` using the parameters we just made in `<parameter-dir>`.
We also specify that we want it to simulate 10 rearrangement events (not 10 sequences) with `<n-max-queries>`.
To get the actual number of sequences, we multiply this by the number of leaves per tree.
There are a couple of levers available for setting the number of leaves per tree.
At the start of a simulation run, we use TreeSim to generate a set of `--n-trees` trees.
Throughout the run, we sample a tree at random from this set for each rearrangement event.
If `--random-number-of-leaves` is false, all the trees in this set will have the same number of leaves (`--n-leaves`).
Otherwise, we choose a number of leaves at random for each tree, from some distribution (subject to change TODO decide what to say about this).

##### `run-viterbi`: find most likely annotations

If you already have parameters and HMM files cached from a previous run, you can just run the Viterbi algorithm by itself.
As an example invocation, we could find the 5 most likely annotations for the first sequence in the previous data set TODO make sure this works without a plotdir specified

```./bin/partis.py --action run-viterbi --seqfile test/A-every-100-subset-0.tsv.bz2 --is-data --parameter-dir _output/example/data/hmm --n-best-events 5 --n-max-queries 1 --debug 1```

Here I've kicked the debug level up to 1, so it prints out a colored summary of the candidate rearrangement events.

##### `run-forward`: find total probability of sequences

Exactly the same as `run-viterbi`, except with the forward algorithm, i.e. it sums over all possible rearrangement events to get the total log probability of the sequence.

### Parallelization

Happily enough, sequence annotation lends itself quite readily to independent parallelization.
You specify the number of processes on your local machine with `--n-procs`.
If you also throw in the `--slurm` flag, subsidiary processes will be run with slurm (under the hood this just adds 'srun' to the front of the commands).
Now, partis writes a lot of temporary files to a working directory, which is by default under `/tmp/$USER`.
If you're running with slurm, though, you need the working directory to be a network mount everybody can see, so if you're slurming you'll need to set `--workdir` to something visible by your batch nodes.
A suitable choice on our system is `_tmp/$RANDOM`.

### Higher Abstractions

The script `bin/run-driver.py` can help to automate some of these steps.
The command

```./bin/run-driver.py --label example --datafname test/A-every-100-subset-0.tsv.bz2 --plotdir _plots/example```

will cache data parameters, run simulation, cache simulation parameters, and then run annotation a final time in order to plot performance.
