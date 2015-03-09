### Introduction

Partis is an HMM-based framework for B-cell receptor annotation, simulation, and partitioning.
It is built on top of the [ham](https://github.com/psathyrella/ham) HMM compiler, and also uses the [ighutil](https://github.com/cmccoy/ighutil) set of Smith-Waterman annotation tools.
Partis is free software under the GPL v3.

This manual is organized into the following sections: Quick Start, if you want to start doing things without understanding what they do; Subcommands, a rundown on each of the actions `partis.py` can run; and Higher Abstractions, for a description of scripts that automate a number of `partis.py` actions.
In general, we will assume that the reader is familiar with the [paper](TODO add a link) describing partis.

### Quick Start

Partis has a lot of dependencies.
We promise that these are all totally necessary (except ROOT! That &*@! is getting exorcised real soon).
But the upshot is that partis is kind of a pain to install from scratch.
So, unless you need to do a lot of mucking about under the hood, you'll have an easier time of it with just the [Docker image](https://registry.hub.docker.com/u/psathyrella/partis/).
Docker images are kind of like lightweight virtual machines, and as such all the dependencies are taken care of automatically.

To use Docker, you'll first want to follow the installation instructions for your particular system.
Once you have Docker installed, pull the partis image from dockerhub

``` sudo docker pull psathyrella/partis ```.

Then enter it with

``` sudo docker run -t -i psathyrella/partis /bin/bash ```.

This will drop you into the main `partis/` directory.
From here, you can run the full analysis chain, if you will, with

```scons validate```.

To find out what the "full analysis chain" is, look elsewhere in this manual.
In short, though, it infers a set of parameters from a test data set, then makes a simulated data set based on these parameters, and finally annotates these simulated sequences.

### Subcommands

`partis.py`, in `bin/`, is the script that drives everything you can get partis to do.
Every time you invoke it, you choose from a list of actions you want it to perform

```./bin/partis.py --action cache-parameters|run-viterbi|run-forward|partition|simulate|build-hmms|generate-trees ```.

Each of these subcommands is described in detail below.

#### cache-parameters

When presented with a new data set, the first thing we need to do is infer a set of parameters, a task for which we need a preliminary annotation.
As such, 


### Higher Abstractions

Starting with a small data set, this runs a preliminary annotation step with ighutil, then uses these parameters to build a preliminary set of HMM model files.
These files are then used for an HMM annotation run which builds a more accurate set of parameters.
These new parameters are then passed to 

This does parameter inference and model building on a sample data file, then makes a sample of simulated sequences, then does parameter inference and model building on the simulation sample, and finally evaluates performance on this simulation.

Under the hood, this just runs `./bin/run-driver.py` once with a few command line options; this in turn runs `./bin/partis.py` several times with different options.
The full command line is printed at each step, so if you want to redo some steps with different data or options you can just rerun those steps by copying and modifying the printed command lines.
Both of these commands will also print usage messages if you pass them `--help`.





Note that if you only want to run things (not edit the source), it is probably easier to just use the [docker image](https://registry.hub.docker.com/u/psathyrella/partis/).

dependecies
==============
required
--------------
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

required only for simulator
--------------
  - dendropy
  - TreeSim

included
--------------
  - tclap
  - yaml-cpp
  - ham
  - samtools
  - bppseqgen

installation
==============

```
git clone git@github.com:psathyrella/partis
cd partis/packages/ham/
scons
cd ../samtools/
make
ln -s $PWD/samtools ~/bin/
export PATH=~/bin:$PATH
cd ../ighutil/
make -C clj
pip install --user ./python
cd ..
```
