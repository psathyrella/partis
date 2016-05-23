### Introduction

Ham is a general-purpose HMM compiler with an emphasis on minimizing both development and run time.
Ham is very young: while we think its virtues will prove useful to others, we've so far focused mostly on implementing features that we needed.
If you'd like something added, submit a [new issue](https://github.com/psathyrella/ham/issues/new) on github -- it might be really easy and we'd love to help.

In writing ham, We were greatly inspired by the [HMMoC](http://genserv.anat.ox.ac.uk/downloads/software/hmmoc/) and [StochHMM](https://github.com/KorfLab/StochHMM) tools, but neither of them fully satisfied our needs.
We found the XML input format in HMMoC complicated to script, and the paradigm of auto-generated C++ code to be extremely difficult to test and debug.
StochHMM, meanwhile, had a custom configuration file format, lacked a pair-HMM implementation, and possessed an extraordinarily extensive but not entirely functional feature set.
However, StochHMM's overall structure, and its basic idea of reading HMM structure from a text file rather than generating code, were similar to what we required, so we used it as a starting point for a complete rewrite.
(As an aside, we note that the excellent HMMER tool only implements profile HMMs, and thus was not appropriate for our needs.)

From a usability standpoint, ham is distinguished by the use of [YAML](http://yaml.org) config files.
These are concise plain text (for example, the CPG island XML config in HMMoC is 5961 characters, while examples/cpg.yaml is 408 characters) files which are also eminently scriptable using libraries for a variety of languages.
Changing anything in the HMMs just requires editing these text files, then rerunning an existing C++ binary.

In performance tests, we've found ham to be a bit faster and somewhat more memory-efficient than HMMoC.

`ham` is free software under the GPL v3.

This manual is divided into sections on installation, input specification and running.


### Installation

As an alternative to installation, consider using the latest [docker image](https://registry.hub.docker.com/u/psathyrella/ham/).

You will need to have recent versions of the following dependencies on your system:
  - [boost](http://www.boost.org/) headers
  - [GCC](https://gcc.gnu.org/)
  - [scons](http://scons.org/)

Ham is guaranteed to build with the [matsengrp/cpp](https://github.com/matsengrp/dockerfiles/blob/master/cpp/Dockerfile) Docker container.

To compile the library and `hample` (the example binary, which is suitable for general-purpose inference) run `scons` in the top-level directory.

### Input specification

YAML is super easy to write by hand for simple cases.
For instance this

```
name: cpg
tracks:
  nucs: [A,C,G,T]
states:
- name: init
  transitions:
    island: 0.3
    ocean: 0.7
- name: island
  emissions:
    probs:
      A: 0.1
      C: 0.4
      G: 0.4
      T: 0.1
  transitions:
    island: 0.8
    ocean: 0.05
    end: 0.15
- name: ocean
  emissions:
    probs:
      A: 0.25
      C: 0.25
      G: 0.25
      T: 0.25
  transitions:
    island: 0.3
    ocean: 0.6
    end: 0.1
```

is the entirety of the canonical CPG island example.

The input file needs to have three things: the HMM name, a list of tracks (which has to have length one at the moment), and a list of states.
The first state has to be `init`.
Besides a name, each state has to have a list of transitions; each transition is just the to-state name and the transition probability.
Most states will also have emissions, which consist of the track on which to emit, and a list of probabilities for each letter in that track's alphabet.
Some of your states also have to go to the `end` state.

For complicated HMMs, you'll want to script the yaml creation.
There exist yaml modules for most popular languages; we like those for [python](http://pyyaml.org/) and [C++](https://code.google.com/p/yaml-cpp/).

### Running

To run the cpg island example in `examples/cpg.yaml`, issue the command

```./hample --hmmfname examples/cpg.yaml --seqs ACTTTTACCGTCAGTGCAGTGCGCGCGCGCGCGCGCCGTTTTAAAAAACCAATT```

This prints out the most likely (Viterbi) path for the sequence we've passed on the command line, as well as the Forward (total) probability.
To simultaneously emit as many sequences as you like (for two sequences we would call this a pair HMM), just tack on more colon-separated sequences:

    ./hample --hmmfname examples/cpg.yaml --seqs CGCCGCACTTTTACCGTCAGTGCAGTGCGC:TGCTTGCGCGCGCGAGCCGTTTGCATTAAC:GCGGCGCAAAAAACCGTCAGTGCAGTGCTT

The Viterbi path should now be thought of as the single most probable path that all of the sequences could have taken together.
Similarly, the Forward probability is now the total probability, summed over all paths, of these sequences having taken the same path.
The example directory also has the occasionally dishonest casino:

    ./hample --hmmfname examples/casino.yaml --seqs 666655666613423414513666666666666

If you want to do more complicated things with the output, start from the source file for this example binary, `src/hample.cc`.
