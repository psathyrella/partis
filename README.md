      ;`.                       ,'/
      |`.`-.      _____      ,-;,'|
      |  `-.\__,-'     `-.__//'   |
      |     `|               \ ,  |
      `.  ```                 ,  .'
        \_`      .     ,   ,  `_/
          \    ^  `   ,   ^ ` /
           | '  |  ____  | , |
           |     ,'    `.    |
           |    (  O' O  )   |
           `.    \__,.__/   ,'
             `-._  `--'  _,'
                 `------'

# ham -- the fantabulous HMM compiler

`ham` began from [StochHMM](https://github.com/KorfLab/StochHMM), but has been completely rewritten.

`ham` is free software under the GPL v3.


## installation

As an alternative to installation, consider using the latest [docker image](https://registry.hub.docker.com/u/psathyrella/ham/).

#### dependencies

  - [boost](http://www.boost.org/) headers
  - [GCC](https://gcc.gnu.org/)
  - [scons](http://scons.org/)

The package is guaranteed to build with the [matsengrp/cpp](https://github.com/matsengrp/dockerfiles/blob/master/cpp/Dockerfile) Docker container.


#### building

Running `scons` in the top-level directory will build the example binary called `hample`, which can be used for general-purpose inference.


## examples

An example of the occasionally dishonest casino is in `examples/casino.yaml`. Run it with:

    ./hample --hmmfname examples/casino.yaml --seqs 666655666613423414513666666666666

This runs the Viterbi algorithm, prints the Viterbi (most probable) path, and its log probability.
It then runs the Forward algorithm and prints the log Forward probability (the probability of taking any path through the HMM).

The similarly-canonical CPG island example may be found in `examples/cpg.yaml`:

    ./hample --hmmfname examples/cpg.yaml --seqs ACTTTTACCGTCAGTGCAGTGCGCGCGCGCGCGCGCCGTTTTAAAAAACCAATT

To simultaneously emit onto as many sequences as you like (for two sequences we would call this a pair HMM), just tack on more colon-separated sequences:

    ./hample --hmmfname examples/cpg.yaml --seqs CGCCGCACTTTTACCGTCAGTGCAGTGCGC:TGCTTGCGCGCGCGAGCCGTTTGCATTAAC:GCGGCGCAAAAAACCGTCAGTGCAGTGCTT

The Viterbi path should now be thought of as the single most probable path that all of the sequences could have taken together.
Similarly, the Forward probability is now the total probability, summed over all paths, of these sequences having taken the same path.


## why ham?

We were greatly inspired by the [HMMoC](http://genserv.anat.ox.ac.uk/downloads/software/hmmoc/) and [StochHMM](https://github.com/KorfLab/StochHMM) tools, but neither of them fully satisfied our needs.
We found the XML input format in HMMoC complicated to script, and the paradigm of auto-generated C++ code to be extremely difficult to test and debug.
StochHMM, meanwhile, had a custom configuration file format, lacked a pair-HMM implementation, and possessed an extraordinarily extensive but not entirely functional feature set.
However, StochHMM's overall structure, and its basic idea of reading HMM structure from a text file rather than generating code, were similar to what we required, so we used it as a starting point for a complete rewrite.
(As an aside, we note that the excellent HMMER tool only implements profile HMMs, and thus is not appropriate for our needs.)

From a usability standpoint, ham is distinguished by the use of [YAML](http://yaml.org) config files.
These are concise plain text (for example, the CPG island XML config in HMMoC is 5961 characters, while examples/cpg.yaml is 440 characters) files which are also eminently scriptable using libraries for a variety of languages.
