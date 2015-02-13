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

Ham is free software under the GPL v3.

# dependencies

  - boost headers
  - g++
  - astyle
  - scons
  - or look at
    - [https://github.com/matsengrp/dockerfiles/blob/master/cpp/Dockerfile]
    - [http://dockerfile.github.io/#/python]
    - [http://dockerfile.github.io/#/ubuntu]

## building

Running `scons` in the top-level directory will build the example binary called `hample`.

    scons

## examples

An example of the occasionally dishonest casino is in `examples/casino.yaml`. Run it with:

    ./hample --hmmfname examples/casino.yaml --seqs 666655666613423414513666666666666

This runs the Viterbi algorithm, prints the Viterbi (most probable) path, and its log probability.
It then runs the Forward algorithm and prints the log Forward probability (the probability of taking any path through the HMM).

The similarly-canonical CPG island example may be found in `examples/cpg.yaml`:

    ./hample --hmmfname examples/cpg.yaml --seqs ACTTTTACCGTCAGTGCAGTGCGCGCGCGCGCGCGCCGTTTTAAAAAACCAATT

To simultaneously emit onto as many sequences as you like (for two sequences we would call this a pair HMM),
just tack on more colon-separated sequences:

    ./hample --hmmfname examples/cpg.yaml --seqs CGCCGCACTTTTACCGTCAGTGCAGTGCGC:TGCTTGCGCGCGCGAGCCGTTTGCATTAAC:GCGGCGCAAAAAACCGTCAGTGCAGTGCTT

The Viterbi path should now be thought of as the single most probable path that all of the sequences could
have taken together. Similarly, the Forward probability is now the total probability, summed over all paths,
of these sequences having taken the same path.

## why ham?

There are already lots of HMM compilers out there. Why do we need a new one?

**HMMoC** Works well once you get it doing what you want, but the XML config files are very difficult
to use. You are in effect using an XML file as a text editor to write c++ code. We found this confusing.
You then generate c++ source from this XML, and compile it. In practice auto-generated
c++ is quite difficult to debug if anything goes wrong. Implementing the really big
HMMs that we need turned out to be infeasible.

**HMMER** works really well, and is super user friendly, but only does profile HMMs.

**StochHMM** has nice, concise plain text config files (although they are a custom format). Unfortunately,
while the overall code structure is great and there is a huge feature set, a large fraction of the features
are not completely implemented, and a (different) large fraction are undocumented, so it is
difficult to know what you are telling it to do. It also does not have pair HMMS.

While trying to implement pair HMMS in StochHMM, it became clear that it was going to be a complete rewrite.

From a usability standpoint, ham is distinguished by the use of [YAML](http://yaml.org) config files. These are concise
plain text (for example, the CPG island xml config in HMMoC is 5961 characters, while examples/cpg.yaml
is 440 characters). Yaml is also emminently scriptable with existing python modules.
