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

ham -- the fantabulous HMM compiler
--------

`ham` is based on [StochHMM](https://github.com/KorfLab/StochHMM).

building
--------
Install [SCons](http://www.scons.org/) and then run

    scons

examples
--------

The occasionally dishonest casino example is in `examples/casino.yaml`. Run it with:

    ./hample --hmmfname examples/casino.yaml --seq 666655666613423414513666666666666

The similarly-canonical CPG island example may be found in `examples/cpg.yaml`. Run it with:

    ./hample --hmmfname examples/cpg.yaml --seq ACTTTTACCGTCAGTGCAGTGCGCGCGCGCGCGCGCCGTTTTAAAAAACCAATT

why ham?
--------

There are already lots of HMM compilers out there. Why do we need a new one?

**HMMOC** Works well once you get it doing what you want, but the xml config files are very difficult
to use. You need to include blocks of actual c++ code in the xml, So
you are in essence using an xml file as a text editor to write c++ code. This is suboptimal.
You then generate c++ code from this xml and compile it, and in practice auto-generated
c++ is extremely hard to debug and understand if anything goes wrong. In practice, something always
goes wrong. In sum, implementing the huge and complex HMM that I needed, correctly, proved to be infeasible.

**HMMER** works really well, and is super user friendly, but only does profile HMMS.

**StochHMM** has more concise plain text config files (although they are a custom format). Unfortunately, 
while the overall code structure is great and there is a huge feature set, a large fraction of the features
are not completely implemented, and a (different) large fraction are undocumented, so it is
very difficult to know what you are telling it to do. It also does not have pair HMMS.

While trying to implement pair HMMS in StochHMM, it became clear that it was going to be a complete rewrite.
From a useability standpoint, the biggest difference is yaml config files. These are plain text files that are
incredibly concise (the CPG island xml config in HMMOC is 5961 characters, while examples/cpg.yaml
is 480 characters). Yaml is also emminently scriptable from within python.
.