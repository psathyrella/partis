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
