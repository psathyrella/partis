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

based on StochHMM: https://github.com/KorfLab/StochHMM

installation
--------
```
scons
```

examples
--------

an example of the occasionally dishonest casino is in ```examples/casino.yaml```. Run it with:

```./hample --hmmfname examples/casino.yaml --seq 666655666613423414513666666666666```

while the similarly-canonical CPG island may be found in ```examples/cpg.yaml```. Run it with:

```./hample --hmmfname examples/cpg.yaml --seq ACTTTTACCGTCAGTGCAGTGCGCGCGCGCGCGCGCCGTTTTAAAAAACCAATT```

