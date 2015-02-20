#### partis

Partis is an HMM-based framework for B-cell receptor annotation, simulation, and partitioning.
It is built on top of the [ham](https://github.com/psathyrella/ham) HMM compiler, and also uses the [ighutil](https://github.com/cmccoy/ighutil) set of Smith-Waterman annotation tools.
A docker image can be found [here](https://registry.hub.docker.com/u/psathyrella/partis/), or if you prefer there also exist [installation instructions](https://github.com/psathyrella/partis/blob/master/install.txt).

Partis is free software under the GPL v3.

#### basic usage

For full analysis chain, run
```
scons validate
```
This does parameter inference and model building on a sample data file, then makes a sample of simulated sequences, then does parameter inference and model building on the simulation sample, and finally evaluates performance on simulation.

This really just runs `./bin/run-driver.py` with a few command line options, which in turn runs `./bin/partis.py` a number of times with different options.
The full command line is printed at each step, so if you want to redo some steps with different data or options you can just rerun those steps by copying and modifying the printed command lines.
Both of these commands will also print usage messages if you pass them `--help`.
