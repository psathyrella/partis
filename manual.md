Partis is an HMM-based framework for B- and T-cell receptor sequence annotation, simulation, clonal family, and germline inference. It uses the [ham](https://github.com/psathyrella/ham) HMM compiler and [ig-sw](https://github.com/matsengrp/ig-sw) set of Smith-Waterman annotation tools.

Partis is free software under the GPL v3, and is described in the following papers:

  * [HMM framework and BCR annotation](http://doi.org/10.1371/journal.pcbi.1004409) Ralph, DK, & Matsen IV, FA (2016). _Consistency of \[...\] Enables Accurate B Cell Receptor Sequence Annotation._ PLOS Computational Biology, 12(1), e1004409.
  * [Clonal family inference](http://dx.doi.org/10.1371/journal.pcbi.1005086) Ralph, DK, & Matsen IV, FA (2016). _Likelihood-based Inference of B-cell Clonal Families._ PLOS Computational Biology, 12(10), e1005086.
  * [Germline inference](https://arxiv.org/abs/1711.05843) Ralph, DK, & Matsen IV, FA (in review at PLOS Computational Biology) _Per-sample immunoglobulin germline inference \[...\]_ 

This manual is divided into the following pages:

  * [Installation](docs/install.md)
  * [Quick start](docs/quick-start.md)
  * [Subcommands](docs/subcommands.md)
  * [Output formats](docs/output-formats.md)
  * [Parallelization](docs/parallel.md)

Details of the multitude of flags not vital enough to document here, but which can be used to control many aspects of behavior can be found by running each subcommand's help, for instance `partis cache-parameters --help` or `partis simulate --help`.

To ask questions, or search through past discussions, please use the [google group](https://groups.google.com/forum/#!forum/partis).
For specific issues with the software, e.g. bug reports or feature requests, on the other hand, [submit an issue](https://github.com/psathyrella/partis/issues?utf8=%E2%9C%93&q=) on github.
