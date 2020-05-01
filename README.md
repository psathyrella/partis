### partis

Partis is an HMM-based framework for B- and T-cell receptor sequence annotation, simulation, clonal family and germline inference, and affinity prediction.
It is built on top of the [ham](https://github.com/psathyrella/ham) HMM compiler and [ig-sw](https://github.com/matsengrp/ig-sw) set of Smith-Waterman annotation tools.
Partis is free software under the GPL v3.

The various components are described in the following papers.
Since they do quite different things, it's best if you can cite the specific paper(s) that describe the components that you're using.

  * [Selection metrics](https://arxiv.org/abs/2004.11868) Ralph, DK, & Matsen IV, FA (2020). Submitted to PLOS Computational Biology. _Using B cell receptor lineage structures to predict affinity._
  * [Germline inference](https://doi.org/10.1371/journal.pcbi.1007133) Ralph, DK, & Matsen IV, FA (2019). _Per-sample immunoglobulin germline inference \[...\]._ PLOS Computational Biology, 15(7), e1007133.
  * [Clonal family inference](http://dx.doi.org/10.1371/journal.pcbi.1005086) Ralph, DK, & Matsen IV, FA (2016). _Likelihood-based Inference of B-cell Clonal Families._ PLOS Computational Biology, 12(10), e1005086.
  * [HMM framework and BCR annotation](http://doi.org/10.1371/journal.pcbi.1004409) Ralph, DK, & Matsen IV, FA (2016). _Consistency of \[...\] Enables Accurate B Cell Receptor Sequence Annotation._ PLOS Computational Biology, 12(1), e1004409.

The manual is divided into the following pages:

  * [Installation](docs/install.md)
  * [Quick start](docs/quick-start.md)
  * [Subcommands](docs/subcommands.md)
  * [Output formats](docs/output-formats.md)
  * [Plotting](docs/plotting.md)
  * [Germline inference confidence](docs/germline-inference.md)
  * [Parallelization](docs/parallel.md)

Details of the multitude of flags not vital enough to document here, but which can be used to control many aspects of behavior can be found by running each subcommand's help, for instance `partis cache-parameters --help` or `partis simulate --help`.

To ask questions, or search through past discussions, please use the [google group](https://groups.google.com/forum/#!forum/partis).
For specific issues with the software, e.g. bug reports or feature requests, on the other hand, [submit an issue](https://github.com/psathyrella/partis/issues?utf8=%E2%9C%93&q=) on github.

* [docker image](https://hub.docker.com/r/psathyrella/partis/)
* [support forum](https://groups.google.com/forum/#!forum/partis) subscribe to this for (very occasional) updates when something important gets changed
