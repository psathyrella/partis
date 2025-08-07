[![](https://img.shields.io/static/v1?label=AIRR-C%20sw-tools%20v1&message=compliant&color=008AFF&labelColor=000000&style=plastic)](https://docs.airr-community.org/en/stable/swtools/airr_swtools_standard.html)
### partis

Partis is an HMM-based framework for B- and T-cell receptor sequence annotation, simulation, clonal family and germline inference, and affinity prediction.
It is built on top of the [ham](https://github.com/psathyrella/ham) HMM compiler and [ig-sw](https://github.com/matsengrp/ig-sw) set of Smith-Waterman annotation tools.
Partis is free software under the GPL v3.

The various components are described in the following papers.
Since they do quite different things, it's best if you can cite the specific paper(s) that describe the components that you're using.

  * [Paired clustering](https://doi.org/10.1371/journal.pcbi.1010723) Ralph, DK, & Matsen IV, FA (2022). _Inference of B cell clonal families using heavy/light chain pairing information._ PLoS Computational Biology, 18(11), e1010723.
  * [Selection metrics](https://doi.org/10.1371/journal.pcbi.1008391) Ralph, DK, & Matsen IV, FA (2020). _Using B cell receptor lineage structures to predict affinity._ PLoS Computational Biology, 16(11), e1008391.
  * [Germline inference](https://doi.org/10.1371/journal.pcbi.1007133) Ralph, DK, & Matsen IV, FA (2019). _Per-sample immunoglobulin germline inference \[...\]._ PLoS Computational Biology, 15(7), e1007133.
  * [Clonal family inference](http://dx.doi.org/10.1371/journal.pcbi.1005086) Ralph, DK, & Matsen IV, FA (2016). _Likelihood-based Inference of B-cell Clonal Families._ PLoS Computational Biology, 12(10), e1005086.
  * [HMM framework and BCR annotation](http://doi.org/10.1371/journal.pcbi.1004409) Ralph, DK, & Matsen IV, FA (2016). _Consistency of \[...\] Enables Accurate B Cell Receptor Sequence Annotation._ PLoS Computational Biology, 12(1), e1004409.

The best place to start reading the manual is either [installation](docs/install.md) or [quick start](docs/quick-start.md), but there's also a [table of contents](docs/contents.md).
Details on the many options not documented there may be found by running each subcommand's help, for instance `partis cache-parameters --help` or `partis simulate --help`.

To ask questions about how to run or how things work, please use the [google group](https://groups.google.com/g/partis).
For specific issues with the software, e.g. bug reports or feature requests, on the other hand, [submit an issue](https://github.com/psathyrella/partis/issues/new) on github.
You can also search through past discussion both on the google group and in [closed issues](https://github.com/psathyrella/partis/issues?q=is%3Aissue+is%3Aclosed).
If you'd like to be notified when something important changes, please subscribe to the google group.
