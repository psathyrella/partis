# ig-sw

Smith-Waterman alignment for immunoglobulin sequences.

## Introduction

The [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) (SW) alignment algorithm is a simple but effective way to align DNA sequences.
Immunoglobulin sequences derive from [VDJ recombination](https://en.wikipedia.org/wiki/V(D)J_recombination), and an important step in analyzing these sequences is to identify the recombination events that led to their creation.
Duncan's [partis](https://github.com/psathyrella/partis) package does this work, and does it well, using a [hidden Markov model](https://en.wikipedia.org/wiki/Hidden_Markov_model) approach.

It depends on a package called [ighutil](https://github.com/matsengrp/ighutil) that uses SW alignment to identify promising candidate alignments to investigate more closely.
Thus far in partis' development we have been treating ighutil as a black box, and just calling it via the command line [here](https://github.com/psathyrella/partis/blob/a6737936ee9fd191271485843ca2014b7fc109a7/python/waterer.py#L227).
However, we are currently reworking some internals of partis, and we need a more flexible SW implementation.

The long-term goal of this project is to build that more flexible implementation.
The first step towards that goal is to put together a C++ implementation of the single command that is needed for our purposes: `vdjalign align-fastq`.
*ighutil does many things, but we aren't concerned about anything else than this single command.*
This takes a DNA sequence in [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format), does alignment, and then writes the output in BAM format, which is a binary version of the [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).

"Unwinding, the stack", this call makes this sequence of calls for alignment:
* [align_fastq.py](https://github.com/matsengrp/ighutil/blob/master/python/vdjalign/subcommands/align_fastq.py), which by using the function `sw.ig_align` calls the Cython script...
* [sw.pyx](https://github.com/matsengrp/ighutil/blob/master/python/vdjalign/sw.pyx), which calls the `ig_align_reads` function of...
* [ig_align.c](https://github.com/matsengrp/ighutil/blob/master/python/src/ig_align.c) which uses much from the [klib](https://github.com/attractivechaos/klib) library.

At this point you don't need to dig any deeper for the moment-- we can just use that C.

The other part of the call encodes all of this alignment information into the BAM format using a command-line tool, `samtools view`.
It would be nicer to just use the [htslib](https://github.com/samtools/htslib) C library directly, which is what [is used to write samtools](https://github.com/samtools/samtools/blob/develop/sam_view.c#L231).
