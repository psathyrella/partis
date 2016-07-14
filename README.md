# ig-sw

Smith-Waterman alignment for immunoglobulin sequences.

## Introduction

This is a C++ implementation of the  `vdjalign align-fastq` command in [ighutil](https://github.com/cmccoy/ighutil). This implementation uses the [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) (SW) alignment algorithm, a simple but effective way to align DNA sequences, as in ighutil.

[Docker Hub Automated Build](https://hub.docker.com/r/matsengrp/ig-sw/builds/)

## Prerequisites

- [zlib](http://www.zlib.net/)
- [scons](http://scons.org/)
- [pthread (POSIX Thread)](https://computing.llnl.gov/tutorials/pthreads/)
- C++03-compatible compiler (tested against [g++](https://gcc.gnu.org/) 4.8)
- a smile


## Running the alignment

Clone the repository and run `scons` from the  `ig-sw/src/ig_align/` directory.

To run the sw-alignment, run `./ig-sw` with the appropriate arguments in the same directory or add a path to the executable.

### Example

`./ig-sw -p test_data/ query_file.fasta output_file.sam`

### Arguments

ig-sw uses the [tclap](http://tclap.sourceforge.net/) command line interface.

##### Required

| Flag         | Description                              |
| ------------ | ---------------------------------------- |
| (none)       | The path to the query file               |
| (none)       | The path to the output file              |
| -p --vdj-dir | Directory from which to read germline genes (ending with /) |

##### Optional

| Flag            | Description                              | Default Value |
| :-------------- | ---------------------------------------- | ------------- |
| -l --locus      | Locus to align reads against (options: IGH, IGK, IGL, DJ) | IGH           |
| -m --match      | Match score                              | 2             |
| -u --mismatch   | Mismatch score                           | 2             |
| -o --gap-open   | Gap opening penalty                      | 3             |
| -e --gap-extend | Gap extension penalty                    | 1             |
| -d --max-drop   | Max drop                                 | 1000          |
| -s --min-score  | Min score                                | 0             |
| -b --bandwidth  | Bandwidth                                | 150           |
| -j --threads    | Number of threads                        | 1             |


## Workflow

[ig_align_main.cpp](https://github.com/matsengrp/ig-sw/blob/master/src/ig_align/ig_align_main.cpp) takes in a DNA sequence in [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format), does alignment, and outputs a file in SAM format.

ig_align_main.cpp calls the `ig_align_reads` function of [ig_align.c](https://github.com/matsengrp/ig-sw/blob/master/src/ig_align/ig_align.c) which heavily utilizes the [klib](https://github.com/attractivechaos/klib) library.
