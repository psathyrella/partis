# ighutil

Tools for working with immunoglobulin heavy chain sequences at high throughput

The project is divided into two sections: a python module for aligning sequences, and a JVM-based tool for post-processing.

# Prerequisites

* [samtools](http://samtools.sourceforge.net/)
* Java 7
* Python 2.7
* zlib
* A C++11-compatible compiler (tested against gcc 4.8)
* [Cython](http://cython.org) (tested against v0.20)

# Installation

* Run `make -C clj` - creates `clj/bin/ighutil`
* Run `pip install ./python` - installs `vdjalign` Python tool.

Alternatively, Debian packages `ighutil` and `vdjalign` can be installed from this repository:

    deb     https://s3.amazonaws.com/cmccoy-debian-repo/ stable main

Packages are signed with key `3AA09EC1`. Add with:

    sudo apt-key adv --keyserver pgp.mit.edu --recv-keys 3AA09EC1

Then:

    sudo apt-get update
    sudo apt-get install -y ighutil vdjalign

# Workflow

## Alignment

IGH reads are aligned via the included python package, via the command `ighutil align-fastq`.
The alignment process is simple: each sequence is aligned to all IGHV genes, keeping all alignments within a user-specified value of the maximum alignment score.

Alignment is multithreaded, and SSE3 accelerated using the excellent Smith-Waterman implementation from [klib](https://github.com/attractivechaos/klib).

The remainder of the sequence is aligned against IGHJ genes; then the remainder from the best IGHJ hit against IGHD genes.

The resulting alignments are stored in BAM format.

## Choosing a V/D/J assignment

Remaining post-processing uses subcommands of the JVM interface.

* `ighutil reset-primary` randomly assigns a primary alignment among ties for each segment, and attempts to filter out alleles not present in the data.
* `ighutil calculate-match-probability` calculates the expectation that each aligned base matches germline, treating all maximally scoring alignmnets as equiprobable.
* `samtools view -F 260` will subset to primary alignments only.

## Annotating sequences

Annotate sequences with V/D/J gene assignment, productivity status, and more using:

    ighutil annotate-vdj -i <BAM_FILE> -o <CSV_FILE>
