#!/bin/bash

VERSION="1.17.1"
ncbdir=packages/ncbi-igblast-${VERSION}
export PATH=$PATH:$ncbdir/bin

AssignGenes.py igblast -s test/mishmash.fa -b $PWD/packages/ncbi-igblast-1.17.1/igbl-db --organism human --loci ig --format blast
