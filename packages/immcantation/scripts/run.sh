#!/bin/bash

VERSION="1.17.1"
ncbdir=packages/ncbi-igblast-${VERSION}
export PATH=$PATH:$ncbdir/bin

AssignGenes.py igblast -s test/mishmash.fa -b $PWD/packages/ncbi-igblast-1.17.1/igbl-db --organism human --loci ig --format airr

# airr: --outfmt 19
export IGDATA=~/share/igblast
igblastn \
    -germline_db_V ~/share/igblast/database/imgt_human_ig_v\
    -germline_db_D ~/share/igblast/database/imgt_human_ig_d \
    -germline_db_J ~/share/igblast/database/imgt_human_ig_v \
    -auxiliary_data ~/share/igblast/optional_file/human_gl.aux \
    -domain_system imgt -ig_seqtype Ig -organism human \
    -outfmt '7 std qseq sseq btop' \
    -query HD13M.fasta \
    -out HD13M.fmt7

