#!/bin/bash

bd=/home/dralph/work/partis/packages/ncbi-igblast-1.17.1/igbl-db/fasta

for locus in igh igk igl; do
    mkdir -p $bd/partis-format/$locus
    for rgn in v d j; do
	grep -iA1 $locus$rgn $bd/imgt_human_ig_$rgn.fasta >$bd/partis-format/$locus/$locus$rgn.fasta
    done
done
