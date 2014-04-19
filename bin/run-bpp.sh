#!/bin/bash

infile=$1
outfile=$2

bpp_dir=$HOME/work/bpp-master-20140414
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$bpp_dir/lib
export PATH=$PATH:$bpp_dir/bin

# count the sequence length
n_sites=`wc -l $infile | cut -f1 -d' '`
((n_sites--))

bppseqgen \
    input.tree.file=$PWD/data/test.tre \
    output.sequence.file=$outfile \
    number_of_sites=$n_sites \
    input.tree.format=Newick \
    output.sequence.format=Fasta\(\) \
    alphabet=DNA \
    --seed=$$ \
    model=JC69 \
    rate_distribution='Constant()' \
    input.infos.states=state \
    input.infos=$infile \
    input.infos.rates=none
