#!/bin/bash

input_seq_file=$1
treefile=$2
outfile=$3

bpp_dir=$HOME/Dropbox/work/bpp-master-20140414
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$bpp_dir/lib
export PATH=$PATH:$bpp_dir/bin

# count the sequence length
n_sites=`wc -l $input_seq_file | cut -f1 -d' '`
((n_sites--))

bppseqgen \
    input.tree.file=$treefile \
    output.sequence.file=$outfile \
    number_of_sites=$n_sites \
    input.tree.format=Newick \
    output.sequence.format=Fasta\(\) \
    alphabet=DNA \
    --seed=$$ \
    model=JC69 \
    rate_distribution='Constant()' \
    input.infos.states=state \
    input.infos=$input_seq_file \
    input.infos.rates=none #| grep 'Random seed' | cut -f5 -d' '  # grep for the random seed -- a.t.m. it's the only thing I need in the output

if [[ $? != 0 ]]; then  # if it failed, return the pid (random seed)
    exit $$
else
    exit 0
fi
