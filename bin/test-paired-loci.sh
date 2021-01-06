#!/bin/bash

label=fix-bad-merge-v0 # pairclean-v2  # re-refactor-v0 #refactor-v4  #paired-clustering-output-v2
# NOTE use tmp.sh for refactor-v{3,4}

nprocs=10
common="--n-sim-events 30 --n-leaves 5 --n-procs $nprocs --no-per-base-mutation --allowed-cdr3-lengths 33 --mutation-multiplier 1" # --debug 1" #  --constant-number-of-leaves
in_param_dir=_output/paired-simulation/parameters
outdir=_output/$label
out_param_dir=$outdir/params

echo ./bin/partis simulate --paired-loci --seed 1 --parameter-dir $in_param_dir --paired-outdir $outdir $common --mean-cells-per-droplet 1.1
for action in cache-parameters partition; do  # merge-paired-partitions
    common="--n-procs $nprocs --is-simu"
    echo ./bin/partis $action --paired-loci --paired-indir $outdir --input-metafname $outdir/meta.yaml --paired-outdir $outdir/inferred $common # >no-auto-cache-$action.log
done
