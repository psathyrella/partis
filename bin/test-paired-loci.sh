#!/bin/bash

label=fix-bad-merge-v3 #2 # pairclean-v2  # re-refactor-v0 #refactor-v4  #paired-clustering-output-v2
# NOTE use tmp.sh for refactor-v{3,4}

nprocs=10
common="--n-sim-events 10 --n-leaves 5 --n-procs $nprocs --no-per-base-mutation --allowed-cdr3-lengths 30:33:36:42:45:48 --mutation-multiplier 1" # --debug 1" #  --constant-number-of-leaves
in_param_dir=_output/paired-simulation/parameters
outdir=_output/$label
out_param_dir=$outdir/params

mkdir -p $outdir  # for the .log files
# ./bin/partis simulate --paired-loci --seed 1 --parameter-dir $in_param_dir --paired-outdir $outdir/simu $common --mean-cells-per-droplet 1.1 >$outdir/simu.log
for action in cache-parameters partition; do  # merge-paired-partitions
    common="--n-procs $nprocs --is-simu"
    echo "./bin/partis $action --paired-loci --paired-indir $outdir/simu --input-metafname $outdir/simu/meta.yaml --paired-outdir $outdir/inferred $common >$outdir/$action.log"
done
