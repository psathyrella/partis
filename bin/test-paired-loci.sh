#!/bin/bash

ploci="--paired-loci igh:igk"

common="$ploci --n-sim-events 100 --n-leaves 15 --n-procs 10 --seed 1 --no-per-base-mutation --allowed-cdr3-lengths 21:27:33" # --debug 1"
outdir=_output/paired-simulation
param_dir=$outdir/parameters
# ./bin/partis simulate --parameter-dir $param_dir/igh --light-chain-parameter-dir $param_dir/igk $common
echo ./bin/partis simulate --parameter-dir $param_dir/igh --light-chain-parameter-dir $param_dir/igk $common --outfname $outdir/simu-igh.yaml --light-chain-outfname $outdir/simu-igk.yaml
# ./bin/partis simulate --simulate-from-scratch $common
# ./bin/partis simulate --simulate-from-scratch $common --outfname $outdir/simu-igh.yaml --light-chain-outfname $outdir/simu-igk.yaml

common="$ploci --n-procs 10"
echo ./bin/partis partition --is-simu --parameter-dir $param_dir/igh --light-chain-parameter-dir $param_dir/igk $common --infname $outdir/simu-igh.yaml --light-chain-infname $outdir/simu-igk.yaml --outfname $outdir/partitions-igh.yaml --light-chain-outfname $outdir/partitions-igk.yaml # --abbreviate
