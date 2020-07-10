#!/bin/bash

common="--paired-loci igh:igk --n-sim-events 1500 --n-procs 10 --seed 1" # --debug 1"
outdir=_output/paired-simulation
param_dir=$outdir/parameters

echo ./bin/partis simulate --parameter-dir $param_dir/igh --light-chain-parameter-dir $param_dir/igk $common
echo ./bin/partis simulate --parameter-dir $param_dir/igh --light-chain-parameter-dir $param_dir/igk $common --outfname $outdir/simu-igh.yaml --light-chain-outfname $outdir/simu-igk.yaml
echo ./bin/partis simulate --simulate-from-scratch $common
echo ./bin/partis simulate --simulate-from-scratch $common --outfname $outdir/simu-igh.yaml --light-chain-outfname $outdir/simu-igk.yaml
