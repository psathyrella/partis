#!/bin/bash

common="--n-sim-events 10 --n-leaves 15 --constant-number-of-leaves --n-procs 1 --no-per-base-mutation --allowed-cdr3-lengths 33" # --debug 1"
param_dir=_output/paired-simulation/parameters
outdir=_output/paired-clustering-test
# these are just to test the main ways of running simulation:
# ./bin/partis simulate --parameter-dir $param_dir/igh --light-chain-parameter-dir $param_dir/igk $common
# ./bin/partis simulate --parameter-dir $param_dir/igh --light-chain-parameter-dir $param_dir/igk $common --outfname $outdir/simu-igh.yaml --light-chain-outfname $outdir/simu-igk.yaml
# ./bin/partis simulate --simulate-from-scratch $common
# ./bin/partis simulate --simulate-from-scratch $common --outfname $outdir/simu-igh.yaml --light-chain-outfname $outdir/simu-igk.yaml
# actually using these for partitioning:
# for lc in k l; do
#     echo ./bin/partis simulate --paired-loci igh:ig$lc --seed 1 --parameter-dir $param_dir/igh --light-chain-parameter-dir $param_dir/ig$lc $common --outfname $outdir/h$lc/simu-igh.yaml --light-chain-outfname $outdir/h$lc/simu-ig$lc.yaml
#     echo ./bin/extract-fasta.py --input-file $outdir/h$lc/simu-igh.yaml --fasta-output-file $outdir/h$lc/simu-igh.fa
#     echo ./bin/extract-fasta.py --input-file $outdir/h$lc/simu-ig$lc.yaml --fasta-output-file $outdir/h$lc/simu-ig$lc.fa
# done
# echo cat $outdir/h?/*.fa >$outdir/simu.fa

common="--paired-loci igh:igk --n-procs 10"
# echo ./bin/partis partition --is-simu --parameter-dir $param_dir/igh --light-chain-parameter-dir $param_dir/igk $common --infname $outdir/simu-igh.yaml --light-chain-infname $outdir/simu-igk.yaml --outfname $outdir/partitions-igh.yaml --light-chain-outfname $outdir/partitions-igk.yaml # --abbreviate
# echo ./bin/partis merge-paired-partitions --is-simu $common --infname $outdir/simu-igh.yaml --light-chain-infname $outdir/simu-igk.yaml --outfname $outdir/partitions-igh.yaml --light-chain-outfname $outdir/partitions-igk.yaml

for action in cache-parameters partition; do
    echo ./bin/partis $action --split-loci --infname $outdir/simu.fa --paired-loci-output-dir $outdir/split-loci-test/no-auto-cache # >no-auto-cache-$action.log
done
# ./bin/partis partition --split-loci --infname $outdir/simu.fa --paired-loci-output-dir $outdir/split-loci-test/auto-cache-no-pdir >auto-cache-no-pdir.log
# ./bin/partis partition --split-loci --infname $outdir/simu.fa --paired-loci-output-dir $outdir/split-loci-test/auto-cache-with-pdir --parameter-dir $outdir/split-loci-test/auto-cache-with-pdir/x/y/z >auto-cache-with-pdir.log
