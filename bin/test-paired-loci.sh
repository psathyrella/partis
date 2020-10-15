#!/bin/bash

label=paired-clustering-output-v2

nprocs=10
common="--n-sim-events 30 --n-leaves 5 --n-procs $nprocs --no-per-base-mutation --allowed-cdr3-lengths 33 --mutation-multiplier 3" # --debug 1" #  --constant-number-of-leaves
in_param_dir=_output/paired-simulation/parameters
outdir=_output/$label
out_param_dir=$outdir/params

# these are just to test the main ways of running simulation:
# ./bin/partis simulate --parameter-dir $in_param_dir/igh --light-chain-parameter-dir $in_param_dir/igk $common
# ./bin/partis simulate --parameter-dir $in_param_dir/igh --light-chain-parameter-dir $in_param_dir/igk $common --outfname $outdir/simu-igh.yaml --light-chain-outfname $outdir/simu-igk.yaml
# ./bin/partis simulate --simulate-from-scratch $common
# ./bin/partis simulate --simulate-from-scratch $common --outfname $outdir/simu-igh.yaml --light-chain-outfname $outdir/simu-igk.yaml
# actually using these for partitioning:
for lc in k l; do
    echo ./bin/partis simulate --paired-loci igh:ig$lc --seed 1 --parameter-dir $in_param_dir/igh --light-chain-parameter-dir $in_param_dir/ig$lc $common --outfname $outdir/h$lc/simu-igh.yaml --light-chain-outfname $outdir/h$lc/simu-ig$lc.yaml
    # echo ./bin/extract-fasta.py --input-file $outdir/h$lc/simu-igh.yaml --fasta-output-file $outdir/h$lc/simu-igh.fa
    # echo ./bin/extract-fasta.py --input-file $outdir/h$lc/simu-ig$lc.yaml --fasta-output-file $outdir/h$lc/simu-ig$lc.fa
done
# echo cat $outdir/h?/*.fa >$outdir/simu.fa

for lc in k l; do
    common="--is-simu --paired-loci igh:ig$lc --n-procs $nprocs"
    ifnstr="--infname $outdir/h$lc/simu-igh.yaml --light-chain-infname $outdir/h$lc/simu-ig$lc.yaml"
    pdstr="--parameter-dir $out_param_dir/igh --light-chain-parameter-dir $out_param_dir/ig$lc"
    ofnstr="--outfname $outdir/h$lc/partitions-igh.yaml --light-chain-outfname $outdir/h$lc/partitions-ig$lc.yaml"
    echo ./bin/partis cache-parameters $pdstr $common $ifnstr
    echo ./bin/partis partition $pdstr $common $ifnstr $ofnstr # --abbreviate
    echo ./bin/partis merge-paired-partitions $pdstr $common $ifnstr $ofnstr
done

# testing split-loci stuff:
# for action in cache-parameters partition; do
#     echo ./bin/partis $action --split-loci --infname $outdir/simu.fa --paired-loci-output-dir $outdir/split-loci-test/no-auto-cache # >no-auto-cache-$action.log
# done
# ./bin/partis partition --split-loci --infname $outdir/simu.fa --paired-loci-output-dir $outdir/split-loci-test/auto-cache-no-pdir >auto-cache-no-pdir.log
# ./bin/partis partition --split-loci --infname $outdir/simu.fa --paired-loci-output-dir $outdir/split-loci-test/auto-cache-with-pdir --parameter-dir $outdir/split-loci-test/auto-cache-with-pdir/x/y/z >auto-cache-with-pdir.log
