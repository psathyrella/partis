#!/bin/bash

bin=./test/cf-paired-loci.py

for act in simulate cache-parameters partition; do
    actstr="--actions simulate"

    largs="--label n-final-test"
    echo $bin $largs $actstr --version v1 --n-sim-events 10  --n-leaf-list 30  --allowed-cdr3-lengths 33:36:42:45 --mutation-multiplier 7
    $bin $largs $actstr --version v2 --n-sim-events 100 --n-leaf-list 30  --allowed-cdr3-lengths 33:36:42:45 --mutation-multiplier 7
    $bin $largs $actstr --version v3 --n-sim-events 3   --n-leaf-list 300 --allowed-cdr3-lengths 33:36:42:45 --mutation-multiplier 7
    $bin $largs $actstr --version v4 --n-sim-events 10  --n-leaf-list 30  --allowed-cdr3-lengths 33:36:42:45 --mutation-multiplier 7 --cells-per-drop-list 1
done
./test/cf-paired-loci.py --n-leaves-list 5:10 --label vs-shm --version v0 --n-replicates 5 --n-sim-events-list 1000 --scratch-mute-freq-list 0.01:0.05:0.10:0.30:0.50:0.80 --n-sub-procs 7 --n-max-procs 5 --single-light-locus igk --base-outdir /fh/local/dralph/partis/paired-loci --actions cache-parameters:partition:vsearch-partition:vjcdr3-0.9:synth-distance-0.03:synth-singletons-0.20 --dry  --final-plot-xvar scratch-mute-freq --combo-extra-str scratch-x-var --pvks-to-plot 10
./test/cf-paired-loci.py --label vs-n-leaves --version v0 --n-replicates 3 --n-leaves-list 1:5:10:50 --simu-extra-args="--constant-number-of-leaves" --n-sim-events-list 1000 --n-sub-procs 5 --n-max-procs 5 --single-light-locus igk --base-outdir /fh/local/dralph/partis/paired-loci --actions simu:cache-parameters:partition:vsearch-partition:vjcdr3-0.9:synth-distance-0.03:synth-singletons-0.20 --dry
