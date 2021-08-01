#!/bin/bash

bin=./test/cf-paired-loci.py

for act in simulate cache-parameters partition; do
    actstr="--actions simulate"

    largs="--label n-final-test"
    echo $bin $largs $actstr --version v1 --n-sim-events 10  --n-leaf-list 30  --allowed-cdr3-lengths 33:36:42:45 --mutation-multiplier 7
    break
    # $bin $largs $actstr --version v2 --n-sim-events 100 --n-leaf-list 30  --allowed-cdr3-lengths 33:36:42:45 --mutation-multiplier 7
    # $bin $largs $actstr --version v3 --n-sim-events 3   --n-leaf-list 300 --allowed-cdr3-lengths 33:36:42:45 --mutation-multiplier 7
    # $bin $largs $actstr --version v4 --n-sim-events 10  --n-leaf-list 30  --allowed-cdr3-lengths 33:36:42:45 --mutation-multiplier 7 --cells-per-drop-list 1
done
