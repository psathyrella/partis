#!/bin/bash

bin=./test/cf-tree-metrics.py

# $bin get-max-lbi --label test-2d #--make-plots
# exit 0

# $bin run-bcr-phylo --label xxx
# exit 0

for action in run-bcr-phylo partition plot; do
    $bin $action --label v --carry-cap 1000 --obs-times 125,150 --n-sim-seqs-per-gen 30:50:75:100:150:200 --n-replicates 3
done

# $bin run-bcr-phylo --label b --carry-cap 300:400 --n-sim-seqs-per-gen 10,20:30,40 --obs-times 40,50:100,150 --lb-tau-list 0.001:0.002
