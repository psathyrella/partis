#!/bin/bash

bin=./test/cf-tree-metrics.py

label=XXX
testlabel=test-v2
testargs="--n-sim-seqs-per-gen-list 50:125 --lb-tau-list 0.002:0.003 --obs-times 100 --carry-cap 1000 --n-generations-list 4:5"

# $bin get-lb-bounds --label v2  #--make-plots
# $bin get-lb-bounds --label $testlabel $testargs --make-plots
# exit 0

for action in run-bcr-phylo partition plot; do
    # $bin $action --label v0 --n-replicates 3 $common
    $bin $action --label $testlabel $testargs
done
