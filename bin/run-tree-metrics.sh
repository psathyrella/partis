#!/bin/bash

bin=./test/cf-tree-metrics.py

# lbtl=0.0005:0.001:0.002:0.0025:0.003:0.004:0.005  # remove the larger ones since they diverge too much (but for tau optimization we do want them, since it extends the x range axis to the right)
# $bin get-lb-bounds --label v2 --lb-tau-list $lbtl  #--make-plots
# exit 0

testargs="--n-sim-seqs-per-gen-list 50:125 --lb-tau-list 0.002:0.003 --obs-times 100 --carry-cap 1000"
common="--carry-cap 1000 --obs-times 125,150 --n-sim-seqs-per-gen-list 30:50:75:100:150:200 --lb-tau-list 0.0005:0.001:0.002:0.003:0.005:0.008:0.012"

for action in run-bcr-phylo partition plot; do
    # $bin $action --label v0 --n-replicates 3 $common
    $bin $action --label test $testargs
done
