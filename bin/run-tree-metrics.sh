#!/bin/bash

bin=./test/cf-tree-metrics.py

lbtl=0.0005:0.001:0.002:0.0025:0.003:0.004:0.005  # remove the larger ones since they diverge too much (but for tau optimization we do want them, since it extends the x range axis to the right)
$bin get-lbi-bounds --label v1 --lb-tau-list $lbtl  #--make-plots
exit 0

# $bin run-bcr-phylo --label xxx
# exit 0

common="--carry-cap 1000 --obs-times 125,150 --n-sim-seqs-per-gen 30:50:75:100:150:200 --lb-tau-list 0.0005:0.001:0.002:0.003:0.005:0.008:0.012"
for action in run-bcr-phylo partition plot; do
    $bin $action --label test-v0 --n-replicates 3 $common
done

# $bin run-bcr-phylo --label b --carry-cap 300:400 --n-sim-seqs-per-gen 10,20:30,40 --obs-times 40,50:100,150 --lb-tau-list 0.001:0.002
