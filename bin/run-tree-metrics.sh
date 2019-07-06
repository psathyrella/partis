#!/bin/bash

bin=./test/cf-tree-metrics.py

label=v0
testlabel=test-old
testargs="--n-sim-seqs-per-gen-list 50:125 --lb-tau-list 0.002:0.003 --obs-times 100 --carry-cap 1000 --n-generations-list 4:5"

# $bin get-lb-bounds --label $label  #--make-plots
# $bin get-lb-bounds --label $testlabel $testargs --make-plots
# exit 0

echo $bin --label $testlabel $testargs --only-csv-plots
echo $bin --label $label --n-replicates 3 $common --only-csv-plots
