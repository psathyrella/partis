#!/bin/bash

bin=./test/cf-tree-metrics.py

label=v0
testlabel=test-old
testargs="--n-sim-seqs-per-gen-list 50:125 --lb-tau-list 0.002:0.003 --obs-times 100 --carry-cap 1000 --n-generations-list 4:5"

# $bin get-lb-bounds --label $label  #--make-plots
# $bin get-lb-bounds --label $testlabel $testargs --make-plots
# exit 0

# echo $bin --label $testlabel $testargs --only-csv-plots
# echo $bin --label $label --n-replicates 3 --only-csv-plots

# label=vary-carry-cap-v0
# echo $bin --label $label --n-replicates 10 --n-sim-events-per-proc 10 --only-csv-plots --carry-cap-list 500:750:1000:2000:5000 --obs-times-list 100,200 --n-sim-seqs-per-gen-list 75
# label=vary-obs-times-v0
# echo $bin --label $label --n-replicates 10 --n-sim-events-per-proc 10 --only-csv-plots --carry-cap-list 1000 --obs-times-list 100:200:300:100,150:200,250:100,200,300 --n-sim-seqs-per-gen-list 100:100:100:50:50:33 --zip-vars obs-times:n-sim-seqs-per-gen
# label=vary-obs-times-v1
# echo $bin --label $label --n-replicates 10 --n-sim-events-per-proc 10 --only-csv-plots --carry-cap-list 1000 --obs-times-list 300:100,200,300:200,250,300 --n-sim-seqs-per-gen-list 100:33:33 --zip-vars obs-times:n-sim-seqs-per-gen
# label=vary-obs-frac-v0
# echo $bin --label $label --n-replicates 30 --n-sim-events-per-proc 10 --only-csv-plots --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 30:50:75:100:150:200
label=vary-metric-v0
echo $bin --label $label --n-replicates 30 --n-sim-events-per-proc 10 --only-csv-plots --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 100 --metric-for-target-distance-list nuc:aa:aa-sim-ascii:aa-sim-blosum
