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

common="--only-csv-plots --slurm --n-sim-events-per-proc 10"
# TODO figure out v2 command, which now needs --legend-var obs_frac
# label=vary-carry-cap-v0
# echo $bin --label $label --n-replicates 10 $common --carry-cap-list 500:750:1000:2000:5000 --obs-times-list 100,200 --n-sim-seqs-per-gen-list 75
# label=vary-obs-times-v0
# echo $bin --label $label --n-replicates 10 $common --carry-cap-list 1000 --obs-times-list 100:200:300:100,150:200,250:100,200,300 --n-sim-seqs-per-gen-list 100:100:100:50:50:33 --zip-vars obs-times:n-sim-seqs-per-gen
# label=vary-obs-times-v1
# echo $bin --label $label --n-replicates 10 $common --carry-cap-list 1000 --obs-times-list 300:100,200,300:200,250,300 --n-sim-seqs-per-gen-list 100:33:33 --zip-vars obs-times:n-sim-seqs-per-gen
# label=vary-obs-frac-v0
# echo $bin --label $label --n-replicates 30 $common --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 30:50:75:100:150:200
# label=vary-metric-v0
# echo $bin --label $label --n-replicates 30 $common --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 100 --metric-for-target-distance-list nuc:aa:aa-sim-ascii:aa-sim-blosum
# label=vary-selection-strength-v0
# echo $bin --label $label --n-replicates 30 $common --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 100 --selection-strength-list 0.1:0.4:0.7:0.8:0.9:1.0
# label=carry-cap-vs-n-obs-v0
# echo $bin --label $label --n-replicates 30 $common --carry-cap-list 260:260:500:500:700:700:1500:1500:3000:3000 --obs-times-list 150 --n-sim-seqs-per-gen-list 13:26:25:50:35:70:75:150:150:300 --lb-tau-list 0.0025 --zip-vars carry-cap:n-sim-seqs-per-gen --final-plot-xvar carry-cap --legend-var obs_frac
# label=carry-cap-vs-n-obs-v1
# echo $bin --label $label --n-replicates 30 $common --carry-cap-list 250:500:1000:3000 --obs-times-list 150 --n-sim-seqs-per-gen-list 15:30:75:150:500 --lb-tau-list 0.0025 --final-plot-xvar carry-cap
# TODO everything before here was using relative affinity, but now I've added --use-relative-affy, for which the default is False. So either add the arg to the previous commands, or decide you'd rather have it turned off for them
label=carry-cap-vs-n-obs-only-leaves-v0
echo $bin --label $label --n-replicates 30 $common --carry-cap-list 250:500:1000:3000 --obs-times-list 150 --n-sim-seqs-per-gen-list 15:75:500 --lb-tau-list 0.0025 --dont-observe-common-ancestors --final-plot-xvar carry-cap
