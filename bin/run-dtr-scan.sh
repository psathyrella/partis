#!/bin/bash

action=test  # plot
trlabel=dtr-train-v1; trseed=0
# trlabel=dtr-train-v2; trseed=0 #1  # NOTE dtr-scan.py doesn't actually check what seed was used for training, it relies on us to pass in the correct one here
common=" --training-label $trlabel --training-seed $trseed --n-max-queries 250000" # --cgroup among-families --tvar affinity"
# ./test/dtr-scan.py $action --label dtr-train-v0 $common
# ./test/dtr-scan.py $action --label dtr-train-v1 $common
./test/dtr-scan.py $action --label dtr-train-v2 $common
