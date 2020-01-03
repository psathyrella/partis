#!/bin/bash

action=test #train  # plot
trlabel=dtr-train-v2
common=" --training-label $trlabel --n-max-queries 150000 --n-max-procs 7" # --cgroup among-families --tvar affinity"# ; trseed=0 --training-seed $trseed
./test/dtr-scan.py $action --label dtr-train-v0 $common
./test/dtr-scan.py $action --label dtr-train-v1 $common
./test/dtr-scan.py $action --label dtr-train-v2 $common
