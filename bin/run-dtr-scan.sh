#!/bin/bash

# action=test #train  # plot
# for trv in v0 v1 v2 v3; do
#     trlabel=dtr-train-$trv
#     common=" --training-label $trlabel --n-max-queries 150000 --n-max-procs 7" # --cgroup among-families --tvar affinity"# ; trseed=0 --training-seed $trseed
#     ./test/dtr-scan.py $action --label dtr-train-v3 $common
# done
# exit 0

action=plot #test  # train
for cg in within-families among-families; do
    for tv in affinity delta-affinity; do
	lfn=_output/dtr-scan/$cg-$tv.txt
	echo "" >$lfn
	for trv in v0 v1 v2 v3; do
	    trlabel=dtr-train-$trv
	    common=" --training-label $trlabel --n-max-queries 150000 --n-max-procs 7 --cgroup $cg --tvar $tv"  # ; trseed=0 --training-seed $trseed
	    ./test/dtr-scan.py $action --label dtr-train-v0 $common >>$lfn
	    ./test/dtr-scan.py $action --label dtr-train-v1 $common >>$lfn
	    ./test/dtr-scan.py $action --label dtr-train-v2 $common >>$lfn
	    ./test/dtr-scan.py $action --label dtr-train-v3 $common >>$lfn
	done
    done
done
