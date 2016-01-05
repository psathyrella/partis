#!/bin/bash

# quick test command
# ./bin/compare-partition-methods.py --actions write-plots --subset 0 --n-subsets 10 --no-mixcr --no-changeo --no-similarity-matrices --humans A

# leaves=1:2:5:10:25:50:100:200  #:500
# leaf_mut_hum="--n-leaf-list $leaves --mutation-multipliers 1:4 --humans A"
# ./bin/compare-partition-methods.py --actions simulate $leaf_mut_hum  #  --n-sim-seqs 100000
# ./bin/compare-partition-methods.py --actions cache-simu-parameters $leaf_mut_hum
# for isub in 0 1 2; do
#     # ./bin/compare-partition-methods.py --actions run-viterbi:vsearch-partition:naive-hamming-partition:partition --subset $isub --n-subsets 10 $leaf_mut_hum &
#     # ./bin/compare-partition-methods.py --actions run-changeo --subset $isub --n-subsets 10 $leaf_mut_hum
#     # sleep 60
#     ./bin/compare-partition-methods.py --actions write-plots --subset $isub --n-subsets 10 --no-mixcr --no-changeo $leaf_mut_hum &  # --no-similarity-matrices &
#     # break
# done
# ./bin/compare-partition-methods.py --actions compare-subsets --plot-mean-of-subsets --n-subsets 3 $leaf_mut_hum --no-mixcr --no-changeo
# ./bin/compare-partition-methods.py --actions compare-subsets --plot-mean-of-subsets --humans A --n-subsets 3 --n-leaf-list 10:25 --mutation-multipliers 1 --no-mixcr --no-changeo

leaf_mut_hum="--n-leaf-list 7 --mutation-multipliers 1 --humans A"  # 100:200:500
# istartstoplist="0:100 100:600 600:1350 1350:2350 2350:3850 3850:5850 5850:8850 8850:12850 12850:17850 17850:25350 25350:35350 35350:50350"  # 50350:70350 70350:100350"
istartstoplist="5850:8850"
istartstopstr=`echo $istartstoplist | sed -e 's/:/,/g' -e 's/ /:/g'`
for istartstop in $istartstoplist; do  # see code below to generate these
    # ./bin/compare-partition-methods.py --actions run-viterbi:vsearch-partition:naive-hamming-partition:partition --istartstop $istartstop $leaf_mut_hum --count-distances & # --overwrite &
    # sleep 60
    ./bin/compare-partition-methods.py --actions write-plots --istartstop $istartstop --no-mixcr --no-changeo $leaf_mut_hum  --no-similarity-matrices --count-distances &
    # break
done
# ./bin/compare-partition-methods.py --actions compare-subsets --istartstoplist $istartstopstr $leaf_mut_hum --no-mixcr --no-changeo

# istart=0
# for s in 100 500 750 1000 1500 2000 3000 4000 5000 7500 10000 15000 20000 30000; do
#     ((istop = istart + s))
#     echo -n " $istart:$istop"
#     ((istart += s))
# done

# rerun with three different 750-seq subsets
# run once:
# 0.965 0.96 0.90
# 0.991 0.98 0.96
# 0.986 0.96 0.95
# run again (basically different random seed):
# 0.964 0.97 0.91
# 0.989 0.97 0.98
# 0.953 0.86 0.97
