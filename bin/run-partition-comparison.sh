#!/bin/bash

# ----------------------------------------------------------------------------------------
# thresholds optimization

humans=A  # 021-018:021-019  #:B:C:021-018:021-044"
# ./bin/compare-partition-methods.py --actions cache-data-parameters --humans $humans  &
# leaf_mut_hum="--n-leaf-list 3:11:51 --mutation-multipliers 0.5:1:2:3 --humans $humans"  # NOTE I accidentally overwrote the new 10-leaf sample (with an older one), so there's threshold optimization plots for a 10-leaf-1.0-mutate sim file that no longer exists. If you absolutely have to remake the plots, use this 11-leaf sample (note that the 0.5, 2, and 3x 10-leaf samples are still here, and are from threshold optimization)
# mimic stuff # leaf_mut_hum="--n-leaf-list 3:10 --mutation-multipliers 1 --humans $humans"  # mimic stuff
leaf_mut_hum="--n-leaf-list 10 --mutation-multipliers 0.5:3 --humans $humans"
# ./bin/compare-partition-methods.py --actions simulate $leaf_mut_hum &
# ./bin/compare-partition-methods.py --actions cache-simu-parameters $leaf_mut_hum &
# action=naive-hamming-partition
# tholds="0.015 0.02 0.025 0.03 0.04 0.05 0.07 0.09 0.11 0.14 0.18"
action=vsearch-partition
tholds="0.015 0.03 0.05 0.07 0.11"
# action=partition
# tholds="0 5 10 12.5 15 16 17.5 19 20 22.5 25 30 35"
for th in $tholds; do
    if [ "$bounds" == "" ]; then
	bounds="$th,$th"
    else
	bounds="$bounds:$th,$th"
    fi
done
./bin/compare-partition-methods.py --actions $action $leaf_mut_hum --hfrac-bound-list $bounds --n-to-partition 3000 &  # --istartstop 0:10000 &  #0:1000 1000:8000 &
# ./bin/compare-partition-methods.py --actions write-plots $leaf_mut_hum --hfrac-bound-list $bounds --expected-methods $action &  # --istartstop 0:10000  &  #1000:8000 &
# ./bin/compare-partition-methods.py --actions compare-subsets $leaf_mut_hum --hfrac-bound-list $bounds --expected-methods $action &  # --istartstop 0:10000 &  #1000:8000 &

# # ----------------------------------------------------------------------------------------
# # leaves=1:2:5:10:25:50:100:200  #:500
# leaves=100:200
# leaf_mut_hum="--n-leaf-list $leaves --mutation-multipliers 1:4 --humans A"
# # echo "dont forget you copied these sim files from the old dir";  # ./bin/compare-partition-methods.py --actions simulate $leaf_mut_hum  #  --n-sim-seqs 100000
# ./bin/compare-partition-methods.py --actions cache-simu-parameters $leaf_mut_hum &
# # for isub in 0 1 2; do
# #     # ./bin/compare-partition-methods.py --actions vsearch-partition:naive-hamming-partition:partition --subset $isub --n-subsets 10 $leaf_mut_hum --overwrite &  #run-viterbi:
# #     # ./bin/compare-partition-methods.py --actions run-changeo --subset $isub --n-subsets 10 $leaf_mut_hum
# #     # sleep 60
# #     ./bin/compare-partition-methods.py --actions write-plots --subset $isub --n-subsets 10 --no-mixcr --no-changeo $leaf_mut_hum --no-similarity-matrices &
# #     # break
# # done
# # ./bin/compare-partition-methods.py --actions compare-subsets --plot-mean-of-subsets --n-subsets 3 $leaf_mut_hum --no-mixcr --no-changeo
# # ./bin/compare-partition-methods.py --actions compare-subsets --plot-mean-of-subsets --humans A --n-subsets 3 --n-leaf-list 10:25 --mutation-multipliers 1 --no-mixcr --no-changeo

# # ----------------------------------------------------------------------------------------
# leaf_mut_hum="--n-leaf-list 7 --mutation-multipliers 1 --humans A"  # 100:200:500
# # # istartstoplist="0:100 100:600 600:1350 1350:2350 2350:3850 3850:5850 5850:8850 8850:12850 12850:17850 17850:25350 25350:35350 35350:50350"  # 50350:70350 70350:100350"
# istartstoplist="5850:8850"
# istartstopstr=`echo $istartstoplist | sed -e 's/:/,/g' -e 's/ /:/g'`
# for istartstop in $istartstoplist; do  # see code below to generate these
#     # ./bin/compare-partition-methods.py --actions vsearch-partition:naive-hamming-partition:partition --istartstop $istartstop $leaf_mut_hum --overwrite &  # run-viterbi: --count-distances
#     # sleep 60
#     ./bin/compare-partition-methods.py --actions write-plots --istartstop $istartstop --no-mixcr --no-changeo $leaf_mut_hum --no-similarity-matrices --count-distances &
#     # break
# done
# # ./bin/compare-partition-methods.py --actions compare-subsets --istartstoplist $istartstopstr $leaf_mut_hum --no-mixcr --no-changeo

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

