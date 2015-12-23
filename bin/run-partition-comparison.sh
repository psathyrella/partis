#!/bin/bash
# leaves=1:2:5:10:25:50:100:200  #:500
# leaf_mut_hum="--n-leaf-list $leaves --mutation-multipliers 4 --humans A"

# ./bin/compare-partition-methods.py --actions simulate $leaf_mut_hum  #  --n-sim-seqs 100000
# ./bin/compare-partition-methods.py --actions cache-simu-parameters $leaf_mut_hum
# for isub in 0 1 2; do
#     # ./bin/compare-partition-methods.py --actions run-viterbi:vsearch-partition:naive-hamming-partition:partition --subset $isub --n-subsets 10 $leaf_mut_hum &
#     ./bin/compare-partition-methods.py --actions run-changeo --subset $isub --n-subsets 10 $leaf_mut_hum
#     # sleep 60
#     # ./bin/compare-partition-methods.py --actions write-plots --subset $isub --n-subsets 10 --no-mixcr --no-changeo $leaf_mut_hum
#     break
# done
# ./bin/compare-partition-methods.py --actions compare-subsets --plot-mean-of-subsets --n-subsets 3 $leaf_mut_hum --no-mixcr --no-changeo
# ./bin/compare-partition-methods.py --actions compare-subsets --plot-mean-of-subsets --humans A --n-subsets 3 --n-leaf-list 10:25 --mutation-multipliers 1 --no-mixcr --no-changeo  # :500

leaf_mut_hum="--n-leaf-list 7 --mutation-multipliers 1 --humans A"  # 100:200:500
istartstoplist="0:100 100:600 600:1350 1350:2350 2350:3850 3850:5850 5850:8850 8850:12850 12850:17850 17850:25350 25350:35350 35350:50350"
istartstopstr=`echo $istartstoplist | sed -e 's/:/,/g' -e 's/ /:/g'`
# for istartstop in $istartstoplist; do  # see code below to generate these
for istartstop in 1000:1750 2000:2750 3000:3750; do
    ./bin/compare-partition-methods.py --actions partition --istartstop $istartstop $leaf_mut_hum &
    # sleep 600
    # ./bin/compare-partition-methods.py --actions write-plots --istartstop $istartstop --no-mixcr --no-changeo $leaf_mut_hum
    # break
done
# ./bin/compare-partition-methods.py --actions compare-subsets --istartstoplist $istartstopstr $leaf_mut_hum --no-mixcr --no-changeo

# istart=0
# for s in 100 500 750 1000 1500 2000 3000 4000 5000 7500 10000 15000; do
#     ((istop = istart + s))
#     echo -n " $istart:$istop"
#     ((istart += s))
# done
