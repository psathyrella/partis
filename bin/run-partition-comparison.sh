#!/bin/bash
# leaf_mut_hum="--n-leaf-list 1:2:5:10:25:50:100:200 --mutation-multipliers 1:4 --humans A"  # 100:200:500
leaf_mut_hum="--n-leaf-list 7 --mutation-multipliers 1 --humans A"  # 100:200:500

# ./bin/compare-partition-methods.py --actions simulate $leaf_mut_hum  #  --n-sim-seqs 100000
# ./bin/compare-partition-methods.py --actions cache-simu-parameters $leaf_mut_hum
# for isub in 0 1 2; do
#     ./bin/compare-partition-methods.py --actions run-viterbi:vsearch-partition:naive-hamming-partition:partition --subset $isub --n-subsets 10 $leaf_mut_hum &
#     sleep 60
#     # ./bin/compare-partition-methods.py --actions write-plots --subset $isub --n-subsets 10 --no-mixcr --no-changeo $leaf_mut_hum
#     # break
# done
# ./bin/compare-partition-methods.py --actions compare-subsets --n-subsets 3 $leaf_mut_hum --no-mixcr --no-changeo  # 
# ./bin/compare-partition-methods.py --actions compare-subsets --humans A --n-subsets 3 --n-leaf-list 25 --mutation-multipliers 1 --no-mixcr --no-changeo  # :500

# for istartstop in 0:100 100:600 600:1350 1350:2350 2350:3850 3850:5850 5850:8850 8850:12850 12850:17850 17850:25350; do  # see code below to generate these
for istartstop in 12850:17850 17850:25350; do
    ./bin/compare-partition-methods.py --actions run-viterbi:vsearch-partition:naive-hamming-partition:partition --istartstop $istartstop $leaf_mut_hum &
    sleep 600
    # ./bin/compare-partition-methods.py --actions write-plots --subset $isub --n-subsets 10 --no-mixcr --no-changeo $leaf_mut_hum
    # break
done

# istart=0
# for s in 100 500 750 1000 1500 2000 3000 4000 5000 7500; do
#     ((istop = istart + s))
#     echo -n " $istart:$istop"
#     ((istart += s))
# done
