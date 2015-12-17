#!/bin/bash
leaf_mut_hum="--n-leaf-list 100:200:500 --mutation-multipliers 4 --humans A"

./bin/compare-partition-methods.py --actions simulate $leaf_mut_hum
./bin/compare-partition-methods.py --actions cache-simu-parameters $leaf_mut_hum
for isub in 0 1 2; do
    ./bin/compare-partition-methods.py --actions run-viterbi:vsearch-partition:naive-hamming-partition:partition --subset $isub --n-subsets 10 $leaf_mut_hum
    ./bin/compare-partition-methods.py --actions write-plots --subset $isub --n-subsets 10 --no-mixcr --no-changeo $leaf_mut_hum
done
./bin/compare-partition-methods.py --actions compare-subsets --humans A --n-subsets 3 --n-leaf-list 1:2:5:10:25:50:100:200:500 --mutation-multipliers 1:4 --no-mixcr --no-changeo
