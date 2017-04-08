#!/bin/bash

# # ----------------------------------------------------------------------------------------
# # weird
# lf=10:50 #:3  #30:150
# leaf_mut_hum="--n-leaf-list $lf --mutation-multipliers 1 --humans A"
# for loc in v cdr3; do
# xtra="--indels --indel-location $loc"  # "--mimic"  # --box
# # ./bin/compare-partition-methods.py --actions simulate $leaf_mut_hum $xtra &
# # ./bin/compare-partition-methods.py --actions cache-simu-parameters $leaf_mut_hum $xtra &
# # for isub in 0 1 2; do
# #     xtra="$xtra --subset $isub --n-subsets 10"
# #     # ./bin/compare-partition-methods.py --actions run-viterbi:vsearch-partition:naive-hamming-partition:partition $leaf_mut_hum $xtra   # --n-to-partition 1000 &
# #     # ./bin/compare-partition-methods.py --actions run-mixcr $leaf_mut_hum $xtra --n-to-partition 1000 &
# #     # ./bin/compare-partition-methods.py --actions write-plots --expected-methods vollmers-0.9:partition:naive-hamming-partition:vsearch-partition $leaf_mut_hum $xtra --no-similarity-matrices & #  --n-to-partition 1000 --no
# #     # ./bin/compare-partition-methods.py --actions write-plots --expected-methods vollmers-0.9:partition:naive-hamming-partition:vsearch-partition $leaf_mut_hum $xtra & # --no-similarity-matrices & #  --n-to-partition 1000 --no
# # done
# ./bin/compare-partition-methods.py --actions compare-subsets --plot-mean-of-subsets --n-subsets 3 $leaf_mut_hum $xtra --expected-methods vollmers-0.9:partition:naive-hamming-partition:vsearch-partition --no-similarity-matrices &
# done
# exit 0

# # ----------------------------------------------------------------------------------------
# # distance plots
# # nl=3 #:13:51
# mm=3  #0.5:1:2:3
# for nl in 3 51; do
#     leaf_mut_hum="--n-leaf-list $nl --mutation-multipliers $mm --humans A"
#     ./bin/compare-partition-methods.py --actions write-plots $leaf_mut_hum --count-distances --bak --no-similarity-matrices &
# done

# ----------------------------------------------------------------------------------------
# thresholds optimization

# humans=A  # 021-018:021-019  #:B:C:021-018:021-044"
# # ./bin/compare-partition-methods.py --actions cache-data-parameters --humans $humans  &
# # leaf_mut_hum="--n-leaf-list 3:11:51 --mutation-multipliers 0.5:1:2:3 --humans $humans"  # NOTE I accidentally overwrote the new 10-leaf sample (with an older one), so there's threshold optimization plots for a 10-leaf-1.0-mutate sim file that no longer exists. If you absolutely have to remake the plots, use this 11-leaf sample (note that the 0.5, 2, and 3x 10-leaf samples are still here, and are from threshold optimization)
# # mimic stuff # leaf_mut_hum="--n-leaf-list 3:10 --mutation-multipliers 1 --humans $humans"  # mimic stuff
# leaf_mut_hum="--n-leaf-list 3:10:51 --mutation-multipliers 0.5:3 --humans $humans"
# # ./bin/compare-partition-methods.py --actions simulate $leaf_mut_hum &
# # ./bin/compare-partition-methods.py --actions cache-simu-parameters $leaf_mut_hum &
# # action=naive-hamming-partition
# # tholds="0.015 0.02 0.025 0.03 0.04 0.05 0.07 0.09 0.11 0.14 0.18"
# action=vsearch-partition
# tholds="0.005 0.01 0.015 0.03 0.05 0.07 0.11 0.13"
# # action=partition
# # tholds="0 5 10 12.5 15 16 17.5 19 20 22.5 25 30 35"
# for th in $tholds; do
#     if [ "$bounds" == "" ]; then
# 	bounds="$th,$th"
#     else
# 	bounds="$bounds:$th,$th"
#     fi
# done
# # ./bin/compare-partition-methods.py --actions $action $leaf_mut_hum --hfrac-bound-list $bounds --n-to-partition 3000 &  # --istartstop 0:10000 &  #0:1000 1000:8000 &
# # ./bin/compare-partition-methods.py --actions write-plots $leaf_mut_hum --hfrac-bound-list $bounds --expected-methods $action &  # --istartstop 0:10000  &  #1000:8000 &
# ./bin/compare-partition-methods.py --actions compare-subsets $leaf_mut_hum --hfrac-bound-list $bounds --expected-methods $action &  # --istartstop 0:10000 &  #1000:8000 &

# # ----------------------------------------------------------------------------------------
# # leaves=1:2:5:10:25:50:100:200  #:500
# leaves=2:50:200  #:100
# mm=1:4
# leaf_mut_hum="--n-leaf-list $leaves --mutation-multipliers $mm --humans A"
# # echo "dont forget you copied these sim files from the old dir";  # ./bin/compare-partition-methods.py --actions simulate $leaf_mut_hum  #  --n-sim-seqs 100000
# # ./bin/compare-partition-methods.py --actions cache-simu-parameters $leaf_mut_hum &
# for isub in 0; do  # 1 2; do
#     # ./bin/compare-partition-methods.py --actions run-viterbi:vsearch-partition:naive-hamming-partition:partition --subset $isub --n-subsets 10 $leaf_mut_hum &  #run-viterbi:
#     # ./bin/compare-partition-methods.py --actions seed-partition --subset $isub --n-subsets 10 $leaf_mut_hum &  # --overwrite &
#     # sleep 60
#     # ./bin/compare-partition-methods.py --actions write-plots --subset $isub --n-subsets 10 --expected-methods vollmers-0.9:changeo:partition:naive-hamming-partition:vsearch-partition:synthetic $leaf_mut_hum --no-similarity-matrices &
#     ./bin/compare-partition-methods.py --actions write-plots --subset $isub --n-subsets 10 --expected-methods seed-partition $leaf_mut_hum & # --no-similarity-matrices &
#     # break
# done
# # ./bin/compare-partition-methods.py --actions compare-subsets --plot-mean-of-subsets --n-subsets 3 $leaf_mut_hum --expected-methods vollmers-0.9:changeo:mixcr:partition:naive-hamming-partition:vsearch-partition:misassign-0.60-singletons:misassign-distance-0.03 --no-similarity-matrices &
# # ./bin/compare-partition-methods.py --actions compare-subsets --plot-mean-of-subsets --n-subsets 3 $leaf_mut_hum --expected-methods vollmers-0.9:changeo:mixcr:partition:naive-hamming-partition:vsearch-partition --no-similarity-matrices &

# # ----------------------------------------------------------------------------------------
# # data
# hum="--humans 021-018"  # A:B:021-018:021-019"
# # ./bin/compare-partition-methods.py --actions cache-data-parameters $hum &
# istartstop=0:20000 #0:10000  # 
# # ./bin/compare-partition-methods.py --actions run-viterbi:vsearch-partition:naive-hamming-partition:partition --istartstop $istartstop $hum &  # --count-distances
# # ./bin/compare-partition-methods.py --actions run-changeo --istartstop $istartstop $hum &  # --count-distances
# ./bin/compare-partition-methods.py --old-output-structure --actions write-plots --istartstop $istartstop $hum --expected-methods vollmers-0.9:mixcr:changeo:partition:naive-hamming-partition:vsearch-partition #--no-similarity-matrices &  # --count-distances &

# ----------------------------------------------------------------------------------------
# different sample sizes
# ----------------------------------------------------------------------------------------
nl=7; xtra="--seed-cluster-bounds 10:15"
# nl=2.3; xtra="--zipf" # --count-distances --hfrac-bound-list 0,0"
leaf_mut_hum="--is-simu --n-leaf-list $nl --mutation-multipliers 1 --humans 021-018"
# ./bin/compare-partition-methods.py --actions simulate --n-sim-seqs 2000000 $leaf_mut_hum $xtra &
# ./bin/compare-partition-methods.py --actions cache-simu-parameters $leaf_mut_hum $xtra --n-simu-to-cache 200000 &

# istartstoplist="0:250 250:750 750:1500 1500:2500 2500:4000 4000:6500 6500:9500 9500:13500 13500:18500 18500:26000 26000:36000 36000:51000 51000:71000 71000:101000 101000:141000 141000:191000 191000:266000 266000:366000 366000:516000 516000:816000 816000:1316000 1316000:2066000"
istartstoplist="26000:36000 36000:51000 51000:71000 71000:101000 101000:141000 141000:191000 191000:266000 266000:366000 366000:516000 516000:816000 816000:1316000 1316000:2066000"
# extra: 141001:191001 141002:191002

# istartstoplist="51000:71000 71000:101000 101000:141000 141000:191000 191000:266000 266000:366000 366000:516000 516000:816000 816000:1316000 1316000:2066000"

# not_for_plotting="7:1000007 0:500000 0:1500 1500:4500 4500:8500 8500:13500 13500:21000 21000:31000 31000:81000 141000:216000 216000:316000 316000:466000 500000:800000 2:100000 1000007:1500007"
# new_seed="81000:141000 466000:666000 666000:966000 966000:1366000 1366000:1866000"
# e6="0:1000000 1:1000000 2:1000000"
# five_hundred_k="7:500007 500007:1000007"
# e5="0:100000 1:100000"
# reverse="800000:950000 950000:1050000 1050000:1125000 1125000:1175000 1175000:1215000 1215000:1245000 1245000:1265000 1265000:1280000 1280000:1290000 1290000:1297500 1297500:1302500 1302500:1306500 1306500:1309500 1309500:1312000 1312000:1313500 1313500:1314500 1314500:1315250 1315250:1315750 1315750:1316000"
# istartstoplist="$reverse $misc $e5 $e6 $five_hundred_k $new_seed"


# 966000:1366000 7:500007 1366000:1866000

# istartstoplist="81000:141000"

# iseed=2
# xtra="$xtra --iseed $iseed"  # --seed-cluster-bounds 3:5"
# istartstoplist="$iseed:1000000"

istartstopstr=`echo $istartstoplist | sed -e 's/:/,/g' -e 's/ /:/g'`

# for istartstop in $istartstoplist; do  # see code below to generate these
for istartstop in 141002:191002; do
    # ./bin/compare-partition-methods.py --actions naive-hamming-partition:partition:run-viterbi:vsearch-partition --istartstop $istartstop $leaf_mut_hum $xtra &
    # ./bin/compare-partition-methods.py --actions seed-partition --istartstop $istartstop $leaf_mut_hum $xtra &  #  --no-slurm
    ./bin/compare-partition-methods.py --actions cache-parameters --istartstop $istartstop $leaf_mut_hum $xtra
    ./bin/compare-partition-methods.py --actions partition --istartstop $istartstop $leaf_mut_hum $xtra
    # break
    # sleep 0.5s
    # ./bin/compare-partition-methods.py --actions write-plots --istartstop $istartstop $leaf_mut_hum --expected-methods run-viterbi:partition:naive-hamming-partition:vsearch-partition:synthetic --no-similarity-matrices $xtra &  # --count-distances &
    # ./bin/compare-partition-methods.py --actions write-plots --istartstop $istartstop $leaf_mut_hum --expected-methods seed-partition $xtra --no-similarity-matrices &  # --count-distances & --no-similarity-matrices 
    # break
    # ./bin/compare-partition-methods.py --actions annotate-seed-clusters --istartstop $istartstop $leaf_mut_hum $xtra  # --count-distances & --no-similarity-matrices 
done
# ./bin/compare-partition-methods.py --actions compare-subsets --istartstoplist $istartstopstr $leaf_mut_hum --expected-methods vollmers-0.9:mixcr:partition:naive-hamming-partition:vsearch-partition:misassign-0.60-singletons:misassign-distance-0.03 --no-similarity-matrices $xtra &
# ./bin/compare-partition-methods.py --actions compare-subsets --istartstoplist $istartstopstr $leaf_mut_hum --expected-methods seed-partition --no-similarity-matrices $xtra &

# istart=0
# # for s in 250 500 750 1000 1500 2500 3000 4000 5000 7500 10000 15000 20000 30000 40000 50000 75000 100000 150000 300000 500000; do
# for s in 1500 3000 4000 5000 7500 10000 50000 60000 75000 100000 150000 200000 300000 400000 500000; do
#     ((istop = istart + s))
#     echo -n " $istart:$istop"
#     ((istart += s))
# done
# exit
