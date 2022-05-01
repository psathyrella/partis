#!/bin/bash

bin=./test/cf-tree-metrics.py

label=v0
testlabel=test-old
testargs="--n-sim-seqs-per-gen-list 50:125 --lb-tau-list 0.002:0.003 --obs-times 100 --carry-cap 1000 --n-generations-list 4:5"

# # $bin --actions get-lb-bounds --label $label  #--make-plots
# # $bin --actions get-lb-bounds --label $testlabel $testargs --make-plots
# # echo $bin --actions get-lb-bounds --seq-len 133 --label aa-lb-bounds-v0 --make-plots
# # echo $bin --actions get-lb-bounds --seq-len 133 --label $testlabel $testargs --make-plots
# for slen in 1200; do #300 400 500 600 700 900; do  # redo normalization for different seq lengths
#     # 1200 doesn't work, some %.xf print things i think need to be changed
#     $bin --actions get-lb-bounds --seq-len $slen --label generalize-lb-bounds-$slen --lb-tau-list auto --make-plots >$slen.log #--make-plots
# done
# exit 0

# echo $bin --label $testlabel $testargs --only-csv-plots
# echo $bin --label $label --n-replicates 3 --only-csv-plots

dtv=3; nest=100; depth=10  # dtv=3; nest=30; depth=10  # dtv=2; nest=100; depth=5
# mm=dtr; dtr_args="--dtr-path /fh/fast/matsen_e/dralph/partis/tree-metrics/dtr-train-v$dtv/seed-0/dtr/train_n-estimators_${nest}_max-depth_${depth}-dtr-models --extra-plotstr v$dtv-$nest-$depth"
# dtr_args="--actions plot --plot-metrics dtr --plot-metric-extra-strs v3-100-10"  # :dtr:dtr  # :v2-100-5:v3-30-10

# for mm in lbr aa-lbr; do #aa-lbi aa-lbr shm delta-lbi cons-dist-aa cons-dist-nuc; do
common="--only-csv-plots --n-max-procs 25 $dtr_args"
# common="--actions get-tree-metrics --only-csv-plots --n-max-procs 25 --metric-method $mm $dtr_args" # # --no-tree-plots --slurm
# common="--actions plot --plot-metrics lbr:aa-lbr --distr-hist-limit N:3"  # cons-dist-nuc:cons-dist-aa"
# common="--actions plot --plot-metrics shm:delta-lbi:lbi:lbr:cons-dist-aa:cons-dist-nuc:aa-lbi:aa-lbr" # --plot-metric-extra-strs ::::"
# common="--actions combine-plots --plot-metrics shm:delta-lbi:lbi:lbr:cons-dist-aa:cons-dist-nuc:aa-lbi:aa-lbr:dtr --combo-extra-str with-dtr-and-aa-lb --plot-metric-extra-strs ::::::::v3-100-10 --dont-plot-extra-str"
# common="--actions combine-plots --plot-metrics shm:delta-lbi:lbi:lbr:cons-dist-aa:cons-dist-nuc:aa-lbi:aa-lbr --combo-extra-str with-aa-lb" #-pvk-$pvk --pvks-to-plot $pvk"  # :dtr --plot-metric-extra-strs :::::v3-100-10 --dont-plot-extra-str
common="$common --base-outdir /fh/local/dralph/partis/tree-metrics"
# ----------------------------------------------------------------------------------------
# $bin --label vary-metric-v2 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 500 --obs-times-list 50:100:500 --n-sim-seqs-per-gen-list 100 --metric-for-target-distance-list aa:aa-sim-blosum --lb-tau-list 0.0025 --final-plot-xvar obs-times --include-relative-affy-plots $common # --pvks-to-plot aa-sim-blosum
# echo $bin --label tau-vs-obs-frac-v1 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 30:50:100:200 --legend-var obs_frac $common  # rerun of "v2" (at top), but redoing things since plots look a bit weird on v2 when I remake them
echo $bin --label vary-selection-strength-v1 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 100 --selection-strength-list 0.1:0.4:0.7:0.8:0.9:1.0 --final-plot-xvar selection-strength $common  # all the other tau values are also there, but I decided I wanted selection strength on the x axis NOTE the other tau values will have lbr-tau-factor 1 and non-normalized lbi, so you *really* need to not mix them with 0.0025 ##  --lb-tau-list 0.0025
# $bin --label carry-cap-vs-n-obs-v1 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 500:1000:3000 --obs-times-list 150 --n-sim-seqs-per-gen-list 30:75:150:500 --lb-tau-list 0.0025 --final-plot-xvar carry-cap $common --pvks-to-plot 30  # full carry cap: 250:500:1000:3000 and n/gen: 15:30:75:150:500 lists (not plotting them all)
# $bin --label vary-obs-times-v2 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 350:1000:2000 --obs-times-list 50:100:250:500:1000:3000 --n-sim-seqs-per-gen-list 100 --lb-tau-list 0.0025 --final-plot-xvar obs-times $common --pvks-to-plot 350  # full carry-cap-list 250:350:500:1000:2000
# $bin --label vary-obs-times-v3 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 350:1000:2000 --obs-times-list 50,100,150,200,250:100,200,300,400,500:200,400,600,800,1000:600,1200,1800,2400,3000 --n-sim-seqs-per-gen-list 20 --lb-tau-list 0.0025 --final-plot-xvar obs-times --include-relative-affy-plots $common --pvks-to-plot 1000  # full carry-cap-list 250:350:500:1000:2000
# $bin --label vary-obs-times-v4 --n-replicates 30 --n-sim-events-per-proc 10 --min-target-distance-list 2 --carry-cap-list 350:1000:2000 --obs-times-list 50,100,150,200,250:100,200,300,400,500:200,400,600,800,1000:600,1200,1800,2400,3000 --n-sim-seqs-per-gen-list 20 --lb-tau-list 0.0025 --final-plot-xvar obs-times --include-relative-affy-plots $common --pvks-to-plot 1000
# $bin --label vary-obs-times-v5 --n-replicates 30 --n-sim-events-per-proc 10 --min-target-distance-list 2 --carry-cap-list 350:1000:2000 --obs-times-list 50:100:250:500:1000:3000 --n-sim-seqs-per-gen-list 100 --lb-tau-list 0.0025 --final-plot-xvar obs-times $common --pvks-to-plot 350
# $bin --label vary-sampling-scheme-v1 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 30:50:100:200 --leaf-sampling-scheme-list uniform-random:affinity-biased:high-affinity --lb-tau-list 0.0025 --final-plot-xvar n-sim-seqs-per-gen $common --pvks-to-plot high-affinity
# $bin --label vary-n-targets-v0 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 50:100:250:500 --n-sim-seqs-per-gen-list 100 --target-count-list 1:2:4 --lb-tau-list 0.0025 --final-plot-xvar obs-times $common # --pvks-to-plot 4
# $bin --label vary-n-targets-v1 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 50:100:250:500 --n-sim-seqs-per-gen-list 100 --target-count-list 4:8:16 --n-target-clusters-list 1:2:4 --zip-vars target-count:n-target-clusters --lb-tau-list 0.0025 --final-plot-xvar obs-times $common --pvks-to-plot="4; 16"
# $bin --label vary-context-dependence-v1 --n-replicates 5 --n-sim-events-per-proc 50 --carry-cap-list 350 --obs-times-list 100:350 --n-sim-seqs-per-gen-list 100 --context-depend-list 0:1 --lb-tau-list 0.0025 $common --final-plot-xvar obs-times  --n-sub-procs 10  # --pvks-to-plot 0 (and 1), do 'em separately
# $bin --label among-family-variation-v0 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 250 --n-sim-seqs-per-gen-list 0 --lb-tau-list 0.0025 --parameter-variances-list n-sim-seqs-per-generation,15..15:n-sim-seqs-per-generation,50..50:n-sim-seqs-per-generation,250..250:n-sim-seqs-per-generation,25..50..75:n-sim-seqs-per-generation,15..50..150..500 --final-plot-xvar parameter-variances $common
# $bin --label among-family-variation-v1 --n-replicates 15 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 250 --n-sim-seqs-per-gen-list 150 --lb-tau-list 0.0025 --parameter-variances-list selection-strength,0.25..0.25:selection-strength,0.75..0.75:selection-strength,1.0..1.0:selection-strength,0.5..0.75..1.0:selection-strength,0.25..0.67..1.0:selection-strength,0.1..0.25..0.75..1.0 --final-plot-xvar parameter-variances $common
# $bin --label among-family-variation-v2 --n-replicates 15 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list=-1 --n-sim-seqs-per-gen-list 150 --lb-tau-list 0.0025 --parameter-variances-list obs-times,50..50:obs-times,150..150:obs-times,500..500:obs-times,50..100..200:obs-times,50..150..500..1000 --final-plot-xvar parameter-variances $common
# $bin --label among-family-variation-v4 --n-replicates 15 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list=-1 --n-sim-seqs-per-gen-list=-1 --lb-tau-list 0.0025 --parameter-variances-list n-sim-seqs-per-generation,150..150_c_selection-strength,0.75..0.75_c_obs-times,250..250:n-sim-seqs-per-generation,50..150..300_c_selection-strength,0.5..0.75..1.0_c_obs-times,100..250..500:n-sim-seqs-per-generation,25..150..500_c_selection-strength,0.25..0.75..1.0_c_obs-times,50..250..1000 --final-plot-xvar parameter-variances $common
# echo $bin --label true-vs-inferred-v0 --n-replicates 2 --n-sim-events-per-proc 30 --carry-cap-list 500 --obs-times-list 150:1500 --n-sim-seqs-per-gen-list 100 --lb-tau-list 0.0025 $common  # NOTE to get the true vs inferred  plots you have to turn off --only-csv-plots for get-tree-metrics (and probably turn on --no-tree-plots) and uncomment two fcn calls in treeutils.plot_tree_metrics(): lbplotting.plot_true_vs_inferred_lb() and lbplotting.plot_cons_seq_accuracy()
# echo $bin --label true-vs-inferred-v1 --n-replicates 2 --n-sim-events-per-proc 50 --carry-cap-list 500:2000 --obs-times-list 150:1500 --n-sim-seqs-per-gen-list 100 --lb-tau-list 0.0025 $common
# echo $bin --label true-vs-inferred-v2 --n-replicates 2 --n-sim-events-per-proc 50 --carry-cap-list 500:2000 --obs-times-list 150:1500 --n-sim-seqs-per-gen-list 100 --min-target-distance-list 0:2:5 --lb-tau-list 0.0025 $common
# $bin --label fix-lbr-tau-sample-time-v0 --n-replicates 10 --n-sim-events-per-proc 50 --carry-cap-list 1000 --obs-times-list 100:350:50,75,100,150:100,150,300,350 --n-sim-seqs-per-gen-list 100:100:25:25 --zip-vars obs-times:n-sim-seqs-per-gen $common
# $bin --label fix-lbr-tau-fraction-sampled-v0 --n-replicates 10 --n-sim-events-per-proc 50 --carry-cap-list 1000 --obs-times-list 50,100,150 --n-sim-seqs-per-gen-list 10:25:50:100:300 $common
# echo $bin --label fix-lbr-tau-obs-time-v0 --n-replicates 5 --n-sim-events-per-proc 30 --carry-cap-list 1000 --obs-times-list 100:250:500 --n-sim-seqs-per-gen-list 100:200 $common
# echo $bin --label fix-lbr-tau-fraction-sampled-v1 --n-replicates 10 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 100,150,200,250,300 --n-sim-seqs-per-gen-list 10:50 $common
 # --final-plot-xvar n-sim-seqs-per-gen --pvks-to-plot 0.002

# some of these may have had different n max queries for different steps (but n max queries was only set for the steps where its value isn't very important, i.e. it was never set for training)
# common="--actions bcr-phylo --bcr-phylo-actions cache-parameters:partition --only-csv-plots" # --base-outdir /fh/local/dralph/partis/tree-metrics" # --sub-slurm"  #  /loc/scratch/dralph/partis/tree-metrics
common="$common --n-max-queries 150000"  # $common  --iseed 1"  for the dtr training samples we usually only want to test on 150k or so queries (and also sometimes need to restrict to certain seeds to avoid training samples)
# echo $bin --label dtr-train-v0 --n-replicates 5 --n-sim-events-per-proc 1000 --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 150 --selection-strength 0.75 --lb-tau-list 0.0025 --no-scan-parameter-variances carry-cap,2000:obs-times,150:n-sim-seqs-per-generation,200:selection-strength,0.5 $common &
# echo $bin --label dtr-train-v1 --n-replicates 5 --n-sub-procs 30 --n-sim-events-per-proc 50000 --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 30 --selection-strength 0.75 --lb-tau-list 0.0025 --no-scan-parameter-variances carry-cap,2000:obs-times,150:n-sim-seqs-per-generation,15:selection-strength,0.5 $common &
# echo $bin --label dtr-train-v2 --n-replicates 2 --n-sub-procs 15 --n-sim-events-per-proc 300000 --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 20 --selection-strength 0.75 --lb-tau-list 0.0025 --no-scan-parameter-variances carry-cap,2000:obs-times,150:n-sim-seqs-per-generation,15:selection-strength,0.5 $common &
# echo $bin --label dtr-train-v3 --n-replicates 2 --n-sub-procs 25 --n-sim-events-per-proc 50000 --carry-cap-list=-1 --obs-times-list=-1 --n-sim-seqs-per-gen-list=-1 --selection-strength=-1. --lb-tau-list 0.0025 --no-scan-parameter-variances carry-cap,250..500..900..1000..1100..1500..5000:obs-times,75..100..150..200..1000:n-sim-seqs-per-generation,15..30..75..150..500:selection-strength,0.5..0.9..0.95..1.0 $common  # NOTE made a second replicate (iseed 1) with only 1000 events, just for testing
# done


# ----------------------------------------------------------------------------------------
# old ones:
# echo $bin --label vary-carry-cap-v0 --n-replicates 10 --n-sim-events-per-proc 10 --carry-cap-list 500:750:1000:2000:5000 --obs-times-list 100,200 --n-sim-seqs-per-gen-list 75 --include-relative-affy-plots $common
# echo $bin --label vary-obs-times-v0 --n-replicates 10 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 100:200:300:100,150:200,250:100,200,300 --n-sim-seqs-per-gen-list 100:100:100:50:50:33 --zip-vars obs-times:n-sim-seqs-per-gen --include-relative-affy-plots $common
# echo $bin --label vary-obs-times-v1 --n-replicates 10 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 300:100,200,300:200,250,300 --n-sim-seqs-per-gen-list 100:33:33 --zip-vars obs-times:n-sim-seqs-per-gen --include-relative-affy-plots $common
# echo $bin --label vary-obs-frac-v0 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 30:50:75:100:150:200 $common
# echo $bin --label vary-metric-v0 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 100 --metric-for-target-distance-list aa:aa-sim-blosum --include-relative-affy-plots $common
# echo $bin --label vary-metric-v1 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 500 --n-sim-seqs-per-gen-list 100 --metric-for-target-distance-list aa:aa-sim-blosum --include-relative-affy-plots $common --pvks-to-plot aa-sim-blosum  # seeing if different parameters will change the fact that lbi does better than cons-dist-aa (as in vary-metric-v0)
# echo $bin --label vary-selection-strength-v0 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 100 --selection-strength-list 0.1:0.4:0.7:0.8:0.9:1.0 $common
# echo $bin --label carry-cap-vs-n-obs-v0 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 260:260:500:500:700:700:1500:1500:3000:3000 --obs-times-list 150 --n-sim-seqs-per-gen-list 13:26:25:50:35:70:75:150:150:300 --lb-tau-list 0.0025 --zip-vars carry-cap:n-sim-seqs-per-gen --final-plot-xvar carry-cap --legend-var obs_frac $common
# echo $bin --label carry-cap-vs-n-obs-only-leaves-v0 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 250:500:1000:3000 --obs-times-list 150 --n-sim-seqs-per-gen-list 15:75:500 --lb-tau-list 0.0025 --dont-observe-common-ancestors --final-plot-xvar carry-cap $common
# echo $bin --label choose-among-families-v1 --n-replicates 30 --n-sim-events-per-proc 30  --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 150 --selection-strength 0.75 --lb-tau-list 0.0025 --dont-observe-common-ancestors --no-scan-parameter-variances carry-cap,2000:obs-times,150:n-sim-seqs-per-generation,200:selection-strength,0.5 $common
# echo $bin --label choose-among-families-v2 --n-replicates 10 --n-sim-events-per-proc 30  --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 150 --selection-strength 0.75 --lb-tau-list 0.0025 --dont-observe-common-ancestors --no-scan-parameter-variances selection-strength,0.5 $common
# echo $bin --label choose-among-families-v3 --n-replicates 10 --n-sim-events-per-proc 30  --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 150 --lb-tau-list 0.0025 --dont-observe-common-ancestors $common
# echo $bin --label choose-among-families-v4 --n-replicates 10 --n-sim-events-per-proc 150 --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 150 --lb-tau-list 0.0025 --dont-observe-common-ancestors $common
# echo $bin --label choose-among-families-v5 --n-replicates 10 --n-sim-events-per-proc 150 --carry-cap-list 1500 --obs-times-list 150 --n-sim-seqs-per-gen-list 150 --lb-tau-list 0.0025 $common
# echo $bin --label vary-sampling-scheme-v0 --n-replicates 30 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list 150 --n-sim-seqs-per-gen-list 100 --leaf-sampling-scheme-list uniform-random:affinity-biased:high-affinity --lb-tau-list 0.0025 --include-relative-affy-plots $common
# echo $bin --label vary-context-dependence-v0 --n-replicates 5 --n-sim-events-per-proc 10 --carry-cap-list 350 --obs-times-list 100 --n-sim-seqs-per-gen-list 30 --context-depend-list 0:1 --lb-tau-list 0.0025 $common --n-sub-procs 10
# $bin --label among-family-variation-v3 --n-replicates 15 --n-sim-events-per-proc 10 --carry-cap-list 1000 --obs-times-list=-1 --n-sim-seqs-per-gen-list 150 --lb-tau-list 0.0025 --parameter-variances-list obs-times,50..50:obs-times,150..150:obs-times,500..500:obs-times,100..150..250:obs-times,50..150..250..500 --final-plot-xvar parameter-variances $common
