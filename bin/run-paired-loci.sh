#!/bin/bash

bin=./test/cf-paired-loci.py

methods=synth-distance-0.03:synth-singletons-0.20:vjcdr3-0.9:scoper:vsearch-partition:partition
astr="--actions scoper" #:$methods"
# astr="--actions plot --plot-metrics $methods" # --perf-metrics f1"
# astr="--actions combine-plots --plot-metrics $methods"
common="--n-sub-procs 15 --n-max-procs 5 --single-light-locus igk --base-outdir /fh/local/dralph/partis/paired-loci $astr --dry"
# echo $bin --label vs-shm          --version v0 --n-replicates 5 --n-leaves-list 5:10 --n-sim-events-list 1000 --scratch-mute-freq-list 0.01:0.05:0.10:0.30:0.50:0.80 --final-plot-xvar scratch-mute-freq --x-legend-var mfreq --combo-extra-str scratch-x-var --pvks-to-plot 10 $common
# echo $bin --label vs-shm          --version v1 --n-replicates 3 --n-leaves-list 10 --n-sim-events-list 1000 --scratch-mute-freq-list 0.01:0.05:0.10:0.30 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs\" --final-plot-xvar scratch-mute-freq $common  # with these simu args, the scratch mute freq is almost identical to the final mean mfreq, so can just use the scratch mute freq on x axis
echo $bin --label vs-shm          --version v2 --n-replicates 3 --n-leaves-list 3 --n-sim-events-list 10000 --scratch-mute-freq-list 0.01:0.05:0.10:0.20:0.30 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs\" --final-plot-xvar scratch-mute-freq $common  # with these simu args, the scratch mute freq is almost identical to the final mean mfreq, so can just use the scratch mute freq on x axis
# echo $bin --label vs-n-leaves     --version v0 --n-replicates 3 --n-leaves-list 1:5:10:50 --n-sim-events-list 1000 --simu-extra-args="--constant-number-of-leaves" $common
# echo $bin --label vs-n-sim-events --version v0 --n-replicates 3 --n-leaves-list 1 --n-sim-events-list 100:1000:10000:50000 --simu-extra-args=\"--constant-number-of-leaves\" $common
echo $bin --label vs-n-sim-events --version v1 --n-replicates 3 --n-leaves-list 1 --n-sim-events-list 100:1000:10000:50000 --scratch-mute-freq-list 0.07:0.15 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs --constant-number-of-leaves\" --final-plot-xvar n-sim-events --pvks-to-plot 0.15 $common
echo $bin --label time-reqd --version v0 --n-replicates 3 --n-leaves-list 10 --n-sim-events-list 100:1000:10000:50000 --perf-metrics time-reqd --x-legend-var n-seqs $common --n-sub-procs 28 --n-max-procs 1  # NOTE duplicate parallelization args also NOTE made a 4th replicate after adding some time printing, but i think don't need it
echo $bin --label pairclean --version v0 --n-replicates 3 --n-leaves-list 3:10 --n-sim-events-list 3000 --scratch-mute-freq-list 0.15 --mean-cells-per-droplet-list None:1:2:5:10 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs\" --final-plot-xvar mean-cells-per-droplet --perf-metrics precision:sensitivity:f1:pcfrac-corr:pcfrac-mis:pcfrac-un $common  #  --fraction-of-reads-to-remove-list 0.05
