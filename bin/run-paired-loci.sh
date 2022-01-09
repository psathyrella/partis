#!/bin/bash

bin=./test/cf-paired-loci.py

methods=synth-distance-0.03:synth-singletons-0.20:single-chain-partis:vjcdr3-0.9:mobille:scoper:vsearch-partition:partition #annotate:star-partition:partition:linearham # 
# astr="--actions $methods"
astr="--actions plot --plot-metrics $methods" # --perf-metrics f1"
# astr="--actions combine-plots --plot-metrics $methods"
common="--n-sub-procs 15 --n-max-procs 5 --single-light-locus igk --base-outdir /fh/fast/matsen_e/dralph/partis/paired-loci $astr" # --dry" # /fh/local/dralph
echo $bin --label vs-shm          --version v2 --n-replicates 3 --n-leaves-list 3 --n-sim-events-list 10000 --scratch-mute-freq-list 0.01:0.05:0.10:0.20:0.30 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs\" --final-plot-xvar scratch-mute-freq $common  # with these simu args, the scratch mute freq is almost identical to the final mean mfreq, so can just use the scratch mute freq on x axis
# # echo $bin --label vs-n-leaves     --version v0 --n-replicates 3 --n-leaves-list 1:5:10:50 --n-sim-events-list 1000 --simu-extra-args="--constant-number-of-leaves" $common
# echo $bin --label vs-n-sim-events --version v1 --n-replicates 3 --n-leaves-list 1 --n-sim-events-list 100:1000:10000:50000 --scratch-mute-freq-list 0.07:0.15 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs --constant-number-of-leaves\" --final-plot-xvar n-sim-events --pvks-to-plot 0.15 $common
# echo $bin --label time-reqd --version v0 --n-replicates 3 --n-leaves-list 10 --n-sim-events-list 100:1000:10000 --perf-metrics time-reqd --x-legend-var n-seqs $common --n-sub-procs 28 --n-max-procs 1  # NOTE duplicate parallelization args also NOTE there's also a 4th replicate, and 50000 samples for at least 3 replicates
# echo $bin --label pairclean --version v0 --n-replicates 3 --n-leaves-list 3:10 --constant-number-of-leaves-list 0:1 --n-sim-events-list 3000 --scratch-mute-freq-list 0.15 --mean-cells-per-droplet-list None:1:2:5:10 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs\" --final-plot-xvar mean-cells-per-droplet --perf-metrics precision:sensitivity:f1:pcfrac-corr:pcfrac-mis:pcfrac-un --use-val-cfgs $common  #  --fraction-of-reads-to-remove-list 0.05
# echo $bin --label test-antn --version imbal-v3   --n-replicates 3 --tree-imbalance-list None:0.04:0.07 --scratch-mute-freq-list 0.15 --n-leaves-list 50 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs\" --n-sim-events-list 50 --antn-perf --perf-metrics naive-hdist $common  # NOTE also made :0.13:0.14:0.16
exit 0

dvsn=v1 #test-ctnt
simvsn=v0

# bin=./datascripts/run.py
# common="--study 10x-examples --version $dvsn --paired-loci --n-procs 20 --no-slurm" # --only-csv-plots"
# infxtra='--extra-args=\"--make-per-gene-plots\"'
# # $bin cache-parameters $common --no-simu $infxtra
# # $bin partition $common --no-simu $infxtra
# # $bin simulate $common --logstr $simvsn # --extra-args="--paired-correlation-values v_gene.d_gene,0.5:v_gene.j_gene,0.5"
# # $bin cache-parameters $common $infxtra
# # $bin partition $common $infxtra
# exit 0

bin=./bin/compare-plotdirs.py
bidir=/fh/fast/matsen_e/processed-data/partis/10x-examples/$dvsn
bodir=$bidir/comparisons

for subj in hs-1-postvax hs-1-prevax hs-2-pbmc mm-balbc; do #hs-1-prevax; do #mm-balbc
    # for ltmp in igh igk; do #  igl
    # 	# ddir=$bidir/$subj/single-chain/plots/$ltmp/parameters/hmm
    # 	ddir=$bidir/$subj/igh+igk/plots/$ltmp
    # 	sdir=$bidir/$subj-simu-$simvsn/single-chain/plots/$ltmp/parameters/true
    # 	# for subd in overall mute-freqs; do
    # 	for subd in mute-freqs/per-gene-per-position/v mute-freqs/per-gene-per-position/d mute-freqs/per-gene-per-position/j; do
    # 	    echo $bin --outdir $bodir/$subj/$subd-$ltmp --plotdirs $ddir/$subd:$sdir/$subd --names data:simu --normalize # --extra-stats mean
    # 	done
    # 	# break
    # done
    for ltmp in igh igk igl; do
	# single-chain/plots/igh/partitions/sizes
	# plots/igh/partitions/sizes
	substr=plots/$ltmp/partitions/sizes
	ddir=$bidir/$subj #/single-chain/$substr
	sdir=$bidir/$subj-simu-$simvsn #/single-chain/$substr
	subd=cluster-sizes
	$bin --outdir $bodir/$subj/$subd-$ltmp --plotdirs $sdir/single-chain/$substr:$ddir/single-chain/$substr:$sdir/$substr:$ddir/$substr --names simu@single:data@single:simu@joint:data@joint \
	     --add-to-title $ltmp --log xy --colors '#006600:#2b65ec:#990012:black' --linewidths 5:2:5:2 &
    done
    # exit 0
done
