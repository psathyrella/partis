#!/bin/bash

bin=./test/cf-paired-loci.py

# methods=synth-distance-0.03:synth-singletons-0.20:vjcdr3-0.8:enclone:mobille:scoper:vsearch-partition:partition  # this is for vs-shm; for time-reqd: enclone:mobille:scoper:vsearch-partition:partition NOTE enclone needs fixing tho (for missing uids)
# methods=igblast:annotate:star-partition:partition:linearham # for test-antn imbal-v3
methods=partition:single-chain-partis; xstr="--combo-extra-str single-vs-joint-partis"
# methods=scoper:single-chain-scoper; xstr="--combo-extra-str single-vs-joint-scoper"  # NOTE this is only for vs-shm (comparing single vs joint); for time-reqd you only need scoper
astr="--actions $methods" #partition --merge-paired-partitions" #$methods"
# astr="--actions plot --plot-metrics $methods" # --perf-metrics f1"
# astr="--actions combine-plots --plot-metrics $methods $xstr"
common="--n-sub-procs 15 --n-max-procs 5 --single-light-locus igk --base-outdir /fh/fast/matsen_e/dralph/partis/paired-loci $astr" # --dry" # /fh/local/dralph
# echo $bin --label vs-shm          --version v3 --n-replicates 3 --n-leaves-list 3 --n-sim-events-list 10000 --scratch-mute-freq-list 0.01:0.05:0.10:0.20:0.30 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs --mutate-stop-codons\" --final-plot-xvar scratch-mute-freq $common  # with these simu args, the scratch mute freq is almost identical to the final mean mfreq, so can just use the scratch mute freq on x axis
# # echo $bin --label vs-n-leaves     --version v0 --n-replicates 3 --n-leaves-list 1:5:10:50 --n-sim-events-list 1000 --simu-extra-args="--constant-number-of-leaves" $common
# echo $bin --label vs-n-leaves     --version v1 --n-replicates 3 --n-leaves-list 1:2:3:5:7:10:25:50 --n-sim-events-list 50 --antn-perf --perf-metrics naive-hdist --simu-extra-args="--constant-number-of-leaves" $common
# # echo $bin --label vs-n-leaves     --version v2 --n-replicates 3 --n-leaves-list 1:2:3:5:7:25:50 --n-sim-events-list 50 --antn-perf --perf-metrics naive-hdist --scratch-mute-freq-list 0.15 --simu-extra-args="--constant-number-of-leaves --flat-mute-freq --same-mute-freq-for-all-seqs" $common
# echo $bin --label vs-n-sim-events --version v1 --n-replicates 3 --n-leaves-list 1 --n-sim-events-list 100:1000:10000:50000 --scratch-mute-freq-list 0.07:0.15 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs --constant-number-of-leaves\" --final-plot-xvar n-sim-events --pvks-to-plot 0.15 $common
# echo $bin --label time-reqd --version v0 --n-replicates 3 --n-leaves-list 10 --n-sim-events-list 100:1000:10000 --perf-metrics time-reqd --x-legend-var n-seqs $common --n-sub-procs 28 --n-max-procs 1  # NOTE duplicate parallelization args also NOTE there's also a 4th replicate, and 50000 samples for at least 3 replicates NOTE also will have to turn up fail % by hand for enclone
# echo $bin --label pairclean --version v0 --n-replicates 3 --n-leaves-list 3:10 --constant-number-of-leaves-list 0:1 --n-sim-events-list 3000 --scratch-mute-freq-list 0.15 --mean-cells-per-droplet-list None:1:2:5:10 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs\" --final-plot-xvar mean-cells-per-droplet --perf-metrics all-pcfrac:f1:precision:sensitivity --make-hist-plots --use-val-cfgs $common  #  --fraction-of-reads-to-remove-list 0.05
# echo $bin --label pairclean --version v1 --n-replicates 3 --n-leaves-list 1:2:3:10:hist --constant-number-of-leaves-list 1:1:1:1:0 --zip-vars n-leaves:constant-number-of-leaves --n-sim-events-list 3000 --scratch-mute-freq-list 0.15 --mean-cells-per-droplet-list None:1:2:5:10 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs\" --final-plot-xvar mean-cells-per-droplet --perf-metrics all-pcfrac:f1:precision:sensitivity --make-hist-plots --use-val-cfgs $common  #  --fraction-of-reads-to-remove-list 0.05
# echo $bin --label pairclean --version v2 --n-replicates 3 --n-leaves-list 1:2:3:10:hist --constant-number-of-leaves-list 1:1:1:1:0 --zip-vars n-leaves:constant-number-of-leaves --n-sim-events-list 3000 --scratch-mute-freq-list 0.07 --mean-cells-per-droplet-list 1:2:5:10 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs --constant-cells-per-droplet\" --final-plot-xvar mean-cells-per-droplet --perf-metrics all-pcfrac:f1:precision:sensitivity --make-hist-plots --use-val-cfgs $common  #  --fraction-of-reads-to-remove-list 0.05
# echo $bin --label pairfix --version v1 --n-replicates 3 --n-leaves-list hist --n-sim-events-list 3000 --scratch-mute-freq-list 0.07 --bulk-data-fraction-list 0:0.5:0.8:0.9:0.95 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs\" --inference-extra-args=\"--pair-unpaired-seqs-with-paired-family\" --final-plot-xvar bulk-data-fraction --perf-metrics all-pcfrac:f1:precision:sensitivity --make-hist-plots --use-val-cfgs --empty-bin-range 0:200 $common
# echo $bin --label test-antn --version imbal-v3   --n-replicates 2 --tree-imbalance-list None:0.04:0.07 --scratch-mute-freq-list 0.15 --n-leaves-list 50 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs\" --n-sim-events-list 50 --antn-perf --perf-metrics naive-hdist $common  # NOTE also made :0.13:0.14:0.16
# NOTE have to set --n-sub-procs to 1 for partition step, and re-set --n-sim-events-list for each --n-leaves value (500 leaves: 10 events, 100:50, 50:100):
# echo $bin --label key-trans --version v0   --n-replicates 2  --biggest-naive-seq-cluster-to-calculate 5:15:10000 --biggest-logprob-cluster-to-calculate 5:15:10000 --zip-vars biggest-naive-seq-cluster-to-calculate:biggest-logprob-cluster-to-calculate --scratch-mute-freq-list 0.15 --n-leaves-list 50:100:500 --simu-extra-args=\"--flat-mute-freq --same-mute-freq-for-all-seqs --constant-number-of-leaves --only-genes IGHV1-2*01:IGHD2-15*01:IGHJ6*02:IGKV1-12*01:IGKJ1*01:IGKDx-x*x --force-dont-generate-germline-set --allowed-cdr3-lengths 66:33\" --n-sim-events-list 10 --inference-extra-args=\"--debug 1 --sw-debug 0\" --perf-metrics time-reqd --final-plot-xvar biggest-logprob-cluster-to-calculate --bcrham-time $common
echo $bin --label lcdr3 --version v0 --n-replicates 2 --n-sim-events-list 1000 --single-light-locus igk --base-outdir /fh/fast/matsen_e/dralph/partis/paired-loci --allowed-cdr3-lengths ighM9-22:ighM15-34:ighM36-43:ighM45-52:ighM54-61:ighM63-73 $common
# echo $bin --label lcdr3 --version v1 --n-replicates 2 --n-sim-events-list 500 --single-light-locus igk --base-outdir /fh/fast/matsen_e/dralph/partis/paired-loci --allowed-cdr3-lengths ighM12-13:ighM15-16:ighM18-19:ighM21-22:ighM24-25:ighM27-28:ighM30-31:ighM33-34:ighM36-37 $common
exit 0

# ----------------------------------------------------------------------------------------
# data
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
fsddir=/fh/fast/matsen_e/processed-data/partis
bidir=$fsddir/10x-examples/$dvsn
bodir=$bidir/comparisons

tsubj=hs-1-postvax
tpd=$bidir/comparisons/$tsubj
$bin --outdir $bodir/ins-del-lens-all/$tsubj --plotdirs $bidir/$tsubj/igh+igk/plots/igh/overall:$bidir/$tsubj/igh+igk/plots/igk/overall:$bidir/$tsubj/igh+igl/plots/igl/overall --names igh:igk:igl --normalize & # --extra-stats mean
$bin --outdir $bodir/ins-del-lens-igh/$tsubj --plotdirs $bidir/$tsubj/igh+igk/plots/igh/overall --names igh --normalize --single-plotdir & # --extra-stats mean
exit 0

for subj in hs-1-postvax hs-1-prevax hs-2-pbmc mm-balbc; do #hs-1-prevax; do #mm-balbc
    for ltmp in igh igk; do #  igl
    	# ddir=$bidir/$subj/single-chain/plots/$ltmp/parameters/hmm
    	ddir=$bidir/$subj/igh+igk/plots/$ltmp
    	sdir=$bidir/$subj-simu-$simvsn/single-chain/plots/$ltmp/parameters/true
    	for subd in overall mute-freqs; do
    	# for subd in mute-freqs/per-gene-per-position/v mute-freqs/per-gene-per-position/d mute-freqs/per-gene-per-position/j; do
	    $bin --outdir $bodir/$subj/$subd-$ltmp --plotdirs $ddir/$subd:$sdir/$subd --names data:simu --normalize & # --extra-stats mean
    	done
    	# break
    done
    # for ltmp in igh igk igl; do
    # 	# single-chain/plots/igh/partitions/sizes
    # 	# plots/igh/partitions/sizes
    # 	substr=plots/$ltmp/partitions/sizes
    # 	ddir=$bidir/$subj #/single-chain/$substr
    # 	sdir=$bidir/$subj-simu-$simvsn #/single-chain/$substr
    # 	subd=cluster-sizes
    # 	echo $bin --outdir $bodir/$subj/$subd-$ltmp --plotdirs $sdir/single-chain/$substr:$ddir/single-chain/$substr:$sdir/$substr:$ddir/$substr --names simu@single:data@single:simu@joint:data@joint \
    # 	     --add-to-title $ltmp --log xy --colors '#006600:#2b65ec:#990012:black' --linewidths 5:2:5:2 &
    # done
done
exit 0

for subjdir in $bidir/hs-1-postvax $fsddir/goo-dengue-10x/v12/d-14; do
    pstr=paired-seqs-per-seq
    cfd=$bodir/`basename $subjdir`/$pstr
    # for tstr in before after; do
    # 	mkdir -p $cfd/$tstr
    # 	cp -v $subjdir/plots/$pstr-$tstr.csv $cfd/$tstr/$pstr.csv
    # done
    $bin --outdir $cfd/cf --plotdirs $cfd/after:$cfd/before --names after:before --normalize --colors='#006600:#990012' --xbounds=-0.5:9.5 --ybounds=0:0.9 --xticks 0:1:2:3:4:5:6:7:8:9 --no-errors --square-bins --ytitle="fraction of seqs" &  #  --add-to-title `basename $subjdir`
done

