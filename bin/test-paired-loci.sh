#!/bin/bash

label=pair-params-v2  # smetrics-v0  # cf-ccfs-v0 #cells-per-drop-$cpd # pair-plots-v0 #fix-bad-merge-v2 #2 # pairclean-v2  # re-refactor-v0 #refactor-v4  #paired-clustering-output-v2
# NOTE use tmp.sh for refactor-v{3,4}

nprocs=10
common="--n-sim-events 50 --n-leaves 100 --n-procs $nprocs --no-per-base-mutation --allowed-cdr3-lengths 30:33:36:42:45:48 --mutation-multiplier 1" # --debug 1" #  --constant-number-of-leaves
in_param_dir=_output/paired-simulation/parameters
outdir=$fs/partis/tmp/$label #_output/$label
out_param_dir=$outdir/params

mkdir -p $outdir  # for the .log files
./bin/partis simulate --paired-loci --seed 1 --parameter-dir $in_param_dir --paired-outdir $outdir/simu $common --mean-cells-per-droplet 1.1 >$outdir/simu.log
for action in cache-parameters partition; do  # merge-paired-partitions plot-partitions
    common="--n-procs $nprocs --is-simu --no-mds-plots"
    ./bin/partis $action --paired-loci --paired-indir $outdir/simu --input-metafname $outdir/simu/meta.yaml --paired-outdir $outdir/inferred $common >$outdir/$action.log
    # --plot-partitions  # don't need to add this to either partition or merge-paired-partitions if we run plot-partitions afterward
done

# bd=_output/cells-per-drop
# subd=inferred/plots
# ./bin/compare-plotdirs.py --outdir ~/Dropbox/tmp-plots/cells-per-drop \
#      --normalize --translegend=-0.2:-0.2 \
#      --plotdirs $bd-1.0/$subd:$bd-1.2/$subd:$bd-1.7/$subd:$bd-2.0/$subd:$bd-3.0/$subd \
#      --names 1.0:1.2:1.7:2.0:3.0
