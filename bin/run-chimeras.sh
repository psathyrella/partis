#!/bin/bash

# ----------------------------------------------------------------------------------------
# dual barcode data from chaim/sai(?)
label=chimera-dual-barcode
datadir=/fh/fast/matsen_e/data/$label

outdir=$fs/partis/$label/fixed-d-synthesis
for fn in good chimeras; do
    mkdir -p $outdir/$fn
    # ./bin/partis cache-parameters --only-sm --infname $datadir/$fn.fa --parameter-dir $outdir/$fn/parameters --plotdir $outdir/$fn/plots --sw-cachefname $outdir/$fn/sw-cache.yaml --n-procs 5 --only-overall-plots --only-csv-plots >$outdir/$fn/log
    ./bin/chimera-plot.py $outdir/$fn/sw-cache.yaml $outdir/$fn/chimera-plots --title $fn --chunk-len 75

    # ./bin/partis simulate --outfname $outdir/$fn/simu.yaml --parameter-dir $outdir/$fn/parameters --parameter-type sw --n-sim-events 1000 --n-leaves 1 --constant-number-of-leaves
    # ./bin/partis cache-parameters --is-simu --only-sm --infname $outdir/$fn/simu.yaml --parameter-dir $outdir/$fn/simu-parameters --plotdir $outdir/$fn/simu-plots --sw-cachefname $outdir/$fn/simu-parameters/sw-cache.yaml --n-procs 5 --only-overall-plots --only-csv-plots >$outdir/$fn/simu-log
    # ./bin/add-chimeras.py $outdir/$fn/simu.yaml $outdir/$fn/simu-all-chimeras.fa --min-chunk-len 1
    # ./bin/partis cache-parameters --only-sm --infname $outdir/$fn/simu-all-chimeras.fa --parameter-dir $outdir/$fn/simu-all-chimeras-parameters --plotdir $outdir/$fn/simu-all-chimeras-plots --sw-cachefname $outdir/$fn/simu-all-chimeras-parameters/sw-cache.yaml --n-procs 5 --only-overall-plots --only-csv-plots >$outdir/$fn/simu-all-chimeras-log
    # ./bin/chimera-plot.py $outdir/$fn/simu-parameters/sw-cache.yaml $outdir/$fn/chimera-plots-simu
    # ./bin/chimera-plot.py $outdir/$fn/simu-all-chimeras-parameters/sw-cache.yaml $outdir/$fn/chimera-plots-simu-all-chimeras
    # ./bin/compare-plotdirs.py --outdir $outdir/comparison-plots/simu-all-chimeras-vs-none/$fn --log y --translegend 0:-0.2 --names no@chimeras:all@chimeras --plotdirs $outdir/$fn/chimera-plots-simu:$outdir/$fn/chimera-plots-simu-all-chimeras
done
./bin/compare-plotdirs.py --outdir $outdir/comparison-plots-good-vs-chimeras --log y --translegend 0:-0.2 --names none:all@chimeras --plotdirs $outdir/good/chimera-plots:$outdir/chimeras/chimera-plots
# ./bin/compare-plotdirs.py --outdir $outdir/comparison-plots-no-indels --log y --translegend 0:-0.2 --names good:chimeras --plotdirs $outdir/good/chimera-plots:$outdir/chimeras/chimera-plots
# ./bin/compare-plotdirs.py --outdir $outdir/comparison-plots XXX --log y --translegend 0:-0.2 --names good:chimeras --plotdirs $outdir/good/chimera-plots:$outdir/chimeras/chimera-plots

# subd=sw/mute-freqs/overall
# ./bin/compare-plotdirs.py --outdir $outdir/comparison-plots/data-vs-simu/$subd --names good:good-simu:chimeras:chimeras-simu --plotdirs $outdir/good/plots/$subd:$outdir/good/simu-plots/$subd:$outdir/chimeras/plots/$subd:$outdir/chimeras/simu-plots/$subd

exit 0

# ----------------------------------------------------------------------------------------
# initial simulation study
pd=/fh/fast/matsen_e/dralph/partis/chimera-testing

simfn=test/reference-results/test/simu.csv
chfn=simu-all-chimeras

label=typical-simu
# ./bin/partis cache-parameters --only-sm --infname $simfn --parameter-dir $pd/parameters/$label --sw-cachefname $pd/parameters/$label/sw-cache.csv --n-procs 8
# ./bin/chimera-plot.py $pd/parameters/$label/sw-cache.csv $pd/plots/$label

# ./bin/add-chimeras.py $simfn $chfn.fa
# ./bin/partis cache-parameters --only-sm --infname $chfn.fa --parameter-dir $pd/parameters/$chfn --sw-cachefname $pd/parameters/$chfn/sw-cache.csv --n-procs 8
# ./bin/chimera-plot.py $pd/parameters/$chfn/sw-cache.csv $pd/plots/$chfn

# ./bin/compare-plotdirs.py --outdir $pd/plots/comparisons/simu --names typical:chimeras --plotdirs $pd/plots/$label:$pd/plots/$chfn
# ./bin/compare-plotdirs.py --outdir $pd/plots/comparisons/simu-log --log y --translegend 0:0.15 --names typical:chimeras --plotdirs $pd/plots/$label:$pd/plots/$chfn

# data
study=jason-mg
samples="MK02-igh MK03-igh MK08-igh AR03-igh AR04-igh HD10-igh HD13-igh"
# for sample in $samples; do
#     echo $sample
#     ./bin/chimera-plot.py /fh/fast/matsen_e/processed-data/partis/$study/gls-gen-paper-v10/$sample/sw-cache.csv $pd/plots/$study/$sample
# done
./bin/compare-plotdirs.py --outdir $pd/plots/comparisons/$study-log --log y --translegend 0:-0.2 --names `echo $samples|sed 's/ /:/g'` --plotdirs $pd/plots/$study

study=jason-influenza
samples="FV-igh IB-igh GMC-igh"
# for sample in $samples; do
#     echo $sample
#     ./bin/chimera-plot.py /fh/fast/matsen_e/processed-data/partis/$study/gls-gen-paper-v10/$sample/sw-cache.csv $pd/plots/$study/$sample
# done
./bin/compare-plotdirs.py --outdir $pd/plots/comparisons/$study-log --log y --translegend 0:0 --names `echo $samples|sed 's/ /:/g'` --plotdirs $pd/plots/$study
