#!/bin/bash

pd=/fh/fast/matsen_e/dralph/partis/chimera-testing

simfn=test/reference-results/test/simu.csv
chfn=simu-all-chimeras

label=typical-simu
# ./bin/partis cache-parameters --only-sm --infname $simfn --parameter-dir $pd/parameters/$label --sw-cachefname $pd/parameters/$label/sw-cache.csv --n-procs 8
# ./bin/chimera-plot.py $pd/parameters/$label/sw-cache.csv $pd/plots/$label

# ./bin/add-chimeras.py $simfn $chfn.fa
# ./bin/partis cache-parameters --only-sm --infname $chfn.fa --parameter-dir $pd/parameters/$chfn --sw-cachefname $pd/parameters/$chfn/sw-cache.csv --n-procs 8
# ./bin/chimera-plot.py $pd/parameters/$chfn/sw-cache.csv $pd/plots/$chfn

./bin/compare-plotdirs.py --outdir $pd/plots/comparisons/simu --names typical:chimeras --plotdirs $pd/plots/$label:$pd/plots/$chfn
./bin/compare-plotdirs.py --outdir $pd/plots/comparisons/simu-log --log y --translegend 0:0.15 --names typical:chimeras --plotdirs $pd/plots/$label:$pd/plots/$chfn

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
