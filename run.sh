#!/bin/bash

# for n_min in 5 10 40; do
#     ./runpart.py --cache_parameters \
# 	--seqfile caches/recombinator/longer-reads/simu.csv \
# 	--parameter_dir caches/$n_min-min-occur \
# 	--n_procs 10 \
# 	--min_observations_to_write $n_min \
# 	--n_max_queries 10000 \
# 	--plotdir /tmp/dralph/partis/$n_min &
# done

cmd=./runpart.py
tests='few many'
for tst in $tests; do
    common_args="--n_procs 10 --datadir data/$tst"
    $cmd --cache_parameters --seqfile test/every-hundredth-data.tsv.bz2 --is_data --skip_unproductive --parameter_dir caches/$tst --plotdir $www/partis/$tst --n_max_queries 1000 $common_args
    simu_file=caches/recombinator/$tst/simu.csv 
    $cmd --simulate --parameter_dir caches/$tst/hmm_parameters --n_max_queries 200 --outfname $simu_file $common_args
    $cmd --cache_parameters --seqfile $simu_file --parameter_dir caches/simu-$tst --plotdir $www/partis/simu-$tst $common_args
    $cmd --point_estimate --seqfile $simu_file --parameter_dir caches/simu-$tst/hmm_parameters --plot_performance --plotdir $www/partis/$tst $common_args
done
