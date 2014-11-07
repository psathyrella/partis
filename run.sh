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

for iproc in {1..8}; do
    ./runpart.py --simulate --parameter_dir caches/test-parameters/hmm_parameters --n_max_queries 1250 --outfname test/recombinator/$iproc/simu.csv &
    limit_procs runpart
done
