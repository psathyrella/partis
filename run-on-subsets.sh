#!/bin/bash

for isubset in 0 1 2 3 4 5; do
    ./runpart.py --cache-parameters --is-data --skip-unproductive --n-procs 10 --datadir data/imgt \
		 --seqfile /tmp/dralph/subsets/every-hundredth-subset-$isubset.tsv.bz2 \
		 --parameter-dir caches/subsets/$isubset --plotdir caches/subsets/$isubset
done
