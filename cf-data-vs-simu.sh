#!/bin/bash

parameters='v_3p_del d_5p_del d_3p_del j_5p_del vd_insertion dj_insertion vd_insertion_content dj_insertion_content'
mute_parameters='all-mean-freq v-mean-freq d-mean-freq j-mean-freq'

label=mimic
basedir=$www/partis/$label/params  #$_output/$label/plots/params
simudir=$basedir/simu/hmm_parameters
datadir=$basedir/data/hmm_parameters

outdir=/var/www/sharing/dralph/partis/$label/data-vs-simu

for param in . $parameters $mute_parameters; do
    if [ -d $simudir/$param ]; then
	echo "comparing $param"
	./python/compare.py --plotdirs $simudir/$param:$datadir/$param --names simu:data --outdir $outdir/$param
	# exit 1
    fi
done



