#!/bin/bash

h1=A #021-018
h2=B #021-081
h3=C #021-019
label1=mimic-$h1
label2=mimic-$h2
label3=mimic-$h3
basedir=/var/www/sharing/dralph/partis  #/$label/params
subdir=dj_insertion

# data vs simu:
./python/compare.py --plotdirs \
		    $basedir/$label1/params/data/hmm_parameters/$subdir:$basedir/$label1/params/simu/hmm_parameters/true/$subdir:$basedir/$label2/params/data/hmm_parameters/$subdir:$basedir/$label2/params/simu/hmm_parameters/true/$subdir:$basedir/$label3/params/data/hmm_parameters/$subdir:$basedir/$label3/params/simu/hmm_parameters/true/$subdir \
		    --names data-$h1:simu-$h1:data-$h2:simu-$h2:data-$h3:simu-$h3 --outdir /var/www/sharing/dralph/partis/$label1/data-vs-true-simu/$subdir \
		    --colors 632:632:596:596:418:418 --linestyles 1:2:1:2:1:2 \
		    --scale-errors 2.24
    # --rebin 3

# # data:
# ./python/compare.py --plotdirs \
# 		    $basedir/$label1/params/data/hmm_parameters/$subdir:$basedir/$label2/params/data/hmm_parameters/$subdir:$basedir/$label3/params/data/hmm_parameters/$subdir \
# 		    --names $h1:$h2:$h3 --outdir /var/www/sharing/dralph/partis/$label1/stanford-data/$subdir \
# # --colors 632:632:596:596:418:418 --linestyles 1:2:1:2:1:2

