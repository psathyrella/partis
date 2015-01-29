#!/bin/bash

dtype=adaptive  #stanford, fiveprace
if [ "$dtype" == "stanford" ]; then
    datadir=/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers
    files=`ls $datadir`
elif [ "$dtype" == "adaptive" ]; then
    datadir=/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw
    files=`ls $datadir/?/*-M_merged.tsv.bz2`  # if you switch to naive (N), be careful 'cause A is split in pieces
    # prefix=03-
    # suffix=-M_merged.tsv.bz2
fi    
# fivepracedir=/shared/silo_researcher/Matsen_F/MatsenGrp/data/adaptive-five-prime-race
# for fname in $fivepracedir/*/*/*; do
# human=`echo $fname | xargs basename | sed 's/[nex.fq]//g'`

for fname in $files; do
    if [ "$dtype" == "stanford" ]; then
	human=`echo $fname | sed 's/_Lineages\.fasta//' | xargs basename`
    elif [ "$dtype" == "adaptive" ]; then
	human=`echo $fname | sed 's@/.*\([ABC]\)/.*@\1@'`
    fi
    echo $human

    # ./python/subset-data.py --infname $fname --outdir test/stanford/$human --modulo 5 --start-indices 0 &
    # ./python/subset-data.py --infname $fname --outdir test/fiveprace/$human --modulo 50 --start-indices 0 &
    # sleep 1
    # limit_procs python
    # continue

    subfname=test/$dtype/$human/*.bz2
    # bzgrep . $subfname | wc
    # continue

    label=mimic-$human
    ./plotperformance.py --label $label --plotdir $www/partis/$label \
			 --datafname $subfname --action cache-data-parameters \
			 --extra-args " --mimic-data-read-length:--slurm:--workdir:tmp/$RANDOM --n-procs 20" &
			 # --extra-args " --mimic-data-read-length" --n-queries 1000
    break
    sleep 30
done
