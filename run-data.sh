#!/bin/bash

dtype=adaptive  #stanford, fiveprace

if [ "$dtype" == "stanford" ]; then
    datadir=/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers
    files=`ls $datadir/*`
    modulo=5
elif [ "$dtype" == "adaptive" ]; then
    datadir=/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw
    files=`ls $datadir/?/*-M_merged.tsv.bz2`  # if you switch to naive (N), be careful 'cause A is split in pieces
    modulo=10
    # prefix=03-
    # suffix=-M_merged.tsv.bz2
fi    
# fivepracedir=/shared/silo_researcher/Matsen_F/MatsenGrp/data/adaptive-five-prime-race
# for fname in $fivepracedir/*/*/*; do
# human=`echo $fname | xargs basename | sed 's/[nex.fq]//g'`
iproc=0
for subset in 4 5 6 7 8 9; do
    for fname in $files; do
        if [ "$dtype" == "stanford" ]; then
    	    human=`echo $fname | sed 's/_Lineages\.fasta//' | xargs basename`
        elif [ "$dtype" == "adaptive" ]; then
    	    human=`echo $fname | sed 's@/.*\([ABC]\)/.*@\1@'`
        fi
        echo $human
	
        # ./python/subset-data.py --infname $fname --outdir test/$dtype/$human --modulo 5 --start-indices 0:1:2:3:4:5:6:7:8:9 &
        # sleep 1
        # limit_procs python
        # continue
        subfname=test/$dtype/$human/every-$modulo-subset-$subset.csv.bz2

        label=every-$modulo-$human-subset-$subset
        ./plotperformance.py --label $label --plotdir $www/partis/$label \
        		     --datafname $subfname --action cache-data-parameters \
        		     --extra-args " --slurm:--workdir:tmp/$RANDOM --n-procs 50" &>>$label.log  &
        # --extra-args " --mimic-data-read-length" --n-queries 1000
        # if (( iproc >= 5 )); then
    	#     break
        # fi
        (( iproc++ ))
        sleep 30
    done
done
