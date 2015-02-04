#!/bin/bash

dtype=stanford  #adaptive fiveprace

modulo=10
if [ "$dtype" == "stanford" ]; then
    datadir=/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers
    files=`ls $datadir/*`
elif [ "$dtype" == "adaptive" ]; then
    datadir=/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw
    files=`ls $datadir/?/*-M_merged.tsv.bz2`  # if you switch to naive (N), be careful 'cause A is split in pieces
    # prefix=03-
    # suffix=-M_merged.tsv.bz2
fi    
# fivepracedir=/shared/silo_researcher/Matsen_F/MatsenGrp/data/adaptive-five-prime-race
# for fname in $fivepracedir/*/*/*; do
# human=`echo $fname | xargs basename | sed 's/[nex.fq]//g'`
iproc=0
missing =
('10-021-019_9',
 '10-021-048_1',
 '10-021-055_1',
 '10-021-063_0',
 '10-021-071_0',
 '10-021-084_0',
 '10-021-018_1')

for subset in {0..9}; do
    for fname in $files; do
        if [ "$dtype" == "stanford" ]; then
    	    human=`echo $fname | sed 's/_Lineages\.fasta//' | xargs basename`
        elif [ "$dtype" == "adaptive" ]; then
    	    human=`echo $fname | sed 's@/.*\([ABC]\)/.*@\1@'`
        fi
        echo $human
	
        # oh wait this loops over subsets DOH fix it ./python/subset-data.py --infname $fname --outdir test/$dtype/$human --modulo $modulo --start-indices 0:1:2:3:4:5:6:7:8:9 &
        # sleep 1
        # limit_procs python
        # continue
        subfname=test/$dtype/$human/every-$modulo-subset-$subset.csv.bz2

        label=every-$modulo-$human-subset-$subset
	./plotperformance.py --label $label --plotdir $www/partis/$label \
        		     --datafname $subfname --action cache-data-parameters \
        		     --extra-args " --slurm:--workdir:tmp/$RANDOM --n-procs 20" &>>$label.log  &
        # if (( iproc >= 5 )); then
    	#     break
        # fi
        (( iproc++ ))
	limit_procs python
        sleep 10
    done
done
