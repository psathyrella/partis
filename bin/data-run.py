#!/usr/bin/env python
import sys
import os
from subprocess import check_output
sys.path.insert(1, './python')
import humans

fs = '/fh/fast/matsen_e/dralph/work/partis-dev'


dset = 'stern'
human = 'SRR1383472' #SRR1383451 #SRR1383326
# outdir = fs + '/_output/' + datadir + '/' + human
nprocs = 30
# human=B
# seqfname=$fs/data/adaptive/$human/shuffled.csv
# outdir=$fs/_output/$human-all
# nprocs=100

outname = 'annotations'
outfname = os.path.dirname(humans.get_fname(dset, human)) + '/_output/' + outname + '.csv'
logfbase = os.path.dirname(humans.get_fname(dset, human)) + '/_logs/' + outname

n_queries = humans.get_nseqs(dset, human)

proc = Popen(cmd.split(), stdout=open(logfbase + '.out', 'w'), stderr=open(logfbase + '.err', 'w'))
# ./bin/partis.py --action cache-parameters --seqfile $seqfname \
# 		--workdir $fs/_tmp/$RANDOM --n-procs $nprocs --slurm \
# 		--parameter-dir $outdir/parameters/data &> _tmp/$human.log &
# cmd./bin/partis.py --action run-viterbi --seqfile $seqfname \
# 		--workdir $fs/_tmp/$RANDOM --n-procs $nprocs --slurm \
# 		--parameter-dir $outdir/parameters/data/hmm --outfname $outdir/annotations.csv &> _tmp/$human.log &

# ----------------------------------------------------------------------------------------
# echo "" >newmash.csv

# lines=200
# for hum in A B C; do
#     head -n$lines /fh/fast/matsen_e/dralph/work/partis-dev/data/adaptive/$hum/shuffled.csv |grep -v name|awk -F, '{print $2,$1}' >>newmash.csv
# done

# # # sort newmash.csv|uniq

# for hum in 18 19 59; do
#     bzgrep . -m$lines /fh/fast/matsen_e/dralph/work/partis-dev/data/stanford/021-0$hum/every-10-subset-0.csv.bz2 |grep -v unique_id|awk -F, '{print $1,$2}'>>newmash.csv
# done

# grep -v unique_id tmp.csv |awk -F, '{print $17,$1}' >>newmash.csv
# grep -v unique_id tmp-high.csv |awk -F, '{print $17,$1}' >>newmash.csv
