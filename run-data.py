#!/usr/bin/env python
import os
import sys
import random
import argparse
import re
import time
from subprocess import check_call, Popen
sys.path.insert(1, './python')

from humans import humans

parser = argparse.ArgumentParser()
parser.add_argument('--dtype', required=True, choices=['stanford', 'adaptive'])
args = parser.parse_args()

# dtype = 'adaptive'  #stanford fiveprace

modulo = 10

if args.dtype == 'stanford':
    datadir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers'
    files = os.listdir(datadir)
elif args.dtype == 'adaptive':
    datadirs = [ '/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/' + h for h in humans['adaptive'] ]
    files = []
    for datadir in datadirs:
        files += [ fname for fname in os.listdir(datadir) if '-M_merged.tsv.bz2' in fname ]  # if you switch to naive (N), be careful 'cause A is split in pieces

# iproc = 0
missing = None  #('10-021-084_2',)
# ('10-021-019_9',
#  '10-021-048_1',
#  '10-021-055_1',
#  '10-021-063_0',
#  '10-021-071_0',
#  '10-021-084_0',
#  '10-021-018_1')

action = 'cache-simu-parameters'

procs = []
for subset in range(1):  #modulo):
    for fname in files:
        if args.dtype == 'stanford':
    	    human = os.path.basename(fname).replace('_Lineages.fasta', '')
        elif args.dtype == 'adaptive':
    	    human = re.findall('[ABC]', fname)[0]

        label = 'every-' + str(modulo) + '-' + human + '-subset-' + str(subset)

        # if missing is not None and str(modulo) + '-' + human + '_' + str(subset) not in missing:
        #     continue
        # if os.path.exists('_output/' + label + '/simu.csv'):
        #     print 'exists:', '_output/' + label + '/simu.csv'
        #     continue
        if os.path.exists('_output/' + label + '/simu'):
            print 'ok', label
            continue

        print human, subset

        # oh wait this loops over subsets DOH fix it ./python/subset-data.py --infname $fname --outdir test/$dtype/$human --modulo $modulo --start-indices 0:1:2:3:4:5:6:7:8:9 &
        # sleep 1
        # check_call(['limit_procs', 'python'])
        # continue
        subfname = 'test/' + args.dtype + '/' + human + '/every-' + str(modulo) + '-subset-' + str(subset) + '.csv.bz2'

        cmd = './plotperformance.py --label ' + label + ' --plotdir ' + os.getenv('www') + '/partis/' + label + ' --action ' + action
        cmd += ' --n-procs 20 --n-sim-events 10000 --datafname ' + subfname + ' --extra-args __slurm:__workdir:tmp/' + str(random.randint(0,99999))
        # cmd += ' --n-queries 1000'
        # + ' --action cache-data-parameters \
        check_call(cmd.split())
        # procs.append(Popen(cmd.split(), stdout=open('_logs/' + label + '.out', 'w'), stderr=open('_logs/' + label + '.err', 'w')))

        # check_call(['limit_procs', 'python'])
        # time.sleep(10)
        # sys.exit()
