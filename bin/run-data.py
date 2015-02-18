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
import utils

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', choices=['stanford', 'adaptive'], default='adaptive')
parser.add_argument('--only-run')  # colon-separated list of human,subset pairs to run, e.g. A,3:C,8
parser.add_argument('--action', required=True)
args = parser.parse_args()
args.only_run = utils.get_arg_list(args.only_run)
if args.only_run is not None:
    tmp_items = []
    for item in args.only_run:
        tmp_items.append(item.split(','))
    args.only_run = tmp_items

if args.dataset == 'stanford':
    datadir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers'
    files = os.listdir(datadir)
elif args.dataset == 'adaptive':
    datadirs = [ '/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/' + h for h in humans['adaptive'] ]
    files = []
    for datadir in datadirs:
        files += [ fname for fname in os.listdir(datadir) if '-M_merged.tsv.bz2' in fname ]  # if you switch to naive (N), be careful 'cause A is split in pieces

modulo = 10

procs = []
for subset in range(modulo):
    for fname in files:
        if args.dataset == 'stanford':
    	    human = os.path.basename(fname).replace('_Lineages.fasta', '')
        elif args.dataset == 'adaptive':
    	    human = re.findall('[ABC]', fname)[0]

        label = 'every-' + str(modulo) + '-' + human + '-subset-' + str(subset)

        if args.only_run is not None:
            run = False
            for h, s in args.only_run:
                if human == h and subset == int(s):
                    run = True
        if not run:
            continue
        # if os.path.exists('_output/' + label + '/simu.csv'):
        #     print 'exists:', '_output/' + label + '/simu.csv'
        #     continue
        # if os.path.exists('_output/' + label + '/simu'):
        #     print 'ok', label
        #     continue

        print 'run', human, subset

        # oh wait this loops over subsets DOH fix it ./python/subset-data.py --infname $fname --outdir test/$dtype/$human --modulo $modulo --start-indices 0:1:2:3:4:5:6:7:8:9 &
        # sleep 1
        # check_call(['limit_procs', 'python'])
        # continue
        subfname = 'test/' + args.dataset + '/' + human + '/every-' + str(modulo) + '-subset-' + str(subset) + '.csv.bz2'

        cmd = './bin/run-driver.py --label ' + label + ' --plotdir ' + os.getenv('www') + '/partis/' + label + ' --action ' + args.action
        cmd += ' --n-procs 20 --n-sim-events 10000 --datafname ' + subfname + ' --extra-args __slurm:__workdir:tmp/' + str(random.randint(0,99999))
        # check_call(cmd.split())
        # procs.append(Popen(cmd.split(), stdout=open('_logs/' + label + '.out', 'w'), stderr=open('_logs/' + label + '.err', 'w')))

        # check_call(['limit_procs', 'python'])
        # time.sleep(10)
        # sys.exit()
