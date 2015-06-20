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

fsdir = '/fh/fast/matsen_e/' + os.getenv('USER') + '/work/partis-dev'

parser = argparse.ArgumentParser()
parser.add_argument('--dataset', choices=['stanford', 'adaptive'], default='adaptive')
parser.add_argument('--only-run')  # colon-separated list of human,subset pairs to run, e.g. A,3:C,8
parser.add_argument('--n-procs', type=int, required=True)
parser.add_argument('--action', required=True)
args = parser.parse_args()
args.only_run = utils.get_arg_list(args.only_run)

if args.dataset == 'stanford':
    datadir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers'
    files = [ datadir + '/' + f for f in os.listdir(datadir)]
elif args.dataset == 'adaptive':
    datadirs = [ '/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/' + h for h in humans['adaptive'] ]
    files = []
    for datadir in datadirs:
        files += [ datadir + '/' + fname for fname in os.listdir(datadir) if '-M_merged.tsv.bz2' in fname ]  # if you switch to naive (N), be careful 'cause A is split in pieces

for fname in files:
    if args.dataset == 'stanford':
	    human = os.path.basename(fname).replace('_Lineages.fasta', '')
    elif args.dataset == 'adaptive':
	    human = re.findall('[ABC]', fname)[0]
    print 'run', human
    label = human

    if args.only_run is not None and human not in args.only_run:
        continue

    # procs = []
    for n_leaves in [3, 5, 10, 25, 50]:
    
        # if os.path.exists('_output/' + label + '/simu.csv'):
        #     print 'exists:', '_output/' + label + '/simu.csv'
        #     continue

        simfname = fsdir + '/_output/' + label + '/simu-' + str(n_leaves) + '-leaves.csv'

        cmd = './bin/run-driver.py --label ' + label + ' --plotdir ' + os.getenv('www') + '/partis/ --action ' + args.action + ' --n-procs ' + str(args.n_procs)
        extras = []
        if args.action == 'cache-data-parameters':
            cmd += ' --datafname ' + fname
        elif args.action == 'simulate':
            cmd += ' --simfname ' + simfname
            extras += ['--n-sim-events', str(int(float(10000) / n_leaves))]
            extras += ['--n-leaves', str(n_leaves)]
        elif args.action == 'cache-simu-parameters':
            cmd += ' --simfname ' + simfname
        elif args.action == 'partition':
            cmd += ' --simfname ' + simfname
            cmd += ' --n-queries 100'
    
        if args.n_procs > 10:
            extras += ['--slurm', '--workdir', '_tmp/' + str(random.randint(0,99999))]
    
        cmd += utils.get_extra_str(extras)
        print '   ' + cmd
        check_call(cmd.split())
        # sys.exit()
    
        # procs.append(Popen(cmd.split(), stdout=open('_logs/' + label + '.out', 'w'), stderr=open('_logs/' + label + '.err', 'w')))
    
        # check_call(['limit_procs', 'python'])
        # time.sleep(10)
    sys.exit()
