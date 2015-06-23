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
all_actions = ['cache-data-parameters', 'simulate', 'cache-simu-parameters', 'partition', 'run-viterbi']
parser.add_argument('--actions', default=':'.join(all_actions))
args = parser.parse_args()
args.only_run = utils.get_arg_list(args.only_run)
args.actions = utils.get_arg_list(args.actions)

if args.dataset == 'stanford':
    datadir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers'
    files = [ datadir + '/' + f for f in os.listdir(datadir)]
elif args.dataset == 'adaptive':
    datadirs = [ '/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/' + h for h in humans['adaptive'] ]
    files = []
    for datadir in datadirs:
        files += [ datadir + '/' + fname for fname in os.listdir(datadir) if '-M_merged.tsv.bz2' in fname ]  # if you switch to naive (N), be careful 'cause A is split in pieces
    
# ----------------------------------------------------------------------------------------
def execute(action, label, n_leaves, fname, mute_mult):
    simfname = fsdir + '/_output/' + label + '/simu-' + str(n_leaves) + '-leaves-' + str(mute_mult) + '-mutate.csv'
    cmd = './bin/run-driver.py --label ' + label + ' --action ' + action
    extras = []

    if action == 'cache-data-parameters':
        cmd += ' --datafname ' + fname
        extras += ['--n-max-queries', + n_data_to_cache]
        n_procs = n_data_to_cache / 500
    elif action == 'simulate':
        cmd += ' --simfname ' + simfname
        extras += ['--n-sim-events', int(float(n_sim_seqs) / n_leaves)]
        extras += ['--n-leaves', n_leaves, '--mutation-multiplier', mute_mult]
        n_procs = 10
    elif action == 'cache-simu-parameters':
        cmd += ' --simfname ' + simfname
        n_procs = 20
    elif action == 'partition':
        cmd += ' --simfname ' + simfname
        extras += ['--n-max-queries', n_to_partition]
        n_procs = n_to_partition / 15  # something like 15 seqs/process to start with
    elif action == 'run-viterbi':
        extras += ['--annotation-clustering', 'vollmers', '--annotation-clustering-thresholds', '0.5:0.7:0.8:0.9:0.95']
        cmd += ' --simfname ' + simfname
        extras += ['--n-max-queries', n_to_partition]
        n_procs = n_to_partition / 100

    cmd +=  ' --plotdir ' + os.getenv('www') + '/partis --n-procs ' + str(n_procs)

    if n_procs > 10:
        extras += ['--slurm', '--workdir', '_tmp/' + str(random.randint(0,99999))]

    # extras += ['--print-partitions', ]

    cmd += utils.get_extra_str(extras)
    print '   ' + cmd
    check_call(cmd.split())
    # sys.exit()
        
# ----------------------------------------------------------------------------------------
n_to_partition = 1000
n_data_to_cache = 50000
mutation_multipliers = ['1', ]  #'2', '4']
n_sim_seqs = 10000
# ----------------------------------------------------------------------------------------
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
    for n_leaves in [3, 5]:  #, 10, 25, 50]:
        # if os.path.exists('_output/' + label + '/simu.csv'):
        #     print 'exists:', '_output/' + label + '/simu.csv'
        #     continue
        for action in args.actions:
            print '  ----> ', action
            for mute_mult in mutation_multipliers:
                print '    ----> ', mute_mult
                execute(action, label, n_leaves, fname, mute_mult)
                if action == 'cache-data-parameters':
                    continue

        # procs.append(Popen(cmd.split(), stdout=open('_logs/' + label + '.out', 'w'), stderr=open('_logs/' + label + '.err', 'w')))
    
        # check_call(['limit_procs', 'python'])
        # time.sleep(10)
        sys.exit()
