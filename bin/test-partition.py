#!/usr/bin/env python
from subprocess import check_call, Popen, PIPE, check_output, call
import random
import time
import sys
import os

cmd = './bin/partis.py'
parameterlabel = 'every-10-A-subset-0'
baselabel = 'vary-n-leaves'
baseoutdir = '_output'
queries_per_proc = 1000
leaves_per_tree = 2  # closeish to the mean

# ----------------------------------------------------------------------------------------
def simulate(label):
    if os.path.exists(baseoutdir + '/' + label + '/simu.csv'):
        print '  simu.csv exists for', label
        return
    cmd_str = cmd + ' --action simulate'
    cmd_str += ' --outfname ' + baseoutdir + '/' + label + '/simu.csv'
    cmd_str += ' --parameter-dir ' + baseoutdir + '/' + parameterlabel + '/simu/hmm'
    cmd_str += ' --n-max-queries ' + str(int(float(queries_per_proc) / leaves_per_tree))
    cmd_str += ' --random-number-of-leaves --n-leaves ' + str(leaves_per_tree)
    cmd_str += ' --n-procs 15'
    print cmd_str
    check_call(cmd_str.split())

# ----------------------------------------------------------------------------------------
def partition(label):
    if os.path.exists('_logs/' + label + '.out'):
        print '  .out exists for', label
        return
    cmd_str = cmd + ' --action partition'
    cmd_str += ' --randomize-input-order'
    cmd_str += ' --seqfile ' + baseoutdir + '/' + label + '/simu.csv'
    cmd_str += ' --parameter-dir ' + baseoutdir + '/' + parameterlabel + '/simu/hmm'
    cmd_str += ' --no-clean --n-procs 100:5 --slurm --workdir tmp/' + str(random.randint(0, 99999))
    cmd_str += ' --n-max-queries ' + str(queries_per_proc)
    # cmd_str += ' --no-clean --n-procs 10:3 --slurm --workdir tmp/' + str(random.randint(0, 99999))
    # cmd_str += ' --n-max-queries 100'
    print cmd_str
    proc = Popen(cmd_str + ' 1>_logs/' + label + '.out 2>_logs/' + label + '.err', shell=True, stdout=PIPE)
    return proc

procs = []
for iproc in range(5, 10):
    label = baselabel + '-' + str(iproc)
    # simulate(label)
    proc = partition(label)
    procs.append(proc)
    time.sleep(30)
    # call('sacct | grep \'RUNNING\|PENDING\'', shell=True)
    n_slurm_jobs = int(check_output('sacct | grep \'RUNNING\|PENDING\' | wc -l', shell=True))
    # print n_slurm_jobs
    while n_slurm_jobs > 35:
        # call('sacct | grep \'RUNNING\|PENDING\'', shell=True)
        print '  %d jobs' % n_slurm_jobs
        time.sleep(10)
        n_slurm_jobs = int(check_output('sacct | grep \'RUNNING\|PENDING\' | wc -l', shell=True))

# cmd +  --action simulate --outfname _output/vary-n-leaves/simu.csv --parameter-dir _output/every-10-A-subset-0/simu/hmm --n-max-queries 5000 --n-procs 15 --n-leaves 2 --random-number-of-leaves
# cmd +  --action partition --randomize-input-order --seqfile _output/vary-n-leaves/simu.csv --parameter-dir _output/every-10-A-subset-0/simu/hmm --no-clean --n-procs 100:5 --slurm --workdir tmp/$RANDOM --n-max-queries 1000
# cmd +  --action run-viterbi --seqfile _output/vary-n-leaves/simu.csv --parameter-dir _output/every-10-A-subset-0/simu/hmm --n-max-queries 1000 --debug 0 --n-procs 10 --pants-seated-clustering --no-clean
