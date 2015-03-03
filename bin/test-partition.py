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
def runstuff(label, voll=False):
    logfbase = '_logs/' + label
    if voll:
        logfbase = logfbase.replace(label, 'voll-' + label)
    if os.path.exists(logfbase + '.out'):
        print '  log exists for', label
        return None

    cmd_str = cmd
    if voll:
        cmd_str += ' --action run-viterbi --pants-seated-clustering'
        cmd_str += ' --n-procs 10'
    else:
        cmd_str += ' --action partition'
        cmd_str += ' --n-procs 100:5 --slurm --workdir tmp/' + str(random.randint(0, 99999))
        # cmd_str += ' --n-procs 10:3 --slurm --workdir tmp/' + str(random.randint(0, 99999))

    cmd_str += ' --randomize-input-order'
    cmd_str += ' --no-clean'
    cmd_str += ' --seqfile ' + baseoutdir + '/' + label + '/simu.csv'
    cmd_str += ' --parameter-dir ' + baseoutdir + '/' + parameterlabel + '/simu/hmm'
    cmd_str += ' --n-max-queries ' + str(queries_per_proc)
    # cmd_str += ' --n-max-queries 100'
    print cmd_str
    proc = Popen(cmd_str + ' 1>' + logfbase + '.out 2>' + logfbase + '.err', shell=True)
    return proc

# ----------------------------------------------------------------------------------------
def plot(n_points):
    import matplotlib as mpl
    mpl.use('Agg')
    from pandas import DataFrame
    import matplotlib.pyplot as plt
    import pandas as pd

    fsize = 20
    mpl.rcParams.update({
        # 'font.size': fsize,
        'legend.fontsize': fsize,
        'axes.titlesize': fsize,
        # 'axes.labelsize': fsize,
        'xtick.labelsize': fsize,
        'ytick.labelsize': fsize,
        'axes.labelsize': fsize
    })

    def get_adj_mi(fname):
        with open(fname) as infile:
            for line in infile.readlines():
                if 'adjusted mi' in line:
                    return float(line.split()[2])
        return -1.0

    # get data from files
    # partis_data, voll_data = [], []
    # for ipt in range(n_points):
    #     label = baselabel + '-' + str(ipt)
    #     logfbase = '_logs/' + label
    #     partis_adj_mi = get_adj_mi(logfbase + '.out')
    #     voll_adj_mi = get_adj_mi(logfbase.replace(label, 'voll-' + label) + '.out')
    #     print partis_adj_mi, voll_adj_mi
    #     partis_data.append(partis_adj_mi)
    #     voll_data.append(voll_adj_mi)

    # ...or copy and paste (after the first time)
    partis_data = [0.975846, 0.994257, 0.995104, 0.991261, 0.993894, 0.981661, 0.989986, 0.977794, 0.982451, 0.993858]
    voll_data = [0.516429,  0.554832,  0.528608,  0.483539,  0.476215,  0.554915,  0.523752,  0.546284,  0.569974,  0.491566]

    df = DataFrame({'partis':partis_data, 'annotation-based':voll_data})
    fig = df.plot(kind='hist', bins=50)
    plt.xlim(0, 1.1)
    plt.xlabel('adjusted MI')
    plt.ylabel('entries')
    plotdir = os.getenv('www') + '/clustering'
    plt.savefig(plotdir + '/adj-mi.svg')
    check_call(['./bin/permissify-www', plotdir])

plot(10)
sys.exit()

procs = []
voll=True
for iproc in range(0, 10):
    label = baselabel + '-' + str(iproc)
    # simulate(label)
    proc = runstuff(label, voll=voll)
    procs.append(proc)
    if voll:
        time.sleep(15)
        continue
    time.sleep(30)
    n_slurm_jobs = int(check_output('sacct | grep \'RUNNING\|PENDING\' | wc -l', shell=True))
    while n_slurm_jobs > 35:
        print '  %d jobs' % n_slurm_jobs
        time.sleep(10)
        n_slurm_jobs = int(check_output('sacct | grep \'RUNNING\|PENDING\' | wc -l', shell=True))

# cmd +  --action simulate --outfname _output/vary-n-leaves/simu.csv --parameter-dir _output/every-10-A-subset-0/simu/hmm --n-max-queries 5000 --n-procs 15 --n-leaves 2 --random-number-of-leaves
# cmd +  --action partition --randomize-input-order --seqfile _output/vary-n-leaves/simu.csv --parameter-dir _output/every-10-A-subset-0/simu/hmm --no-clean --n-procs 100:5 --slurm --workdir tmp/$RANDOM --n-max-queries 1000
# cmd +  --action run-viterbi --seqfile _output/vary-n-leaves/simu.csv --parameter-dir _output/every-10-A-subset-0/simu/hmm --n-max-queries 1000 --debug 0 --n-procs 10 --pants-seated-clustering --no-clean
