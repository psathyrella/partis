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
    # cmd_str += ' --no-clean'
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
    import seaborn as sns
    import numpy
    sns.set_style("ticks")

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
        if not os.path.exists(fname):
            print '  missing', fname
            return -1.0
        with open(fname) as infile:
            for line in infile.readlines():
                if 'adjusted mi' in line:
                    return float(line.split()[2])
        return -1.0

    # # get data from files
    # partis_data, voll_data = [], []
    # for ipt in range(n_points):
    #     label = baselabel + '-' + str(ipt)
    #     logfbase = '_logs/' + label
    #     partis_adj_mi = get_adj_mi(logfbase + '.out')
    #     voll_adj_mi = get_adj_mi(logfbase.replace(label, 'voll-' + label) + '.out')
    #     if partis_adj_mi == -1.0 or voll_adj_mi == -1.0:
    #         print 'missing', label
    #         continue
    #     print partis_adj_mi, voll_adj_mi
    #     partis_data.append(partis_adj_mi)
    #     voll_data.append(voll_adj_mi)

    # print partis_data
    # print voll_data
    # ...or copy and paste (after the first time)
    partis_data = [0.975846, 0.994257, 0.995104, 0.991261, 0.993894, 0.981661, 0.989986, 0.977794, 0.982451, 0.993858, 0.987053, 0.992395, 0.990745, 0.984822, 0.983463, 0.993787, 0.994101, 0.9859, 0.994276, 0.986123, 0.99385, 0.989804, 0.984593, 0.989288, 0.973727, 0.976963, 0.985062, 0.985637, 0.984954, 0.977783]
    voll_data = [0.516429, 0.554832, 0.528608, 0.483539, 0.476215, 0.554915, 0.523752, 0.546284, 0.569974, 0.491566, 0.558272, 0.520357, 0.573231, 0.527581, 0.565255, 0.502774, 0.506896, 0.550683, 0.523568, 0.486133, 0.520546, 0.535586, 0.526053, 0.516613, 0.578367, 0.490506, 0.536998, 0.532947, 0.494028, 0.530365]

    df = DataFrame({'partis':partis_data, 'annotation-based':voll_data})
    fig = df.plot(kind='hist', bins=50)  #, figsize=(6,4))
    fig.grid(False)
    # plt.locator_params(nbins=16)
    
    axes = plt.gca()
    ylimits = axes.get_ylim()
    # axes.set_ylim(ylimits[0], 1.2*ylimits[1])
    xmin, xmax = 0.3, 1.02
    plt.xlim(xmin, xmax)
    # axes.set_xticks([x for x in numpy.arange(xmin, xmax, 0.1)])
    plt.xlabel('adjusted mutual information')
    plt.ylabel('entries')
    plt.subplots_adjust(bottom=0.14)   # <--
    axes.spines["right"].set_visible(False)
    axes.spines["top"].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')

    plotdir = os.getenv('www') + '/clustering'
    plt.savefig(plotdir + '/adj-mi.svg')
    check_call(['./bin/permissify-www', plotdir])

plot(30)
sys.exit()

procs = []
voll=True
for iproc in range(15, 30):
# for iproc in range(15, 16):
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
