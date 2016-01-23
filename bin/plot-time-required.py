#!/usr/bin/env python
from collections import OrderedDict
import sys
sys.path.insert(1, './python')
import time
import os
from subprocess import Popen, PIPE, check_call, check_output, CalledProcessError
import argparse
import random

import utils

parser = argparse.ArgumentParser()
parser.add_argument('--actions', required=True)
parser.add_argument('--timegrep', action='store_true')
parser.add_argument('--plot', action='store_true')
args = parser.parse_args()
args.actions = utils.get_arg_list(args.actions)

fsdir = '/fh/fast/matsen_e/dralph/work/partis-dev/_output'
human = '021-018'
istartstopstr_list = '0:250 250:750 750:1500 1500:2500 2500:4000 4000:6500 6500:9500 9500:13500 13500:18500 18500:26000 26000:36000 36000:51000 51000:71000 71000:101000 101000:141000 141000:191000'.split(' ')
istartstoplist = []
for istartstopstr in istartstopstr_list:
    istartstoplist.append([int(iss) for iss in istartstopstr.split(':')])
n_query_list = [istartstop[1] - istartstop[0] for istartstop in istartstoplist]

timeinfo = OrderedDict()
# NOTE edited some of these values for current code version Dec 22 2015
# n_query_list =                        [  100,   200,   500,  1000,  1500,  2000,  3000,  5000,  7000, 10000, 12000, 15000, 20000, 30000,  50000, 75000, 100000]
# timeinfo['vollmers-0.9'] =            [   30,    34,    43,   385,   396,   217,   398,   477,  None,  1241,  None,  1681,  2247,  None,   4153,  5590,   None]
# timeinfo['mixcr'] =                   [    7,     7,     8,     9,     9,    10,    10,    11,  None,    13,  None,  None,    16,  None,     20,    20,   None]
# timeinfo['changeo'] =                 [    4,     4,     4,     6,  None,     7,  None,     9,  None,    15,  None,  None,  None,  None,   None,  None,    457]
# timeinfo['vsearch-partition'] =       [   52,    53,    62,    70,   303,   408,   460,   477,  None,  1208,  None,  1586,  2561,  None,  11209, 11413,   None]
# timeinfo['naive-hamming-partition'] = [   33,    39,   208,    93,   170,   216,   649,   947,  None,  2372,  None,  4655,  None,  None,   None,  None,   None]
# timeinfo['partition'] =               [   53,    67,   138,   217,  1026,  1623,  2847,  3052,  None,  8228,  None, 11542, 21391, 31418,   None,  None,   None]
def make_plot():
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    from plotting import legends, colors, linewidths, interpolate_values
    fsize = 20
    mpl.rcParams.update({
        # 'font.size': fsize,
        'legend.fontsize': 15,
        'axes.titlesize': fsize,
        # 'axes.labelsize': fsize,
        'xtick.labelsize': fsize,
        'ytick.labelsize': fsize,
        'axes.labelsize': fsize
    })
    # sns.set_style('ticks')  # hm, it actually works here
    fig, ax = plt.subplots()
    fig.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.95)
    
    plots = {}
    for meth, vals in timeinfo.items():
        interpolate_values(n_query_list, vals)

        linestyle = '-'
        alpha = 1.
        if 'vollmers' in meth:
            alpha = 0.5
        if 'vsearch' in meth:
            linestyle = '-.'
        elif 'naive-hamming-partition' in meth:
            linestyle = '--'
        elif 'true' in meth:
            linestyle = '--'
            alpha = 0.5
        plots[meth] = ax.plot(n_query_list, vals, linewidth=linewidths.get(meth, 4), label=legends.get(meth, meth), color=colors.get(meth, 'grey'), linestyle=linestyle, alpha=alpha)  #, markersize=1000)
    
    legend = ax.legend(loc='upper left')
    sns.despine()  #trim=True, bottom=True)
    plt.xlabel('sample size')
    plt.ylabel('time required')
    plt.subplots_adjust(bottom=0.14, left=0.18)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(True)
    yticks = [1, 60, 3600, 86400, 604800]  # seconds
    yticklabels = ['1 sec', '1 min', '1 hour', '1 day', '1 week']
    plt.yticks(yticks, yticklabels)
    plt.savefig(os.getenv('www') + '/partis/clustering/time-required.svg')
    sys.exit()

def get_clock_time(istart, istop, action):
    logfname = fsdir + '/' + human + '/istartstop-' + str(istart) + '-' + str(istop) + '/_logs/simu-7-leaves-1.0-mutate-' + action + '.out'
    if args.timegrep:
        # check_call(['ls', '-ltrh', logfname])
        try:
            outstr = check_output('grep \'total time\|mixcr time\' ' + logfname, shell=True)
            secs = float(outstr.split()[2])
            return secs
            # print '%5.0f,' % secs,
        except CalledProcessError:
            return None
            # print '%5s,' % 'None',
        # print '  %5d\n' % n_queries,

# ----------------------------------------------------------------------------------------
for action in args.actions:
    aname = action
    if action == 'run-viterbi':
        aname = 'vollmers-0.9'
    timeinfo[aname] =  []
    for istart, istop in istartstoplist:
        timeinfo[aname].append(get_clock_time(istart, istop, action))

make_plot()
