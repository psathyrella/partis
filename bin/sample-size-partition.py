#!/usr/bin/env python
from collections import OrderedDict
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.insert(1, './python')
import time
import os
from subprocess import Popen, PIPE, check_call, check_output
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument('--action', required=True)
args = parser.parse_args()

# ----------------------------------------------------------------------------------------
from plotting import legends, colors, linewidths
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
fig, ax = plt.subplots()
fig.tight_layout()
plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.95)

n_query_list = [100, 200, 500, 1000, 2000, 5000]
timeinfo = OrderedDict()
timeinfo['vollmers-0.9'] = [22, 22, 24, 310, 377, 571]
timeinfo['mixcr'] = [5, 5, 6, 7, 7, 9]
timeinfo['changeo'] = [4, 4, 4, 6, 7, 9]
timeinfo['vsearch-partition'] = [29, 33, 83, 83, 332, 595]
timeinfo['naive-hamming-partition'] = [62, 108, 1404, 1589, 1887, 4110]
timeinfo['partition'] = [230, 1482, 2357, 9071, 16763, 91098]
plots = {}
for meth, vals in timeinfo.items():
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
    plots[meth] = ax.plot(n_query_list, vals, linewidth=linewidths.get(meth, 4), label=legends.get(meth, meth), color=colors.get(meth, 'grey'), linestyle=linestyle, alpha=alpha)

legend = ax.legend(loc='upper left')
sns.despine(trim=True, bottom=True)
sns.set_style('ticks')
plt.xlabel('sample size')
plt.ylabel('time required')
plt.subplots_adjust(bottom=0.14, left=0.18)
ax.set_xscale('log')
ax.set_yscale('log')
yticks = [1, 60, 3600, 86400, 604800]  # seconds
yticklabels = ['1 sec', '1 min', '1 hour', '1 day', '1 week']
plt.yticks(yticks, yticklabels)
plt.savefig(os.getenv('www') + '/partis/tmp/time-required.svg')
sys.exit()

# ----------------------------------------------------------------------------------------
istart = 0
# for n_queries in [4, 5, 8, 15, 20, 25, 40, 50, 60, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]:
for n_queries in [100, 200, 500, 1000, 2000]:  #, 5000, ]:  # NOTE you've already used pretty much the whole sim file
    cmd = './bin/compare-partition-methods.py --actions ' + args.action
    istop = istart + n_queries
    # print '  %d queries from %d --> %d' % (n_queries, istart, istop)

    # grep for timing info:
    logfname = '/fh/fast/matsen_e/dralph/work/partis-dev/_output/A/_logs/istartstop-' + str(istart) + '-' + str(istop) + '/simu-10-leaves-1-mutate-' + args.action + '.out'
    outstr = check_output('grep \'total time\' ' + logfname, shell=True)
    secs = float(outstr.split()[2])
    print n_queries
    print '     %.0f' % secs,

    # actually run stuff:
    # cmd += ' --istartstop ' + str(istart) + ':' + str(istop)
    # cmd += ' --mutation-multiplier 1 --n-leaf-list 10'
    # print cmd
    # check_call(cmd.split())

    istart = istop


