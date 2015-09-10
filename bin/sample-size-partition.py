#!/usr/bin/env python
from collections import OrderedDict
import sys
sys.path.insert(1, './python')
import time
import os
from subprocess import Popen, PIPE, check_call, check_output, CalledProcessError
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument('--action', required=True)
parser.add_argument('--timegrep', action='store_true')
parser.add_argument('--plot', action='store_true')
args = parser.parse_args()

fsdir = '/fh/fast/matsen_e/dralph/work/partis-dev/_output'
human = 'A'

timeinfo = OrderedDict()
n_query_list =                        [  100,   200,   500,  1000,  1500,  2000,  3000,  5000,  7000, 10000, 12000, 15000, 20000, 50000, 75000, 100000]
timeinfo['vollmers-0.9'] =            [   30,    34,    43,   385,   396,   217,   398,   534,  None,   635,  None,  None,  2247,  4153,  5590,   None]  #   30,    34,    43,   218,   396,   217,   398,   534,  None,   635,  None,  None,  2247,  4152,  5589]
timeinfo['mixcr'] =                   [    7,     7,     8,     9,     9,    10,    10,    11,  None,    13,  None,  None,    16,    20,    20,   None]
timeinfo['changeo'] =                 [    4,     4,     4,     6,  None,     7,  None,     9,  None,    15,  None,  None,  None,  None,  None,    457]
timeinfo['vsearch-partition'] =       [   52,    53,    62,    70,   303,   408,   460,   498,  None,   893,  None,  None,  2561, 11209, 11413,   None]
timeinfo['naive-hamming-partition'] = [   33,    39,   208,    93,   170,   216,   649,  1938,  None, 13137,  None, 23316,  None,  None,  None,   None]  #,42,     52,   258,   138,   294,   277,   795,  2325,  None, 13137,  None, 23316,  None,  None,  None]
timeinfo['partition'] =               [   53,    67,   138,   217,  1026,  1623,  2847,  5419,  None, 24248, 38905, 44530,  None,  None,  None,   None]  #   80,    87,   147,   214,  1005,  1191,  2644,  7165,  None, 24248, 38904, 44530,  None,  None,  None]  # new 15k: 63965

if args.action == 'plot':
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
    sns.despine(trim=True, bottom=True)
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
# for n_queries in [100, 200, 500, 1000, 1500, 2000, 3000, 5000, 10000, 20000, 50000]:
# for n_queries in [7000, 10000, 20000, 50000]:
# for n_queries in n_query_list:
for n_queries in [10000]:
    cmd = './bin/compare-partition-methods.py --actions ' + args.action
    istop = istart + n_queries
    # print '  %d queries from %d --> %d' % (n_queries, istart, istop)

    logfname = fsdir + '/' + human + '/_logs/istartstop-' + str(istart) + '-' + str(istop) + '/data-' + args.action + '.out'
    if args.timegrep:  # grep for timing info:
        # check_call(['ls', '-ltrh', logfname])
        try:
            outstr = check_output('grep \'total time\|mixcr time\' ' + logfname, shell=True)
            secs = float(outstr.split()[2])
            print '%5.0f,' % secs,
        except CalledProcessError:
            print '%5s,' % 'None',
        # print '  %5d\n' % n_queries,

        continue

    # actually run stuff:
    cmd += ' --istartstop ' + str(istart) + ':' + str(istop)
    cmd += ' --data'
    cmd += ' --humans 021-018'  #021-050:021-055:021-057:021-059:021-060:021-061:021-063:021-068:021-071:021-081:021-084'
    cmd += ' --dataset stanford'  #adaptive'
    # cmd += ' --overwrite'

    # uncomment this if you want to be able to grep for the timing info
    # if args.action == 'run-mixcr':
    #     cmd += ' >' + logfname

    print cmd
    check_call(cmd, shell=True)
    # time.sleep(900)
    
    # istart = istop
    break
