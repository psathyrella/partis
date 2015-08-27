#!/usr/bin/env python
from collections import OrderedDict
import sys
sys.path.insert(1, './python')
import time
import os
from subprocess import Popen, PIPE, check_call, check_output
import argparse
import random

parser = argparse.ArgumentParser()
parser.add_argument('--action', required=True)
parser.add_argument('--timegrep', action='store_true')
parser.add_argument('--plot', action='store_true')
args = parser.parse_args()

fsdir = '/fh/fast/matsen_e/dralph/work/partis-dev/_output'
human = 'A'

if args.action == 'plot':
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
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
    # sns.set_style('ticks')  # hm, it actually works here
    fig, ax = plt.subplots()
    fig.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.95)
    
    timeinfo = OrderedDict()
    n_query_list =                        [100, 200, 500, 1000, 1500, 2000, 3000, 5000, 10000, 12000, 20000, 50000, 75000]
    timeinfo['vollmers-0.9'] =            [30,   34,  43,  218,  396,  217,  398,  534,   635,  None,  2247,  4152,  5589]
    timeinfo['mixcr'] =                   [7,     7,   8,    9,    9,   10,   10,   11,    13,  None,    16,    20,    20]
    timeinfo['changeo'] =                 [4,     4,   4,    6, None,    7, None,    9,  None,  None,  None,  None,  None]
    timeinfo['vsearch-partition'] =       [52,   53,  62,   70,  303,  408,  460,  498,   893,  None,  2561, 11209, 11413]
    timeinfo['naive-hamming-partition'] = [42,   52, 258,  138,  294,  277,  795, 2325, 13137,  None,  None,  None,  None]
    timeinfo['partition'] =               [80,   87, 147,  544, 1005, 1191, 2644, 7165, 24248, 38904,  None,  None,  None]

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
for n_queries in [12000, 15000]:  #[75000, 50000]:  #, 20000, 10000]:
    cmd = './bin/compare-partition-methods.py --actions ' + args.action
    istop = istart + n_queries
    # print '  %d queries from %d --> %d' % (n_queries, istart, istop)

    logfname = fsdir + '/' + human + '/_logs/istartstop-' + str(istart) + '-' + str(istop) + '/data-' + args.action + '.out'
    if args.timegrep:  # grep for timing info:
        # check_call(['ls', '-ltrh', logfname])
        outstr = check_output('grep \'total time\|mixcr time\' ' + logfname, shell=True)
        secs = float(outstr.split()[2])
        # print ' %9.0f' % secs,
        print '%.0f,' % secs,
        # print '  %5d\n' % n_queries,
        continue

    # actually run stuff:
    cmd += ' --istartstop ' + str(istart) + ':' + str(istop)
    cmd += ' --data'
    cmd += ' --humans A'
    # cmd += ' --overwrite'
    if args.action == 'run-mixcr':
        cmd += ' >' + logfname
    print cmd
    check_call(cmd, shell=True)
    # time.sleep(900)
    
    # istart = istop
    # break
