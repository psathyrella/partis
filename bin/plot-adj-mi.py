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

n_leaf_list = [5, 10, 25, 50]
adj_mis = OrderedDict()
adj_mis[1] = OrderedDict()
adj_mis[4] = OrderedDict()
adj_mis[1]['vollmers-0.9'] =               [(0.31, 0.02), (0.28, 0.05), (0.31, 0.05), (0.28, 0.06)]
# adj_mis[1]['mixcr'] =                    [(None, None), (None, None), (None, None), (None, None)]
adj_mis[1]['changeo'] =                    [(0.43, 0.06), (0.38, 0.08), (0.43, 0.09), (0.34, 0.13)]
adj_mis[1]['vsearch-partition'] =          [(0.90, 0.01), (0.90, 0.03), (0.90, 0.03), (0.89, 0.05)]
adj_mis[1]['naive-hamming-partition'] =    [(0.96, 0.01), (0.95, 0.02), (0.97, 0.02), (0.91, 0.05)]
adj_mis[1]['partition'] =                  [(0.98, 0.01), (0.94, 0.03), (0.95, 0.02), (0.84, 0.05)]
adj_mis[4]['vollmers-0.9'] =               [(0.05, 0.01), (0.05, 0.01), (0.04, 0.01), (0.04, 0.02)]
# adj_mis[4]['mixcr'] =                    [(None, None), (None, None), (None, None), (None, None)]
adj_mis[4]['changeo'] =                    [(0.20, 0.05), (0.14, 0.07), (0.13, 0.08), (0.07, 0.03)]
adj_mis[4]['vsearch-partition'] =          [(0.67, 0.03), (0.70, 0.02), (0.70, 0.05), (0.67, 0.08)]
adj_mis[4]['naive-hamming-partition'] =    [(0.76, 0.03), (0.85, 0.02), (0.88, 0.02), (0.87, 0.04)]
adj_mis[4]['partition'] =                  [(0.86, 0.02), (0.93, 0.01), (0.86, 0.03), (0.82, 0.03)]

sns.set_style('ticks')
def make_plot(mut_mult):
    fig, ax = plt.subplots()
    fig.tight_layout()
    plots = {}
    for meth, valerrs in adj_mis[mut_mult].items():
        vals = [ve[0] for ve in valerrs]
        errs = [ve[1] for ve in valerrs]
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
        if meth is 'mixcr':
            plots[meth] = ax.plot(n_leaf_list, vals, linewidth=linewidths.get(meth, 4), label=legends.get(meth, meth), color=colors.get(meth, 'grey'), linestyle=linestyle, alpha=alpha)
        else:
            plots[meth] = ax.errorbar(n_leaf_list, vals, yerr=errs, linewidth=linewidths.get(meth, 4), label=legends.get(meth, meth), color=colors.get(meth, 'grey'), linestyle=linestyle, alpha=alpha, fmt='-o')
        # plots[meth][-1][0].set_linewidth(linewidths.get(meth, 4))
    
    # legend = ax.legend(loc='center left')
    if mut_mult == 1:
        ly = 0.85
    else:
        ly = 0.62
    legend = ax.legend(bbox_to_anchor=(0.5, ly))
    ax.set_xlim(3, 55)
    ax.set_ylim(0, 1)
    sns.despine(trim=True, bottom=True)
    sns.set_style('ticks')
    plt.title('%dx mutation' % mut_mult)
    plt.xlabel('mean N leaves')
    plt.ylabel('adjusted MI')
    plt.gcf().subplots_adjust(bottom=0.14, left=0.18, right=0.78, top=0.95)
    xticks = n_leaf_list
    xticklabels = [str(xt) for xt in xticks]
    plt.xticks(xticks, xticklabels)
    yticks = [0., .2, .4, .6, .8, 1.]
    yticklabels = [str(yt) for yt in yticks]
    plt.yticks(yticks, yticklabels)
    plt.savefig(os.getenv('www') + '/partis/tmp/adj-mi-%d-mutation.svg' % mut_mult)

for mut_mult in [1, 4]:
    make_plot(mut_mult)
