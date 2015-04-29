#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import argparse
import sys
import os
import glob
import csv
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from subprocess import check_call
from pandas import DataFrame
sns.set_style("ticks")

parser = argparse.ArgumentParser()
parser.add_argument('--scatter', action='store_true')
parser.add_argument('--zoom', action='store_true')
args = parser.parse_args()

fsize = 26
mpl.rcParams.update({
    'font.size': 26,
    'axes.labelsize': 26,
    'xtick.labelsize':20,
    'ytick.labelsize':20,
    'font.family': 'Lato',
    'font.weight': 600,
    'axes.labelweight': 600,
    # 'font.size': fsize,
    'legend.fontsize': fsize,
    'axes.titlesize': fsize
    # 'axes.labelsize': fsize,
})

cmap = 'BuGn'  #mpl.cm.jet

logprobs, adj_mis = {}, {}
fnames = glob.glob('24[0-9]*.csv')
for fname in fnames:
    with open(fname) as infile:
        for line in csv.DictReader(infile):
            ipath = int(line['path_index'])
            if ipath not in logprobs:
                logprobs[ipath] = []
                adj_mis[ipath] = []
            # if len(logprobs[ipath]) > 0 and float(line['score']) <= logprobs[ipath][-1]:
            #     continue
            logprobs[ipath].append(float(line['score']))
            adj_mis[ipath].append(float(line['adj_mi']))

# n_skip = 350
# for ipath in logprobs.keys():
#     logprobs[ipath] = logprobs[ipath][n_skip:]
#     adj_mis[ipath] = adj_mis[ipath][n_skip:]
    
max_length = -1
for ipath in logprobs.keys():
    if len(logprobs[ipath]) > max_length:
        max_length = len(logprobs[ipath])

min_logprob, max_logprob = None, None
print '%d paths' % len(logprobs.keys())
for ipath in logprobs.keys():
    while len(logprobs[ipath]) < max_length:
        logprobs[ipath].append(None)
        adj_mis[ipath].append(None)
    for il in range(len(logprobs[ipath])):
        if min_logprob is None or logprobs[ipath][il] < min_logprob:
            min_logprob = logprobs[ipath][il]
        if max_logprob is None or logprobs[ipath][il] > max_logprob:
            max_logprob = logprobs[ipath][il]
    # for il in range(len(logprobs[ipath])):
    #     logprobs[ipath][il] /= max_logprob
    # for il in range(len(logprobs[ipath])):
    #     print logprobs[ipath][il]

fig = plt.figure(1)
fig.clf()
fig, ax = plt.subplots()
adj_mi_color = '#980000'
logprob_color = '#1947A3'

if args.scatter:
    # linex = (0, max_length)
    # liney = (1, 1)
    # ax.plot(linex, liney)

    ax2 = ax.twinx()
    nxbins = 8
    nybins = 5

    yextrafactor = 1.
    if args.zoom:
        xmin = 350
        ymin = 0.7
        markersize = 25
    else:
        xmin = 0
        ymin = 0.
        markersize = 4
    ax.set_xlim(xmin, max_length)
    ax.set_ylim(ymin, yextrafactor * 1.)
    ax2.set_ylim(min_logprob, yextrafactor * max_logprob)
    # ax.set_ylim(min_logprob, yextrafactor * max_logprob)
    ax.set_xlabel('agglomeration step', fontweight='bold')
    ax.set_ylabel('adjusted MI', color=adj_mi_color, fontweight='bold')
    ax2.set_ylabel('log prob', color=logprob_color, fontweight='bold')
    # ax.set_ylabel('log prob', color=logprob_color, fontweight='bold')
    fig.tight_layout()
    plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.95)

    ax.locator_params(nbins=nxbins, axis='x')
    ax.locator_params(nbins=nybins, axis='y')
    # ax2.locator_params(nbins=nybins, axis='y')

    one_trace = False
    if one_trace:
        npoints = 10
        xposes = [30 for i in range(npoints)]
        yposes = [0.4*i/npoints for i in range(npoints)]
        legvals = [min_logprob + (max_logprob - min_logprob)*float(i)/npoints for i in range(npoints)]
        ax.scatter(xposes, yposes, c=legvals, cmap=cmap)
        ax.annotate('log likelihood', (30, yposes[-1]+.04))
        for i, txt in enumerate(legvals):
            ax.annotate(txt, (xposes[i] + 1, yposes[i] - .01))
        
        for ipath in range(len(logprobs)):
            ax.scatter([i for i in range(len(logprobs[ipath]))], adj_mis[ipath], c=logprobs[ipath], cmap=cmap)
    else:
        for ipath in logprobs.keys():
            steps = [i for i in range(len(logprobs[ipath]))]
            sizes = [markersize for i in range(len(logprobs[ipath]))]
            fig_adj_mi = ax.scatter(steps, adj_mis[ipath], color=adj_mi_color, alpha=1, s=sizes)
            fig_logprob = ax2.scatter(steps, logprobs[ipath], color=logprob_color, alpha=1, s=sizes)
        
else:
    # iplot = 0
    for iplot in logprobs.keys():
        # istart = 140
        # logprobs[iplot] = logprobs[iplot][istart:]
        # adj_mis[iplot] = adj_mis[iplot][istart:]
    
        plt.scatter(logprobs[iplot], adj_mis[iplot], c=[i for i in range(len(logprobs[iplot]))], cmap=cmap)

    npoints = 10
    xposes = [max_logprob for i in range(npoints)]
    yposes = [0.5*i/npoints for i in range(npoints)]
    legvals = [i*len(logprobs[iplot])/npoints for i in range(npoints)]
    ax.annotate('agglom. step', (max_logprob, yposes[-1]+.04))
    ax.scatter(xposes, yposes, c=legvals, cmap=cmap)
    for i, txt in enumerate(legvals):
        ax.annotate(txt, (xposes[i] + 75, yposes[i] - .02))
    
    linex = (min_logprob, max_logprob)
    liney = (1, 1)
    ax.plot(linex, liney)

    ax.set_xlabel('log likelihood')
    ax.set_ylabel('adjusted MI')

plotdir = os.getenv('www') + '/tmp'
if args.zoom:
    plotname = 'foo-zoom'
else:
    plotname = 'foo'
plt.savefig(plotdir + '/' + plotname + '.png')
check_call(['permissify-www', plotdir])
