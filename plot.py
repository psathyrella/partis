#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import os
import csv
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from subprocess import check_call
from pandas import DataFrame
sns.set_style("ticks")

fsize = 20
mpl.rcParams.update({
    # 'font.size': fsize,
    'font.family': 'serif',
    'legend.fontsize': fsize,
    'axes.titlesize': fsize,
    # 'axes.labelsize': fsize,
    'xtick.labelsize': fsize,
    'ytick.labelsize': fsize,
    'axes.labelsize': fsize
})

cmap = mpl.cm.jet

logprobs, adj_mis = {}, {}
max_length = -1
with open('new.csv') as infile:
    for line in csv.DictReader(infile):
        ipath = int(line['path_index'])
        if ipath not in logprobs:
            logprobs[ipath] = []
            adj_mis[ipath] = []
        # if len(logprobs[ipath]) > 0 and float(line['score']) <= logprobs[ipath][-1]:
        #     continue
        logprobs[ipath].append(float(line['score']))
        adj_mis[ipath].append(float(line['adj_mi']))
        if len(logprobs[ipath]) > max_length:
            max_length = len(logprobs[ipath])

min_logprob, max_logprob = None, None
for ipath in range(len(logprobs)):
    while len(logprobs[ipath]) < max_length:
        logprobs[ipath].append(logprobs[ipath][-1])
        adj_mis[ipath].append(adj_mis[ipath][-1])
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
ax = fig.add_subplot(1, 1, 1)

scatter = False
if scatter:
    linex = (0, max_length)
    liney = (1, 1)
    ax.plot(linex, liney)
    ax.set_xlim(0, max_length)
    ax.set_ylim(-0.1, 1.1)
    ax.set_xlabel('step')
    ax.set_ylabel('adjusted MI')
    fig.tight_layout()

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
    # iplot = 0
    for iplot in range(len(logprobs)):
        # istart = 140
        # logprobs[iplot] = logprobs[iplot][istart:]
        # adj_mis[iplot] = adj_mis[iplot][istart:]
    
        plt.scatter(logprobs[iplot], adj_mis[iplot], c=[i for i in range(len(logprobs[iplot]))], cmap=cmap)
        plt.plot(logprobs[iplot], adj_mis[iplot], linewidth=0.2, color='black')

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
plt.savefig(plotdir + '/foo.png')
check_call(['permissify-www', plotdir])
