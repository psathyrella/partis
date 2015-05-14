#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import argparse
import sys
import os
import glob
import csv
import seaborn as sns
import matplotlib.pyplot as plt
from subprocess import check_call
sns.set_style("ticks")
sys.path.insert(1, './python')

from utils import get_arg_list

parser = argparse.ArgumentParser()
parser.add_argument('--zoom', action='store_true')
parser.add_argument('--normalize-logprobs', action='store_true')
parser.add_argument('--xbounds')
parser.add_argument('--logprob-bounds')
parser.add_argument('--adjmi-bounds')
args = parser.parse_args()
if args.xbounds is not None:
    args.xbounds = get_arg_list(args.xbounds, floatify=True)
if args.logprob_bounds is not None:
    args.logprob_bounds = get_arg_list(args.logprob_bounds, floatify=True)
if args.adjmi_bounds is not None:
    args.adjmi_bounds = get_arg_list(args.adjmi_bounds, floatify=True)

fsize = 26
mpl.rcParams.update({
    'font.size': 26,
    'axes.labelsize': 26,
    'xtick.labelsize':20,
    'ytick.labelsize':20,
    'font.family': 'Lato',
    'font.weight': 600,
    'axes.labelweight': 600,
    'legend.fontsize': fsize,
    'axes.titlesize': fsize
})

cmap = 'BuGn'  #mpl.cm.jet

logprobs, adj_mis = {}, {}
final_logweights = []
# fnames = ['/fh/fast/matsen_e/dralph/work/partis-dev/_output/sample-sizes/' + str(n) + '-queries.csv' for n in (10, 20, 50, 100)]
fnames = ['/fh/fast/matsen_e/dralph/work/partis-dev/_output/sample-sizes/' + str(n) + '-queries.csv' for n in (3000, )]
# fnames = glob.glob('[23][0-9]*.csv')
for fname in fnames:
    print fname
    with open(fname) as infile:
        for line in csv.DictReader(infile):
            ipath = int(line['path_index'])
            if ipath not in logprobs:
                logprobs[ipath] = []
                adj_mis[ipath] = []
            logprobs[ipath].append(float(line['score']))
            adj_mis[ipath].append(float(line['adj_mi']))
            if float(line['adj_mi']) == 1.:
                final_logweights.append(float(line['logweight']))

# for lw in final_logweights:
#     print lw
# sys.exit()
max_length = -1
for ipath in logprobs.keys():
    if len(logprobs[ipath]) > max_length:
        max_length = len(logprobs[ipath])

min_logprob, max_logprob = None, None
min_logprobs, max_logprobs = {}, {}  #[None for _ in range(len(logprobs.keys()))], [None for _ in range(len(logprobs.keys()))]
print '%d paths' % len(logprobs.keys())
for ipath in logprobs.keys():
    while len(logprobs[ipath]) < max_length:
        logprobs[ipath].append(None)
        adj_mis[ipath].append(None)
    min_logprobs[ipath], max_logprobs[ipath] = None, None
    for il in range(len(logprobs[ipath])):
        # min/max for this path
        if min_logprobs[ipath] is None or logprobs[ipath][il] < min_logprobs[ipath]:
            min_logprobs[ipath] = logprobs[ipath][il]
        if max_logprobs[ipath] is None or logprobs[ipath][il] > max_logprobs[ipath]:
            max_logprobs[ipath] = logprobs[ipath][il]
        # global min/max
        if min_logprob is None or logprobs[ipath][il] < min_logprob:
            min_logprob = logprobs[ipath][il]
        if max_logprob is None or logprobs[ipath][il] > max_logprob:
            max_logprob = logprobs[ipath][il]

if args.normalize_logprobs:
    for ipath in logprobs.keys():
        print ipath, min_logprobs[ipath], max_logprobs[ipath]
        for il in range(len(logprobs[ipath])):
            if logprobs[ipath][il] is not None:
                logprobs[ipath][il] = max_logprobs[ipath] / logprobs[ipath][il]
                # print logprobs[ipath][il]

    min_logprob = 0.
    max_logprob = 1.

fig = plt.figure(1)
fig.clf()
fig, ax = plt.subplots()
adj_mi_color = '#980000'
logprob_color = '#1947A3'

ax2 = ax.twinx()
nxbins = 8
nybins = 5

yextrafactor = 1.
if args.zoom:
    xmin = 340
    min_adjmi = 0.92
    min_logprob = -14400
    markersize = 32
else:
    xmin = 0
    min_adjmi = 0.
    markersize = 10

max_adjmi = 1.
if args.xbounds is not None:
    xmin = args.xbounds[0]
    max_length = args.xbounds[1]
if args.adjmi_bounds is not None:
    min_adjmi = args.adjmi_bounds[0]
    max_adjmi = args.adjmi_bounds[1]
if args.logprob_bounds is not None:
    min_logprob = args.logprob_bounds[0]
    max_logprob = args.logprob_bounds[1]

ax.set_xlim(xmin, max_length)
ax.set_ylim(min_adjmi, yextrafactor * max_adjmi)
ax2.set_ylim(min_logprob, yextrafactor * max_logprob)
ax.set_xlabel('agglomeration step', fontweight='bold')
ax.set_ylabel('adjusted MI', color=adj_mi_color, fontweight='bold')
ax2.set_ylabel('log prob', color=logprob_color, fontweight='bold')
fig.tight_layout()
plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.95)

ax.locator_params(nbins=nxbins, axis='x')
ax.locator_params(nbins=nybins, axis='y')

for ipath in logprobs.keys():
    steps = [i for i in range(len(logprobs[ipath]))]
    sizes = [markersize for i in range(len(logprobs[ipath]))]
    fig_logprob = ax2.plot(steps, logprobs[ipath], color=logprob_color, alpha=.5, linewidth=1)
    fig_adj_mi = ax.plot(steps, adj_mis[ipath], color=adj_mi_color, alpha=1, linewidth=1)
    if True:  # args.zoom:
        fig_logprob_sc = ax2.scatter(steps, logprobs[ipath], color=logprob_color, alpha=.5, s=sizes)
        fig_adj_mi_sc = ax.scatter(steps, adj_mis[ipath], color=adj_mi_color, alpha=1, s=sizes)

plotdir = os.getenv('www') + '/tmp'
if args.zoom:
    plotname = 'foo-zoom'
else:
    plotname = 'foo'
plt.savefig(plotdir + '/' + plotname + '.png')
check_call(['permissify-www', plotdir])
