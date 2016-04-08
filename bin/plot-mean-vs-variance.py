#!/usr/bin/env python
import csv
import glob
import argparse
import math
import os
from subprocess import check_call
import numpy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
sys.path.insert(1, './python')

import utils
from humans import humans, colors, all_subdirs

sns.set_style("ticks")
parser = argparse.ArgumentParser()
parser.add_argument('--subdirs', default='v_3p_del:d_5p_del:d_3p_del:j_5p_del:vd_insertion:dj_insertion', help='Which variable categories?')
parser.add_argument('--dataset', choices=('adaptive', 'vollmers', 'both'), default='adaptive')
args = parser.parse_args()
if args.subdirs == 'all':
    args.subdirs = all_subdirs
else:
    args.subdirs = utils.get_arg_list(args.subdirs)

modulo = '10'

assert os.getenv('www') is not None
baseplotdir = os.getenv('www') + '/partis'

max_cols = 3
n_humans = len(humans[args.dataset])
n_cols = min(max_cols, n_humans)
n_rows = max(1, int(math.ceil(float(n_humans)/n_cols)))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 3*n_rows), sharey=True, sharex=True)
plt.locator_params(nbins=3)

fsize = 25
fdict = {'fontsize': fsize}
mpl.rcParams.update({
    # 'font.size': fsize,
    'legend.fontsize': fsize,
    'axes.titlesize': fsize,
    # 'axes.labelsize': fsize,
    'xtick.labelsize': fsize,
    'ytick.labelsize': fsize,
    'axes.labelsize': fsize
})

print n_cols, n_rows
for indiv_n in range(len(humans[args.dataset])):
    indiv = humans[args.dataset][indiv_n]
    meanlist, variancelist = [], []
    for meanfname in [ baseplotdir + '/every-' + modulo + '-' + indiv + '/cf-subsets/' + sd + '/plots/means.csv' for sd in args.subdirs ]:
        if not os.path.exists(meanfname):
            raise Exception('ERROR means files d.n.e.: ' + meanfname)
        with open(meanfname, 'r') as meanfile:
            reader = csv.DictReader(meanfile)
            for line in reader:
                means = [ float(m) for m in line['means'].split(':') ]
                meanlist.append(numpy.mean(means))
                variancelist.append(numpy.var(means))
    irow = indiv_n / n_cols
    icol = indiv_n % n_cols
    # print indiv_n, irow, '/', len(axes),
    # print '   ', icol, '/', len(axes[irow])
    if n_rows > 1:
        ax = axes[irow][icol]
    else:
        ax = axes[icol]
    df = pd.DataFrame({'mean': meanlist, 'variance': variancelist})
    df.plot(kind='hist', ax=ax, bins=50, title=indiv)
    # ax.set_title(indiv)
    ax.set_ylabel('count', fontdict=fdict)  #, labelpad=0)
    if irow != n_rows-1:
        ax.get_xaxis().set_visible(False)
    if icol != 0:
        ax.get_yaxis().set_visible(False)
    if irow != 0 or icol != n_cols-1:
        ax.legend().set_visible(False)
    ax.grid(False)
    ax.set_xlim([0, 10])

fig.tight_layout()
if n_rows > 1:
    fig.subplots_adjust(left=0.1, bottom=0.1, wspace=0.25)
else:
    fig.subplots_adjust(left=0.1, bottom=0.2, wspace=0.25)
# fig.xlim(0,15)
sns.despine(fig=fig, offset=5, trim=True);

outdir = baseplotdir + '/every-' + modulo + '-' + args.dataset
outfname = outdir + '/mean-variance-hist.svg'
if not os.path.exists(outdir):
    os.makedirs(outdir)
if os.path.exists(outfname):
    os.remove(outfname)
fig.savefig(outfname)
check_call(['./bin/permissify-www', outdir])
