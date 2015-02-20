#!/usr/bin/env python
import csv
import glob
import argparse
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
mpl.rcParams['font.size'] = 48.0

parser = argparse.ArgumentParser()
parser.add_argument('--subdirs', default='v_3p_del:d_5p_del:d_3p_del:j_5p_del:vd_insertion:dj_insertion', help='Which variable categories?')
parser.add_argument('--dataset', choices=('adaptive', 'stanford', 'both'), default='adaptive')
args = parser.parse_args()
if args.subdirs == 'all':
    args.subdirs = all_subdirs
else:
    args.subdirs = utils.get_arg_list(args.subdirs)

modulo = '10'

assert os.getenv('www') is not None
baseplotdir = os.getenv('www') + '/partis'

fig, axes = plt.subplots(1, len(humans[args.dataset]), figsize=(5*len(humans[args.dataset]), 5), sharey=True,)

for indiv_n in range(len(humans[args.dataset])):
    indiv = humans[args.dataset][indiv_n]
    meanlist, variancelist = [], []
    for meanfname in [ baseplotdir + '/every-' + modulo + '-' + indiv + '/cf-subsets/' + sd + '/plots/means.csv' for sd in args.subdirs ]:
        if not os.path.exists(meanfname):
            raise Exception('ERROR means files d.n.e.: ' + meanfname)
        # print meanfname
        with open(meanfname, 'r') as meanfile:
            reader = csv.DictReader(meanfile)
            for line in reader:
                means = [ float(m) for m in line['means'].split(':') ]
                meanlist.append(numpy.mean(means))
                variancelist.append(numpy.var(means))
                # print ' ', meanlist[-1], variancelist[-1]
    ax = axes[indiv_n]
    df = pd.DataFrame({'mean': meanlist, 'variance': variancelist})
    df.plot(kind='hist', ax=ax, bins=50)
    ax.set_title(indiv)
    ax.set_ylabel('count')
    ax.grid(False)

sns.despine(fig=fig,offset=10, trim=True);

outdir = baseplotdir + '/every-' + modulo + '-' + args.dataset
os.makedirs(outdir)
fig.savefig(outdir + '/mean-variance-hist.svg')
check_call(['./bin/permissify-www', outdir])
