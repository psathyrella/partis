#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
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
    'legend.fontsize': fsize,
    'axes.titlesize': fsize,
    # 'axes.labelsize': fsize,
    'xtick.labelsize': fsize,
    'ytick.labelsize': fsize,
    'axes.labelsize': fsize
})

values = {'obs':[], 'mean':[]}
with open('out.csv') as infile:
    for line in csv.DictReader(infile):
        values['obs'].append(float(line['xobs']))
        values['mean'].append(float(line['xmean']))

obsdf = DataFrame(values)
obsfig = obsdf.plot(color='b')  #, figsize=(6,4))
obsfig.grid(False)
plotdir = os.getenv('www') + '/tmp'
plt.savefig(plotdir + '/foo.png')
check_call(['permissify-www', plotdir])
