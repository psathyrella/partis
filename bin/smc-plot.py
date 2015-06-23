#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import argparse
import sys
import math
import numpy
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
parser.add_argument('--infnames')
parser.add_argument('--normalize-axes')
parser.add_argument('--only-last-step', action='store_true', default=True)
parser.add_argument('--is-data', action='store_true')
parser.add_argument('--xbounds')
parser.add_argument('--logprob-bounds')
parser.add_argument('--adjmi-bounds')
args = parser.parse_args()
args.infnames = get_arg_list(args.infnames)
if args.xbounds is not None:
    args.xbounds = get_arg_list(args.xbounds, floatify=True)
if args.logprob_bounds is not None:
    args.logprob_bounds = get_arg_list(args.logprob_bounds, floatify=True)
if args.adjmi_bounds is not None:
    args.adjmi_bounds = get_arg_list(args.adjmi_bounds, floatify=True)
if args.normalize_axes is not None:
    args.normalize_axes = get_arg_list(args.normalize_axes)

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
print 'need to add stuff to only use maximum when n_procs = 1'
max_logprob = {'logprobs' : [], 'adj_mis' : [], 'imaxes' : [], 'weights' : []}  # values at the point of maximum logprob
max_adj_mi = {'logprobs' : [], 'adj_mis' : [], 'imaxes' : [], 'weights' : []}  # values at the point of maximum adj_mi
imaxes = {'logprob' : None, 'adj_mi' : None}
def process_paths(logprobs, adj_mis, logweights):
    for ipath in logprobs.keys():
        imax_logprob = logprobs[ipath].index(max(logprobs[ipath]))
        max_logprob['logprobs'].append(logprobs[ipath][imax_logprob])
        max_logprob['adj_mis'].append(adj_mis[ipath][imax_logprob])
        max_logprob['imaxes'].append(imax_logprob)
        max_logprob['weights'].append(logweights[ipath][imax_logprob])

        imax_adj_mi = adj_mis[ipath].index(max(adj_mis[ipath]))
        max_adj_mi['logprobs'].append(logprobs[ipath][imax_adj_mi])
        max_adj_mi['adj_mis'].append(adj_mis[ipath][imax_adj_mi])
        max_adj_mi['imaxes'].append(imax_adj_mi)
        max_adj_mi['weights'].append(logweights[ipath][imax_adj_mi])

    def weighted_mean_rootvariance(vals, wgts):
        mean = numpy.average(vals, weights=wgts)
        variance = numpy.average((vals - mean)**2, weights=wgts)
        return (mean, math.sqrt(variance))

    imaxes['logprob'] = weighted_mean_rootvariance(max_logprob['imaxes'], max_logprob['weights'])
    imaxes['adj_mi'] = weighted_mean_rootvariance(max_adj_mi['imaxes'], max_adj_mi['weights'])
    delta_imaxes = [i_adj_mi - i_logprob for i_adj_mi, i_logprob in zip(max_adj_mi['imaxes'], max_logprob['imaxes'])]
    delta_wgts = [0.5*(adj_mi_wgt + logprob_wgt) for adj_mi_wgt, logprob_wgt in zip(max_adj_mi['weights'], max_logprob['weights'])]
    mean_delta_imax = weighted_mean_rootvariance(delta_imaxes, delta_wgts)

    max_logprob_means = {'logprobs' : weighted_mean_rootvariance(max_logprob['logprobs'], max_logprob['weights']),
                         'adj_mis' : weighted_mean_rootvariance(max_logprob['adj_mis'], max_logprob['weights'])}
    max_adj_mi_means = {'logprobs' : weighted_mean_rootvariance(max_adj_mi['logprobs'], max_adj_mi['weights']),
                        'adj_mis' : weighted_mean_rootvariance(max_adj_mi['adj_mis'], max_adj_mi['weights'])}
    print 'max logprob step %.1f:' % imaxes['logprob'][0]
    print '     logprob: %.1f +/- %.1e' % max_logprob_means['logprobs']
    print '      adj_mi: %.4f +/- %.1e' % max_logprob_means['adj_mis']
    print 'max adj_mi step %.1f:' % imaxes['adj_mi'][0]
    print '     logprob: %.1f +/- %.1e' % max_adj_mi_means['logprobs']
    print '      adj_mi: %.4f +/- %.1e' % max_adj_mi_means['adj_mis']
    print 'change: %.1f steps' % mean_delta_imax[0]
    logprob_at_adj_mi = max_adj_mi_means['logprobs'][0]
    logprob_at_logprob = max_logprob_means['logprobs'][0]
    print '     logprob: %.1f (delta %.3f)' % (logprob_at_adj_mi - logprob_at_logprob, (logprob_at_adj_mi - logprob_at_logprob) / logprob_at_logprob)
    print '      adj_mi: %.4f' % (max_adj_mi_means['adj_mis'][0] - max_logprob_means['adj_mis'][0])

if args.infnames is None:
    samples = (4000, )
    args.infnames = ['/fh/fast/matsen_e/dralph/work/partis-dev/_output/sample-sizes/' + str(n) + '-*-queries.csv' for n in samples]
    for ifn in range(len(args.infnames)):
        args.infnames[ifn] = glob.glob(args.infnames[ifn])[0]

logprobs, adj_mis, logweights = {}, {}, {}
final_logweights = []
for fname in args.infnames:
    print fname
    with open(fname) as infile:
        for line in csv.DictReader(infile):
            if args.only_last_step and int(line['n_procs']) != 1:  # skip any preliminary steps with more than one process
                continue
            ipath = int(line.get('path_index', 0))
            if ipath not in logprobs:
                logprobs[ipath], adj_mis[ipath], logweights[ipath] = [], [], []
            logprobs[ipath].append(float(line['logprob']))
            adj_mis[ipath].append(float(line.get('adj_mi', -1)))
            logweights[ipath].append(float(line.get('logweight', 1)))

process_paths(logprobs, adj_mis, logweights)

max_length = -1
for ipath in logprobs.keys():
    if len(logprobs[ipath]) > max_length:
        max_length = len(logprobs[ipath])

min_logprob, max_logprob = None, None
min_logprobs, max_logprobs = {}, {}  #[None for _ in range(len(logprobs.keys()))], [None for _ in range(len(logprobs.keys()))]
print '%d paths' % len(logprobs.keys())
for ipath in logprobs.keys():
    # while len(logprobs[ipath]) < max_length:
    #     logprobs[ipath].append(None)
    #     adj_mis[ipath].append(None)
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

if args.normalize_axes is not None and 'logprob' in args.normalize_axes:
    for ipath in logprobs.keys():
        # print ipath, min_logprobs[ipath], max_logprobs[ipath]
        for il in range(len(logprobs[ipath])):
            if logprobs[ipath][il] is not None:
                logprobs[ipath][il] = max_logprobs[ipath] / logprobs[ipath][il]
                # print logprobs[ipath][il]

    min_logprob = 0.
    max_logprob = 1.

# fig = plt.figure(1)
# fig.clf()
fig, ax = plt.subplots()
# sns.despine(fig=fig,offset=10, trim=True)
adj_mi_color = '#980000'
logprob_color = '#1947A3'

ax2 = ax.twinx()
nxbins = 8
nybins = 5

yextrafactor = 1.
xmin = 0
min_adj_mi = 0.
markersize = 30

max_adj_mi = 1.
if args.xbounds is None:
    args.xbounds = (0, 1 if 'step' in args.normalize_axes else max_length)
if args.adjmi_bounds is not None:
    min_adj_mi = args.adjmi_bounds[0]
    max_adj_mi = args.adjmi_bounds[1]
if args.logprob_bounds is not None:
    min_logprob = args.logprob_bounds[0]
    max_logprob = args.logprob_bounds[1]

ax.set_xlim(args.xbounds[0], args.xbounds[1])
ax.set_ylim(min_adj_mi, yextrafactor * max_adj_mi)
ax2.set_ylim(min_logprob, yextrafactor * max_logprob)
ax.set_xlabel('agglomeration step', fontweight='bold')
ax.set_ylabel('adjusted MI', color=adj_mi_color, fontweight='bold')
ax2.set_ylabel('log prob', color=logprob_color, fontweight='bold')
fig.tight_layout()
plt.gcf().subplots_adjust(bottom=0.16, left=0.2, right=0.78, top=0.95)

ax.locator_params(nbins=nxbins, axis='x')
ax.locator_params(nbins=nybins, axis='y')

for ipath in logprobs.keys():
    if 'step' in args.normalize_axes:
        steps = [float(i) / len(logprobs[ipath]) for i in range(len(logprobs[ipath]))]
    else:
        steps = [i for i in range(len(logprobs[ipath]))]
    sizes = [markersize for i in range(len(logprobs[ipath]))]
    fig_logprob = ax2.plot(steps, logprobs[ipath], color=logprob_color, alpha=.5, linewidth=1)
    if not args.is_data:
        fig_adj_mi = ax.plot(steps, adj_mis[ipath], color=adj_mi_color, alpha=1, linewidth=1)

    fig_logprob_sc = ax2.scatter(steps, logprobs[ipath], color=logprob_color, alpha=.5, s=sizes)
    if not args.is_data:
        fig_adj_mi_sc = ax.scatter(steps, adj_mis[ipath], color=adj_mi_color, alpha=1, s=sizes)


xlinepos_logprob = [imaxes['logprob'][0] for _ in range(2)]
if not args.is_data:
    xlinepos_adj_mi = [imaxes['adj_mi'][0] for _ in range(2)]
if 'step' in args.normalize_axes:
    assert False
    # xlinepos_logprob = [ float(x) / len(logprobs xlinepos_logprob
plt.plot(xlinepos_logprob, [min_logprob, max_logprob], color=logprob_color, linestyle='--')  #, color='k', linestyle='-', linewidth=2)
if not args.is_data:
    plt.plot(xlinepos_adj_mi, [min_adj_mi, max_adj_mi], color=adj_mi_color, linestyle='--')  #, color='k', linestyle='-', linewidth=2)

plotdir = os.getenv('www') + '/tmp'
plotname = 'foo'
plt.savefig(plotdir + '/' + plotname + '.png')
check_call(['permissify-www', plotdir])
