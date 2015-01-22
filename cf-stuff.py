#!/usr/bin/env python
from subprocess import check_call
import sys

h_vs_h = True  # human vs human, or data vs simu?

humans = ['A', 'B', 'C']
labels = ['mimic-' + h for h in humans]
basedir = '/var/www/sharing/dralph/partis'  #/$label/params
subdir = 'mute-freqs/j'
if h_vs_h:
    outlabel = 'cf-adaptive-data'
else:
    outlabel = 'data-vs-true-simu'

plotdirs = []
plotdirs.append('data/hmm_parameters')
plotdirs.append('data/hmm_parameters')
plotdirs.append('data/hmm_parameters')
if not h_vs_h:
    plotdirs.insert(1, 'simu/hmm_parameters/true')
    plotdirs.insert(3, 'simu/hmm_parameters/true')
    plotdirs.insert(5, 'simu/hmm_parameters/true')

if h_vs_h:
    plotdirs = [ basedir + '/' + labels[ipd] + '/params/' + plotdirs[ipd] + '/' + subdir for ipd in range(len(plotdirs))]
else:
    plotdirs = [ basedir + '/' + labels[ipd/2] + '/params/' + plotdirs[ipd] + '/' + subdir for ipd in range(len(plotdirs))]

if h_vs_h:
    names = [ humans[i%3] for i in range(len(plotdirs)) ]
else:
    names = []
    for i in range(len(plotdirs)):
        if i%2:
            names.append('simu@')
        else:
            names.append('data@')
    names = [ names[i] + humans[i/2] for i in range(len(plotdirs)) ]

# data vs simu:
cmd = './python/compare.py --plotdirs ' + ':'.join(plotdirs) + ' --names ' + ':'.join(names) + ' --outdir /var/www/sharing/dralph/partis/' + labels[0] + '/' + outlabel + '/' + subdir
if 'mute-freqs' not in subdir:
    cmd += ' --scale-errors 2.24'
if h_vs_h:
    cmd += ' --colors 632:596:418'
else:
    cmd += ' --colors 632:596:418:632:596:418 --linestyles 1:1:1:2:2:2'
# + ' --rebin 3'

check_call(cmd.split(' '))
