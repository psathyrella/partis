#!/usr/bin/env python
from subprocess import check_call
import sys
import os
sys.path.insert(1, './python')

import utils
import plotting

h_vs_h = True  # human vs human, or data vs simu?
dataset = 'both'

if dataset == 'stanford':
    humans = ['021-019', '021-044', '021-048', '021-050', '021-055', '021-057', '021-059', '021-060', '021-061', '021-063', '021-068', '021-071', '021-084', '021-018']
elif dataset == 'adaptive':
    humans = ['A', 'B', 'C']
    # humans = ['A', 'B', 'C', '021-019', '021-044', '021-048', '021-050', '021-055', '021-057']
    # humans = ['021-019', '021-044', '021-048', '021-050', '021-055', '021-057']
elif dataset == 'both':
    humans = ['A', 'B', 'C', '021-019', '021-044', '021-048', '021-050', '021-055', '021-057', '021-059', '021-060', '021-061', '021-063', '021-068', '021-071', '021-084', '021-018']
else:
    assert False

labels = ['mimic-' + h for h in humans]
basedir = '/var/www/sharing/dralph/partis'  #/$label/params
subdirs = [ e + '_del' for e in utils.real_erosions ] + [ i + '_insertion' for i in utils.boundaries]
# subdirs = [ 'dj_insertion', ]
if h_vs_h:
    outlabel = 'cf-' + dataset + '-data'
else:
    outlabel = 'data-vs-true-simu'

baseoutdir = '/var/www/sharing/dralph/partis/' + labels[0] + '/' + outlabel
plotting.make_mean_plots(baseoutdir, subdirs,  baseoutdir + '/hexmean')
sys.exit()

for subdir in subdirs:
    print subdir
    plotdirs = ['data/hmm_parameters' for _ in range(len(humans)) ]
    if not h_vs_h:
        assert not stanford  # would need to update for stanford humans
        plotdirs.insert(1, 'simu/hmm_parameters/true')
        plotdirs.insert(3, 'simu/hmm_parameters/true')
        plotdirs.insert(5, 'simu/hmm_parameters/true')

    if h_vs_h:
        plotdirs = [ basedir + '/' + labels[ipd] + '/params/' + plotdirs[ipd] + '/' + subdir for ipd in range(len(plotdirs))]
    else:
        plotdirs = [ basedir + '/' + labels[ipd/2] + '/params/' + plotdirs[ipd] + '/' + subdir for ipd in range(len(plotdirs))]
    
    if h_vs_h:
        names = [ humans[i] for i in range(len(plotdirs)) ]
    else:
        names = []
        for i in range(len(plotdirs)):
            if i%2:
                names.append('simu@')
            else:
                names.append('data@')
        names = [ names[i] + humans[i/2] for i in range(len(plotdirs)) ]
    
    # data vs simu:
    cmd = './python/compare.py --plotdirs ' + ':'.join(plotdirs) + ' --names ' + ':'.join(names) + ' --outdir ' + baseoutdir + '/' + subdir
    # if 'mute-freqs' not in subdir:
    #     print 'NOTE need to decide what to use for final uncertainties'
    #     cmd += ' --scale-errors 2.24'
    leaves_per_trees = [ '5' if 'simu' in pd else '1' for pd in plotdirs ]
    cmd += ' --leaves-per-tree ' + ':'.join(leaves_per_trees)
    
    colorlist = ['632', '596', '418', '632', '596', '418', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1']
    cmd += ' --colors ' + ':'.join(colorlist)
    if not h_vs_h:
        cmd += ' --linestyles 1:1:1:2:2:2'
    # + ' --rebin 3'

    check_call(cmd.split(' '))
