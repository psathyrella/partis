#!/usr/bin/env python
from subprocess import check_call
import sys
import os
sys.path.insert(1, './python')

import utils
import plotting

basedir = '/var/www/sharing/dralph/partis'

h_vs_h = True  # human vs human, or data vs simu?
dataset = 'adaptive'

if dataset == 'stanford':
    humans = ['021-019', '021-044', '021-048', '021-050', '021-055', '021-057', '021-059', '021-060', '021-061', '021-063', '021-068', '021-071', '021-084', '021-018']
elif dataset == 'adaptive':
    # humans = ['A', 'B', 'C']
    humans = ['C',]
elif dataset == 'both':
    humans = ['A', 'B', 'C', '021-019', '021-044', '021-048', '021-050', '021-055', '021-057', '021-059', '021-060', '021-061', '021-063', '021-068', '021-071', '021-084', '021-018']
else:
    assert False

subsets = [ str(i) for i in range(8) ]

labels = ['every-10-' + h + '-subset-' + s for h in humans for s in subsets]

# subdirs = [ e + '_del' for e in utils.real_erosions ] + [ i + '_insertion' for i in utils.boundaries] + [ '.', ]
subdirs = [ '.', ]

if h_vs_h:
    # outlabel = 'cf-' + dataset + '-data'
    outlabel = 'cf-subsets'
else:
    outlabel = 'data-vs-true-simu'

baseoutdir = '/var/www/sharing/dralph/partis/' + labels[0] + '/' + outlabel
# plotting.make_mean_plots(baseoutdir, subdirs,  baseoutdir + '/hexmean')
# sys.exit()

for subdir in subdirs:
    print subdir
    plotdirs = ['data/hmm' for _ in range(len(humans)) for _ in range(len(subsets)) ]
    if not h_vs_h:
        assert not stanford  # would need to update for stanford humans
        assert subsets is None
        plotdirs.insert(1, 'simu/hmm/true')
        plotdirs.insert(3, 'simu/hmm/true')
        plotdirs.insert(5, 'simu/hmm/true')

    if h_vs_h:
        plotdirs = [ basedir + '/' + labels[ipd] + '/params/' + plotdirs[ipd] + '/' + subdir for ipd in range(len(plotdirs))]
    else:
        plotdirs = [ basedir + '/' + labels[ipd/2] + '/params/' + plotdirs[ipd] + '/' + subdir for ipd in range(len(plotdirs))]
    
    if h_vs_h:
        names = [ h + '@' + s for h in humans for s in subsets ]
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
    
    # colorlist = ['632', '596', '418', '632', '596', '418', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1']
    colors = ('632', '596', '418')
    colorlist = [ colors[i] for i in range(len(humans)) for _ in range(len(subsets)) ]
    cmd += ' --colors ' + ':'.join(colorlist)
    if not h_vs_h:
        cmd += ' --linestyles 1:1:1:2:2:2'
    # + ' --rebin 3'

    check_call(cmd.split(' '))
