#!/usr/bin/env python
from subprocess import check_call
import argparse
import sys
import os
sys.path.insert(1, './python')

import utils
import plotting
from humans import humans, colors

dataset = 'adaptive'  #stanford  #both

modulo = '10'

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
# parser.add_argument('--write-means-over-subsets', action='store_true')
args = parser.parse_args()

# subdirs = [ e + '_del' for e in utils.real_erosions ] \
#           + [ i + '_insertion' for i in utils.boundaries]
# subdirs = [ '.', ] \
#           + [ e + '_del' for e in utils.real_erosions ] \
#           + [ i + '_insertion' for i in utils.boundaries] \
#           + [ 'mute-freqs', ] \
#           + [ 'mute-freqs/' + r for r in utils.regions ]
subdirs = [ 'mute-freqs', ]

webdir = '/var/www/sharing/dralph/partis'
subset = 0

for subdir in subdirs:
    print subdir, '-----------------'
    plotdirs, names, colorlist, linestyles = [], [], [], []
    for human in humans[dataset]:
        if human != 'A':
            continue
        print '  ', human
        baselabel = 'every-' + modulo + '-' + human
        datadir = webdir + '/' + baselabel + '/cf-subsets/' + subdir
        simudir = webdir + '/' + baselabel + '-subset-' + str(subset) + '/params/simu/hmm/true/' + subdir
        plotdirs += [datadir, simudir]
        names += [ human + '@data', human + '@simu' ]
        colorlist += [colors[human], colors[human]]
        linestyles += ['1', '2']

    cmd = './python/compare.py --dont-calculate-mean-info --colors ' + ':'.join(colorlist)
    cmd += ' --linestyles ' + ':'.join(linestyles)
    cmd += ' --plotdirs ' + ':'.join(plotdirs)
    cmd += ' --names ' + ':'.join(names)
    cmd += ' --outdir ' + webdir + '/data-vs-simu-' + dataset + '/' + subdir
    cmd += ' --strings-to-ignore _mean-bins:'
    cmd += ' --leaves-per-tree 1:5'

    if dataset == 'both':
        cmd += ' --linewidth 2'  # --linestyles ' + ':'.join()
    if subdir == 'mute-freqs':
        cmd += ' --no-errors --graphify --rebin 5'
    if 'mute-freqs/v' not in subdir and 'mute-freqs/d' not in subdir and 'mute-freqs/j' not in subdir:
        cmd += ' --normalize'
    if 'mute-freqs/v' in subdir or 'mute-freqs/d' in subdir:
        cmd += ' --markersize 1 --linewidth 1'
    else:
        cmd += ' --markersize 0'
    check_call(cmd.split(' '))
