#!/usr/bin/env python
from subprocess import check_call
import argparse
import sys
import os
sys.path.insert(1, './python')

import utils
import plotting
from humans import humans, colors

dataset = 'adaptive'  #both  #stanford

modulo = '10'

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
# parser.add_argument('--write-means-over-subsets', action='store_true')
args = parser.parse_args()

# subdirs = [ e + '_del' for e in utils.real_erosions ] \
#           + [ i + '_insertion' for i in utils.boundaries]
subdirs = [ '.', ] \
          + [ e + '_del' for e in utils.real_erosions ] \
          + [ i + '_insertion' for i in utils.boundaries] \
          + [ 'mute-freqs', ] \
          + [ 'mute-freqs/' + r for r in utils.regions ]
# subdirs = [ 'mute-freqs/j', ]

webdir = '/var/www/sharing/dralph/partis'
subset = 0

for subdir in subdirs:
    print subdir, '-----------------'
    plotdirs, names, colorlist, linestyles, linewidths, markersizes, scale_errors, strings_to_ignore = [], [], [], [], [], [], [], []
    for human in humans[dataset]:
        print '  ', human
        baselabel = 'every-' + modulo + '-' + human
        datadir = webdir + '/' + baselabel + '/cf-subsets/' + subdir
        simudir = webdir + '/' + baselabel + '-subset-' + str(subset) + '/params/simu/hmm/true/' + subdir
        # print datadir
        # print simudir
        # sys.exit()
        plotdirs += [datadir, simudir]
        names += [ human + '@data', human + '@simu' ]
        colorlist += [colors[human], colors[human]]
        linestyles += ['1', '2']
        linewidths += ['1', '2']
        markersizes += ['1', '0']
        scale_errors += ['1.414', '2.24']
        strings_to_ignore += ['_mean-bins', '']

    cmd = './python/compare.py --dont-calculate-mean-info --colors ' + ':'.join(colorlist)
    cmd += ' --linestyles ' + ':'.join(linestyles)
    cmd += ' --plotdirs ' + ':'.join(plotdirs)
    cmd += ' --names ' + ':'.join(names)
    cmd += ' --outdir ' + webdir + '/data-vs-simu-' + dataset + '/' + subdir
    cmd += ' --strings-to-ignore ' + ':'.join(strings_to_ignore)
    cmd += ' --scale-errors ' + ':'.join(scale_errors)
    # cmd += ' --no-errors'

    if dataset == 'both':
        cmd += ' --linewidth 2'  # --linestyles ' + ':'.join()
    if subdir == 'mute-freqs':
        # cmd += ' --no-errors --rebin 5'
        cmd += '  --rebin 5'
    if 'mute-freqs/v' not in subdir and 'mute-freqs/d' not in subdir and 'mute-freqs/j' not in subdir:
        cmd += ' --normalize'
    if 'mute-freqs/v' in subdir or 'mute-freqs/d' in subdir or 'mute-freqs/j' in subdir:
        markersizes = ' '.join(markersizes).replace('1', '0').split()  # switch data marker sizes to zero
    # else:
    #     cmd += ' --markersize 0'

    cmd += ' --graphify --linewidths ' + ':'.join(linewidths) + ' --markersizes ' + ':'.join(markersizes)
    check_call(cmd.split(' '))
