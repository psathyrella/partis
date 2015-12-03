#!/usr/bin/env python
from subprocess import check_call
import argparse
import sys
import os
sys.path.insert(1, './python')

import utils
import plotting
from humans import humans, colors, all_subdirs
raise Exception('needs to be tested since root removal')

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
parser.add_argument('--subdirs', default='all', help='Which variable categories?')
parser.add_argument('--dataset', choices=('adaptive', 'stanford', 'both'), default='adaptive')
args = parser.parse_args()
if args.subdirs == 'all':
    args.subdirs = all_subdirs
else:
    args.subdirs = utils.get_arg_list(args.subdirs)

modulo = '10'
webdir = '/var/www/sharing/dralph/partis'
subset = 0

for subdir in args.subdirs:
    print subdir, '-----------------'
    if subdir not in all_subdirs:
        raise Exception('ERROR bad subdir: ' + str(subdir))
    plotdirs, names, colorlist, linestyles, linewidths, markersizes, scale_errors, strings_to_ignore = [], [], [], [], [], [], [], []
    for human in humans[args.dataset]:
        print '  ', human
        baselabel = 'every-' + modulo + '-' + human
        datadir = webdir + '/' + baselabel + '/cf-subsets/' + subdir
        simudir = webdir + '/' + baselabel + '-subset-' + str(subset) + '/params/simu/hmm/true/' + subdir
        plotdirs += [datadir, simudir]
        names += [ human + '@data', human + '@simu' ]
        colorlist += [colors[human], colors[human]]
        linestyles += ['1', '2']
        linewidths += ['1', '2']
        markersizes += ['1', '0']
        scale_errors += ['1.414', '2.24']
        strings_to_ignore += ['_mean-bins', '']

    cmd = './bin/compare.py --colors ' + ':'.join(colorlist)
    cmd += ' --linestyles ' + ':'.join(linestyles)
    cmd += ' --plotdirs ' + ':'.join(plotdirs)
    cmd += ' --names ' + ':'.join(names)
    cmd += ' --outdir ' + webdir + '/data-vs-simu-' + args.dataset + '/' + subdir
    cmd += ' --strings-to-ignore ' + ':'.join(strings_to_ignore)
    cmd += ' --scale-errors ' + ':'.join(scale_errors)

    if subdir == 'mute-freqs':
        cmd += ' --rebin 5 --no-errors'
    if 'mute-freqs/v' not in subdir and 'mute-freqs/d' not in subdir and 'mute-freqs/j' not in subdir:
        cmd += ' --normalize'
    if 'mute-freqs/v' in subdir or 'mute-freqs/d' in subdir or 'mute-freqs/j' in subdir:
        markersizes = ' '.join(markersizes).replace('1', '0').split()  # switch data marker sizes to zero

    cmd += ' --linewidths ' + ':'.join(linewidths) + ' --markersizes ' + ':'.join(markersizes)
    check_call(cmd.split(' '))
