#!/usr/bin/env python
from subprocess import check_call
import argparse
import sys
import os
sys.path.insert(1, './python')

import utils
import plotting
from humans import humans, colors, all_subdirs

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
outlabel = 'cf-subsets'
webdir = os.getenv('www') + '/partis' if os.getenv('www') is not None else '_plots'
subsets = [ str(i) for i in range(int(modulo)) ]

for subdir in args.subdirs:
    if subdir not in all_subdirs:
        raise Exception('ERROR bad subdir: ' + subdir)
    print subdir, '-----------------'
    final_plotdirs = []
    for human in humans[args.dataset]:  # loop over each human, making mean/variance over subset plots
        print '  ', human
        baselabel = 'every-' + modulo + '-' + human
        baseoutdir = webdir + '/' + baselabel + '/' + outlabel
        labels = [ baselabel + '-subset-' + s for s in subsets]

        plotdirs = ['data/hmm' for _ in range(len(subsets)) ]
        plotdirs = [ webdir + '/' + labels[ipd] + '/params/' + plotdirs[ipd] + '/' + subdir for ipd in range(len(plotdirs))]
        
        names = [ human + '@' + s for s in subsets ]

        cmd = './bin/compare.py --no-errors --linewidth 1 --plotdirs ' + ':'.join(plotdirs) + ' --names ' + ':'.join(names) + ' --outdir ' + baseoutdir + '/' + subdir
        final_plotdirs.append(baseoutdir + '/' + subdir)
        
        colorlist = [ colors[human] for _ in range(len(subsets)) ]
        cmd += ' --colors ' + ':'.join(colorlist)

        if 'mute-freqs/v' not in subdir and 'mute-freqs/d' not in subdir and 'mute-freqs/j' not in subdir:
            cmd += ' --normalize'

        check_call(cmd.split(' '))
        # print 'SKIPPING'

    # then make the human-comparison plots
    final_cmd = './bin/compare.py --colors ' + ':'.join([ colors[human] for human in humans[args.dataset] ]) \
                + ' --plotdirs ' + ':'.join(final_plotdirs) + ' --names ' + ':'.join(humans[args.dataset]) + ' --outdir ' + webdir + '/cf-data-' + args.dataset + '/' + subdir
    final_cmd += ' --graphify --linewidth 1'
    final_cmd += ' --scale-errors 1.414'
    if subdir == 'mute-freqs':
        final_cmd += ' --rebin 5'
    if 'mute-freqs/v' not in subdir and 'mute-freqs/d' not in subdir and 'mute-freqs/j' not in subdir:
        final_cmd += ' --normalize'
    if 'mute-freqs/v' in subdir or 'mute-freqs/d' in subdir:
        final_cmd += ' --markersize 1 --linewidth 1'

    check_call(final_cmd.split(' '))
