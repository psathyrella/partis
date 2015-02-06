#!/usr/bin/env python
from subprocess import check_call
import argparse
import sys
import os
sys.path.insert(1, './python')

import utils
import plotting
from humans import humans, colors

dataset = 'both'  #stanford  #adaptive

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

outlabel = 'cf-subsets'
webdir = '/var/www/sharing/dralph/partis'
subsets = [ str(i) for i in range(int(modulo)) ]

#('632', '596', '418')
colors = {'A':'595', 'B':'807', 'C':'834',
          '021-019': '1', '021-044': '1', '021-048': '1', '021-050': '1', '021-055': '1', '021-057': '1', '021-059': '1', '021-060': '1', '021-061': '1', '021-063': '1', '021-068': '1', '021-071': '1', '021-084': '1', '021-018': '1'}

for subdir in subdirs:
    print subdir, '-----------------'
    final_plotdirs = []
    for human in humans[dataset]:
        print '  ', human
        baselabel = 'every-' + modulo + '-' + human
        baseoutdir = webdir + '/' + baselabel + '/' + outlabel
        labels = [ baselabel + '-subset-' + s for s in subsets]

        # plotting.make_mean_plots(baseoutdir , subdirs,  baseoutdir + '/hexmean')
        # sys.exit()
    
        plotdirs = ['data/hmm' for _ in range(len(subsets)) ]
        plotdirs = [ webdir + '/' + labels[ipd] + '/params/' + plotdirs[ipd] + '/' + subdir for ipd in range(len(plotdirs))]
        # for p in plotdirs:
        #     print p
        
        names = [ human + '@' + s for s in subsets ]
        # for n in names:
        #     print n

        cmd = './python/compare.py --no-errors --linewidth 1 --plotdirs ' + ':'.join(plotdirs) + ' --names ' + ':'.join(names) + ' --outdir ' + baseoutdir + '/' + subdir
        final_plotdirs.append(baseoutdir + '/' + subdir)
        # leaves_per_trees = [ '5' if 'simu' in pd else '1' for pd in plotdirs ]
        # cmd += ' --leaves-per-tree ' + ':'.join(leaves_per_trees)
        
        colorlist = [ colors[human] for _ in range(len(subsets)) ]
        cmd += ' --colors ' + ':'.join(colorlist)

        if 'mute-freqs/v' not in subdir and 'mute-freqs/d' not in subdir and 'mute-freqs/j' not in subdir:
            cmd += ' --normalize'

        # check_call(cmd.split(' '))

    print '  YOOOOO skipping subset plots'

    final_cmd = './python/compare.py --dont-calculate-mean-info --colors ' + ':'.join([ colors[human] for human in humans[dataset] ]) + ' --plotdirs ' + ':'.join(final_plotdirs) + ' --names ' + ':'.join(humans[dataset]) + ' --outdir ' + webdir + '/cf-data-' + dataset + '/' + subdir
    if dataset == 'both':
        final_cmd += ' --linewidth 2'  # --linestyles ' + ':'.join()
    if 'mute-freqs/v' not in subdir and 'mute-freqs/d' not in subdir and 'mute-freqs/j' not in subdir:
        final_cmd += ' --normalize'
    if 'mute-freqs/v' in subdir or 'mute-freqs/d' in subdir:
        final_cmd += ' --markersize 1 --linewidth 1'
    else:
        final_cmd += ' --markersize 0'
    check_call(final_cmd.split(' '))
