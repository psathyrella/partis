#!/usr/bin/env python
from subprocess import check_call
import argparse
import sys
import os
sys.path.insert(1, './python')

import utils
import plotting

basedir = '/var/www/sharing/dralph/partis'

dataset = 'adaptive'

if dataset == 'stanford':
    humans = ['021-019', '021-044', '021-048', '021-050', '021-055', '021-057', '021-059', '021-060', '021-061', '021-063', '021-068', '021-071', '021-084', '021-018']
    modulo = '5'
elif dataset == 'adaptive':
    humans = ['A', 'B', 'C']
    modulo = '10'
elif dataset == 'both':
    humans = ['A', 'B', 'C', '021-019', '021-044', '021-048', '021-050', '021-055', '021-057', '021-059', '021-060', '021-061', '021-063', '021-068', '021-071', '021-084', '021-018']
else:
    assert False


parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
# parser.add_argument('--write-means-over-subsets', action='store_true')
args = parser.parse_args()

subdirs = [ e + '_del' for e in utils.real_erosions ] \
          + [ i + '_insertion' for i in utils.boundaries]
# subdirs = [ '.', ] \
#           + [ e + '_del' for e in utils.real_erosions ] \
#           + [ i + '_insertion' for i in utils.boundaries] \
#           + [ 'mute-freqs', ] \
#           + [ 'mute-freqs/' + r for r in utils.regions ]
# subdirs = [ 'd_3p_del', ]

outlabel = 'cf-subsets'
webdir = '/var/www/sharing/dralph/partis'
subsets = [ str(i) for i in range(int(modulo)-1) ]
print '\nNOT USING ALL SUBSETS\n'

colors = {'A':'595', 'B':'807', 'C':'834'}  #('632', '596', '418')

for subdir in subdirs:
    print subdir, '-----------------'
    final_plotdirs = []
    for human in humans:
        baselabel = 'every-' + modulo + '-' + human
        baseoutdir = webdir + '/' + baselabel + '/' + outlabel
        labels = [ baselabel + '-subset-' + s for s in subsets]

        plotting.make_mean_plots(baseoutdir , subdirs,  baseoutdir + '/hexmean')
        sys.exit()
    
        plotdirs = ['data/hmm' for _ in range(len(subsets)) ]
        plotdirs = [ basedir + '/' + labels[ipd] + '/params/' + plotdirs[ipd] + '/' + subdir for ipd in range(len(plotdirs))]
        
        names = [ human + '@' + s for s in subsets ]

        cmd = './python/compare.py --no-errors --linewidth 1 --plotdirs ' + ':'.join(plotdirs) + ' --names ' + ':'.join(names) + ' --outdir ' + baseoutdir + '/' + subdir
        final_plotdirs.append(baseoutdir + '/' + subdir)
        # leaves_per_trees = [ '5' if 'simu' in pd else '1' for pd in plotdirs ]
        # cmd += ' --leaves-per-tree ' + ':'.join(leaves_per_trees)
        
        colorlist = [ colors[human] for _ in range(len(subsets)) ]
        cmd += ' --colors ' + ':'.join(colorlist)
        check_call(cmd.split(' '))

    final_cmd = './python/compare.py --dont-calculate-mean-info --colors ' + ':'.join(colors.values()) + ' --plotdirs ' + ':'.join(final_plotdirs) + ' --names ' + ':'.join(humans) + ' --outdir ' + webdir + '/cf-data/' + subdir
    check_call(final_cmd.split(' '))
