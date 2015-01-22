#!/usr/bin/env python

import argparse

import sys
sys.path.append('python')

import plotting
import utils

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
parser.add_argument('--outdir', required=True)
parser.add_argument('--plotdirs', required=True)
parser.add_argument('--names', required=True)
parser.add_argument('--stats', default='')
parser.add_argument('--no-errors', action='store_true')
parser.add_argument('--plot-performance', action='store_true')
parser.add_argument('--scale-errors', type=float)
parser.add_argument('--rebin', type=int)
parser.add_argument('--colors')
parser.add_argument('--linestyles')

args = parser.parse_args()
args.plotdirs = utils.get_arg_list(args.plotdirs)
args.colors = utils.get_arg_list(args.colors, intify=True)
args.linestyles = utils.get_arg_list(args.linestyles, intify=True)
args.names = utils.get_arg_list(args.names)
for iname in range(len(args.names)):
    args.names[iname] = args.names[iname].replace('@', ' ')

assert len(args.plotdirs) == len(args.names)

plotting.compare_directories(args.outdir,
                             dirs = args.plotdirs,
                             names = args.names, stats=args.stats, errors=(not args.no_errors), scale_errors=args.scale_errors, rebin=args.rebin,
                             colors=args.colors, linestyles=args.linestyles, plot_performance=args.plot_performance)

# label = 'check-new-imgt'
# plotdir = '/var/www/sharing/dralph/partis/performance/'
# plotting.compare_directories(plotdir + '/' + label + '/igblast-vs-partis-vs-imgt',
#                              dirs = [plotdir + '/' + label + '/hmm/plots',
#                                      plotdir + '/' + label + '/sw/plots',
#                                      plotdir + 'igblast/' + label + '/plots',
#                                      plotdir + 'imgt/' + label + '/plots'],
#                              names = ['partis', 'sw', 'igblast', 'imgt'])

# plotdir = '/var/www/sharing/dralph/partis/performance'
# plotting.compare_directories(plotdir + '/adaptive-vs-vollmers',
#                              dirs = [plotdir + '/compare-to-vollmers/params/data/hmm_parameters/plots',
#                                      plotdir + '/vollmers/params/data/hmm_parameters/plots'],
#                              names = ['adaptive', 'vollmers'])
