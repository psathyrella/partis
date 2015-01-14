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

args = parser.parse_args()
args.plotdirs = utils.get_arg_list(args.plotdirs)
args.names = utils.get_arg_list(args.names)

assert len(args.plotdirs) == len(args.names)
assert len(args.names) < 5  # need to change plotting function to allow more

plotting.compare_directories(args.outdir,
                             dirs = args.plotdirs,
                             names = args.names, stats=args.stats, errors=(not args.no_errors))

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
