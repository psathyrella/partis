#!/usr/bin/env python
import argparse
import json
import csv
import os
import sys
current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '/python')
if not os.path.exists(current_script_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % current_script_dir
sys.path.insert(1, current_script_dir)

import plotting
import utils
import glutils
from opener import opener

parser = argparse.ArgumentParser()
parser.add_argument('--outdir', required=True)
parser.add_argument('--plotdirs', required=True)
parser.add_argument('--names', required=True)
parser.add_argument('--plot-performance', action='store_true')
parser.add_argument('--rebin', type=int)
parser.add_argument('--colors', default='#006600:#cc0000:#990012:#3333ff:#3399ff:#2b65ec:#2b65ec:#808080')
parser.add_argument('--linestyles')
parser.add_argument('--datadir', default='data/imgt')
parser.add_argument('--leaves-per-tree')
parser.add_argument('--linewidths')
parser.add_argument('--alphas')
parser.add_argument('--markersizes')
parser.add_argument('--normalize', action='store_true')
parser.add_argument('--only-csv-plots', action='store_true')
parser.add_argument('--strings-to-ignore')  # remove this string from the plot names in each dir (e.g. '-mean-bins') NOTE replaces '_' with '-'

args = parser.parse_args()
if args.strings_to_ignore is not None:
    args.strings_to_ignore = args.strings_to_ignore.replace('_', '-')
args.plotdirs = utils.get_arg_list(args.plotdirs)
args.colors = utils.get_arg_list(args.colors)
args.linestyles = utils.get_arg_list(args.linestyles, intify=True, translation={1 : '-',2 : '--'})
args.names = utils.get_arg_list(args.names)
args.leaves_per_tree = utils.get_arg_list(args.leaves_per_tree, intify=True)
args.strings_to_ignore = utils.get_arg_list(args.strings_to_ignore)
args.markersizes = utils.get_arg_list(args.markersizes, intify=True)
args.linewidths = utils.get_arg_list(args.linewidths, intify=True)
args.alphas = utils.get_arg_list(args.alphas, floatify=True)
for iname in range(len(args.names)):
    args.names[iname] = args.names[iname].replace('@', ' ')

assert len(args.plotdirs) == len(args.names)

# glfo = glutils.read_glfo(args.datadir)
# args.cyst_positions, args.tryp_positions = glfo['cyst-positions'], glfo['tryp-positions']
args.tryp_positions, args.tryp_positions = None, None

plotting.compare_directories(args)
