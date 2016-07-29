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
parser.add_argument('--colors', default='#006600:#cc0000:#990012:#3333ff:#3399ff:#2b65ec:#2b65ec:#808080')
parser.add_argument('--linewidths', default='5:3:2:2:2')
parser.add_argument('--gldir', default='data/imgt')
parser.add_argument('--chain', default='h')
parser.add_argument('--normalize', action='store_true')

args = parser.parse_args()
args.plotdirs = utils.get_arg_list(args.plotdirs)
args.names = utils.get_arg_list(args.names)
args.colors = utils.get_arg_list(args.colors)
args.linewidths = utils.get_arg_list(args.linewidths)
for iname in range(len(args.names)):
    args.names[iname] = args.names[iname].replace('@', ' ')

if len(args.plotdirs) != len(args.names):
    raise Exception('poorly formatted args:\n  %s\n  %s' % (' '.join(args.plotdirs), ' '.join(args.names)))

args.glfo = glutils.read_glfo(args.gldir, args.chain)

plotting.compare_directories(args)
