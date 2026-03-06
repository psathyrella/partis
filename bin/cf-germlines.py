#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import sys
import os
import copy
import collections
import colored_traceback.always

import partis.utils as utils
import partis.glutils as glutils

parser = argparse.ArgumentParser()
parser.add_argument('gldir1')
parser.add_argument('gldir2')
parser.add_argument('--names', default='+gl-1:+gl-2', help='colon-separated list of length 2 with labels for gldir1 and gldir2, which will be appended to each gene name in the ascii output')
parser.add_argument('--locus', default='igh')
args = parser.parse_args()
args.names = utils.get_arg_list(args.names)

# ----------------------------------------------------------------------------------------
def clrname(name):
    return utils.color('blue', name)

# ----------------------------------------------------------------------------------------
glfos = []
for name, gldir in zip(args.names, [args.gldir1, args.gldir2]):
    print('%s:' % clrname(name))
    glfos.append(glutils.read_glfo(gldir, args.locus, debug=True))

glutils.compare_glfos(glfos, args.names, args.locus)
