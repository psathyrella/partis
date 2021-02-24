#!/usr/bin/env python
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import argparse
import colored_traceback.always

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
from clusterpath import ClusterPath

import gex
feature_matrix_fname = "/fh/fast/matsen_e/data/goo-dengue-10x/data/gex/filtered_feature_bc_matrix"
outdir = "/home/dralph/Dropbox/tmp-gex" #"/fh/fast/matsen_e/dralph/partis/tmp/gex"
# gex.run_msigdbr('datascripts/meta/goo-dengue-10x')
for mname in ['waick', 'msigdb']: # ['hvg', 'fabio', 'waick', 'msigdb']:
    gex.run_gex(feature_matrix_fname, mname, outdir + '/' + mname)
    # gex.read_gex(outdir + '/' + mname)
