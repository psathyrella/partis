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
for subj in ['JL029-001-03', 'JL029-002-02', 'JL029-012-03', 'd14']:
    feature_matrix_fname = "/fh/fast/matsen_e/data/goo-dengue-10x/data/%s/filtered_feature_bc_matrix" % subj
    outdir = '/fh/fast/matsen_e/data/goo-dengue-10x/gex-output/%s' % subj
    assert False  # need to double check this <outdir> is right
    # gex.run_msigdbr('datascripts/meta/goo-dengue-10x')
    for mname in ['waick', 'fabio']: # ['hvg', 'fabio', 'waick', 'msigdb']:
        gex.run_gex(feature_matrix_fname, mname, outdir + '/' + mname)
        # gex.read_gex(outdir + '/' + mname)
