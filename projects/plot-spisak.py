#!/usr/bin/env python
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import argparse
import colored_traceback.always
import glob
import itertools
import json
import collections
import numpy
import math

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/projects', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import glutils
import hutils

# ----------------------------------------------------------------------------------------
def getfam(line):
    if 'true_FAMILYs' in line:
        return line['true_FAMILYs']
    assert len(line['unique_ids']) == 1
    return line['unique_ids']

# ----------------------------------------------------------------------------------------
def get_true_partition(sample):
    cstr, clen, sstr, sval = sample.split('_')
    assert cstr == 'cdr' and sstr == 'set'
    simfn = '%s/families_cdr3l_%s_1e5_set_%s-simu.yaml' % (simdir, clen, sval)
    _, _, tcpath = utils.read_output(simfn, skip_annotations=True)
    return tcpath.best()

# ----------------------------------------------------------------------------------------
vsn = 'v0'
logstr = '' #synth-dist-0.02' # '5k'  # ''
ptndir = '/fh/fast/matsen_e/processed-data/partis/spisak-simu/%s' % vsn
simdir = '/fh/fast/matsen_e/data/spisak-simu/processed-data/simu'

outdirs = sorted(glob.glob('%s/partitions/cdr_*set_?%s'%(ptndir, '' if logstr=='' else '-%s'%logstr)))
print '  found %d output dirs (e.g. %s)' % (len(outdirs), outdirs[0])
for iod, odir in enumerate(outdirs):
    sample = os.path.basename(odir)

    # print sample
    # true_ptn = get_true_partition(sample)
    _, _, cpath = utils.read_output('%s/partition.yaml' % odir, skip_annotations=True)
    cpath.print_partitions(dont_print_clusters=True, print_header=iod==0, extrastr='    %s  '%sample)
    # inf_ptn = cpath.best()
    # utils.per_seq_correct_cluster_fractions(inf_ptn, true_ptn, debug=True)
