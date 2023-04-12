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
def get_true_partition(antn_list, inf_ptn=None):
    assert False  # don't need this any more
    bad_antns = [l for l in antn_list if l['invalid']]
    if inf_ptn is not None and len(bad_antns) > 0:
        print '  %s adding %d failed annotations as singletons to inf_ptn' % (utils.wrnstr(), len(bad_antns))
        inf_ptn += [l['unique_ids'] for l in bad_antns]
    all_fids = {u : fid for l in antn_list for u, fid in zip(l['unique_ids'], getfam(l))}
    def kfcn(u): return all_fids[u]
    true_partition = [list(group) for _, group in itertools.groupby(sorted(all_fids, key=kfcn), key=kfcn)]
    return true_partition, inf_ptn

# ----------------------------------------------------------------------------------------
vsn = 'test'
logstr = '' #synth-dist-0.02' # '5k'  # ''
bd = '/fh/fast/matsen_e/processed-data/partis/spisak-simu/%s' % vsn

outdirs = sorted(glob.glob('%s/partitions/cdr3l_*set_?%s'%(bd, '' if logstr=='' else '-%s'%logstr)))
print '  found %d output dirs (e.g. %s)' % (len(outdirs), outdirs[0])
for odir in outdirs:
    sample = os.path.basename(odir)
    ofn = '%s/partition.yaml' % odir
    print sample
    _, antn_list, cpath = utils.read_output(ofn)
    true_ptn, inf_ptn = get_true_partition(antn_list, inf_ptn=cpath.best())
    utils.per_seq_correct_cluster_fractions(inf_ptn, true_ptn, debug=True)
