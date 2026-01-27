#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import argparse
import colored_traceback.always
import glob
import itertools
import json
import yaml
import collections
import numpy
import math

import partis.utils as utils
import partis.glutils as glutils
import partis.hutils as hutils

# NOTE see also reviews/spisak-2023/notes.txt in github papers repo

# ----------------------------------------------------------------------------------------
def getsimfn(sample, csv=False):
    cstr, clen, sstr, sval = sample.split('_')
    assert cstr == 'cdr' and sstr == 'set'
    return '%s/%s/families_cdr3l_%s_1e5_set_%s%s' % (simdir, 'data' if csv else 'processed-data/simu', clen, sval, '.csv' if csv else '-simu.yaml')

# ----------------------------------------------------------------------------------------
def get_true_partition(sample):
    _, _, tcpath = utils.read_output(getsimfn(sample), skip_annotations=True)
    return tcpath.best()

# ----------------------------------------------------------------------------------------
def get_true_naive_seqs(sample, cluster):
    # clines = utils.csvlines(getsimfn(sample, csv=True), n_max_queries=100)
    glfo, true_antn_list, _ = utils.read_output(getsimfn(sample), dont_add_implicit_info=True)
    uid_lines = [l for l in true_antn_list if len(set(l['unique_ids']) & set(cluster)) > 0]
    reco_ids, naive_seqfos = [], []
    for tln in uid_lines:
        if tln['reco_id'] in reco_ids:
            continue
        reco_ids.append(tln['reco_id'])
        naive_seqfos.append({'name' : tln['reco_id'], 'seq' : tln['naive_seq']})
    utils.align_many_seqs(naive_seqfos, debug=True)

# ----------------------------------------------------------------------------------------
vsn = 'v0'
logstr = '5k' #simu-v0' #synth-dist-0.005' # '5k'  # ''
ptndir = '/fh/fast/matsen_e/processed-data/partis/spisak-simu/%s' % vsn
simdir = '/fh/fast/matsen_e/data/spisak-simu'

outdirs = sorted(glob.glob('%s/partitions/cdr_*set_?%s'%(ptndir, '' if logstr=='' else '-%s'%logstr)))
print('  found %d output dirs (e.g. %s)' % (len(outdirs), outdirs[0]))
for iod, odir in enumerate(outdirs):
    sample = os.path.basename(odir).replace('-'+logstr, '')

    # true_ptn = get_true_partition(sample)
    _, _, cpath = utils.read_output('%s/partition.yaml' % odir, skip_annotations=True)
    largest_cluster = sorted(cpath.best(), key=len, reverse=True)[0]
    # get_true_naive_seqs(sample, largest_cluster)
    # sys.exit()
    cpath.print_partitions(dont_print_clusters=True, print_header=iod==0, extrastr='    %s  '%sample)
    # inf_ptn = cpath.best()
    # utils.per_seq_correct_cluster_fractions(inf_ptn, true_ptn, debug=True)
    # break
