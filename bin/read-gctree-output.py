#!/usr/bin/env python
import glob
import sys
import csv
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import copy
import argparse
import colored_traceback.always
import json

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils
import paircluster
import glutils
from clusterpath import ClusterPath

helpstr = """
Run partis selection metrics on gctree output dir (gctree docs: https://github.com/matsengrp/gctree/).
Plots are written to <--outdir>/selection-metrics/plots.
Log files are written to <--outdir>; get-selection-metrics.log has most of the interesting information (view with less -RS).
Example usage:
  ./bin/read-gctree-output.py --seqfname <fasta-input-file> --gctreedir <gctree-output-dir> --outdir <dir-for-partis-output>
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('--seqfname', required=True, help='fasta file with all input sequences (if --paired-loci is set, this should include, separately, all heavy and all light sequences, where the two sequences in a pair have identical uids [at least up to the first \'_\'])')
parser.add_argument('--gctreedir', required=True, help='gctree output dir (to get --tree-basename and abundances.csv')
parser.add_argument('--outdir', required=True, help='directory to which to write partis output files')
parser.add_argument('--paired-loci', action='store_true', help='run on paired heavy/light data. <gctreedir> must contain a subdir for each locus with names among %s' % utils.sub_loci('ig'))
parser.add_argument('--locus', default='igh', choices=utils.loci)
parser.add_argument('--kdfname', help='csv file with kd values')
parser.add_argument('--kd-file-headers', default='name:kd', help='colon-separated list of two column names in --kdfname, the first for sequence name (which must match node names in --treefname), the second for kd values.')
parser.add_argument('--dont-invert-kd', action='store_true', help='by default with invert (take 1/kd) to convert to \'affinity\', or at least something monotonically increasing with affinity. This skips that step, e.g. if you\'re passing in affinity.')
parser.add_argument('--species', default='mouse', choices=('human', 'macaque', 'mouse'))
parser.add_argument('--tree-basename', default='gctree.out.inference.1.nk', help='basename of tree file to take from --gctreedir')  # .1 is the most likely one (all trees are also in the pickle file as ete trees: gctree.out.inference.parsimony_forest.p
parser.add_argument('--abundance-basename', default='abundances.csv', help='basename of tree file to take from --gctreedir')  # .1 is the most likely one (all trees are also in the pickle file as ete trees: gctree.out.inference.parsimony_forest.p
parser.add_argument('--dry', action='store_true')
args = parser.parse_args()
args.kd_file_headers = utils.get_arg_list(args.kd_file_headers, key_list=['name', 'kd'])
if not args.paired_loci:
    raise Exception('needs testing')

# ----------------------------------------------------------------------------------------
def metafname():
    return '%s/gctree-meta.yaml' % args.outdir

# ----------------------------------------------------------------------------------------
def run_cmd(action):
    locstr = '--paired-loci' if args.paired_loci else '--locus %s'%args.locus
    cmd = './bin/partis %s %s --species %s --guess-pairing-info' % (action, locstr, args.species)
    if action in ['cache-parameters', 'annotate']:
        cmd += ' --infname %s' % args.seqfname
        if args.paired_loci:
            cmd += ' --paired-outdir %s' % args.outdir
        else:
            cmd += ' --parameter-dir %s/parameters' % args.outdir
        if args.kdfname is not None:
            cmd += ' --input-metafnames %s' % metafname()
    if action == 'annotate':
        cmd += ' --all-seqs-simultaneous'
    if action in ['annotate', 'get-selection-metrics'] and '--paired-outdir' not in cmd:
        cmd += ' --%s %s%s' % ('paired-outdir' if args.paired_loci else 'outfname', args.outdir, '' if args.paired_loci else '/partition.yaml')
    if action == 'get-selection-metrics':
        cmd += ' --treefname %s/%s --plotdir %s --selection-metrics-to-calculate cons-dist-aa:aa-lbi:aa-lbr --queries-to-include-fname %s' % (args.gctreedir, args.tree_basename, 'paired-outdir' if args.paired_loci else '%s/selection-metrics/plots'%args.outdir, paircluster.paired_fn(args.outdir, 'igh'))
    utils.simplerun(cmd, logfname='%s/%s.log'%(args.outdir, action), dryrun=args.dry)

# ----------------------------------------------------------------------------------------
utils.mkdir(args.outdir)
metafos = {}
with open('%s/%s'%(args.gctreedir, args.abundance_basename)) as afile:
    reader = csv.DictReader(afile, fieldnames=('name', 'abundance'))
    for line in reader:
        for ltmp in utils.sub_loci('ig'):
            metafos['%s-%s' % (line['name'], ltmp)] = {'abundance' : max(1, int(line['abundance']))}  # increase 0s (inferred ancestors) to 1
if args.kdfname is not None:
    with open(args.kdfname) as kfile:
        reader = csv.DictReader(kfile)
        for line in reader:
            base_id = line[args.kd_file_headers['name']]
            for ltmp in utils.sub_loci('ig'):
                affy = float(line[args.kd_file_headers['kd']])
                if not args.dont_invert_kd:
                    affy = 1. / affy
                uid = '%s-%s' % (base_id, ltmp)
                if uid not in metafos:
                    metafos[uid] = {}
                metafos[uid]['affinity'] = affy
with open(metafname(), 'w') as mfile:
    json.dump(metafos, mfile)

for action in ['cache-parameters', 'annotate', 'get-selection-metrics']:
    run_cmd(action)
