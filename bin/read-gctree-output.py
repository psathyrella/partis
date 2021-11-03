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
import glutils
from clusterpath import ClusterPath

helpstr = """
Run partis selection metrics on gctree output dir https://github.com/matsengrp/gctree/.
Plots are written to <--outdir>/selection-metrics/plots.
Log files are written to <--outdir>; get-selection-metrics.log has most of the information you probably want (view with less -RS).
Example usage:
  ./bin/read-gctree-output.py --seqfname <fasta-sequence-file> --treefname <gctree-tree-file> --outdir <dir-for-partis-output>
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('--seqfname', required=True, help='fasta file with all input sequences (if --paired-loci is set, this should include, separately, all heavy and all light sequences, where the two sequences in a pair have identical uids [at least up to the first \'_\'])')
parser.add_argument('--treefname', required=True, help='input tree file from gctree output dir, e.g. /path/to/output/gctree.out.inference.1.fasta')  # .1 is the most likely one (all trees are also in the pickle file as ete trees: gctree.out.inference.parsimony_forest.p
parser.add_argument('--outdir', required=True, help='directory to which to write partis output files')
parser.add_argument('--paired-loci', action='store_true', help='run on paired heavy/light data. <gctreedir> must contain a subdir for each locus with names among %s' % utils.sub_loci('ig'))
parser.add_argument('--locus', default='igh', choices=utils.loci)
parser.add_argument('--kdfname', help='csv file with kd values')
parser.add_argument('--kd-file-headers', default='name:kd', help='colon-separated list of two column names in --kdfname, the first for sequence name (which must match node names in --treefname), the second for kd values.')
parser.add_argument('--dont-invert-kd', action='store_true', help='by default with invert (take 1/kd) to convert to \'affinity\', or at least something monotonically increasing with affinity. This skips that step, e.g. if you\'re passing in affinity.')
parser.add_argument('--species', default='mouse', choices=('human', 'macaque', 'mouse'))
parser.add_argument('--add-dj-seqs', action='store_true', help='if input seqs only contain V (up to cdr3), you\'ll need to set this so some dummy D and J bits are appended to the actual input seqs.')
parser.add_argument('--drop-zero-abundance-seqs', action='store_true', help='By default we increase the abundance of any seqs with 0 abundance (typically inferred ancestors). If you set this, we instead ignore them when running partis annotate (although if they\'re in the tree, they still show up in the selection metrics).')
parser.add_argument('--dry', action='store_true')
args = parser.parse_args()
args.kd_file_headers = utils.get_arg_list(args.kd_file_headers, key_list=['name', 'kd'])
if not args.paired_loci:
    raise Exception('needs testing')

# ----------------------------------------------------------------------------------------
def metafname():
    return '%s/kd-meta-info.yaml' % args.outdir

# ----------------------------------------------------------------------------------------
def run_cmd(action):
    locstr = '--paired-loci' if args.paired_loci else '--locus %s'%args.locus
    cmd = './bin/partis %s %s --species %s --read-gctree-output' % (action, locstr, args.species)
    if action in ['cache-parameters', 'annotate']:
        cmd += ' --infname %s' % args.seqfname
        if args.paired_loci:
            cmd += ' --paired-outdir %s' % args.outdir
        else:
            cmd += ' --parameter-dir %s/parameters' % args.outdir
        if args.kdfname is not None:
            cmd += ' --input-metafnames %s' % metafname()
    if args.add_dj_seqs:
        cmd += ' --add-dj-seqs'
    if action == 'annotate':
        cmd += ' --all-seqs-simultaneous'
    if action in ['annotate', 'get-selection-metrics'] and '--paired-outdir' not in cmd:
        cmd += ' --%s %s%s' % ('paired-outdir' if args.paired_loci else 'outfname', args.outdir, '' if args.paired_loci else '/partition.yaml')
    if action == 'get-selection-metrics':
        cmd += ' --treefname %s --plotdir %s/selection-metrics/plots --selection-metrics-to-calculate cons-dist-aa:aa-lbi:aa-lbr --debug 1' % (args.treefname, args.outdir) #  --queries-to-include-fname %s, args.seqfname)# TODO wtf was qti set here?
    utils.simplerun(cmd, logfname='%s/%s.log'%(args.outdir, action), dryrun=args.dry)

# ----------------------------------------------------------------------------------------
utils.mkdir(args.outdir)
if args.kdfname is not None:
    kdvals = {}
    with open(args.kdfname) as kfile:
        reader = csv.DictReader(kfile)
        for line in reader:
            base_id = line[args.kd_file_headers['name']]
            for ltmp in utils.sub_loci('ig'):
                affy = float(line[args.kd_file_headers['kd']])
                if not args.dont_invert_kd:
                    affy = 1. / affy
                kdvals['%s-%s' % (base_id, ltmp)] = {'affinity' : affy}
    with open(metafname(), 'w') as mfile:
        json.dump(kdvals, mfile)
for action in ['cache-parameters', 'annotate', 'get-selection-metrics']:
    run_cmd(action)
