#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import glob
import sys
import csv
from io import open
csv.field_size_limit(sys.maxsize)  # make sure we can write very large csv fields
import os
import copy
import argparse
import colored_traceback.always
import json

import partis.utils as utils
import partis.paircluster as paircluster
import partis.glutils as glutils
from partis.clusterpath import ClusterPath
import partis.treeutils as treeutils

gctree_outstr = 'gctree.out.inference.1'

helpstr = """
Run partis selection metrics on gctree output dir (gctree docs: https://github.com/matsengrp/gctree/).
Plots are written to --outdir; it is probably best to start by looking at the summary html with a browser, e.g.: firefox <--outdir>/selection-metrics/plots/inferred-tree-metrics/overview.html.
Log files are written to <--outdir>/*.log; get-selection-metrics.log has most of the interesting information (view with less -RS).
Example usage:
  single chain:
    ./bin/read-gctree-output.py --locus igh --gctreedir <gctree-output-dir> --outdir <dir-for-partis-output>
  paired:
    ./bin/read-gctree-output.py --paired-loci --seqfname <fasta-input-file> --gctreedir <gctree-output-dir> --outdir <dir-for-partis-output>
  other args, with examples (see datascripts/meta/taraki-gctree-2021-10/partis-run.py):
    --kdfname /fh/fast/matsen_e/data/taraki-gctree-2021-10/processed-data/determistic/gc1/kdvals.csv
    --tree-basename tree.nwk --kd-columns delta_bind_CGG_FVS_additive --dont-invert-kd --multiplicity-column multiplicity --species mouse --no-insertions-or-deletions
    --initial-germline-dir /home/dralph/work/partis/datascripts/meta/taraki-gctree-2021-10/germlines --parameter-plots
    --slice-bin-fname /home/dralph/work/partis/datascripts/meta/taraki-gctree-2021-10/slice-bins.yaml
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
all_actions = ['cache-parameters', 'annotate', 'get-selection-metrics']
parser.add_argument('--actions', default=':'.join(all_actions), help='colon-separated list of actions to run')
parser.add_argument('--seqfname', help='Fasta file with input sequences. If single chain, this defaults to the standard location in --gctreedir. If --paired-loci is set, this should include, separately, all heavy and all light sequences, where the two sequences in a pair have identical uids (at least up to the first \'_\').')
parser.add_argument('--gctreedir', required=True, help='gctree output dir (to get --tree-basename, and maybe abundances.csv, --seqfname).')
parser.add_argument('--outdir', required=True, help='directory to which to write partis output files')
parser.add_argument('--input-partition-fname', help='partis style yaml file with a partition grouping seqs into clonal families; if set, input data is assumed to contain many families (if not set, we assume it\'s only one fmaily).')
parser.add_argument('--paired-loci', action='store_true', help='run on paired heavy/light data')
parser.add_argument('--locus', choices=utils.loci, help='locus of sequences (required for single chain).')
parser.add_argument('--kdfname', help='csv file with kd values (and, optionally, multiplicities), with header names as specified in subsequent args.')
parser.add_argument('--name-column', default='name', help='column name in --kdfname from which to take sequence name')
parser.add_argument('--kd-columns', default='kd', help='colon-separated list of column name[s] in --kdfname from which to take kd values. If more than one, the values are added.')
parser.add_argument('--multiplicity-column', help='If set, column name in --kdfname from which to take multiplicity value (which must be >0, i.e. inferred ancestors should have multiplicity 1). If not set, abundances are read from --abundance-basename in --gctreedir and converted to multiplicities.')
parser.add_argument('--dont-invert-kd', action='store_true', help='by default we invert (take 1/kd) to convert to \'affinity\' (after adding multiple kd columns, if specified), or at least something monotonically increasing with affinity. This skips that step, e.g. if you\'re passing in affinity.')
parser.add_argument('--species', default='mouse', choices=('human', 'macaque', 'mouse'))
parser.add_argument('--tree-basename', default='%s.nk'%gctree_outstr, help='basename of tree file to take from --gctreedir')  # .1 is the most likely one (all trees are also in the pickle file as ete trees: gctree.out.inference.parsimony_forest.p
parser.add_argument('--abundance-basename', default='abundances.csv', help='basename of abundance file in --gctreedir. Abundance of 0 (inferred ancestor) is converted to multiplicity of 1. Not used if multiplicities are read from kdfname')  # .1 is the most likely one (all trees are also in the pickle file as ete trees: gctree.out.inference.parsimony_forest.p
parser.add_argument('--dry', action='store_true')
parser.add_argument('--no-tree-plots', action='store_true')
parser.add_argument('--parameter-plots', action='store_true')
parser.add_argument('--no-insertions-or-deletions', action='store_true', help='see partis help')
parser.add_argument('--n-procs', type=int)
parser.add_argument('--slice-bin-fname')
parser.add_argument('--initial-germline-dir', help='see partis help')
args = parser.parse_args()
args.actions = utils.get_arg_list(args.actions, choices=all_actions)
args.kd_columns = utils.get_arg_list(args.kd_columns)
if args.multiplicity_column is not None and args.kdfname is None:
    raise Exception('have to set --kdfname if --multiplicity-column is set')
if args.paired_loci:
    assert args.seqfname is not None
else:
    if args.seqfname is None:
        args.seqfname = '%s/%s.fasta' % (args.gctreedir, gctree_outstr)
        print('    set --seqfname to default location in --gctreedir: %s' % args.seqfname)
    if args.locus is None:
        raise Exception('have to set --locus for single chain')

# ----------------------------------------------------------------------------------------
def metafname():
    return '%s/gctree-meta.yaml' % args.outdir

# ----------------------------------------------------------------------------------------
def run_cmd(action):
    locstr = '--paired-loci' if args.paired_loci else '--locus %s'%args.locus
    # this doesn't work since all gcs together needs 0:1 but single gc runs need 0, and the guessing functionality is working fine atm --droplet-id-separators - --droplet-id-indices 0:1
    # NOTE would maybe be better to not guess pair info?
    cmd = './bin/partis %s %s --species %s --guess-pairing-info --input-metafnames %s' % (action, locstr, args.species, metafname())
    if args.no_insertions_or_deletions:
        cmd += ' --no-insertions-or-deletions'
    if action in ['cache-parameters', 'annotate']:
        cmd += ' --infname %s' % args.seqfname
        if args.paired_loci:
            cmd += ' --paired-outdir %s' % args.outdir
        else:
            cmd += ' --parameter-dir %s/parameters' % args.outdir
        if args.input_partition_fname is None:  # one gc at a time
            cmd += ' --all-seqs-simultaneous'
        else:  # many gcs together
            cmd += ' --input-partition-fname %s' % args.input_partition_fname
    if action == 'cache-parameters':
        if args.initial_germline_dir is not None:
            cmd += ' --initial-germline-dir %s' % args.initial_germline_dir
        if args.parameter_plots:
            cmd += ' --plotdir %s' % args.outdir
    if action in ['annotate', 'get-selection-metrics'] and '--paired-outdir' not in cmd:
        cmd += ' --%s %s%s' % ('paired-outdir' if args.paired_loci else 'outfname', args.outdir, '' if args.paired_loci else '/partition.yaml')
    if action == 'get-selection-metrics':
        cmd += ' --min-selection-metric-cluster-size 3 --treefname %s/%s --plotdir %s --selection-metrics-to-calculate lbi:aa-lbi:cons-dist-aa:lbr:aa-lbr:lbf:aa-lbf' % (args.gctreedir, args.tree_basename, 'paired-outdir' if args.paired_loci else '%s/selection-metrics/plots'%args.outdir)
        cmd += ' --extra-daffy-metrics lbi:aa-lbi'
        cmd += ' --label-root-node'
        plt_cfg = treeutils.default_plot_cfg + ['distr', 'tree-mut-stats']
        if args.no_tree_plots:
            plt_cfg = [t for t in plt_cfg if t != 'tree']
        cmd += ' --add-selection-metrics-to-outfname --use-droplet-id-for-combo-id --selection-metric-plot-cfg %s' % ':'.join(plt_cfg)
        if args.slice_bin_fname is not None:
            cmd += ' --slice-bin-fname %s' % args.slice_bin_fname
        cmd += ' --choose-all-abs --chosen-ab-fname %s/chosen-abs.csv' % args.outdir  #  --debug 1
    if args.n_procs is not None:
        cmd += ' --n-procs %d' % args.n_procs
    utils.simplerun(cmd, logfname='%s/%s.log'%(args.outdir, action), dryrun=args.dry)

# ----------------------------------------------------------------------------------------
utils.mkdir(args.outdir)
metafos = {}
if args.multiplicity_column is None:  # if not set, read abundances from args.abundance_basename
    abfn = '%s/%s' % (args.gctreedir, args.abundance_basename)
    print('    reading abundance info from %s' % abfn)
    with open(abfn) as afile:
        reader = csv.DictReader(afile, fieldnames=('name', 'abundance'))
        for line in reader:
            if line['name'] not in metafos:
                metafos[line['name']] = {}
            metafos[line['name']]['multiplicity'] = max(1, int(line['abundance']))  # increase 0s (inferred ancestors) to 1
if args.kdfname is not None:
    print('    reading kd info%s from %s' % ('' if args.multiplicity_column is None else ' and multiplicity info', args.kdfname))
    with open(args.kdfname) as kfile:
        reader = csv.DictReader(kfile)
        for line in reader:
            uid = line[args.name_column]
            if uid not in metafos:
                metafos[uid] = {}
            if all(line[k] not in ['None', None, ''] for k in args.kd_columns):
                kdval = sum(float(line[k]) for k in args.kd_columns)
                metafos[uid]['affinity'] = kdval if args.dont_invert_kd else 1. / kdval
            if args.multiplicity_column is not None:
                metafos[uid]['multiplicity'] = int(line[args.multiplicity_column])

if args.paired_loci:  # convert metafos to per-locus names
    for base_id in list(metafos.keys()):
        for ltmp in utils.sub_loci('ig'):
            new_id = '%s-%s' % (base_id, ltmp)
            metafos[new_id] = metafos[base_id]
        del metafos[base_id]

# and write to json/yaml
print('    writing input meta info to %s' % metafname())
utils.jsdump(metafname(), metafos)

for action in args.actions:
    run_cmd(action)
