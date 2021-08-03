#!/usr/bin/env python
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
run partis selection metrics on gctree output dir https://github.com/matsengrp/gctree/
Example usage:
  ./bin/read-gctree-output.py <gctree-output-dir> <dir-for-partis-output>
"""
class MultiplyInheritedFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass
formatter_class = MultiplyInheritedFormatter
parser = argparse.ArgumentParser(formatter_class=MultiplyInheritedFormatter, description=helpstr)
parser.add_argument('gctreedir', help='directory with gctree output files')
parser.add_argument('outdir', help='directory to which to write partis output files')
parser.add_argument('--locus', default='igh', choices=utils.loci)
parser.add_argument('--species', default='mouse', choices=('human', 'macaque', 'mouse'))
parser.add_argument('--add-dj-seqs', action='store_true', help='if input seqs only contain V (up to cdr3), you\'ll need to set this so some dummy D and J bits are appended to the actual input seqs.')
parser.add_argument('--drop-zero-abundance-seqs', action='store_true', help='By default we increase the abundance of any seqs with 0 abundance (typically inferred ancestors). If you set this, we instead ignore them when running partis annotate (although if they\'re in the tree, they still show up in the selection metrics).')
parser.add_argument('--tree-index', type=int, default=1, help='gctree outputs a list of trees in order of decreasing likelihood, and by default we take the first one, but this arg lets you select one of the others.')
args = parser.parse_args()

utils.mkdir(args.outdir)
gcoutstr = 'gctree.out.inference.%d' % args.tree_index
metafname = '%s/meta.yaml' % args.outdir
infname = '%s/input-seqs.fa' % args.outdir
outfname = '%s/partition.yaml' % args.outdir

# ----------------------------------------------------------------------------------------
def write_input():
    seqfos = utils.read_fastx('%s/%s.fasta'%(args.gctreedir, gcoutstr))  # .1 is the most likely one (all trees are also in the pickle file as ete trees: gctree.out.inference.parsimony_forest.p
    metafos = {}
    n_to_skip, n_abundance_increased = 0, 0
    for sfo in seqfos:
        name, abfo_str = sfo['infostrs']
        assert abfo_str.count('=') == 1
        abstr, abval = abfo_str.split('=')
        assert abstr == 'abundance'
        if int(abval) == 0:
            if args.drop_zero_abundance_seqs:
                n_to_skip += 1
            else:
                abval = '1'
                n_abundance_increased += 1
        metafos[sfo['name']] = {'multiplicity' : int(abval)}
    if n_to_skip > 0:
        n_before = len(seqfos)
        seqfos = [s for s in seqfos if metafos[s['name']]['multiplicity'] > 0]
        print '    skipping %d / %d seqs with zero abundance' % (n_to_skip, n_before)
        assert n_before == len(seqfos) + n_to_skip
    if n_abundance_increased > 0:
        print '    increased abundance of %d / %d seqs from 0 to 1 (these are presumably inferred ancestral sequences)' % (n_abundance_increased, len(seqfos))
    with open(metafname, 'w') as mfile:
        json.dump(metafos, mfile)

    if args.add_dj_seqs:
        glfo = glutils.read_glfo('data/germlines/%s'%args.species, args.locus)
        tgenes = {r : sorted(glfo['seqs'][r].keys())[0] for r in 'dj'}
        d_seq, j_seq = [glfo['seqs'][r][tgenes[r]] for r in 'dj']
        print '    appending d and j gene sequences: %s %s' % (utils.color_gene(tgenes['d']), utils.color_gene(tgenes['j']))
        for sfo in seqfos:
            sfo['seq'] = '%s%s%s' % (sfo['seq'], d_seq, j_seq)
    print '    writing input file to %s' % infname
    utils.write_fasta(infname, seqfos)

# ----------------------------------------------------------------------------------------
write_input()
cmd = './bin/partis cache-parameters --infname %s --input-metafname %s --locus %s --species %s --parameter-dir %s/parameters' % (infname, metafname, args.locus, args.species, args.outdir)
utils.simplerun(cmd, logfname='%s/cache-parameters.log'%args.outdir)
cmd = './bin/partis annotate --infname %s --locus %s --species %s --parameter-dir %s/parameters --all-seqs-simultaneous --outfname %s' % (infname, args.locus, args.species, args.outdir, outfname)
utils.simplerun(cmd, logfname='%s/annotate.log'%args.outdir)
cmd = './bin/partis get-selection-metrics --outfname %s --treefname %s/%s.nk --plotdir %s/selection-metrics/plots --queries-to-include-fname %s --debug 1' % (outfname, args.gctreedir, gcoutstr, args.outdir, infname)
utils.simplerun(cmd, logfname='%s/get-selection-metrics.log'%args.outdir)
