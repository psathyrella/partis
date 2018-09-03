#!/usr/bin/env python
import yaml
import sys
import time
import argparse
import os
import colored_traceback.always

# start = time.time()
# print '      time: %.1f' % (time.time() - start)

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')
sys.path.insert(1, partis_dir + '/packages/baltic')
import utils
import treeutils

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
available_metrics = ['lbi', 'lonr']
parser.add_argument('--metrics', default=':'.join(available_metrics), help='colon-separated list of tree metrics to calculate (choose from: %s)' % ' '.join(available_metrics))
parser.add_argument('--lonr-tree-method', default='dnapars', choices=['dnapars', 'neighbor'], help='which phylip method should lonr use to infer the tree (maximum parsimony or neighbor-joining)? (their original defaults were dnapars for less than 100 sequences, neighbor for more)')
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--debug', action='store_true')

# input
parser.add_argument('--naive-seq-name', required=True, help='uid of inferred naive sequence')
parser.add_argument('--treefile', help='input newick tree file, for lbi calculation. If not set, lonr has to be in --metrics so we can get the tree from there to use for lbi. If both --treefile is set and lonr is in --metrics, --treefile takes precedence')
parser.add_argument('--seqfile', help='input fasta file with aligned sequences corresponding to <treefile>')
parser.add_argument('--phylip-treefile', help='newick file with inferred tree from a previous phylip run (i.e. from cft)')
parser.add_argument('--phylip-seqfile', help='fasta file with node sequences (including internal nodes) from a previous phylip run (i.e. from cft)')

# output
parser.add_argument('--outfile', help='output file name in yaml format')

args = parser.parse_args()
args.outfile = os.path.abspath(args.outfile)

if args.phylip_treefile is not None:
    raise Exception('doesn\'t work (and lonr.r is probably poor enough quality code that it isn\'t worthwhile)')
    args.phylip_treefile = os.path.abspath(args.phylip_treefile)
if args.phylip_seqfile is not None:
    args.phylip_seqfile = os.path.abspath(args.phylip_seqfile)

args.metrics = utils.get_arg_list(args.metrics)
if len(set(args.metrics) - set(available_metrics)) > 0:
    raise Exception('unhandled metric(s): %s (choose from: %s)' % (' '.join(set(args.metrics) - set(available_metrics)), ' '.join(available_metrics)))

# ----------------------------------------------------------------------------------------
output_info = {}
if 'lonr' in args.metrics:
    output_info['lonr'] = treeutils.calculate_lonr(utils.read_fastx(args.seqfile), args.naive_seq_name, args.lonr_tree_method, phylip_treefile=args.phylip_treefile, phylip_seqfile=args.phylip_seqfile, seed=args.seed, debug=args.debug)
if 'lbi' in args.metrics:
    if args.treefile is not None:  # convert to nexml
        treestr = treeutils.get_dendro_tree(treefname=args.treefile, schema='newick', set_internal_node_labels=True, ignore_existing_internal_node_labels=('asttree' in args.treefile)).as_string('nexml')  # fasttree output is the only one so far that has shitty node labels
        print '  using --treefile for lbi tree'
    elif 'lonr' in output_info:
        treestr = output_info['lonr']['tree']
        print '  using lonr tree also for lbi'
    else:
        raise Exception('have to either set --treefile, or run lonr so we can get the tree from the lonr output')
    output_info['lbi'] = treeutils.calculate_lbi(args.naive_seq_name, treestr=treestr, debug=args.debug)

if not os.path.exists(os.path.dirname(args.outfile)):
    os.makedirs(os.path.dirname(args.outfile))
with open(args.outfile, 'w') as outfile:
    yaml.dump(output_info, outfile, width=400)
