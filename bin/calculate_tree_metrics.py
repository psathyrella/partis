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
def run_lbi(args):
    if args.overwrite:
        print '%s --overwrite not implemented for lbi' % utils.color('red', 'warning')
    if args.treefile is None:
        raise Exception('need to set --treefile to run lbi (could instead use tree from lonr, maybe I should implement that?)')

    if args.reroot_at_naive:
        print treeutils.get_ascii_tree(treeutils.get_treestr(args.treefile))
        dendro_tree = treeutils.get_dendro_tree(treefname=args.treefile)
        dendro_tree.reroot_at_node(dendro_tree.find_node_with_taxon_label(args.naive_seq_name), update_bipartitions=True)
        # print dendro_tree.as_ascii_plot(width=100)  # why tf does this show them as all the same depth?
        treestr = dendro_tree.as_string(schema='newick')  #, suppress_rooting=True)
        print treeutils.get_ascii_tree(treestr)
        bio_tree = treeutils.get_bio_tree(treestr=treestr)
    else:
        bio_tree = treeutils.get_bio_tree(treefname=args.treefile)

    treeutils.calculate_LBI(bio_tree)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
available_metrics = ['lbi', 'lonr']
parser.add_argument('--metrics', default=':'.join(available_metrics), help='colon-separated list of tree metrics to calculate (choose from: %s)' % ' '.join(available_metrics))
parser.add_argument('--reroot-at-naive', action='store_true')
parser.add_argument('--lonr-tree-method', default='dnapars', choices=['dnapars', 'neighbor'], help='which phylip method should lonr use to infer the tree (maximum parsimony or neighbor-joining)? (their original defaults were dnapars for less than 100 sequences, neighbor for more)')
parser.add_argument('--debug', action='store_true')

# input
parser.add_argument('--treefile', help='input tree file name in newick format')
parser.add_argument('--seqfile', help='input fasta file with aligned sequences corresponding to <treefile>')
parser.add_argument('--naive-seq-name', help='uid of inferred naive sequence')

# output
parser.add_argument('--outfile', help='output file name in yaml format')
parser.add_argument('--lonr-outdir', help='directory for the various lonr output files')
parser.add_argument('--overwrite', action='store_true')

args = parser.parse_args()
args.outfile = os.path.abspath(args.outfile)
args.metrics = utils.get_arg_list(args.metrics)
if len(set(args.metrics) - set(available_metrics)) > 0:
    raise Exception('unhandled metric(s): %s (choose from: %s)' % (' '.join(set(args.metrics) - set(available_metrics)), ' '.join(available_metrics)))
if args.reroot_at_naive and args.naive_seq_name is None:
    raise Exception('have to specify --naive-seq-name if --reroot-at-naive is set')

# ----------------------------------------------------------------------------------------
output_info = {}
if 'lbi' in args.metrics:
    run_lbi(args, debug=args.debug)
if 'lonr' in args.metrics:
    output_info['lonr'] = treeutils.calculate_lonr(args.seqfile, args.naive_seq_name, args.lonr_outdir, args.lonr_tree_method, treefile=args.treefile, overwrite=args.overwrite, reroot_at_naive=args.reroot_at_naive, debug=args.debug)

if not os.path.exists(os.path.dirname(args.outfile)):
    os.makedirs(os.path.dirname(args.outfile))
with open(args.outfile, 'w') as outfile:
    yaml.dump(output_info, outfile, width=400)
