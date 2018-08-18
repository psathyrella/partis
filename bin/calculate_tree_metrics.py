#!/usr/bin/env python
import sys
import time
import argparse
import dendropy
import os

# start = time.time()
# print '      time: %.1f' % (time.time() - start)

partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')
sys.path.insert(1, partis_dir + '/packages/baltic')
import utils
import treeutils

parser = argparse.ArgumentParser()
parser.add_argument('treefile')
parser.add_argument('--naive-seq-name', required=True)
args = parser.parse_args()

tree = treeutils.get_bio_tree(treefname=args.treefile)

treeutils.calculate_LBI(tree, debug=True)
sys.exit()

# print treeutils.get_ascii_tree(treeutils.get_treestr(args.treefile))
# dendro_tree = treeutils.get_dendro_tree(treefname=args.treefile)
# dendro_tree.reroot_at_node(dendro_tree.find_node_with_taxon_label(args.naive_seq_name), update_bipartitions=True)
# # print dendro_tree.as_ascii_plot(width=100)  # why tf does this show them as all the same depth?
# treestr = dendro_tree.as_string(schema='newick')  #, suppress_rooting=True)
# print treeutils.get_ascii_tree(treestr)
