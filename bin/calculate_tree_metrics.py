#!/usr/bin/env python
import sys
import time
from cStringIO import StringIO
import Bio.Phylo
import argparse
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
args = parser.parse_args()

tree = treeutils.get_bio_tree(treefname=args.treefile)

treeutils.calculate_LBI(tree)
sys.exit()
# ----------------------------------------------------------------------------------------

import dendropy
with open(args.treefile) as treefile:
    treestr = '\n'.join(treefile.readlines())
# dendro_tree = dendropy.Tree.get_from_string(treestr, 'newick')
# dendro_depths = {str(n.taxon).strip('\'') : n.distance_from_root() for n in dendro_tree if n.taxon is not None}
# dendro_depths = treeutils.get_depths(treestr, 'dendropy')

# out_dtree.reroot_at_node(out_dtree.find_node_with_taxon_label(naive_seq_name), update_bipartitions=True)
# out_tree = out_dtree.as_string(schema='newick', suppress_rooting=True)

dendro_depths = treeutils.get_depths(treestr, 'Bio')
depths = tree.depths()  # keyed by clade, not clade name
for node in tree.find_clades():
    node.num_date = depths[node]
    if node.name is not None:
        print '%20s  %16.18f   %16.18f %s' % (node.name, depths[node], dendro_depths[node.name], utils.color('red', 'x',) if depths[node] != dendro_depths[node.name] else '')
sys.exit()
