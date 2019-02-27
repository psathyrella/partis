#!/usr/bin/env python
import os
import sys
import argparse
import pickle
import csv
from ete3 import TreeNode, TreeStyle, NodeStyle, SVG_COLORS

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/packages/bcr-phylo-benchmark/bin')
import GCutils

parser = argparse.ArgumentParser()
parser.add_argument('--pickle-tree-file')
parser.add_argument('--kdfile', required=True)
parser.add_argument('--newick-tree-file', required=True)
args = parser.parse_args()

with open(args.pickle_tree_file, 'rb') as lfile:
    tree = pickle.load(lfile)

# print tree
# print tree.name
# print tree.sequence

with open(args.newick_tree_file, 'w') as ntfile:
    treestr = tree.write(format=1)  # default format ignores internal node names (numbers listed here: http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees)
    treestr = treestr.replace(';', '%s;' % tree.name)  # add root node name by hand (none of the format integers seem to add the root node name)
    ntfile.write(treestr)

with open(args.kdfile, 'w') as kdfile:
    writer = csv.DictWriter(kdfile, ('uid', 'kd'))
    writer.writeheader()
    for node in tree.traverse():  # small kd is higher affinity
        writer.writerow({'uid' : node.name, 'kd' : node.Kd})
