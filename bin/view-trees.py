#!/usr/bin/env python
import argparse
import pickle
import csv
from ete3 import TreeNode, TreeStyle, NodeStyle, SVG_COLORS

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
    ntfile.write(tree.write(format=1))  # default format ignores internal node names (numbers listed here: http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees)

with open(args.kdfile, 'w') as kdfile:
    writer = csv.DictWriter(kdfile, ('uid', 'kd'))
    writer.writeheader()
    for node in tree.traverse():  # small kd is higher affinity
        writer.writerow({'uid' : node.name, 'kd' : node.Kd})
