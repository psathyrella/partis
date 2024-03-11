#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
import os
import sys
import argparse
import pickle
import csv
from ete3 import TreeNode, TreeStyle, NodeStyle, SVG_COLORS
from io import open

# NOTE this probably doesn't need to be a separate script any more, it used to be necessary since ete3 requires python 3

parser = argparse.ArgumentParser()
parser.add_argument('--pickle-tree-file', help='bcr phylo output pickle tree file')
parser.add_argument('--kdfile', help='output csv file with kd info')
parser.add_argument('--newick-tree-file', required=True, help='output newick tree file')
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

if args.kdfile is not None:
    with open(args.kdfile, 'wb' if sys.version_info.major < 3 else 'w') as kdfile:
        writer = csv.DictWriter(kdfile, ('uid', 'kd', 'time', 'relative_kd', 'lambda', 'target_index', 'target_distance'))
        writer.writeheader()
        for node in tree.traverse():  # small kd is higher affinity
            if node.name == '':
                continue
            writer.writerow({'uid' : node.name, 'kd' : node.Kd, 'time' : node.time, 'relative_kd' : node.relative_Kd, 'lambda' : node.lambda_, 'target_index' : node.target_index, 'target_distance' : node.target_distance})
