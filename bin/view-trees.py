#!/usr/bin/env python
import argparse
import pickle
import csv
from ete3 import TreeNode, TreeStyle, NodeStyle, SVG_COLORS

parser = argparse.ArgumentParser()
parser.add_argument('treefname')
parser.add_argument('kdfname')
args = parser.parse_args()

with open(args.treefname, 'rb') as lfile:
    tree = pickle.load(lfile)
  
# print tree
# print tree.name
# print tree.sequence

with open(args.kdfname, 'w') as kdfile:
    writer = csv.DictWriter(kdfile, ('uid', 'kd'))
    writer.writeheader()
    for node in tree.traverse():  # small kd is higher affinity
        writer.writerow({'uid' : node.name, 'kd' : node.Kd})
