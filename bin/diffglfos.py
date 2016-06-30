#!/usr/bin/env python
import argparse
import sys
import os
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
if not os.path.exists(partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % partis_dir
sys.path.insert(1, partis_dir + '/python')
import utils

parser = argparse.ArgumentParser()
parser.add_argument('file1')
parser.add_argument('file2')
args = parser.parse_args()

genes = {fn : set() for fn in (args.file1, args.file2)}

def readfile(fname):
    with open(fname) as fastafile:
        for line in fastafile:
            line = line.strip()
            if line == '':
                continue
            if line[0] != '>':
                continue
            gene = line.replace('>', '')
            genes[fname].add(gene)

readfile(args.file1)
readfile(args.file2)

print 'file1: %d' % len(genes[args.file1])
print 'file2: %d' % len(genes[args.file2])
print 'both: %d' % len(genes[args.file1] & genes[args.file2])
only_file1 = genes[args.file1] - genes[args.file2]
print 'only file1: %d  (%s)' % (len(only_file1), ' '.join([utils.color_gene(g) for g in only_file1]))
only_file2 = genes[args.file2] - genes[args.file1]
print 'only file2: %d  (%s)' % (len(only_file2), ' '.join([utils.color_gene(g) for g in only_file2]))
