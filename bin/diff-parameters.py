#!/usr/bin/env python
""" 
Replacement for diff -qr  -x\'*.svg\' -x params -x plots.html <dir1> <dir2>, but with more intelligent
control of yamls (e.g. floating point precision) and csvs (e.g. line order)
"""

import argparse
import os
import sys
sys.path.insert(1, './python')
import csv
import yaml
from subprocess import check_output

parser = argparse.ArgumentParser()
parser.add_argument('--dir1', required=True)
parser.add_argument('--dir2', required=True)
args = parser.parse_args()

def get_lines(fname):
    try:
        with open(fname) as infile:
            return sorted(infile.readlines())
    except IOError:
        raise Exception('ERROR ' + os.path.basename(fname) + ' not found in ' + os.path.dirname(fname))

def check_lines(lines1, lines2):
    difflines = []
    for line in lines1:
        if line not in lines2:
            difflines.append(line)
            # print 'not in', args.dir2, ':\n  ', line
    return difflines

def check_csv(fname):
    print fname
    lines1 = get_lines(args.dir1 + fname)
    lines2 = get_lines(args.dir2 + fname)
    difflines1 = check_lines(lines1, lines2)
    difflines2 = check_lines(lines2, lines1)
    if len(difflines1 + difflines2) > 0:
        print 'differing lines from', args.dir1
        for line in difflines1:
            print ' ', line.strip()
        print 'differing lines from', args.dir2
        for line in difflines2:
            print ' ', line.strip()
        sys.exit(1)

csvlist = check_output(['find', args.dir1, '-name', '*.csv']).replace(args.dir1, '').split()
for fname in csvlist:
    check_csv(fname)
