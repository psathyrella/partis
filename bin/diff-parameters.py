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
import re
import yaml
from subprocess import check_output

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true', help='Passed on to ROOT when plotting')
sys.argv.append('-b')
parser.add_argument('--dir1', required=True)
parser.add_argument('--dir2', required=True)
args = parser.parse_args()

def reduce_float_precision(line):
    n_digits_to_keep = 10
    regex = '([0-9]\.' + n_digits_to_keep*'[0-9]' + ')([0-9][0-9]*)'
    match = re.search(regex, line)
    while match is not None:
        start, end = match.groups()
        line = line.replace(start + end, start)
        match = re.search(regex, line)
    return line

def get_lines(fname, reduce_precision=False):
    try:
        with open(fname) as infile:
             lines = sorted(infile.readlines())
             if reduce_precision:
                 for iline in range(len(lines)):
                     lines[iline] = reduce_float_precision(lines[iline])
             return lines
        
    except IOError:
        raise Exception('ERROR ' + os.path.basename(fname) + ' not found in ' + os.path.dirname(fname))

def check_lines(lines1, lines2):
    difflines = []
    for line in lines1:
        if line not in lines2:
            difflines.append(line)
            # print 'not in', args.dir2, ':\n  ', line
    return difflines

def check_textfile(fname):
    lines1 = get_lines(args.dir1 + fname, reduce_precision=True)
    lines2 = get_lines(args.dir2 + fname, reduce_precision=True)
    difflines1 = check_lines(lines1, lines2)
    difflines2 = check_lines(lines2, lines1)
    if len(difflines1 + difflines2) > 0:
        print 'differing lines from', args.dir1 + fname
        for line in difflines1:
            print ' ', line.strip()
        print 'differing lines from', args.dir2 + fname
        for line in difflines2:
            print ' ', line.strip()
        sys.exit(1)

# def diff_models(model1, model2):
#     for state in model1.states:
#         if state
#         for transition in 
#         sys.exit()

# def check_yaml(fname):
#     print fname
#     with open(args.dir1 + fname) as infile1:
#         with open(args.dir2 + fname) as infile2:
#             model1 = yaml.load(infile1)
#             model2 = yaml.load(infile2)
#             diff_models(model1, model2)
#             diff_models(model2, model1)

def get_file_list(extension=''):
    cmd = ['find', args.dir1]
    if extension == '':
        cmd += ['-type', 'f']
    else:
        cmd += ['-name', '*.' + extension]
    output = check_output(cmd)
    return output.replace(args.dir1, '').split()

# ----------------------------------------------------------------------------------------
# for fname in get_file_list('csv'):
#     check_textfile(fname)
# for fname in get_file_list('yaml'):
#     check_textfile(fname)

for fname in get_file_list():
    if '.csv' in fname or '.yaml' in fname:
        check_textfile(fname)
    elif '.svg' in fname or '.html' in fname:
        continue
    else:
        raise Exception('ERROR ' + fname + ' has an extension I can\'t handle')

sys.exit(0)
