#!/usr/bin/env python3
""" 
Replacement for diff -qr  -x\'*.svg\' -x params -x plots.html <arg1> <arg2>, but with more intelligent
control of yamls (e.g. floating point precision) and csvs (e.g. line order)
"""

from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import os
import sys
import csv
import re
import yaml
from subprocess import check_output
from io import open

# current_script_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '') #'/python')
# if not os.path.exists(current_script_dir):
#     print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % current_script_dir
# sys.path.insert(1, current_script_dir)

parser = argparse.ArgumentParser()
parser.add_argument('arg1')
parser.add_argument('arg2')
parser.add_argument('--keep-going', action='store_true', help='Don\'t fail on differences, instead just keep on chugging')
parser.add_argument('--precision', type=int, default=9, help='number of digits after the decimal place to keep when comparing floating point numbers')
args = parser.parse_args()

if os.path.isdir(args.arg1):  # can either pass arg[12] as directories in which to look
    args.dir1 = args.arg1
    assert os.path.isdir(args.arg2)
    args.dir2 = args.arg2
    args.fname = None
else:  # ...or as single files
    args.dir1 = os.path.dirname(args.arg1)
    assert os.path.exists(args.dir1)
    args.dir2 = os.path.dirname(args.arg2)
    assert os.path.exists(args.dir2)
    args.fname = os.path.basename(args.arg1)
    assert args.fname == os.path.basename(args.arg2)

def reduce_float_precision(line):
    regex = '([0-9]\.' + args.precision*'[0-9]' + ')([0-9][0-9]*)'
    match = re.search(regex, line)
    while match is not None:
        left, right = match.groups()  # left-hand (more significant) and right-hand parts of number
        assert '.' in left
        clipped_left = left  # clip trailing zeros in decimal numbers
        while clipped_left[-1] == '0':
            clipped_left = clipped_left[:-1]
        line = line.replace(left + right, clipped_left)
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
        raise Exception(os.path.basename(fname) + ' not found in ' + os.path.dirname(fname))

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
        print('differing lines from', args.dir1 + fname)
        for line in difflines1:
            print(' ', line.strip())
        print('differing lines from', args.dir2 + fname)
        for line in difflines2:
            print(' ', line.strip())
        if not args.keep_going:
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
    output = check_output(cmd, universal_newlines=True)
    return output.replace(args.dir1, '').split()

# ----------------------------------------------------------------------------------------
# for fname in get_file_list('csv'):
#     check_textfile(fname)
# for fname in get_file_list('yaml'):
#     check_textfile(fname)

if args.fname is None:
    filelist = get_file_list()
else:
    filelist = ['/' + args.fname, ]
for fname in filelist:
    if '.csv' in fname or '.yaml' in fname:
        check_textfile(fname)
    elif '.svg' in fname or '.html' in fname:
        continue
    else:
        raise Exception(fname + ' has an extension I can\'t handle')

sys.exit(0)
