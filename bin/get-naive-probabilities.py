#!/usr/bin/env python
import os
import sys
import csv
import argparse
import operator
import argparse
import yaml
import colored_traceback.always

# if you move this script, you'll need to change this method of getting the imports
partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')

import utils

parser = argparse.ArgumentParser()
parser.add_argument('infname')
parser.add_argument('--non-summed-column', default='v_gene', help='column for which we print counts for each value')
parser.add_argument('--skip-column-fname', help='yaml file with columns : list of str(values) for which we want to specify particular values (and skip others). Default/example set below.')
args = parser.parse_args()

if args.skip_column_fname is None:
    skip_column_vals = {  # to input your own dict on the command line, just convert with str() and quote it
        # 'cdr3_length' : ['33', '36', '39', '42', '45', '48'],  # <value> is list of acceptable values NOTE need to all be strings, otherwise you have to worry about converting the values in the csv file

        # bf520.1:
        'v_gene' : ['IGHV1-2*02+G35A', 'IGHV1-2*02+T147C', 'IGHV1-2*02'],
        # 'd_gene' : ['IGHD3-22*01'],
        'j_gene' : ['IGHJ4*02'],
        'cdr3_length' : ['66',],  #  TGTGCGAGAGGGCCATTCCCGAATTACTATGGTCCGGGGAGTTATTGGGGGGGTTTTGACCACTGG
    }
else:
    with open(args.skip_column_fname) as yamlfile:
        skip_column_vals = yaml.load(yamlfile)

info = {}
lines_skipped, lines_used = 0, 0
counts_skipped, counts_used = 0, 0
with open(args.infname) as csvfile:
    reader = csv.DictReader(csvfile)
    print '  all columns in file: %s' % ' '.join(reader.fieldnames)
    if len(set(skip_column_vals) - set(reader.fieldnames)) > 0:
        raise Exception('keys in --skip-column-fname not in file: %s' % ' '.join(set(skip_column_vals) - set(reader.fieldnames)))
    for line in reader:
        skip_this_line = False
        for scol, acceptable_values in skip_column_vals.items():
            if line[scol] not in acceptable_values:
                skip_this_line = True
                lines_skipped += 1
                counts_skipped += int(line['count'])
                break
        if skip_this_line:
            continue

        if line[args.non_summed_column] not in info:
            info[line[args.non_summed_column]] = 0
        info[line[args.non_summed_column]] += int(line['count'])
        lines_used += 1
        counts_used += int(line['count'])

print '  applied restrictions:'
for scol, acceptable_values in skip_column_vals.items():
    print '      %15s in %s' % (scol, acceptable_values)
print '   used:'
print '     %6d / %-6d = %.3f  lines'  % (lines_used, lines_used + lines_skipped, lines_used / float(lines_used + lines_skipped))
print '     %6d / %-6d = %.3f  counts'  % (counts_used, counts_used + counts_skipped, counts_used / float(counts_used + counts_skipped))

# print '  writing %d lines to %s' % (len(info), outfile)
# with open(outfile, 'w') as csvfile:
#     writer = csv.DictWriter(csvfile, (args.non_summed_column, 'count'))
#     writer.writeheader()
#     for val, count in sorted(info.items(), key=operator.itemgetter(1), reverse=True):
#         writer.writerow({args.non_summed_column : val, 'count' : count})
